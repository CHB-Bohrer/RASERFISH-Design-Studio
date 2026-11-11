"""Exact distinct-subject counting despite bounded duplicate-insert caching."""
import random
import sqlite3
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from specificity_store import HitStore, ScoreStore, DiskHits, SUBJECT_INSERT_CACHE_SIZE


POLICY = SimpleNamespace(human_policy='intended_only', avoid_gene_rna=False)


class SubjectInsertCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.stores = []

    def tearDown(self):
        for store in reversed(self.stores):
            store.close()
        self.temporary.cleanup()

    def store(self, kind):
        store = kind(self.root / (kind.__name__ + str(len(self.stores)) + '.sqlite'),
                     identity_digest='test-subjects')
        self.stores.append(store)
        return store

    @staticmethod
    def connection(store):
        return store._work if isinstance(store, ScoreStore) else store.connection

    def rows(self, store):
        return self.connection(store).execute(
            'SELECT source,candidate,subject FROM subjects ORDER BY source,candidate,subject').fetchall()

    def test_repeated_hsps_skip_only_duplicate_subject_inserts(self):
        for kind in (HitStore, ScoreStore):
            with self.subTest(store=kind.__name__):
                store = self.store(kind)
                statements = []
                self.connection(store).set_trace_callback(statements.append)
                keys = [('human', 'C1', 'chr1'), ('human', 'C2', 'chr1'),
                        ('rna_human', 'C1', 'chr1'), ('human', 'C1', 'chr2')]
                for _ in range(1000):
                    for key in keys:
                        store.record_subject(*key)
                self.assertEqual(self.rows(store), sorted(keys))
                inserts = [s for s in statements if s.startswith('INSERT OR IGNORE INTO subjects')]
                self.assertEqual(len(inserts), 4)
                if isinstance(store, HitStore):
                    self.assertEqual(store._subject_pairs, {key[:2] for key in keys})

    def test_cache_eviction_does_not_cap_or_lose_distinct_subjects(self):
        count = SUBJECT_INSERT_CACHE_SIZE + 17
        for kind in (HitStore, ScoreStore):
            with self.subTest(store=kind.__name__):
                store = self.store(kind)
                for n in range(count):
                    store.record_subject('human', 'C1', 'subject_' + str(n))
                self.assertEqual(len(store._subject_insert_cache), SUBJECT_INSERT_CACHE_SIZE)
                self.assertEqual(len(self.rows(store)), count)
                # The oldest key has been evicted; repeating it still retains
                # exactly one SQL row rather than dropping it or double counting.
                self.assertNotIn(('human', 'C1', 'subject_0'), store._subject_insert_cache)
                store.record_subject('human', 'C1', 'subject_0')
                self.assertEqual(len(self.rows(store)), count)
                if isinstance(store, ScoreStore):
                    store.stage_scores('blast:human:1', ['C1'], POLICY, 'human')
                    store.commit_checkpoint('blast:human:1', {'status': 'COMPLETE'})
                self.assertEqual(store.capped_candidates('human', count), ['C1'])
                self.assertEqual(store.capped_candidates('human', count + 1), [])

    def test_rollback_cannot_leave_cached_uncommitted_subjects(self):
        for kind in (HitStore, ScoreStore):
            with self.subTest(store=kind.__name__):
                store = self.store(kind)
                store.record_subject('human', 'C1', 'chr1')
                store.record_subject('human', 'C1', 'chr1')
                store.rollback()
                self.assertEqual(len(store._subject_insert_cache), 0)
                self.assertEqual(self.rows(store), [])
                store.record_subject('human', 'C1', 'chr1')
                self.assertEqual(self.rows(store), [('human', 'C1', 'chr1')])

    def test_score_scratch_reset_reinserts_and_preserves_counts(self):
        store = self.store(ScoreStore)
        for n in range(2):
            store.record_subject('human', 'C1', 'chr1')
            store.stage_scores('batch_' + str(n), ['C1'], POLICY, 'human')
            self.assertEqual(len(store._subject_insert_cache), 0)
            self.assertEqual(self.rows(store), [])
            store.commit_checkpoint('batch_' + str(n), {'status': 'COMPLETE'})
            row = store.connection.execute(
                'SELECT payload,sha256 FROM scores WHERE scope=?', ('batch_' + str(n),)).fetchone()
            score = store._decode_score(('batch_' + str(n), 'C1', 'human'), *row)
            self.assertEqual(score['subject_count'], 1)

    def test_failed_sql_does_not_poison_cache(self):
        for kind in (HitStore, ScoreStore):
            with self.subTest(store=kind.__name__):
                store = self.store(kind)
                con = self.connection(store)
                con.execute('PRAGMA query_only=ON')
                with self.assertRaises(sqlite3.OperationalError):
                    store.record_subject('human', 'C1', 'chr1')
                self.assertEqual(len(store._subject_insert_cache), 0)
                con.execute('PRAGMA query_only=OFF')
                store.record_subject('human', 'C1', 'chr1')
                self.assertEqual(self.rows(store), [('human', 'C1', 'chr1')])

    def test_random_repeats_match_uncached_database_exactly(self):
        rng = random.Random(20261102)
        keys = [(rng.choice(('human', 'rna_human')), 'C' + str(rng.randrange(3)),
                 'subject_' + str(rng.randrange(1300))) for _ in range(14000)]
        keys += keys[::3]
        for kind in (HitStore, ScoreStore):
            with self.subTest(store=kind.__name__):
                store = self.store(kind)
                baseline = sqlite3.connect(':memory:')
                try:
                    baseline.execute('CREATE TABLE subjects(source TEXT,candidate TEXT,subject TEXT,PRIMARY KEY(source,candidate,subject)) WITHOUT ROWID')
                    for key in keys:
                        baseline.execute('INSERT OR IGNORE INTO subjects VALUES(?,?,?)', key)
                        store.record_subject(*key)
                    expected = baseline.execute('SELECT * FROM subjects ORDER BY source,candidate,subject').fetchall()
                    self.assertEqual(self.rows(store), expected)
                    self.assertLessEqual(len(store._subject_insert_cache), SUBJECT_INSERT_CACHE_SIZE)
                finally:
                    baseline.close()

    def test_checkpoint_subject_evidence_stays_valid_after_reopen(self):
        store = self.store(HitStore)
        for _ in range(100):
            store.record_subject('human', 'C1', 'chr1')
            store.record_subject('human', 'C1', 'chr2')
        store.commit_checkpoint('blast:human:1', {'status': 'COMPLETE'})
        self.assertEqual(len(store._subject_insert_cache), 0)
        expected = store.checkpoint('blast:human:1')['_store_evidence']['subjects']
        path = self.root / 'HitStore0.sqlite'
        store.close()
        self.stores.remove(store)
        resumed = HitStore(path, resume=True, identity_digest='test-subjects')
        self.stores.append(resumed)
        self.assertEqual(resumed.checkpoint('blast:human:1')['_store_evidence']['subjects'], expected)
        self.assertEqual(expected[0]['count'], 2)
        self.assertEqual(resumed.capped_candidates('human', 2), ['C1'])

    def test_bulk_scratch_index_preserves_unsorted_cluster_counts_and_resets(self):
        store = self.store(ScoreStore)
        # Compare the production reduction with the old incremental-index
        # algorithm, including bridging overlaps, flags, allowed rows, two
        # candidates, and distinct strands/sources.
        reference = sqlite3.connect(':memory:')
        try:
            reference.execute('CREATE TABLE hits(candidate TEXT,source TEXT,reference TEXT,strand TEXT,start INTEGER,end INTEGER,allowed INTEGER,review INTEGER,high INTEGER,moderate INTEGER,unknown INTEGER,strong_rna INTEGER,rna_exact INTEGER,expected INTEGER)')
            reference.execute('CREATE INDEX hit_order ON hits(candidate,source,reference,strand,start,end)')
            rows = [
                ('C1', 'human', 'chr1', '+', 150, 180, 0, 1, 0, 1, 0, 0, 0, 0),
                ('C1', 'human', 'chr1', '+', 100, 130, 0, 1, 1, 1, 0, 0, 0, 0),
                ('C1', 'human', 'chr1', '+', 120, 170, 0, 1, 0, 0, 1, 0, 0, 0),
                ('C1', 'human', 'chr1', '-', 120, 170, 0, 1, 0, 0, 0, 0, 0, 0),
                ('C1', 'human', 'chr1', '+', 200, 240, 1, 1, 0, 0, 0, 0, 0, 1),
                ('C1', 'rna_human', 'RNA1', 'RNA', 100, 170, 0, 1, 0, 0, 0, 1, 1, 0),
                ('C2', 'human', 'chr2', '+', 100, 170, 0, 1, 1, 1, 0, 0, 0, 0),
            ]
            reference.executemany('INSERT INTO hits VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)', rows)
            expected = {(candidate, source): DiskHits(SimpleNamespace(connection=reference), candidate).specificity_counts(POLICY, source)
                        for candidate, source in (('C1', 'human'), ('C1', 'rna_human'), ('C2', 'human'))}
            for batch in range(2):
                self.assertEqual(store._work.execute("SELECT COUNT(*) FROM sqlite_master WHERE name='hit_order'").fetchone()[0], 0)
                self.assertEqual(store._work.execute('PRAGMA cache_size').fetchone()[0], -4096)
                self.assertEqual(store._work.execute('PRAGMA temp_store').fetchone()[0], 1)
                store._work.executemany('INSERT INTO hits VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)', rows)
                statements = []
                store._work.set_trace_callback(statements.append)
                scope = 'batch_' + str(batch)
                store.stage_scores(scope, ['C1', 'C2'], POLICY)
                creates = [sql for sql in statements if sql.startswith('CREATE INDEX hit_order')]
                self.assertEqual(len(creates), 1)
                for candidate, source in expected:
                    saved = store.connection.execute('SELECT payload,sha256 FROM scores WHERE scope=? AND candidate=? AND source=?', (scope, candidate, source)).fetchone()
                    value = store._decode_score((scope, candidate, source), *saved)
                    self.assertEqual(value['counts'], expected[candidate, source])
                store.commit_checkpoint(scope, {'status': 'COMPLETE'})
            self.assertEqual(expected['C1', 'human']['review_site_clusters'], 2)
            self.assertEqual(expected['C1', 'human']['high_tm_sites'], 1)
        finally:
            reference.close()


if __name__ == '__main__':
    unittest.main()
