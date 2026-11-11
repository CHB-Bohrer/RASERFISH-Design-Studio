"""Exact specificity counting with concise Standard and optional Full audit.

Standard retains completed per-candidate/source scores in ScoreStore. Only the
current candidate's interval geometry and risk flags occupy temporary SQLite;
sequence alignment payloads are not stored, and scratch is recreated each batch.
Full audit uses HitStore to retain losslessly compressed enriched alignments.
Both stores apply the same overlap-cluster and intended-site scientific rules.
"""
import json
import hashlib
import sqlite3
import zlib
from collections import OrderedDict
from pathlib import Path

from sequence_tools import canonical_ref


# Repeated HSPs commonly report the same subject thousands of times. This is
# only a bounded successful-insert cache: evicted or unseen subjects still go
# through SQLite's authoritative unique key, with no limit on distinct subjects.
SUBJECT_INSERT_CACHE_SIZE = 4096


def _record_subject(connection, cache, source, candidate, subject):
    key = (source, candidate, subject)
    if key in cache:
        return
    connection.execute('INSERT OR IGNORE INTO subjects VALUES(?,?,?)', key)
    # Cache only after SQLite accepts the insertion. SQL failure or rollback
    # must never turn a missing subject into an apparent cache hit.
    if len(cache) >= SUBJECT_INSERT_CACHE_SIZE:
        cache.popitem(last=False)
    cache[key] = None


# Format 1 uses a fixed vocabulary of field names and constant classification
# labels, never learned sequence data. These bytes are part of the file format:
# keep this dictionary unchanged and assign a new prefix for any future format.
_PAYLOAD_PREFIX=b'RFH\x01'
_PAYLOAD_DICTIONARY=(
    b'"RNA_SEQUENCE_ONLY" "GAPPED_DUPLEX_NOT_MODELED" "TOO_SHORT_OR_AMBIGUOUS" '
    b'"SUPPORTED_MISMATCH_DNA_NN" "UNSUPPORTED_MISMATCH_CONTEXT" '
    b'"source":"rna_human" "source":"targets" "source":"construct" '
    b'"bound_strand":"RNA" "bound_strand":"-" "allowed":true '
    b'"full_primary":true "rna":true "high_tm":true "moderate_tm":true '
    b'"unknown_long":true "strong_rna":true "rna_exact_risk":true '
    b'{"reference":"ref|NC_|","source":"human","start":0,"end":0,'
    b'"bound_strand":"+","qseq":"","sseq":"","alignment_length":0,'
    b'"query_start":1,"query_end":0,"query_coverage":100.0,"identity":100.0,'
    b'"mismatches":0,"gap_opens":0,"full_primary":false,"rna":false,'
    b'"evalue":0.0,"bitscore":0.0,"allowed":false,"exact_tract":0,'
    b'"offtarget_tm":null,"tm_model":"EXACT_SEGMENT_DNA_NN",'
    b'"exact_tract_tm":null,"risk_tm":null,"on_minus_off_tm":null,'
    b'"review":true,"high_tm":false,"moderate_tm":false,"unknown_long":false,'
    b'"strong_rna":false,"rna_exact_risk":false}'
)


def _encode_payload(hit):
    raw=json.dumps(hit,separators=(',',':')).encode('utf-8')
    compressor=zlib.compressobj(level=1,zdict=_PAYLOAD_DICTIONARY)
    return _PAYLOAD_PREFIX+compressor.compress(raw)+compressor.flush()


def _decode_payload(payload):
    if isinstance(payload,bytes):
        if payload.startswith(b'RFH'):
            if not payload.startswith(_PAYLOAD_PREFIX):
                raise ValueError('Unsupported retained-hit payload format')
            decompressor=zlib.decompressobj(zdict=_PAYLOAD_DICTIONARY)
            raw=decompressor.decompress(payload[len(_PAYLOAD_PREFIX):])+decompressor.flush()
            if not decompressor.eof or decompressor.unused_data:
                raise ValueError('Incomplete or invalid retained-hit payload')
            payload=raw
        else:
            payload=zlib.decompress(payload)  # Earlier unversioned compressed rows.
    return json.loads(payload)


def _hash_hit_row(digest,row):
    """Hash the exact indexed columns AND compressed payload, one row at a time."""
    columns=json.dumps(list(row[:13]),separators=(',',':'),allow_nan=False).encode('utf-8')
    payload=row[13].encode('utf-8') if isinstance(row[13],str) else bytes(row[13])
    digest.update(len(columns).to_bytes(8,'big'));digest.update(columns)
    digest.update(len(payload).to_bytes(8,'big'));digest.update(payload)


class HitStore:
    def __init__(self, path, *, resume=False, identity_digest=None):
        path=Path(path)
        if resume:
            if not path.is_file() or path.is_symlink():
                raise ValueError('Saved checkpoint database is missing or unsafe; analysis remains incomplete.')
            self.connection=sqlite3.connect(path.resolve().as_uri()+'?mode=rw',uri=True)
        else:
            self.connection = sqlite3.connect(str(path))
        self.connection.execute('PRAGMA cache_size=-4096')
        self.connection.execute('PRAGMA temp_store=FILE')
        self.connection.execute('PRAGMA journal_mode=DELETE')
        self.connection.execute('PRAGMA synchronous=FULL')
        self._pending_hash=hashlib.sha256()
        self._subject_pairs=set()
        self._subject_insert_cache=OrderedDict()
        if resume:
            try:
                checked=self.connection.execute('PRAGMA quick_check').fetchall()
                if checked!=[('ok',)]:raise ValueError('Saved checkpoint database integrity check failed.')
                if self.connection.execute('SELECT version FROM checkpoint_store').fetchall()!=[(1,)]:
                    raise ValueError('Unsupported saved checkpoint database schema.')
                if identity_digest is not None and self.connection.execute('SELECT digest FROM checkpoint_identity').fetchall()!=[(identity_digest,)]:
                    raise ValueError('Saved checkpoint database belongs to different design inputs or tools.')
                # Every durable marker is checked before any saved evidence is
                # trusted. Hashes detect damaged metadata; SQLite transactions
                # bind the marker to the corresponding hit/subject insertions.
                subject_total=0;subject_pairs=set();covered_ranges=[]
                for key,payload,digest in self.connection.execute('SELECT key,payload,sha256 FROM checkpoints'):
                    value=self._decode_checkpoint(key,payload,digest)
                    evidence=value.get('_store_evidence')
                    if not isinstance(evidence,dict):raise ValueError('Saved checkpoint lacks transactional evidence metadata.')
                    if evidence['through_rowid']<evidence['after_rowid']:
                        raise ValueError('Saved checkpoint has an invalid evidence range.')
                    if evidence['through_rowid']>evidence['after_rowid']:
                        covered_ranges.append((evidence['after_rowid'],evidence['through_rowid']))
                    found=0;hit_digest=hashlib.sha256()
                    for row in self.connection.execute('SELECT * FROM hits WHERE rowid>? AND rowid<=? ORDER BY rowid',
                            (evidence['after_rowid'],evidence['through_rowid'])):
                        _hash_hit_row(hit_digest,row);found+=1
                    if found!=evidence['hit_count'] or hit_digest.hexdigest()!=evidence['hit_sha256']:
                        raise ValueError('Saved checkpoint retained-hit count or checksum differs from its committed evidence.')
                    for saved in evidence['subjects']:
                        pair=(saved['source'],saved['candidate'])
                        if pair in subject_pairs:raise ValueError('Duplicated saved subject-checkpoint identity.')
                        subject_pairs.add(pair)
                        actual=self._subjects_digest(*pair)
                        if actual!=saved:raise ValueError('Saved checkpoint BLAST subject evidence differs from its committed checksum.')
                        subject_total+=actual['count']
                if self.connection.execute('SELECT COUNT(*) FROM subjects').fetchone()[0]!=subject_total:
                    raise ValueError('Saved BLAST subject evidence contains uncheckpointed records.')
                maximum=self.connection.execute('SELECT COALESCE(MAX(rowid),0) FROM hits').fetchone()[0]
                covered_through=0
                for begin,end in sorted(covered_ranges):
                    if begin!=covered_through:raise ValueError('Saved checkpoint evidence ranges have a gap or overlap.')
                    covered_through=end
                if covered_through!=maximum:raise ValueError('Saved retained-hit evidence lacks complete checkpoint coverage.')
                committed_max=max((self._decode_checkpoint(*row)['_store_evidence']['through_rowid']
                    for row in self.connection.execute('SELECT key,payload,sha256 FROM checkpoints')),default=0)
                if maximum!=committed_max:raise ValueError('Saved evidence contains data without a completed checkpoint.')
                self._checkpoint_highwater=maximum
            except BaseException:
                self.connection.close()
                raise
            return
        self.connection.execute('''CREATE TABLE hits (
            candidate TEXT, source TEXT, reference TEXT, strand TEXT,
            start INTEGER, end INTEGER, allowed INTEGER, review INTEGER,
            high INTEGER, moderate INTEGER, unknown INTEGER,
            strong_rna INTEGER, rna_exact INTEGER, payload BLOB)''')
        self.connection.execute('CREATE INDEX hit_order ON hits(candidate, source, reference, strand, start, end)')
        self.connection.execute('CREATE TABLE subjects (source TEXT, candidate TEXT, subject TEXT, PRIMARY KEY(source,candidate,subject)) WITHOUT ROWID')
        self.connection.execute('CREATE TABLE checkpoint_store (version INTEGER NOT NULL)')
        self.connection.execute('INSERT INTO checkpoint_store VALUES(1)')
        self.connection.execute('CREATE TABLE checkpoint_identity (digest TEXT NOT NULL)')
        if identity_digest is not None:self.connection.execute('INSERT INTO checkpoint_identity VALUES(?)',(identity_digest,))
        self.connection.execute('CREATE TABLE checkpoints (key TEXT PRIMARY KEY, payload TEXT NOT NULL, sha256 TEXT NOT NULL) WITHOUT ROWID')
        self.connection.commit()
        self._checkpoint_highwater=0

    @staticmethod
    def _decode_checkpoint(key,payload,digest):
        if len(payload)>16*1024*1024 or hashlib.sha256(payload.encode('utf-8')).hexdigest()!=digest:
            raise ValueError('Saved checkpoint metadata is damaged: '+str(key))
        value=json.loads(payload)
        if not isinstance(value,dict):raise ValueError('Invalid saved checkpoint metadata: '+str(key))
        return value

    def checkpoint(self,key):
        row=self.connection.execute('SELECT payload,sha256 FROM checkpoints WHERE key=?',(key,)).fetchone()
        return None if row is None else self._decode_checkpoint(key,*row)

    def commit_checkpoint(self,key,value):
        """Commit completed work and its authoritative marker together.

        Call only after the whole batch has parsed successfully. Pending hits,
        subjects and this record are a single SQLite transaction. A process
        interruption before commit leaves none of that batch reusable.
        """
        end=self.connection.execute('SELECT COALESCE(MAX(rowid),0) FROM hits').fetchone()[0]
        count=self.connection.execute('SELECT COUNT(*) FROM hits WHERE rowid>?',(self._checkpoint_highwater,)).fetchone()[0]
        value=dict(value,_store_evidence=dict(after_rowid=self._checkpoint_highwater,through_rowid=end,hit_count=count,
            hit_sha256=self._pending_hash.hexdigest(),subjects=[self._subjects_digest(*pair) for pair in sorted(self._subject_pairs)]))
        payload=json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)
        if len(payload)>16*1024*1024:raise ValueError('Checkpoint metadata is unexpectedly large.')
        digest=hashlib.sha256(payload.encode('utf-8')).hexdigest()
        self.connection.execute('INSERT INTO checkpoints VALUES(?,?,?)',(key,payload,digest))
        self.connection.commit()
        self._checkpoint_highwater=end
        self._pending_hash=hashlib.sha256();self._subject_pairs=set()
        self._subject_insert_cache.clear()

    def checkpoint_count(self,prefix):
        return self.connection.execute('SELECT COUNT(*) FROM checkpoints WHERE key LIKE ?',(prefix+'%',)).fetchone()[0]

    def rollback(self):
        self.connection.rollback()
        self._pending_hash=hashlib.sha256();self._subject_pairs=set()
        self._subject_insert_cache.clear()

    def _subjects_digest(self,source,candidate):
        digest=hashlib.sha256();count=0
        for (subject,) in self.connection.execute('SELECT subject FROM subjects WHERE source=? AND candidate=? ORDER BY subject',(source,candidate)):
            value=subject.encode('utf-8');digest.update(len(value).to_bytes(8,'big'));digest.update(value);count+=1
        return dict(source=source,candidate=candidate,count=count,sha256=digest.hexdigest())

    def hits(self, candidate_id):
        return DiskHits(self, candidate_id)

    def record_subject(self, source, candidate, subject):
        _record_subject(self.connection,self._subject_insert_cache,source,candidate,subject)
        self._subject_pairs.add((source,candidate))

    def capped_candidates(self, source, cap):
        return [row[0] for row in self.connection.execute(
            'SELECT candidate FROM subjects WHERE source=? GROUP BY candidate HAVING COUNT(*)>=?', (source,cap))]

    def flush(self):
        self.connection.commit()

    def close(self):
        self.connection.close()


class DiskHits:
    def __init__(self, store, candidate):
        self.store, self.candidate = store, candidate
        self.count = store.connection.execute('SELECT COUNT(*) FROM hits WHERE candidate=?',(candidate,)).fetchone()[0]

    def __len__(self):
        return self.count

    def append(self, hit):
        values=(
            self.candidate, hit['source'], canonical_ref(hit['reference']), hit['bound_strand'],
            hit['start'], hit['end'], int(bool(hit['allowed'])), int(bool(hit['review'])),
            int(bool(hit['high_tm'])), int(bool(hit['moderate_tm'])), int(bool(hit['unknown_long'])),
            int(bool(hit['strong_rna'])), int(bool(hit.get('rna_exact_risk',False))),
            _encode_payload(hit))
        self.store.connection.execute('INSERT INTO hits VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)',values)
        _hash_hit_row(self.store._pending_hash,values)
        self.count += 1

    def __iter__(self):
        cursor = self.store.connection.execute('SELECT payload FROM hits WHERE candidate=? ORDER BY rowid', (self.candidate,))
        for row in cursor:
            # Accept legacy TEXT when inspecting a mixed-format store. This is
            # a decoding compatibility path, not permission to resume old runs.
            yield _decode_payload(row[0])

    def specificity_counts(self, policy, source=None):
        """Exactly match cluster_hits + any(flag) without materializing HSP lists."""
        con = self.store.connection
        source_filter=' AND source=?' if source is not None else ''
        parameters=(self.candidate,source) if source is not None else (self.candidate,)
        allowed = con.execute('SELECT COUNT(*) FROM hits WHERE candidate=? AND allowed=1'+source_filter, parameters).fetchone()[0]
        ignored = con.execute("SELECT COUNT(*) FROM hits WHERE candidate=? AND allowed=0 AND review=1 AND source='human'"+source_filter, parameters).fetchone()[0] if policy.human_policy=='report_only' else 0
        predicate = source_filter+(" AND source!='human'" if policy.human_policy=='report_only' else '')
        rows = con.execute('''SELECT source,reference,strand,start,end,high,moderate,unknown,strong_rna,rna_exact
            FROM hits WHERE candidate=? AND allowed=0 AND review=1''' + predicate +
            ' ORDER BY source,reference,strand,start,end', parameters)
        totals = dict(review_site_clusters=0, high_tm_sites=0, moderate_tm_sites=0,
                      unknown_long_sites=0, strong_rna_sites=0, rna_risk_sites=0,
                      allowed_hsps=allowed, ignored_human_review_hsps=ignored)
        current = None
        flags = [False]*5
        end = 0

        def finish():
            if current is None:
                return
            totals['review_site_clusters'] += 1
            for name, value in zip(('high_tm_sites','moderate_tm_sites','unknown_long_sites','strong_rna_sites','rna_risk_sites'),flags):
                totals[name] += bool(value)

        for source,ref,strand,start,stop,high,moderate,unknown,rna,exact in rows:
            key = source,ref,strand
            values = high,moderate,unknown,rna,(rna or (policy.avoid_gene_rna and exact))
            if key!=current or start>end:
                finish()
                current,end,flags = key,stop,list(values)
            else:
                end=max(end,stop)
                flags=[a or b for a,b in zip(flags,values)]
        finish()
        return totals


# Standard output retains only completed QC counts. Per-hit geometry exists only
# in a current-batch work file, never in the durable score checkpoint database.
_SCORE_FIELDS=('review_site_clusters','high_tm_sites','moderate_tm_sites',
               'unknown_long_sites','strong_rna_sites','rna_risk_sites',
               'allowed_hsps','ignored_human_review_hsps')


class ScoreStore:
    """Small transactional checkpoints plus bounded current-batch interval data.

    Each candidate/source is summarized only after every alignment in BOTH its
    arm and full-primary queries has been consumed. No hit-count cutoff is used.
    The existing DiskHits overlap algorithm supplies the exact scientific counts.
    """
    def __init__(self,path,*,resume=False,identity_digest=None,deferred_local_sources=()):
        self.path=Path(path)
        if self.path.is_symlink() or (resume and not self.path.is_file()):
            raise ValueError('Saved score checkpoint database is missing or unsafe.')
        self.connection=sqlite3.connect(self.path.resolve().as_uri()+'?mode=rw',uri=True) if resume else sqlite3.connect(str(self.path))
        self.connection.execute('PRAGMA cache_size=-4096')
        self.connection.execute('PRAGMA temp_store=FILE')
        self.connection.execute('PRAGMA journal_mode=DELETE')
        self.connection.execute('PRAGMA synchronous=FULL')
        self.deferred_local_sources=set(deferred_local_sources)
        self._pending=[]
        self._work=None
        # This file has a separate lifecycle from committed scores: partial
        # intervals may be discarded after any interruption, never reused.
        self.work_path=self.path.with_name(self.path.stem+'.working.sqlite')
        try:
            if resume:
                if self.connection.execute('PRAGMA quick_check').fetchall()!=[('ok',)]:
                    raise ValueError('Saved score checkpoint integrity check failed.')
                if self.connection.execute('SELECT version FROM checkpoint_store').fetchall()!=[(2,)]:
                    raise ValueError('Unsupported score checkpoint schema; raw-hit checkpoints cannot be reused.')
                if identity_digest is not None and self.connection.execute('SELECT digest FROM checkpoint_identity').fetchall()!=[(identity_digest,)]:
                    raise ValueError('Saved score checkpoint belongs to different design inputs or tools.')
                covered=set()
                for key,payload,digest in self.connection.execute('SELECT key,payload,sha256 FROM checkpoints'):
                    value=HitStore._decode_checkpoint(key,payload,digest)
                    for evidence in value.get('_score_evidence',[]):
                        ident=(evidence['scope'],evidence['candidate'],evidence['source'])
                        if ident in covered:raise ValueError('Saved score evidence has duplicated completion coverage.')
                        covered.add(ident)
                        row=self.connection.execute('SELECT payload,sha256 FROM scores WHERE scope=? AND candidate=? AND source=?',ident).fetchone()
                        if row is None or row[1]!=evidence['sha256']:
                            raise ValueError('Saved score evidence is missing or its checksum differs.')
                        self._decode_score(ident,*row)
                    if '_score_evidence' not in value:raise ValueError('Saved score checkpoint lacks committed evidence metadata.')
                actual=set(self.connection.execute('SELECT scope,candidate,source FROM scores'))
                if covered!=actual:raise ValueError('Saved scores lack complete checkpoint coverage.')
            else:
                self.connection.execute('CREATE TABLE checkpoint_store(version INTEGER NOT NULL)')
                self.connection.execute('INSERT INTO checkpoint_store VALUES(2)')
                self.connection.execute('CREATE TABLE checkpoint_identity(digest TEXT NOT NULL)')
                if identity_digest is not None:self.connection.execute('INSERT INTO checkpoint_identity VALUES(?)',(identity_digest,))
                self.connection.execute('CREATE TABLE checkpoints(key TEXT PRIMARY KEY,payload TEXT NOT NULL,sha256 TEXT NOT NULL) WITHOUT ROWID')
                self.connection.execute('CREATE TABLE scores(scope TEXT,candidate TEXT,source TEXT,payload TEXT NOT NULL,sha256 TEXT NOT NULL,PRIMARY KEY(scope,candidate,source)) WITHOUT ROWID')
                self.connection.commit()
            self._reset_work()
        except BaseException:
            self.connection.close()
            raise

    def _reset_work(self):
        self._subject_insert_cache=OrderedDict()
        if self._work is not None:self._work.close()
        self._work=None
        for path in (self.work_path,Path(str(self.work_path)+'-journal'),Path(str(self.work_path)+'-wal'),Path(str(self.work_path)+'-shm')):
            if path.is_symlink():raise ValueError('Current-batch score work file must not be a symbolic link.')
            if path.exists():path.unlink()
        self._work=sqlite3.connect(str(self.work_path))
        self._work.execute('PRAGMA cache_size=-4096')
        self._work.execute('PRAGMA temp_store=FILE')
        # This database is explicitly scratch, never checkpoint evidence.
        self._work.execute('PRAGMA journal_mode=OFF')
        self._work.execute('PRAGMA synchronous=OFF')
        self._work.execute('CREATE TABLE hits(candidate TEXT,source TEXT,reference TEXT,strand TEXT,start INTEGER,end INTEGER,allowed INTEGER,review INTEGER,high INTEGER,moderate INTEGER,unknown INTEGER,strong_rna INTEGER,rna_exact INTEGER,expected INTEGER)')
        self._work.execute('CREATE TABLE subjects(source TEXT,candidate TEXT,subject TEXT,PRIMARY KEY(source,candidate,subject)) WITHOUT ROWID')

    @staticmethod
    def _decode_score(identity,payload,digest):
        value=HitStore._decode_checkpoint(str(identity),payload,digest)
        counts=value.get('counts',{})
        if set(counts)!=set(_SCORE_FIELDS) or any(type(n) is not int or n<0 for n in counts.values()):
            raise ValueError('Saved specificity score has invalid counts.')
        if type(value.get('alignment_count')) is not int or value['alignment_count']<0:
            raise ValueError('Saved specificity score has invalid alignment count.')
        if type(value.get('expected_target')) is not bool or type(value.get('subject_count')) is not int or value['subject_count']<0:
            raise ValueError('Saved specificity score has invalid intended-site or subject metadata.')
        return value

    def checkpoint(self,key):
        row=self.connection.execute('SELECT payload,sha256 FROM checkpoints WHERE key=?',(key,)).fetchone()
        return None if row is None else HitStore._decode_checkpoint(key,*row)

    def commit_checkpoint(self,key,value):
        value=dict(value,_score_evidence=list(self._pending),evidence_format='EXACT_COUNTS_NO_ALIGNMENT_PAYLOADS')
        payload=json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)
        if len(payload)>16*1024*1024:raise ValueError('Score checkpoint metadata is unexpectedly large.')
        self.connection.execute('INSERT INTO checkpoints VALUES(?,?,?)',(key,payload,hashlib.sha256(payload.encode()).hexdigest()))
        self.connection.commit()
        self._pending=[]

    def checkpoint_count(self,prefix):
        return self.connection.execute('SELECT COUNT(*) FROM checkpoints WHERE key LIKE ?',(prefix+'%',)).fetchone()[0]

    def record_subject(self,source,candidate,subject):
        _record_subject(self._work,self._subject_insert_cache,source,candidate,subject)

    def restore_local_source(self,source,candidate):
        row=self.connection.execute("SELECT payload,sha256 FROM scores WHERE scope='local' AND candidate=? AND source=?",(candidate,source)).fetchone()
        if row:
            value=self._decode_score(('local',candidate,source),*row)
            if 'deferred_intervals' not in value:
                raise ValueError('Overlapping local/database sources require verified minimal interval evidence.')
            self._work.executemany('INSERT INTO hits VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)',value['deferred_intervals'])

    def stage_scores(self,scope,candidate_ids,policy,source=None):
        """Reduce completed current-batch intervals before committing a marker.

        Scratch is reset immediately after reduction, so disk use cannot grow
        with the number of processed candidates. Only scores enter the durable
        transaction. An unusual API source collision keeps minimal local rows
        for that source, allowing its later BLAST union to be counted exactly.
        """
        from types import SimpleNamespace
        # No interval query is needed during ingestion. Build the same index
        # once from the completed batch rather than updating its B-tree for
        # every HSP. SQLite's file-backed sort and bounded page cache remain in
        # force; the owning worker supervisor also covers this native SQL step.
        self._work.execute('CREATE INDEX hit_order ON hits(candidate,source,reference,strand,start,end)')
        for candidate in candidate_ids:
            sources=[source] if source is not None else [r[0] for r in self._work.execute('SELECT DISTINCT source FROM hits WHERE candidate=?',(candidate,))]
            for current_source in sources:
                # Use the unchanged cluster union algorithm with a SQL source
                # predicate; no copying or materializing an HSP list is needed.
                counts=DiskHits(SimpleNamespace(connection=self._work),candidate).specificity_counts(policy,current_source)
                alignment_count,expected=self._work.execute('SELECT COUNT(*),COALESCE(MAX(expected),0) FROM hits WHERE candidate=? AND source=?',(candidate,current_source)).fetchone()
                subjects=self._work.execute('SELECT COUNT(*) FROM subjects WHERE source=? AND candidate=?',(current_source,candidate)).fetchone()[0]
                value=dict(counts=counts,alignment_count=alignment_count,expected_target=bool(expected),subject_count=subjects)
                if scope=='local' and current_source in self.deferred_local_sources:
                    intervals=[];interval_bytes=0
                    for row in self._work.execute('SELECT * FROM hits WHERE candidate=? AND source=? ORDER BY rowid',(candidate,current_source)):
                        item=list(row);interval_bytes+=len(json.dumps(item,separators=(',',':')).encode('utf-8'))+1
                        if interval_bytes>8*1024*1024:
                            raise ValueError('Local same-source interval evidence exceeds the bounded score-checkpoint limit; analysis incomplete.')
                        intervals.append(item)
                    value['deferred_intervals']=intervals
                payload=json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)
                if len(payload)>16*1024*1024:raise ValueError('Local same-source interval evidence exceeds the bounded score-checkpoint limit.')
                digest=hashlib.sha256(payload.encode()).hexdigest()
                self.connection.execute('INSERT INTO scores VALUES(?,?,?,?,?)',(scope,candidate,current_source,payload,digest))
                self._pending.append(dict(scope=scope,candidate=candidate,source=current_source,sha256=digest))
        self._reset_work()

    def _components(self,candidate):
        result={}
        for scope,source,payload,digest in self.connection.execute('SELECT scope,source,payload,sha256 FROM scores WHERE candidate=? ORDER BY scope',(candidate,)):
            value=self._decode_score((scope,candidate,source),payload,digest)
            # BLAST for a locally named source includes its deferred local
            # intervals. Its union replaces that local source's subtotal.
            if source not in result or scope!='local':result[source]=value
        return result.values()

    def capped_candidates(self,source,cap):
        found=[]
        for candidate,payload,digest in self.connection.execute("SELECT candidate,payload,sha256 FROM scores WHERE source=? AND scope!='local'",(source,)):
            if self._decode_score((candidate,source),payload,digest)['subject_count']>=cap:found.append(candidate)
        return sorted(set(found))

    def hits(self,candidate):return ScoreHits(self,candidate)
    def flush(self):
        self.connection.commit();self._pending=[]
    def rollback(self):
        self.connection.rollback();self._pending=[];self._reset_work()
    def close(self):
        try:
            if self._work is not None:self._work.close()
            for path in (self.work_path,Path(str(self.work_path)+'-journal')):
                if path.is_file() and not path.is_symlink():path.unlink()
        finally:self.connection.close()


class ScoreHits(DiskHits):
    """DiskHits-compatible scientific summaries; raw hits cannot be iterated."""
    def __init__(self,store,candidate):self.store,self.candidate=store,candidate
    def append(self,hit):
        expected=hit['source']=='human' and hit['allowed'] and not hit['full_primary'] and hit['identity']==100 and hit['query_coverage']==100
        self.store._work.execute('INSERT INTO hits VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(
            self.candidate,hit['source'],canonical_ref(hit['reference']),hit['bound_strand'],hit['start'],hit['end'],
            int(bool(hit['allowed'])),int(bool(hit['review'])),int(bool(hit['high_tm'])),int(bool(hit['moderate_tm'])),
            int(bool(hit['unknown_long'])),int(bool(hit['strong_rna'])),int(bool(hit.get('rna_exact_risk',False))),int(expected)))
    def __len__(self):return sum(v['alignment_count'] for v in self.store._components(self.candidate))
    def __iter__(self):
        raise RuntimeError('Standard output retains exact QC scores, not individual alignments. Choose Full audit before a new run for alignment detail.')
    def specificity_counts(self,policy):
        totals=dict.fromkeys(_SCORE_FIELDS,0)
        for component in self.store._components(self.candidate):
            for key in totals:totals[key]+=component['counts'][key]
        return totals
    def has_expected_genomic_target(self):return any(v['expected_target'] for v in self.store._components(self.candidate))
