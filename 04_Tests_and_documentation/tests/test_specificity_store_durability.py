"""Lossless bounded storage and recovery of committed tiny fixture evidence."""
import copy
import hashlib
import json
import signal
import sqlite3
import subprocess
import sys
import tempfile
import time
import tracemalloc
import unittest
import zlib
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from design_engine import Policy, summarize_specificity
from specificity_store import HitStore, _encode_payload, _decode_payload, _PAYLOAD_DICTIONARY


def hit(**changes):
    row=dict(source='human',reference='chr1',bound_strand='+',start=100,end=139,
             allowed=False,review=True,high_tm=False,moderate_tm=True,
             unknown_long=False,strong_rna=False,rna_exact_risk=False,
             qseq='ACGT'*10,sseq='ACGT'*10,identity=100.0,query_coverage=100.0,
             query_start=1,query_end=40,alignment_length=40,mismatches=0,
             gap_opens=0,full_primary=False,rna=False,evalue=1e-12,bitscore=80.0,
             duplex_tm=48.123456789,notes=['alignment fixture',None,True,'µ'])
    row.update(changes)
    return row


class DurableSpecificityTests(unittest.TestCase):
    def test_format_one_dictionary_remains_backward_readable(self):
        # Changing a dictionary under an existing prefix would make stored
        # evidence unreadable. New dictionaries require a new format version.
        self.assertEqual(hashlib.sha256(_PAYLOAD_DICTIONARY).hexdigest(),
                         'b4a3ad92ff7bd66593f3c986181ef2b448f38504e928aede551df80e00702c04')

    def test_compressed_roundtrip_preserves_every_hit_and_insertion_order(self):
        rows=[hit(start=i*100,end=i*100+39,high_tm=bool(i%3),label=str(i)) for i in range(400)]
        with tempfile.TemporaryDirectory() as tmp:
            store=HitStore(Path(tmp)/'hits.sqlite')
            try:
                disk=store.hits('candidate')
                for row in rows:disk.append(row)
                store.flush()
                self.assertEqual(len(disk),len(rows))
                self.assertEqual(list(disk),rows)
                self.assertEqual(store.connection.execute('SELECT DISTINCT typeof(payload) FROM hits').fetchall(),[('blob',)])
                self.assertEqual(store.connection.execute('PRAGMA journal_mode').fetchone()[0],'delete')
                self.assertEqual(store.connection.execute('PRAGMA synchronous').fetchone()[0],2)
                self.assertEqual(store.connection.execute('PRAGMA cache_size').fetchone()[0],-4096)
                self.assertEqual(store.connection.execute('PRAGMA temp_store').fetchone()[0],1)
                raw=sum(len(json.dumps(row,separators=(',',':')).encode()) for row in rows)
                compressed=store.connection.execute('SELECT SUM(LENGTH(payload)) FROM hits').fetchone()[0]
                self.assertLess(compressed,raw*.8)
            finally:store.close()

    def test_legacy_text_and_compressed_payloads_decode_in_original_order(self):
        rows=[hit(start=1,end=40),hit(start=99,end=138),hit(start=222,end=261)]
        with tempfile.TemporaryDirectory() as tmp:
            store=HitStore(Path(tmp)/'hits.sqlite')
            try:
                disk=store.hits('mixed')
                for row in rows:disk.append(row)
                store.connection.execute('UPDATE hits SET payload=? WHERE rowid=2',(json.dumps(rows[1]),))
                store.connection.execute('UPDATE hits SET payload=? WHERE rowid=3',(zlib.compress(json.dumps(rows[2]).encode()),))
                store.flush()
                self.assertEqual(list(disk),rows)
                self.assertEqual(len(disk),3)
            finally:store.close()

    def test_versioned_payload_rejects_unknown_format_and_truncation(self):
        payload=_encode_payload(hit())
        with self.assertRaisesRegex(ValueError,'Unsupported'):_decode_payload(b'RFH\x02'+payload[4:])
        with self.assertRaises((ValueError,zlib.error)):_decode_payload(payload[:-3])
        with self.assertRaisesRegex(ValueError,'invalid'):_decode_payload(payload+b'extra')

    def test_all_classification_flags_cluster_counts_and_rejections_match_lists(self):
        rows=[hit(allowed=True),hit(start=1000,end=1040,high_tm=True),
              hit(start=1030,end=1070,unknown_long=True),
              hit(start=2000,end=2040,reference='ref|NC_000001.11|',high_tm=True),
              hit(start=1000,end=1040,bound_strand='-',moderate_tm=False),
              hit(source='rna_human',reference='NM_FIXTURE.1',bound_strand='RNA',strong_rna=True),
              hit(source='rna_local',reference='TX1',bound_strand='RNA',rna_exact_risk=True),
              hit(source='construct',reference='construct',moderate_tm=True),
              hit(start=3000,end=3040,review=False,high_tm=True)]
        with tempfile.TemporaryDirectory() as tmp:
            store=HitStore(Path(tmp)/'hits.sqlite')
            try:
                disk=store.hits('candidate')
                for row in rows:disk.append(row)
                for human in ('intended_only','report_only'):
                    for rna in (False,True):
                        policy=Policy(human_policy=human,avoid_gene_rna=rna)
                        baseline=dict(offtarget_hits=copy.deepcopy(rows),reject_reasons=[])
                        tested=dict(offtarget_hits=disk,reject_reasons=[])
                        summarize_specificity(baseline,policy)
                        summarize_specificity(tested,policy)
                        self.assertEqual(baseline,tested|{'offtarget_hits':rows})
                self.assertEqual(len(disk),len(rows))
            finally:store.close()

    def test_streaming_payload_encoding_and_decoding_do_not_accumulate_hits(self):
        with tempfile.TemporaryDirectory() as tmp:
            store=HitStore(Path(tmp)/'hits.sqlite')
            try:
                disk=store.hits('candidate')
                tracemalloc.start()
                baseline=tracemalloc.get_traced_memory()[0]
                for index in range(5000):disk.append(hit(index=index))
                store.flush()
                count=0
                for row in disk:
                    self.assertEqual(row['index'],count)
                    count+=1
                peak=tracemalloc.get_traced_memory()[1]-baseline
                tracemalloc.stop()
                self.assertEqual(count,5000)
                self.assertEqual(len(disk),5000)
                self.assertLess(peak,2*1024*1024)
            finally:
                if tracemalloc.is_tracing():tracemalloc.stop()
                store.close()

    @unittest.skipUnless(hasattr(signal,'SIGKILL'),'Requires an owned-process hard-stop test')
    def test_sigkill_mid_transaction_recovers_only_committed_rows(self):
        # The child owns only a tiny temporary fixture. Force page spills to
        # exercise a hot rollback journal, rather than killing before writes.
        child_code='''
import json,sys,time
from pathlib import Path
sys.path.insert(0,sys.argv[1])
from specificity_store import HitStore
path=Path(sys.argv[2]);ready=Path(sys.argv[3]);row=json.loads(sys.argv[4])
store=HitStore(path);disk=store.hits('committed');disk.append(row)
store.record_subject('human','committed','chr1');store.flush()
store.connection.execute('PRAGMA cache_size=-64')
store.connection.execute('UPDATE hits SET start=999999 WHERE candidate=?',('committed',))
pending=store.hits('uncommitted')
for index in range(2000):
    pending.append(dict(row,start=index*100,end=index*100+39,index=index))
    store.record_subject('human','uncommitted','subject'+str(index))
ready.write_text('transaction is uncommitted')
while True:time.sleep(.1)
'''
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path=root/'hits.sqlite';ready=root/'ready';row=hit()
            proc=subprocess.Popen([sys.executable,'-c',child_code,str(Path(__file__).resolve().parents[1]),
                                   str(path),str(ready),json.dumps(row)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            try:
                deadline=time.monotonic()+15
                while not ready.exists() and proc.poll() is None and time.monotonic()<deadline:time.sleep(.02)
                self.assertTrue(ready.exists(),'Fixture child did not finish its uncommitted writes')
                self.assertGreater(Path(str(path)+'-journal').stat().st_size,512)
                proc.send_signal(signal.SIGKILL);proc.wait(timeout=5)
                self.assertEqual(proc.returncode,-signal.SIGKILL)
                # Normal read/write opening performs SQLite's journal recovery.
                con=sqlite3.connect(str(path))
                try:
                    self.assertEqual(con.execute('PRAGMA integrity_check').fetchone(),('ok',))
                    self.assertEqual(con.execute('SELECT candidate,start FROM hits').fetchall(),[('committed',row['start'])])
                    self.assertEqual(con.execute('SELECT * FROM subjects').fetchall(),[('human','committed','chr1')])
                    payload=con.execute('SELECT payload FROM hits').fetchone()[0]
                    self.assertEqual(_decode_payload(payload),row)
                finally:con.close()
            finally:
                if proc.poll() is None:proc.kill();proc.wait(timeout=5)
                proc.communicate(timeout=5)


if __name__=='__main__':unittest.main()
