"""Exports are streamed without reading whole result files into memory."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import studio

class CountingSink:
    def __init__(self):self.total=0;self.maximum=0;self.chunks=0
    def write(self,data):
        self.total+=len(data);self.maximum=max(self.maximum,len(data));self.chunks+=1

class StreamingDownloadTests(unittest.TestCase):
    def test_completed_zip_streams_in_bounded_chunks_and_never_read_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);archive=root/'fixture.zip'
            with archive.open('wb') as out:
                for _ in range(5):out.write(b'A'*(1024*1024))
                out.write(b'end')
            # Downloads require durable completion, including after server restart.
            run=root/'fixture';(run/'result').mkdir(parents=True)
            (run/'analysis_status.json').write_text(json.dumps(dict(status='complete',analysis_complete=True,job_id='fixture')))
            (run/'job_status.json').write_text(json.dumps(dict(status='complete',id='fixture')))
            (run/'result/design_result.json').write_text(json.dumps(dict(status='DRAFT',counts={'selected':1})))
            handler=studio.Handler.__new__(studio.Handler)
            handler.path='/download?id=fixture';handler.headers={'Host':'127.0.0.1'}
            status=[];headers={};handler.wfile=CountingSink()
            handler.send_response=status.append
            handler.send_header=lambda key,value:headers.update({key:value})
            handler.end_headers=lambda:None
            with patch.object(studio,'JOBS',root),patch.object(studio,'STATE',{'fixture':{'status':'complete'}}),patch.object(Path,'read_bytes',side_effect=AssertionError('Unbounded file read prohibited')):
                handler.do_GET()
            self.assertEqual(status,[200])
            self.assertEqual(int(headers['Content-Length']),archive.stat().st_size)
            self.assertEqual(handler.wfile.total,archive.stat().st_size)
            self.assertLessEqual(handler.wfile.maximum,1024*1024)
            self.assertEqual(handler.wfile.chunks,6)

    def test_incomplete_download_is_blocked_before_stream_open(self):
        handler=studio.Handler.__new__(studio.Handler)
        handler.path='/download?id=fixture';handler.headers={'Host':'127.0.0.1'}
        sent=[];handler.send=lambda data,status=200,**kw:sent.append((data,status))
        with patch.object(studio,'STATE',{'fixture':{'status':'incomplete'}}),patch.object(studio.Handler,'send_file') as stream:
            handler.do_GET()
        stream.assert_not_called()
        self.assertEqual(sent[-1][1],400)
        self.assertIn('incomplete',sent[-1][0]['error'])

if __name__=='__main__':unittest.main()
