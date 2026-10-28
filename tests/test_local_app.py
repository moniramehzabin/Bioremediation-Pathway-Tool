"""Integration checks for app input, job handling and localhost protection."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
import urllib.error
import urllib.request
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
spec = importlib.util.spec_from_file_location('bpt_app', Path(__file__).resolve().parents[1] / 'app.py')
app = importlib.util.module_from_spec(spec)
spec.loader.exec_module(app)

GENOME = b'LOCUS       SAMPLE\nFEATURES             Location/Qualifiers\n'


class AppTests(unittest.TestCase):
    def test_input_validation(self):
        with self.assertRaisesRegex(ValueError, 'email'):
            app.validate_inputs({'mode': 'web', 'email': 'bad'}, {'genome': ('a.gbk', GENOME)})
        with self.assertRaisesRegex(ValueError, 'GenBank'):
            app.validate_inputs({'mode': 'off'}, {'genome': ('a.fasta', GENOME)})
        with self.assertRaisesRegex(ValueError, 'TSV'):
            app.validate_inputs({'mode': 'precomputed'}, {'genome': ('a.gbk', GENOME)})
        self.assertEqual(app.validate_inputs({'mode': 'off'}, {'genome': ('a.gbk', GENOME)}), ('off', ''))

    def test_multipart_preserves_bytes(self):
        body = (b'--BOUND\r\nContent-Disposition: form-data; name="genome"; filename="x.gbk"\r\n'
                b'Content-Type: application/octet-stream\r\n\r\n' + GENOME + b'\r\n--BOUND--\r\n')
        fields, files = app.parse_form('multipart/form-data; boundary=BOUND', body)
        self.assertEqual(files['genome'], ('x.gbk', GENOME))

    def test_job_runs_and_preserves_evidence_warning(self):
        with tempfile.TemporaryDirectory() as folder:
            runner = Path(folder) / 'fake_engine.py'
            runner.write_text('''import sys,time
from pathlib import Path
prefix=Path(sys.argv[sys.argv.index('--out-prefix')+1])
print('engine started',flush=True)
time.sleep(.2)
Path(str(prefix)+'_VISUAL_REPORT.html').write_text('<h1>Report</h1>')
Path(str(prefix)+'_rescue_interpro_failed.tsv').write_text('locus\\tstatus\\nA\\tFAILED\\n')
''')
            engine = app.Analysis(Path(folder) / 'runs')
            with patch.object(app, 'RUNNER', runner), patch.object(app, 'dependencies', return_value={'ready': True, 'diamond': 'diamond'}):
                job = engine.start({'mode': 'off'}, {'genome': ('x.gbk', GENOME)})
                with self.assertRaisesRegex(ValueError, 'already running'):
                    engine.start({'mode': 'off'}, {'genome': ('x.gbk', GENOME)})
                deadline = time.monotonic() + 5
                while engine.active and time.monotonic() < deadline:
                    time.sleep(.03)
                status = engine.status(job)
                self.assertEqual(status['state'], 'finished')
                self.assertTrue(status['report'])
                self.assertEqual(status['interpro_failures'], 1)
                self.assertIn('engine started', status['log'])

    def test_http_token_origin_and_archive(self):
        with tempfile.TemporaryDirectory() as folder:
            server = app.ThreadingHTTPServer(('127.0.0.1', 0), app.Handler)
            server.token = 'test-private-token'
            server.analysis = app.Analysis(folder)
            threading.Thread(target=server.serve_forever, daemon=True).start()
            base = f'http://127.0.0.1:{server.server_port}/'
            try:
                with self.assertRaises(urllib.error.HTTPError) as err:
                    urllib.request.urlopen(base)
                self.assertEqual(err.exception.code, 403)
                with urllib.request.urlopen(base + server.token + '/') as response:
                    self.assertIn(b'Annotated genome', response.read())
                request = urllib.request.Request(base + server.token + '/api/run', data=b'', headers={'Origin': 'https://untrusted.example'})
                with self.assertRaises(urllib.error.HTTPError) as err:
                    urllib.request.urlopen(request)
                self.assertEqual(err.exception.code, 403)
                request = urllib.request.Request(base + server.token + '/api/run', data=b'')
                request.add_header('Content-Length', str(app.MAX_UPLOAD + 1))
                with self.assertRaises(urllib.error.HTTPError) as err:
                    urllib.request.urlopen(request)
                self.assertEqual(err.exception.code, 413)
                jobfolder = Path(folder) / 'abc'
                jobfolder.mkdir()
                (jobfolder / 'console.log').write_text('done')
                (jobfolder / 'analysis_VISUAL_REPORT.html').write_text('<h1>Result</h1>')
                server.analysis.jobs['abc'] = {'folder': jobfolder, 'state': 'finished'}
                with urllib.request.urlopen(base + server.token + '/download/abc') as response:
                    import io, zipfile
                    with zipfile.ZipFile(io.BytesIO(response.read())) as archive:
                        self.assertIn('console.log', archive.namelist())
                with urllib.request.urlopen(base + server.token + '/report/abc') as response:
                    self.assertIn('sandbox', response.headers['Content-Security-Policy'])
            finally:
                server.shutdown()
                server.server_close()


if __name__ == '__main__':
    unittest.main()
