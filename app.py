#!/usr/bin/env python3
"""Local browser interface for BPT. No additional web-framework dependency."""
from __future__ import annotations

import argparse
import importlib
import json
import os
from pathlib import Path
import re
import secrets
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
import webbrowser
import zipfile
from excel_export import export_results
from email.parser import BytesParser
from email.policy import default
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / 'core/runner_v0.8.10/bioremediation_unified_runner_v0810.py'
MAX_UPLOAD = 80 * 1024 * 1024


def dependencies():
    problems = []
    for module in ('Bio.SeqIO', 'requests'):
        try:
            importlib.import_module(module)
        except Exception as exc:
            problems.append(f'{module}: {exc}')
    diamond = shutil.which('diamond') or shutil.which('diamond.exe')
    if not diamond:
        for name in ('diamond', 'diamond.exe'):
            candidate = ROOT / name
            if candidate.is_file():
                diamond = str(candidate)
                break
    if not RUNNER.is_file():
        problems.append('The analysis engine is missing. Extract the complete app folder.')
    return {'ready': not problems and bool(diamond), 'problems': problems,
            'diamond': diamond or '', 'python': sys.version.split()[0]}


def parse_form(content_type, data):
    if not content_type.startswith('multipart/form-data;'):
        raise ValueError('Choose a genome file using the form.')
    message = BytesParser(policy=default).parsebytes(
        b'Content-Type: ' + content_type.encode('ascii') + b'\r\nMIME-Version: 1.0\r\n\r\n' + data)
    if not message.is_multipart():
        raise ValueError('The uploaded form could not be read.')
    fields, files = {}, {}
    for part in message.iter_parts():
        name = part.get_param('name', header='content-disposition')
        if not name:
            continue
        payload = part.get_payload(decode=True) or b''
        if part.get_filename():
            files[name] = (part.get_filename(), payload)
        else:
            fields[name] = payload.decode('utf-8')
    return fields, files


def validate_inputs(fields, files):
    genome = files.get('genome')
    if not genome or not genome[1]:
        raise ValueError('Choose a nonempty annotated GenBank file.')
    if Path(genome[0]).suffix.lower() not in {'.gb', '.gbk', '.gbff', '.genbank'}:
        raise ValueError('Use an annotated GenBank file (.gb, .gbk, .gbff or .genbank).')
    if b'LOCUS ' not in genome[1] or b'FEATURES' not in genome[1]:
        raise ValueError('This does not look like an annotated GenBank file.')
    mode = fields.get('mode', 'web')
    if mode not in {'web', 'precomputed', 'off'}:
        raise ValueError('Choose a valid domain-validation mode.')
    email = fields.get('email', '').strip()
    if mode == 'web' and not re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+', email):
        raise ValueError('Enter a valid email for EMBL-EBI InterProScan.')
    if mode == 'precomputed' and not files.get('interpro', ('', b''))[1]:
        raise ValueError('Choose the existing InterProScan TSV file.')
    return mode, email


class Analysis:
    def __init__(self, directory):
        self.directory = Path(directory).expanduser().resolve()
        self.jobs = {}
        self.lock = threading.RLock()
        self.active = None

    def start(self, fields, files):
        mode, email = validate_inputs(fields, files)
        checks = dependencies()
        if not checks['ready']:
            raise ValueError('Setup is incomplete. Check the setup panel before running.')
        with self.lock:
            if self.active:
                raise ValueError('An analysis is already running. Wait for it or cancel it first.')
            job_id = uuid.uuid4().hex
            folder = self.directory / job_id
            folder.mkdir(parents=True)
            genome = folder / 'input.gbk'
            genome.write_bytes(files['genome'][1])
            prefix = folder / 'analysis'
            command = [sys.executable, '-u', str(RUNNER), str(genome),
                       '--diamond', checks['diamond'], '--interpro-mode', mode,
                       '--out-prefix', str(prefix)]
            if mode == 'web':
                command += ['--email', email]
            elif mode == 'precomputed':
                tsv = folder / 'input_interpro.tsv'
                tsv.write_bytes(files['interpro'][1])
                command += ['--interpro-tsv', str(tsv)]
            job = {'id': job_id, 'state': 'starting', 'folder': folder,
                   'mode': mode, 'cancel': False, 'process': None,
                   'filename': Path(files['genome'][0]).name, 'started': time.time()}
            (folder / 'app_run.json').write_text(json.dumps({
                'input_name': job['filename'], 'mode': mode, 'command': command,
                'started': job['started'], 'ui_version': '0.2.0',
                'note': 'Adds separately documented DyP annotation-seeded domain validation; DIAMOND thresholds unchanged.'}, indent=2))
            self.jobs[job_id] = job
            self.active = job_id
            threading.Thread(target=self._run, args=(job, command), daemon=True).start()
            return job_id

    def _run(self, job, command):
        log = job['folder'] / 'console.log'
        try:
            env = dict(os.environ, PYTHONUNBUFFERED='1')
            kwargs = {'start_new_session': True} if os.name != 'nt' else {
                'creationflags': subprocess.CREATE_NEW_PROCESS_GROUP}
            with log.open('wb') as output:
                with self.lock:
                    if job['cancel']:
                        job['state'] = 'cancelled'
                        return
                    process = subprocess.Popen(command, cwd=ROOT, env=env,
                        stdin=subprocess.DEVNULL, stdout=output, stderr=subprocess.STDOUT, **kwargs)
                    job['process'] = process
                    job['state'] = 'running'
                code = process.wait()
            if code == 0 and not job['cancel']:
                try:
                    export_results(job['folder'])
                except Exception as exc:
                    job['export_error'] = str(exc)
            with self.lock:
                job['exit_code'] = code
                job['state'] = 'cancelled' if job['cancel'] else ('finished' if code == 0 else 'failed')
        except Exception as exc:
            with log.open('ab') as output:
                output.write(('\nCould not run analysis: ' + str(exc)).encode())
            job['state'] = 'failed'
        finally:
            with self.lock:
                job['process'] = None
                self.active = None

    def cancel(self, job_id):
        with self.lock:
            job = self.jobs[job_id]
            if job['state'] not in {'starting', 'running'}:
                return
            job['cancel'] = True
            process = job['process']
            if process and process.poll() is None:
                if os.name == 'nt':
                    subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'],
                                   capture_output=True, check=False)
                else:
                    os.killpg(process.pid, signal.SIGTERM)
                    # Retried web requests must not leave child processes running indefinitely.
                    def force_stop():
                        time.sleep(3)
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                    threading.Thread(target=force_stop, daemon=True).start()

    def status(self, job_id):
        with self.lock:
            job = self.jobs[job_id]
            log = job['folder'] / 'console.log'
            text = ''
            if log.exists():
                with log.open('rb') as handle:
                    handle.seek(max(0, log.stat().st_size - 65000))
                    text = handle.read().decode('utf-8', errors='replace')
            failed = job['folder'] / 'analysis_rescue_interpro_failed.tsv'
            failures = 0
            if failed.exists():
                failures = max(0, len(failed.read_text().splitlines()) - 1)
            return {'id': job_id, 'state': job['state'], 'log': text,
                    'mode': job['mode'], 'interpro_failures': failures,
                    'excel': (job['folder'] / 'BPT_results.xlsx').is_file(),
                    'export_error': job.get('export_error', ''),
                    'report': (job['folder'] / 'analysis_VISUAL_REPORT.html').is_file(),
                    'filename': job['filename']}


PAGE = r'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Bioremediation Pathway Tool</title>
<style>
:root{font:16px/1.5 system-ui,sans-serif;color:#18372f;background:#f4f7f4}*{box-sizing:border-box}
body{margin:0}main{max-width:1000px;margin:48px auto;padding:0 24px}header{margin-bottom:32px}
.eyebrow{color:#557367;font-size:13px;letter-spacing:.1em;text-transform:uppercase}h1{font-size:36px;line-height:1.15;margin:12px 0}
h2{font-size:21px;margin:0 0 18px}p{color:#53685f}section{background:white;border:1px solid #dce5de;border-radius:14px;padding:28px;margin:20px 0}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:24px}label{display:block;font-weight:600;margin-bottom:7px}
input,select{width:100%;padding:12px;border:1px solid #bdcec3;border-radius:7px;background:white;color:#18372f;font:inherit}
small{display:block;color:#61786b;margin-top:7px}button,.button{display:inline-block;background:#1d654f;color:white;border:0;border-radius:7px;padding:12px 20px;font:600 15px system-ui;cursor:pointer;text-decoration:none}
button:disabled{opacity:.45;cursor:not-allowed}.secondary{background:#e7efe9;color:#214c3b}#cancel{background:#f5e7e0;color:#813a24}
.actions{display:flex;gap:12px;flex-wrap:wrap;margin-top:24px}.badge{font-size:13px;background:#e6eee7;border-radius:30px;padding:5px 12px;display:inline-block}
#error{color:#943622;white-space:pre-wrap}#notice{color:#72521c}pre{background:#152f28;color:#e3eee6;padding:20px;border-radius:8px;max-height:340px;overflow:auto;font:12px/1.6 ui-monospace,monospace;white-space:pre-wrap;word-break:break-word}
details{margin-top:20px}summary{cursor:pointer}footer{font-size:13px;color:#63776c;margin:25px 0} [hidden]{display:none!important}
@media(max-width:650px){.grid{grid-template-columns:1fr}h1{font-size:29px}main{margin-top:28px}section{padding:20px}}
</style><main><header><span class="eyebrow">Genome evidence · local analysis</span>
<h1>Bioremediation Pathway Tool</h1><p>Explore supported reactions, pathway gaps and candidates for review.</p>
<span class="badge">Runs on your computer · App preview 0.2</span></header>
<section><h2>New analysis</h2><form id="form"><div class="grid"><div>
<label for="genome">Annotated genome</label><input id="genome" name="genome" type="file" accept=".gb,.gbk,.gbff,.genbank" required>
<small>GenBank format with translated coding sequences. Maximum upload: 80 MB including other files.</small></div><div>
<label for="mode">Domain validation</label><select id="mode" name="mode"><option value="web">InterPro online — recommended</option><option value="precomputed">Use existing InterProScan TSV</option><option value="off">Skip InterPro — preliminary analysis</option></select>
<small id="modeHelp">Selected protein sequences and your email are sent to EMBL-EBI for domain analysis.</small></div>
<div id="emailField"><label for="email">Email for InterPro</label><input id="email" name="email" type="email" placeholder="you@example.com" required autocomplete="email"></div>
<div id="tsvField" hidden><label for="interpro">InterProScan results</label><input id="interpro" name="interpro" type="file" accept=".tsv,.txt"></div></div>
<details><summary>Setup status</summary><p id="setup">Checking Python packages and DIAMOND…</p></details>
<p id="error" role="alert"></p><div class="actions"><button id="run" disabled>Run analysis</button></div></form></section>
<section id="progress" hidden><h2>Analysis progress</h2><span id="state" class="badge">Starting</span><p id="filename"></p><p id="notice"></p>
<pre id="log" aria-label="Analysis log"></pre><div class="actions"><button id="cancel" type="button">Cancel run</button>
<a id="report" class="button" target="_blank" rel="noopener" hidden>Open report</a><a id="excel" class="button secondary" hidden>Download Excel</a><a id="download" class="button secondary" hidden>Download all results</a></div></section>
<footer>Predicted functions are computational evidence, not proof of pollutant degradation. Closing this tab does not stop the analysis. Return to this page to check progress.</footer></main>
<script>
const $=id=>document.getElementById(id);let job=null,ready=false,timer;
async function api(path,opts){const r=await fetch(path,opts);const d=await r.json();if(!r.ok)throw Error(d.error||'Request failed');return d}
$('mode').onchange=()=>{const m=$('mode').value;$('emailField').hidden=m!=='web';$('email').required=m==='web';$('tsvField').hidden=m!=='precomputed';$('interpro').required=m==='precomputed';$('modeHelp').textContent=m==='web'?'Selected protein sequences and your email are sent to EMBL-EBI for domain analysis.':m==='off'?'Domain validation will be skipped. Treat the results as preliminary.':'Use matching protein identifiers from this genome.'};
async function check(){try{const d=await api('api/setup');ready=d.ready;$('setup').textContent=d.ready?'Ready · Python '+d.python+' · DIAMOND found':('Setup needed: '+d.problems.join('; ')+(d.diamond?'':' DIAMOND was not found.'));$('run').disabled=!ready}catch(e){$('error').textContent=e.message}}
async function poll(){if(!job)return;try{const d=await api('api/status/'+job);$('progress').hidden=false;$('state').textContent=d.state;$('filename').textContent=d.filename;const box=$('log');const near=box.scrollHeight-box.scrollTop-box.clientHeight<70;box.textContent=d.log;if(near)box.scrollTop=box.scrollHeight;const done=['finished','failed','cancelled'].includes(d.state);$('cancel').hidden=done;$('run').disabled=!ready||!done;$('report').hidden=!done||!d.report;$('report').href='report/'+job;$('excel').hidden=!done||!d.excel;$('excel').href='excel/'+job;$('download').hidden=!done;$('download').href='download/'+job;$('notice').textContent=d.export_error?'Analysis finished, but Excel export failed: '+d.export_error:d.state==='failed'?'The run stopped with an error. Download the logs for troubleshooting.':d.state==='cancelled'?'Cancelled. Partial outputs remain available.':d.interpro_failures?'Warning: '+d.interpro_failures+' InterPro checks were unavailable. Review these before interpreting results.':d.mode==='off'?'InterPro was skipped; this is a preliminary analysis.':d.state==='finished'?'Analysis finished. Review evidence and any warnings in the report.':'InterPro may take time. Keep the app running.';if(done)clearInterval(timer)}catch(e){$('error').textContent=e.message}}
$('form').onsubmit=async e=>{e.preventDefault();$('error').textContent='';$('run').disabled=true;try{const d=await api('api/run',{method:'POST',body:new FormData($('form'))});job=d.id;sessionStorage.setItem('bptJob',job);clearInterval(timer);await poll();timer=setInterval(poll,2000)}catch(e){$('error').textContent=e.message;$('run').disabled=!ready}};
$('cancel').onclick=async()=>{try{await api('api/cancel/'+job,{method:'POST'});await poll()}catch(e){$('error').textContent=e.message}};
check().then(()=>{job=sessionStorage.getItem('bptJob');if(job){poll();timer=setInterval(poll,2000)}});
</script></html>'''


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def send(self, payload, status=200, kind='application/json', attachment=None):
        if kind == 'application/json':
            payload = json.dumps(payload).encode()
        elif isinstance(payload, str):
            payload = payload.encode()
        self.send_response(status)
        self.send_header('Content-Type', kind)
        self.send_header('Content-Length', str(len(payload)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Referrer-Policy', 'no-referrer')
        if kind == 'text/html; charset=utf-8' and attachment is None and '/report/' in self.path:
            self.send_header('Content-Security-Policy', "sandbox allow-scripts; default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data:")
        if attachment:
            self.send_header('Content-Disposition', f'attachment; filename="{attachment}"')
        self.end_headers()
        self.wfile.write(payload)

    def route(self):
        # Exact Host and per-launch token protect local file upload and job control.
        if self.headers.get('Host') != f'127.0.0.1:{self.server.server_port}':
            raise PermissionError('Invalid host.')
        path = unquote(urlsplit(self.path).path)
        prefix = '/' + self.server.token + '/'
        if not path.startswith(prefix):
            raise PermissionError('Open the private app link printed in the launcher.')
        return path[len(prefix):]

    def do_GET(self):
        try:
            path = self.route()
            if path == '':
                return self.send(PAGE, kind='text/html; charset=utf-8')
            if path == 'api/setup':
                return self.send(dependencies())
            if path.startswith('api/status/'):
                return self.send(self.server.analysis.status(path.split('/')[-1]))
            if path.startswith(('report/', 'download/', 'excel/')):
                job_id = path.split('/')[-1]
                job = self.server.analysis.jobs[job_id]
                if job['state'] in {'starting', 'running'}:
                    raise ValueError('The run is still active.')
                if path.startswith('report/'):
                    return self.send((job['folder'] / 'analysis_VISUAL_REPORT.html').read_bytes(),
                                     kind='text/html; charset=utf-8')
                if path.startswith('excel/'):
                    return self.send((job['folder'] / 'BPT_results.xlsx').read_bytes(),
                        kind='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet', attachment='BPT_results.xlsx')
                archive = job['folder'] / 'results.zip'
                with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as output:
                    for item in sorted(job['folder'].iterdir()):
                        if item.is_file() and item != archive:
                            output.write(item, item.name)
                return self.send(archive.read_bytes(), kind='application/zip', attachment='BPT_results.zip')
            self.send({'error': 'Not found.'}, 404)
        except PermissionError as exc:
            self.send({'error': str(exc)}, 403)
        except (KeyError, FileNotFoundError):
            self.send({'error': 'This result is unavailable. Start a new analysis from the app home page.'}, 404)
        except ValueError as exc:
            self.send({'error': str(exc)}, 400)

    def do_POST(self):
        try:
            path = self.route()
            origin = self.headers.get('Origin')
            if origin and origin != f'http://127.0.0.1:{self.server.server_port}':
                raise PermissionError('Invalid origin.')
            length = int(self.headers.get('Content-Length', '0'))
            if length < 0 or length > MAX_UPLOAD:
                return self.send({'error': 'The total upload must be under 80 MB.'}, 413)
            if path == 'api/run':
                fields, files = parse_form(self.headers.get('Content-Type', ''), self.rfile.read(length))
                return self.send({'id': self.server.analysis.start(fields, files)}, 202)
            if path.startswith('api/cancel/'):
                self.server.analysis.cancel(path.split('/')[-1])
                return self.send({'ok': True})
            self.send({'error': 'Not found.'}, 404)
        except PermissionError as exc:
            self.send({'error': str(exc)}, 403)
        except (ValueError, UnicodeError) as exc:
            self.send({'error': str(exc)}, 400)
        except KeyError:
            self.send({'error': 'Run not found.'}, 404)


def open_browser(url):
    if os.name != 'nt' and 'microsoft' in __import__('platform').release().lower():
        try:
            subprocess.Popen(['cmd.exe', '/c', 'start', '', url], stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)
            return
        except OSError:
            pass
    webbrowser.open(url)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Check required packages and DIAMOND, then exit.')
    parser.add_argument('--no-browser', action='store_true')
    parser.add_argument('--port', type=int, default=0)
    parser.add_argument('--runs-dir', type=Path, default=Path.home() / 'BPT_results')
    args = parser.parse_args()
    if args.check:
        result = dependencies()
        print(json.dumps(result))
        return 0 if result['ready'] else 1
    server = ThreadingHTTPServer(('127.0.0.1', args.port), Handler)
    server.token = secrets.token_urlsafe(24)
    server.analysis = Analysis(args.runs_dir)
    url = f'http://127.0.0.1:{server.server_port}/{server.token}/'
    print('Bioremediation Pathway Tool · local app', flush=True)
    print(url, flush=True)
    print(f'Results: {args.runs_dir}\nKeep this window open. Press Ctrl+C to stop the app.', flush=True)
    if not args.no_browser:
        threading.Timer(0.5, open_browser, args=(url,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        if server.analysis.active:
            server.analysis.cancel(server.analysis.active)
    finally:
        server.server_close()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
