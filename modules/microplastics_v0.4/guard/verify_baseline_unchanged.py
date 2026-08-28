from pathlib import Path
import hashlib,json,sys
here=Path(__file__).resolve().parent
manifest=json.loads((here/'BASELINE_v0.7.6_SHA256.json').read_text())
base=Path(sys.argv[1]) if len(sys.argv)>1 else Path('Bioremediation_DB_v0.7.6_dyes')
if not base.is_dir(): raise SystemExit(f'FAIL: baseline directory not found: {base}')
actual={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in base.iterdir() if p.is_file()}
missing=sorted(set(manifest)-set(actual)); extra=sorted(set(actual)-set(manifest)); changed=sorted(k for k in manifest.keys()&actual.keys() if manifest[k]!=actual[k])
if missing or extra or changed:
 print('BASELINE GUARD: FAIL'); print('missing:',missing); print('extra:',extra); print('changed:',changed); raise SystemExit(2)
print('BASELINE GUARD: PASS')
print(f'Protected files unchanged: {len(manifest)}')
