#!/usr/bin/env python3
from pathlib import Path
import csv, subprocess, sys, tempfile
ROOT=Path(__file__).resolve().parents[1]
scanner=ROOT/'modules/microplastics_v0.4/runner/resolved_annotation_scanner_v0810.py'
gbk=ROOT/'examples/microplastics_annotation_positive_control.gbk'
db=ROOT/'modules/microplastics_v0.4/database'
expected={'BRXN_PET_HYD','BRXN_MHET_HYD','BRXN_TPA_DIOX','BRXN_TPA_DCD_DH','BRXN_PCL_DEP','BRXN_PBAT_DEP','BRXN_PU_DEP'}
with tempfile.TemporaryDirectory() as td:
    out=Path(td)/'pairs.tsv'
    subprocess.run([sys.executable,str(scanner),str(gbk),'--db',str(db),'--out',str(out)],check=True)
    with out.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    got={r['reaction_id'] for r in rows}
    assert got==expected,(got,expected)
    tpa={r['gene'].lower() for r in rows if r['reaction_id']=='BRXN_TPA_DIOX'}
    assert {'tpha1','tpha2','tpha3'} <= tpa,tpa
print('PASS: synthetic annotation positive control resolves all 7 microplastics reactions; TPADO has A1/A2/A3 components')
