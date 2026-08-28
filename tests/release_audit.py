#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, subprocess, sys, re
try:
    from Bio import SeqIO
except ImportError:
    raise SystemExit('ERROR: biopython is required')
ROOT=Path(__file__).resolve().parents[1]
core=ROOT/'core/database_v0.7.6'; micro=ROOT/'modules/microplastics_v0.4'
def tsv(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
assert len(list(core.glob('*')))==26, 'core database must contain exactly 26 protected files'
assert len({r['pathway_id'] for r in tsv(core/'pathways.tsv')})==34
assert len({r['reaction_id'] for r in tsv(core/'reactions.tsv')})==81
assert sum(1 for _ in SeqIO.parse(core/'reference_proteins.faa','fasta'))==379
expected={'BRXN_PET_HYD','BRXN_MHET_HYD','BRXN_TPA_DIOX','BRXN_TPA_DCD_DH','BRXN_PCL_DEP','BRXN_PBAT_DEP','BRXN_PU_DEP'}
assert {r['reaction_id'] for r in tsv(micro/'database/reactions.tsv')}==expected
assert 'BRXN_COP_ATPASE' not in (micro/'database/reactions.tsv').read_text(encoding='utf-8')
# No user-machine absolute paths in release text/code.
pat=re.compile(r'(?:[A-Z]:\\Users\\|E:\\26-08-26|C:/Users/|E:/26-08-26)',re.I)
hits=[]
for p in ROOT.rglob('*'):
    if p.resolve() == Path(__file__).resolve(): continue
    if p.is_file() and p.suffix.lower() in {'.py','.md','.txt','.tsv','.json','.toml','.cff','.gitignore'}:
        try:s=p.read_text(encoding='utf-8',errors='ignore')
        except:continue
        if pat.search(s):hits.append(str(p.relative_to(ROOT)))
assert not hits, f'personal absolute paths found: {hits}'
# Existing branch integrity tests + synthetic annotation positive control.
for test in [micro/'tests/test_branch_integrity_v04.py',micro/'tests/test_pet_v02_regression_lock.py',micro/'tests/test_pcl_strict_gate_v03.py',ROOT/'tests/test_microplastics_annotation_positive.py']:
    subprocess.run([sys.executable,str(test)],check=True)
print('RELEASE AUDIT: PASS')
print('Core: 34 pathways/modules; 81 reactions; 379 reference proteins; 26 protected files')
print('Microplastics: 7 exact reactions; PET/PCL/PBAT/PU; annotation positive control PASS')
print('Personal absolute path scan: PASS')
