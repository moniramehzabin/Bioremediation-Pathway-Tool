#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, sys

ROOT=Path(__file__).resolve().parents[1]
DB=ROOT/'database'
EXPECTED={'BRXN_PET_HYD','BRXN_MHET_HYD','BRXN_TPA_DIOX','BRXN_TPA_DCD_DH','BRXN_PCL_DEP'}

def tsv(name):
    with (DB/name).open(encoding='utf-8') as f: return list(csv.DictReader(f,delimiter='\t'))

def fasta(path):
    seqs={}; key=None
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]; seqs[key]=''
        elif key: seqs[key]+=line.strip()
    return seqs

rx=tsv('reactions.tsv')
assert {r['reaction_id'] for r in rx}==EXPECTED
assert all(r['reaction_id']!='BRXN_COP_ATPASE' for r in rx)
paths=tsv('pathways.tsv')
assert any(r['pathway_id']=='BMOD_PCL_DEPOLYMERIZATION' and r['reaction_id']=='BRXN_PCL_DEP' and r['required']=='yes' for r in paths)
refs=fasta(DB/'reference_proteins.faa')
assert len(refs)==7, len(refs)
assert 'BRT_PCL_DEP|Q6A0I4' in refs
assert len(refs['BRT_PCL_DEP|Q6A0I4'])==301
sha=hashlib.sha256(refs['BRT_PCL_DEP|Q6A0I4'].encode()).hexdigest()
assert sha=='a54bad3ec73c049eb45a846ac4803088b66fd0b38ad407d9090357fa20b68cc9'
meta=tsv('reference_metadata.tsv')
m=next(r for r in meta if r['reference_id']=='BRT_PCL_DEP|Q6A0I4')
assert m['sequence_sha256']==sha and m['reaction_id']=='BRXN_PCL_DEP'
text='\n'.join(p.read_text(encoding='utf-8',errors='ignore') for p in DB.glob('*') if p.is_file())
assert 'BRXN_COP_ATPASE' not in text
unres=(DB/'unresolved_reference_targets.tsv').read_text(encoding='utf-8')
assert 'WP_004373894.1' in unres and 'WP_003239806.1' in unres
print('PASS v0.3 branch integrity')
print('Reactions: 5 (PET/TPA 4 + PCL 1)')
print('Reference proteins: 7 (PET/TPA 6 + PCL 1)')
print('PCL direct seed Q6A0I4 length/hash: PASS')
print('Copper ATPase contamination: NO')
