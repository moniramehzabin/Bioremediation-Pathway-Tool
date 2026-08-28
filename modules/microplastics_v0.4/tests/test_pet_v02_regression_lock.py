#!/usr/bin/env python3
from pathlib import Path
import csv,hashlib,json
ROOT=Path(__file__).resolve().parents[1]; DB=ROOT/'database'
PET={'BRXN_PET_HYD','BRXN_MHET_HYD','BRXN_TPA_DIOX','BRXN_TPA_DCD_DH'}
M=json.loads((ROOT/'guard'/'PET_v0.2_SUBSET_SHA256.json').read_text())
def rows(path):
    with path.open(encoding='utf-8') as f: rr=list(csv.DictReader(f,delimiter='\t'))
    return [r for r in rr if r.get('reaction_id') in PET]
def digest(rr):
    s='\n'.join('\t'.join(f'{k}={r.get(k,"")}' for k in r.keys()) for r in rr)
    return hashlib.sha256(s.encode()).hexdigest(),len(rr)
for fn in ['reactions.tsv','pathways.tsv','enzyme_aliases.tsv','reference_metadata.tsv','reaction_components.tsv']:
    d,n=digest(rows(DB/fn)); assert d==M[fn]['sha256'],(fn,d,M[fn]['sha256']); assert n==M[fn]['rows']
d,n=digest(rows(ROOT/'reaction_interpro_policy_microplastics_v03.tsv')); assert d==M['policy_pet_rows']['sha256']; assert n==M['policy_pet_rows']['rows']
seqs={}; k=None
for line in (DB/'reference_proteins.faa').read_text(encoding='utf-8').splitlines():
    if line.startswith('>'): k=line[1:].split()[0]; seqs[k]=''
    elif k: seqs[k]+=line.strip()
for k,h in M['reference_sequences'].items(): assert hashlib.sha256(seqs[k].encode()).hexdigest()==h,k
print('PASS PET v0.2 regression lock: all PET/TPA rows and 6 reference sequences unchanged')
