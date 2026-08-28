#!/usr/bin/env python3
import csv,re
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DB=ROOT/'core'/'database_v0.7.6'
def read_tsv(p):
    with open(p,encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))
def toks(v):
    return {t.lower() for t in re.split(r'[;,/|\s]+',str(v or '')) if t.strip()}
meta={r['reference_id']:r for r in read_tsv(DB/'reference_metadata.tsv')}
comps=read_tsv(DB/'reaction_components.tsv')
by=defaultdict(list)
for cr in comps:
    if cr.get('reaction_id') and cr.get('component_id'):
        by[cr['reaction_id']].append(cr)
assigned={}
amb=[]
for ref_id,m in meta.items():
    if str(m.get('component_id','')).strip():
        assigned[ref_id]=str(m.get('component_id')).strip(); continue
    defs=by.get(str(m.get('reaction_id','')).strip(),[])
    st=str(m.get('source_target_id','')).strip(); gt=toks(m.get('gene_primary',''))
    matches=[]
    for cr in defs:
        target=str(cr.get('target_id','')).strip(); gf=toks(cr.get('gene_family',''))
        if (st and target and st==target) or (gt and gf and gt & gf):
            matches.append(cr)
    uniq={str(x.get('component_id','')).strip():x for x in matches if str(x.get('component_id','')).strip()}
    if len(uniq)==1:
        assigned[ref_id]=next(iter(uniq))
    elif len(uniq)>1:
        amb.append((ref_id,sorted(uniq)))
assert not amb, f'Ambiguous component mappings: {amb}'
assert assigned['BRT_PCA_I|Q01103']=='A'
assert assigned['BRT_PCA_I|Q43973']=='A'
assert assigned['BRT_PCA_J|P0A102']=='B'
assert assigned['BRT_PCA_J|Q59091']=='B'
assert assigned['BRT_PCA_J|P0A101']=='B'
print('PASS component metadata propagation is unambiguous')
print('PASS PcaI references -> component A')
print('PASS PcaJ references -> component B')
