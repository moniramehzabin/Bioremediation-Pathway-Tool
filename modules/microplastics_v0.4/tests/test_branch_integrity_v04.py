from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parents[1]
with (ROOT/'database/reactions.tsv').open() as f: ids={r['reaction_id'] for r in csv.DictReader(f,delimiter='\t')}
exp={'BRXN_PET_HYD','BRXN_MHET_HYD','BRXN_TPA_DIOX','BRXN_TPA_DCD_DH','BRXN_PCL_DEP','BRXN_PBAT_DEP','BRXN_PU_DEP'}
assert ids==exp,(ids,exp)
text='\n'.join(p.read_text(errors='ignore') for p in ROOT.rglob('*.tsv'))
assert 'BRXN_COP_ATPASE' not in text
print('V0.4 INTEGRITY: PASS')
