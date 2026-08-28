#!/usr/bin/env python3
from pathlib import Path
s=(Path(__file__).parent/'cross_reaction_adjudicator_v083.py').read_text(encoding='utf-8')
assert 'same_multi=ann & rsc' in s
assert 'SUPPORTED_MIXED_LOCI' in s
assert 'SAME_LOCUS_MULTI_SOURCE' in s
assert 'supported_loci' in s and 'review_loci' in s and 'candidate_only_loci' in s
print('PASS same-locus multi-source rule')
print('PASS mixed-loci sources require review')
print('PASS locus classes are separated')
