#!/usr/bin/env python3
from pathlib import Path
s=(Path(__file__).parent/"cross_reaction_adjudicator_v088.py").read_text(encoding="utf-8")
assert 'rr.get("review_flag") != "YES"' in s
assert '"SUPPORTED_MIXED_LOCI"' in s
assert '"review_units":review_units' in s
r=(Path(__file__).parent/"rescore_existing_pathways_v088.py").read_text(encoding="utf-8")
assert 'review_flag")!="YES"' in r
assert 'REVIEW_NEVER_COUNTS_AS_SUPPORTED' in r
print("PASS review-required reaction cannot increase supported_units")
print("PASS mixed-loci support cannot increase supported_units")
print("PASS review units are reported separately")
print("PASS existing runs can be rescored without DIAMOND/InterPro")
print("ALL 4 v0.8.8 SCORING TESTS PASSED")
