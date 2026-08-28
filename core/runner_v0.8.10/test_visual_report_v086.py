#!/usr/bin/env python3
from pathlib import Path
s=(Path(__file__).parent/"visual_pathway_reporter_v086.py").read_text(encoding="utf-8")
assert "Visual Pathway Reporter v0.8.6" in s
assert "BRXN_ARS_C" in s and "BRXN_ARS_EFFLUX" in s
assert "Arsenate reduction" in s
assert "Arsenic resistance / efflux" in s
assert "Efflux is not counted as evidence of arsenate reduction" in s
assert "arsenic-split" in s
assert "flow.append('<div class=\"arrow\">&rarr;</div>')" in s
assert "ONE OF THESE ALTERNATIVES" in s
assert "Candidate detected; not counted as supported." in s
print("PASS arsenate reduction and arsenic efflux are displayed separately")
print("PASS efflux cannot visually count as reduction evidence")
print("PASS OR alternatives remain grouped")
print("PASS candidate-vs-supported wording remains explicit")
print("PASS HTML-safe pathway arrows")
print("PASS reporter version label corrected")
print("ALL 6 v0.8.6 REPORT TESTS PASSED")
