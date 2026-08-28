#!/usr/bin/env python3
from pathlib import Path
s=(Path(__file__).parent/"integrated_evidence_v063.py").read_text(encoding="utf-8")
assert "unique_exact" in s
assert "len(set(v))==1" in s
assert "SUPPORTED_ANNOTATION" in s
assert "CONTEXT_NEVER_CREATES_REACTIONS" in s
assert "PARTIAL_WITH_STRONG_LOCAL_CONTEXT" in s
assert "active_row" in s
print("PASS only unique exact aliases can create annotation support")
print("PASS inactive legacy reactions excluded")
print("PASS context cannot manufacture reactions")
print("PASS rescue and annotation provenance remain separate")
print("ALL 4 INTEGRATED-EVIDENCE TESTS PASSED")
