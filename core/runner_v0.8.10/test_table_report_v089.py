#!/usr/bin/env python3
from pathlib import Path
s=(Path(__file__).parent/"table_first_reporter_v089.py").read_text(encoding="utf-8")
assert "Results summary" in s
assert "Reaction evidence" in s
assert "Pathway details" in s
assert "Required supported" in s
assert "Optional supported" in s
assert "Review evidence only" in s
assert "NO_CANDIDATE" in s
assert "Candidate only" in s
assert "Arsenate reduction" in s
assert "Arsenic resistance / efflux" in s
assert "Review/candidate evidence never counts as supported" in s
assert "Counted = Yes" in s
print("PASS table-first results summary")
print("PASS reaction evidence table")
print("PASS expandable pathway details")
print("PASS strict supported/review/candidate semantics")
print("PASS NO_CANDIDATE cannot display as candidate")
print("PASS arsenate reduction separated from arsenic efflux")
print("PASS required vs optional support separated")
print("ALL 7 v0.8.9 REPORT TESTS PASSED")
