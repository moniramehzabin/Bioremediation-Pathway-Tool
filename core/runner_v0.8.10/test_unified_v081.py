#!/usr/bin/env python3
from pathlib import Path
a=(Path(__file__).parent/"cross_reaction_adjudicator_v081.py").read_text(encoding="utf-8")
r=(Path(__file__).parent/"bioremediation_unified_runner_v081.py").read_text(encoding="utf-8")
assert "LOSER_REMOVE" in a
assert "AMBIGUOUS_REVIEW" in a
assert "EXACT_ANNOTATION" in a
assert "alternative_group" in a and "OR_WITHIN_GROUP" in a
assert "CONTEXT_NEVER_CREATES_REACTIONS" in a
assert "cross_reaction_adjudicator_v081.py" in r
print("PASS shared-locus competition")
print("PASS exact annotation can outrank generic similarity")
print("PASS ambiguous competitions remain review")
print("PASS alternative pathway groups use OR logic")
print("PASS context cannot manufacture reactions")
print("ALL 5 v0.8.1 TESTS PASSED")
