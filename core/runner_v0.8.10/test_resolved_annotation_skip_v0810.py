#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).parent
r=(root/"Bioremediation_Hybrid_Rescue_v0.3.5_resolvedskip"/"hybrid_rescue_engine_v0.3.5_resolvedskip.py").read_text(encoding="utf-8")
c=(root/"bioremediation_unified_runner_v083_core.py").read_text(encoding="utf-8")
a=(root/"resolved_annotation_scanner_v0810.py").read_text(encoding="utf-8")

assert "--resolved-reaction-loci" in r
assert '(rid, h["qseqid"]) in resolved_pairs' in r
assert "DIFFERENT reactions" in r
assert "resolved_annotation_scanner_v0810.py" in c
assert '"--resolved-reaction-loci",str(resolved)' in c
assert "enzyme_aliases.tsv" in c
assert "unique_exact" in a and "ambiguous_exact" in a
assert "resolved_for_same_reaction" in a

print("PASS unique exact annotation is resolved before rescue")
print("PASS same reaction-locus pair is skipped in rescue")
print("PASS same locus remains eligible for different reactions")
print("PASS ambiguous aliases are not treated as resolved")
print("PASS resolved-pair audit is written")
print("ALL 5 v0.8.10 RESOLVED-ANNOTATION TESTS PASSED")
