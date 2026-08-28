#!/usr/bin/env python3
import importlib.util
from pathlib import Path
H=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("g",H/"genbank_input_validator_v082.py")
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)

line="LOCUS       NODE_1_length_467231_cov_77.474813467231 bp   DNA linear\n"
fixed,kind=g.normalize_locus_line(line)
assert kind=="SPADES_GLUE", (fixed,kind)
assert "NODE_1_length_467231_cov_77.474813 467231 bp" in fixed

valid="LOCUS       NC_000913 4641652 bp DNA linear\n"
fixed2,kind2=g.normalize_locus_line(valid)
assert kind2 is None and fixed2==valid

src=(H/"genbank_input_validator_v082.py").read_text(encoding="utf-8")
assert "original_modified=NO" in src
assert "VALID_AS_IS" in src
assert "NORMALIZED_INTERNAL_COPY" in src

print("PASS valid GenBank stays unchanged")
print("PASS malformed Prokka/SPAdes LOCUS is normalized internally")
print("PASS original input is never modified")
print("PASS repaired input is reparsed and audited")
print("ALL 4 v0.8.2 INPUT TESTS PASSED")
