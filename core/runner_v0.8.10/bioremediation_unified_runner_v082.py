#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path

def run(cmd):
    print("RUN:"," ".join(map(str,cmd)))
    p=subprocess.run(cmd)
    if p.returncode: raise SystemExit(p.returncode)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("genbank")
    ap.add_argument("--db",default="Bioremediation_DB_v0.7.6_dyes")
    ap.add_argument("--interpro-tsv",required=True)
    ap.add_argument("--diamond",default=".\\diamond.exe")
    ap.add_argument("--out-prefix",required=True)
    ap.add_argument("--python",default=sys.executable)
    a=ap.parse_args()

    here=Path(__file__).resolve().parent
    norm=Path(a.out_prefix+"_normalized_input.gbk")
    audit=Path(a.out_prefix+"_GENBANK_INPUT_AUDIT.txt")

    run([a.python,str(here/"genbank_input_validator_v082.py"),a.genbank,
         "--out",str(norm),"--audit",str(audit)])

    used=norm if norm.exists() and norm.stat().st_size else Path(a.genbank)
    run([a.python,str(here/"bioremediation_unified_runner_v081.py"),str(used),
         "--db",a.db,"--interpro-tsv",a.interpro_tsv,
         "--diamond",a.diamond,"--out-prefix",a.out_prefix])

if __name__=="__main__":
    main()
