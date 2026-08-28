#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--db",default="Bioremediation_DB_v0.7.6_dyes")
    ap.add_argument("--run-prefix",required=True)
    ap.add_argument("--pathway-summary",default="")
    ap.add_argument("--sample-name",default="")
    a=ap.parse_args()
    here=Path(__file__).resolve().parent
    pre=a.run_prefix
    ps=Path(a.pathway_summary) if a.pathway_summary else Path(pre+"_pathway_summary_STRICT_v088.tsv")
    if not ps.exists(): ps=Path(pre+"_adjudicated_adjudicated_pathway_summary.tsv")
    re=Path(pre+"_adjudicated_adjudicated_reaction_evidence.tsv")
    dc=Path(pre+"_adjudicated_cross_reaction_locus_decisions.tsv")
    for p in (ps,re):
        if not p.exists(): raise SystemExit(f"Missing: {p}")
    out=Path(pre+"_TABLE_REPORT_v089.html")
    cmd=[sys.executable,str(here/"table_first_reporter_v089.py"),"--db",a.db,
         "--pathway-summary",str(ps),"--reaction-evidence",str(re),
         "--out-html",str(out),"--sample-name",a.sample_name or pre]
    if dc.exists(): cmd += ["--cross-reaction-decisions",str(dc)]
    print("RUN:"," ".join(cmd))
    rc=subprocess.run(cmd).returncode
    if rc==0:
        print("Open with:")
        print("start",out)
    raise SystemExit(rc)

if __name__=="__main__":
    main()
