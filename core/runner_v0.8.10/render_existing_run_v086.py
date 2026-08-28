#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path

def main():
    ap=argparse.ArgumentParser(description='Create a visual HTML report from an existing adjudicated run')
    ap.add_argument('--db',default='Bioremediation_DB_v0.7.6_dyes')
    ap.add_argument('--run-prefix',required=True)
    ap.add_argument('--sample-name',default='')
    a=ap.parse_args(); here=Path(__file__).resolve().parent; pre=a.run_prefix
    ps=Path(pre+'_adjudicated_adjudicated_pathway_summary.tsv'); re=Path(pre+'_adjudicated_adjudicated_reaction_evidence.tsv'); dc=Path(pre+'_adjudicated_cross_reaction_locus_decisions.tsv')
    for p in (ps,re):
        if not p.exists(): raise SystemExit(f'Missing: {p}')
    out=Path(pre+'_VISUAL_REPORT.html')
    cmd=[sys.executable,str(here/'visual_pathway_reporter_v086.py'),'--db',a.db,'--pathway-summary',str(ps),'--reaction-evidence',str(re),'--out-html',str(out),'--sample-name',a.sample_name or pre]
    if dc.exists(): cmd += ['--cross-reaction-decisions',str(dc)]
    print('RUN:',' '.join(cmd)); raise SystemExit(subprocess.run(cmd).returncode)
if __name__=='__main__': main()
