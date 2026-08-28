#!/usr/bin/env python3
"""
Bioremediation Unified Runner v0.8.1

GenBank
 -> Hybrid Rescue v0.3.4
 -> Integrated Evidence v0.6.3
 -> Cross-Reaction Locus Adjudication v0.8.1
 -> alternative-aware pathway scoring
 -> final report

Underlying scientific thresholds are preserved. v0.8.1 adds only the competition
layer needed to prevent a single broad homolog from inflating unrelated pathways.
"""
from __future__ import annotations
import argparse,csv,json,subprocess,sys
from pathlib import Path
from datetime import datetime

VERSION="0.8.1"

def check(p,label,dir=False):
    p=Path(p)
    ok=p.is_dir() if dir else p.is_file()
    if not ok: raise SystemExit(f"ERROR: {label} not found: {p}")
    return p

def run(cmd,log):
    print("\nRUN:"," ".join(map(str,cmd)))
    p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                     text=True,encoding="utf-8",errors="replace")
    print(p.stdout,end="")
    with open(log,"a",encoding="utf-8") as f:
        f.write("\nRUN: "+" ".join(map(str,cmd))+"\n"+p.stdout)
    if p.returncode: raise SystemExit(f"ERROR: command failed ({p.returncode}); see {log}")

def tsv(p):
    with open(p,encoding="utf-8") as f:return list(csv.DictReader(f,delimiter="\t"))

def report(prefix,pathway,reaction,decisions):
    P=tsv(pathway); R=tsv(reaction); D=tsv(decisions)
    rank={"SUPPORTED_COMPLETE":0,"COMPLETE_CANDIDATE_REVIEW_REQUIRED":1,
          "PARTIAL_WITH_STRONG_LOCAL_CONTEXT":2,"PARTIAL_WITH_LOCAL_CONTEXT":3,
          "PARTIAL":4,"NOT_DETECTED":5}
    P=sorted(P,key=lambda r:(rank.get(r.get("adjudicated_pathway_state",""),9),r.get("pathway_name","")))
    fp=Path(prefix+"_REPORT.txt")
    with open(fp,"w",encoding="utf-8") as f:
        f.write("BIOREMEDIATION PATHWAY REPORT — v0.8.1\n")
        f.write("="*76+"\n\n")
        f.write("Cross-reaction specificity competition is enabled.\n")
        f.write("One shared protein is not automatically counted for several unrelated reactions.\n")
        f.write("Alternative reaction groups are scored as OR.\n")
        f.write("Genomic context never creates a missing reaction.\n\n")
        for r in P:
            f.write(f"{r.get('pathway_name')} [{r.get('pathway_id')}]\n")
            f.write(f"  State: {r.get('adjudicated_pathway_state')}\n")
            f.write(f"  Candidate units: {r.get('candidate_units')}/{r.get('required_units')} ({r.get('candidate_percent')}%)\n")
            f.write(f"  Supported units: {r.get('supported_units')}/{r.get('required_units')} ({r.get('supported_percent')}%)\n")
            if r.get("cluster_reaction_count") not in {"","0"}:
                f.write(f"  Local cluster: {r.get('cluster_reaction_count')} reactions; {r.get('cluster_contig')}; span {r.get('cluster_span_bp')} bp\n")
            if r.get("missing_required_units"):
                f.write(f"  Missing: {r.get('missing_required_units')}\n")
            f.write("\n")

        f.write("\nINFORMATIVE REACTION CALLS\n"+"="*76+"\n")
        for r in sorted(R,key=lambda x:x.get("reaction_id","")):
            if r.get("integrated_state")=="NO_CANDIDATE": continue
            f.write(f"{r.get('reaction_id')} | {r.get('enzyme_name','')}\n")
            f.write(f"  State: {r.get('integrated_state')}\n")
            f.write(f"  Loci: {r.get('all_evidence_loci') or '-'}\n")
            if r.get("removed_competing_loci"):
                f.write(f"  Competing loci removed: {r.get('removed_competing_loci')}\n")
            if r.get("ambiguous_loci"):
                f.write(f"  Ambiguous loci: {r.get('ambiguous_loci')}\n")
            if r.get("review_flag")=="YES":
                f.write("  REVIEW REQUIRED\n")
            f.write("\n")

        f.write("\nCROSS-REACTION COMPETITIONS\n"+"="*76+"\n")
        for d in D:
            if d.get("decision") in {"WINNER_KEEP","LOSER_REMOVE","AMBIGUOUS_REVIEW"}:
                f.write(f"{d.get('locus_tag')} | {d.get('reaction_id')} | {d.get('decision')} | score={d.get('score')}")
                if d.get("winner_reaction"):
                    f.write(f" | winner={d.get('winner_reaction')}")
                f.write("\n")
    return fp

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("genbank")
    ap.add_argument("--db",default="Bioremediation_DB_v0.7.6_dyes")
    ap.add_argument("--interpro-tsv",required=True)
    ap.add_argument("--diamond",default=".\\diamond.exe")
    ap.add_argument("--rescue-engine",default="Bioremediation_Hybrid_Rescue_v0.3.4_canonical\\hybrid_rescue_engine_v0.3.4.py")
    ap.add_argument("--integrated-engine",default="Bioremediation_Integrated_Evidence_v0.6.3\\integrated_evidence_v063.py")
    ap.add_argument("--adjudicator",default="Bioremediation_Unified_Runner_v0.8.1\\cross_reaction_adjudicator_v081.py")
    ap.add_argument("--policy",default="Bioremediation_Hybrid_Rescue_v0.3.4_canonical\\reaction_interpro_policy_v054.tsv")
    ap.add_argument("--out-prefix",required=True)
    ap.add_argument("--python",default=sys.executable)
    a=ap.parse_args()

    gb=check(a.genbank,"GenBank"); db=check(a.db,"database",True)
    ipr=check(a.interpro_tsv,"InterPro TSV")
    rescue=check(a.rescue_engine,"rescue engine")
    integ=check(a.integrated_engine,"integrated engine")
    adjud=check(a.adjudicator,"adjudicator")
    policy=check(a.policy,"policy")
    for fn in ["reference_proteins.faa","reference_metadata.tsv","reactions.tsv","pathways.tsv","reaction_components.tsv"]:
        check(db/fn,f"database {fn}")

    prefix=a.out_prefix
    rpre=prefix+"_rescue"
    ipre=prefix+"_integrated"
    apre=prefix+"_adjudicated"
    log=prefix+"_RUN.log"

    run([a.python,str(rescue),str(gb),
         "--references",str(db/"reference_proteins.faa"),
         "--reference-metadata",str(db/"reference_metadata.tsv"),
         "--reactions",str(db/"reactions.tsv"),
         "--pathways",str(db/"pathways.tsv"),
         "--reaction-components",str(db/"reaction_components.tsv"),
         "--policy",str(policy),"--diamond",a.diamond,
         "--interpro-mode","precomputed","--interpro-tsv",str(ipr),
         "--out-prefix",rpre],log)

    rsum=Path(rpre+"_reaction_evidence_summary.tsv")
    rall=Path(rpre+"_all_candidates.tsv")
    run([a.python,str(integ),str(gb),"--db",str(db),
         "--reaction-evidence",str(rsum),"--out-prefix",ipre],log)

    ire=Path(ipre+"_integrated_reaction_evidence.tsv")
    run([a.python,str(adjud),str(gb),"--db",str(db),
         "--integrated-reaction-evidence",str(ire),
         "--all-candidates",str(rall),
         "--out-prefix",apre],log)

    are=Path(apre+"_adjudicated_reaction_evidence.tsv")
    aps=Path(apre+"_adjudicated_pathway_summary.tsv")
    dec=Path(apre+"_cross_reaction_locus_decisions.tsv")
    for p in [are,aps,dec]:
        check(p,"final output")

    rep=report(prefix,aps,are,dec)
    manifest={
        "runner_version":VERSION,
        "finished":datetime.now().isoformat(timespec="seconds"),
        "genbank":str(gb.resolve()),"database":str(db.resolve()),
        "outputs":{
            "final_pathway_summary":str(aps),
            "final_reaction_evidence":str(are),
            "cross_reaction_decisions":str(dec),
            "report":str(rep),"log":log
        }
    }
    with open(prefix+"_RUN_MANIFEST.json","w",encoding="utf-8") as f:
        json.dump(manifest,f,indent=2)

    print("\n"+"="*76)
    print("UNIFIED RUN v0.8.1 COMPLETE")
    print("="*76)
    print("Final pathway summary:",aps)
    print("Final reaction evidence:",are)
    print("Cross-reaction decisions:",dec)
    print("Human-readable report:",rep)

if __name__=="__main__":
    main()
