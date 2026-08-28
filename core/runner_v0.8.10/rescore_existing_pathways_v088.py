#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv
from collections import defaultdict
from pathlib import Path

def tsv(p):
    with open(p,encoding="utf-8") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def wt(p,rows,fields):
    with open(p,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def truth(x):
    return str(x or "").strip().lower() in {"true","yes","1","supported"}

def i(x):
    try:return int(float(x))
    except:return 0

def main():
    ap=argparse.ArgumentParser(description="Strictly rescore pathway support from an existing adjudicated run")
    ap.add_argument("--db",required=True)
    ap.add_argument("--reaction-evidence",required=True)
    ap.add_argument("--old-pathway-summary",required=True)
    ap.add_argument("--out",required=True)
    a=ap.parse_args()

    paths=tsv(Path(a.db)/"pathways.tsv")
    ev={r["reaction_id"]:r for r in tsv(a.reaction_evidence)}
    old={r["pathway_id"]:r for r in tsv(a.old_pathway_summary)}

    byp=defaultdict(list)
    for p in paths:
        byp[p.get("pathway_id","")].append(p)

    out=[]
    for pid,rows in byp.items():
        units=[]; seen=set()
        for p in rows:
            required=str(p.get("required","")).lower() in {"yes","true","1","required"}
            if not required: continue
            ag=p.get("alternative_group","").strip()
            if ag:
                if ag in seen: continue
                members=[x.get("reaction_id","") for x in rows if x.get("alternative_group","").strip()==ag]
                units.append(("ALT:"+ag,members)); seen.add(ag)
            else:
                units.append((p.get("reaction_id",""),[p.get("reaction_id","")]))

        cand=supp=review=0; missing=[]
        for uname,members in units:
            cm=sm=rv=False
            for rid in members:
                r=ev.get(rid,{})
                state=r.get("integrated_state","")
                c=truth(r.get("candidate_present")) or state not in {"","NO_CANDIDATE","NOT_DETECTED"}
                # Strict: review-required and mixed-loci evidence never counts as supported.
                sp=(truth(r.get("exact_supported"))
                    and r.get("review_flag")!="YES"
                    and state not in {"CANDIDATE_CONFLICTING","NO_CANDIDATE","SUPPORTED_MIXED_LOCI"})
                cm|=c; sm|=sp; rv|=(r.get("review_flag")=="YES")
            cand+=int(cm); supp+=int(sm); review+=int(cm and not sm and rv)
            if not cm: missing.append(uname if not uname.startswith("ALT:") else uname+"["+"|".join(members)+"]")

        total=len(units)
        cp=100*cand/total if total else 0
        sp=100*supp/total if total else 0
        prev=old.get(pid,{})
        cluster_n=i(prev.get("cluster_reaction_count"))
        if total and supp==total:
            state="SUPPORTED_COMPLETE"
        elif total and cand==total:
            state="COMPLETE_CANDIDATE_REVIEW_REQUIRED"
        elif cluster_n>=3 and supp>=2:
            state="PARTIAL_WITH_STRONG_LOCAL_CONTEXT"
        elif cluster_n>=2 and cand>=2:
            state="PARTIAL_WITH_LOCAL_CONTEXT"
        elif cand:
            state="PARTIAL"
        else:
            state="NOT_DETECTED"

        nr=dict(prev)
        nr.update({
            "pathway_id":pid,
            "pathway_name":prev.get("pathway_name") or rows[0].get("pathway_name",""),
            "required_units":total,
            "candidate_units":cand,
            "candidate_percent":f"{cp:.1f}",
            "supported_units":supp,
            "supported_percent":f"{sp:.1f}",
            "review_units":review,
            "adjudicated_pathway_state":state,
            "missing_required_units":"|".join(missing),
            "support_rule":"REVIEW_NEVER_COUNTS_AS_SUPPORTED",
        })
        out.append(nr)

    fields=[]
    for r in out:
        for k in r:
            if k not in fields: fields.append(k)
    wt(a.out,out,fields)
    print("Strict Pathway Rescorer v0.8.8")
    print("Pathways rescored:",len(out))
    print("Rule: review_flag=YES never increases supported_units")
    print("Wrote:",a.out)

if __name__=="__main__":
    main()
