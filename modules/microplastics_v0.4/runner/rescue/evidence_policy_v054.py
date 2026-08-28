#!/usr/bin/env python3
import argparse,csv
from collections import defaultdict
from pathlib import Path

def read_tsv(p):
    with open(p,encoding="utf-8",errors="replace") as f:return list(csv.DictReader(f,delimiter="\t"))
def parse_ipr(p):
    out=defaultdict(list)
    with open(p,encoding="utf-8",errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            c=line.rstrip("\n").split("\t")
            if len(c)<6: continue
            out[c[0]].append({"analysis":c[3] if len(c)>3 else "","signature_accession":c[4] if len(c)>4 else "","signature_description":c[5] if len(c)>5 else "","interpro_id":c[11] if len(c)>11 and c[11]!="-" else "","interpro_description":c[12] if len(c)>12 and c[12]!="-" else ""})
    return out
def splitv(x): return [z.strip() for z in str(x or "").split("|") if z.strip()]
def load_policy(p): return {r["reaction_id"]:r for r in read_tsv(p)}
def hits(rows,ids,terms):
    ids=set(ids); terms=[t.lower() for t in terms]; out=[]
    for r in rows:
        sig=r.get("signature_accession",""); ipr=r.get("interpro_id","")
        blob=" ".join(r.values()).lower()
        if (sig and sig in ids) or (ipr and ipr in ids) or any(t and t in blob for t in terms):
            txt=f"{r.get('analysis','')}:{sig or ipr}:{r.get('signature_description','') or r.get('interpro_description','')}"
            if txt not in out: out.append(txt)
    return out
def assess(rid,locus,policy,ipr):
    p=policy.get(rid)
    if not p or p.get("policy_status")=="NOT_CURATED": return "NOT_CURATED","No curated reaction-family InterPro policy"
    rr=ipr.get(locus,[])
    if not rr:return "NO_RESULT","No parsed InterPro result"
    sh=hits(rr,splitv(p.get("specific_ids")),splitv(p.get("specific_terms")))
    bh=hits(rr,splitv(p.get("broad_ids")),splitv(p.get("broad_terms")))
    ch=hits(rr,splitv(p.get("conflicting_ids")),splitv(p.get("conflicting_terms")))
    if sh and ch:return "AMBIGUOUS","SPECIFIC=["+"; ".join(sh[:3])+"] CONFLICTING=["+"; ".join(ch[:3])+"]"
    if ch:return "CONFLICTING","; ".join(ch[:5])
    if sh:return "SPECIFIC","; ".join(sh[:5])
    if bh:return "COMPATIBLE","; ".join(bh[:5])
    return "NO_MATCH","InterPro result present, but no curated rule matched"
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--rescued-evidence",required=True)
    ap.add_argument("--interpro-tsv",required=True)
    ap.add_argument("--policy",default=str(Path(__file__).with_name("reaction_interpro_policy_v054.tsv")))
    ap.add_argument("--out-prefix",default="evidence_policy_v054")
    a=ap.parse_args()
    ev=read_tsv(a.rescued_evidence); ipr=parse_ipr(a.interpro_tsv); policy=load_policy(a.policy)
    out=[]
    for r in ev:
        st,detail=assess(r.get("reaction_id",""),r.get("locus_tag",""),policy,ipr)
        x=dict(r); x["policy_interpro_status"]=st; x["policy_interpro_detail"]=detail
        x["policy_recommendation"]={"SPECIFIC":"SUPPORT_EXACT_REACTION","COMPATIBLE":"FAMILY_COMPATIBLE_ONLY","CONFLICTING":"REVIEW_OR_DOWNGRADE","AMBIGUOUS":"MANUAL_REVIEW_REQUIRED","NOT_CURATED":"NO_POLICY_YET","NO_RESULT":"INTERPRO_UNAVAILABLE","NO_MATCH":"DIAMOND_ONLY_NO_POLICY_MATCH"}[st]
        out.append(x)
    of=a.out_prefix+"_evidence_assessed.tsv"
    with open(of,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(out)
    counts=defaultdict(lambda:defaultdict(int)); loci=defaultdict(lambda:defaultdict(set))
    for r in out:
        rid=r["reaction_id"]; st=r["policy_interpro_status"]; counts[rid][st]+=1; loci[rid][st].add(r["locus_tag"])
    sf=a.out_prefix+"_reaction_summary.tsv"
    fields=["reaction_id","SPECIFIC","COMPATIBLE","CONFLICTING","AMBIGUOUS","NOT_CURATED","NO_RESULT","NO_MATCH","specific_loci","conflicting_loci","ambiguous_loci"]
    with open(sf,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader()
        for rid in sorted(counts):
            d=counts[rid]
            w.writerow({"reaction_id":rid,"SPECIFIC":d["SPECIFIC"],"COMPATIBLE":d["COMPATIBLE"],"CONFLICTING":d["CONFLICTING"],"AMBIGUOUS":d["AMBIGUOUS"],"NOT_CURATED":d["NOT_CURATED"],"NO_RESULT":d["NO_RESULT"],"NO_MATCH":d["NO_MATCH"],"specific_loci":"|".join(sorted(loci[rid]["SPECIFIC"])),"conflicting_loci":"|".join(sorted(loci[rid]["CONFLICTING"])),"ambiguous_loci":"|".join(sorted(loci[rid]["AMBIGUOUS"]))})
    print("Evidence policy evaluation complete.")
    print("Wrote:",of); print("Wrote:",sf)
if __name__=="__main__": main()
