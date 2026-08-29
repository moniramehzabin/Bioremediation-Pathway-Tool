#!/usr/bin/env python3
"""
Cross-reaction locus adjudicator v0.8.3

Purpose
-------
Prevent one broadly conserved protein from automatically supporting several
different pollutant-specific reactions.

Evidence hierarchy (general, not reaction-specific):
1. unique exact GenBank annotation support
2. multi-source support
3. rescue support strength / DIAMOND confidence
4. pathway-neighborhood coherence
5. protein-name/enzyme-name specificity

Rules are conservative:
- a strong winner can keep the locus;
- a weaker competing assignment loses that locus;
- close contests stay AMBIGUOUS / REVIEW;
- genomic context never creates a reaction;
- losing one locus does not erase a reaction if other non-conflicting loci remain;
- alternative pathway groups are scored as OR, not AND.
"""
from __future__ import annotations
import argparse, csv, re
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

STOP = {
    "protein","enzyme","family","putative","probable","like","activity","system",
    "subunit","component","associated","metabolic","oxidoreductase","dehydrogenase",
    "dioxygenase","monooxygenase","hydrolase","reductase","transferase"
}

def tsv(p):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def wt(p, rows, fields):
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def split_tags(x):
    if not x: return []
    x = str(x).replace(";", "|").replace(",", "|")
    return [q.strip() for q in x.split("|") if q.strip()]

def truth(x):
    return str(x or "").strip().lower() in {"true","yes","1","supported"}

def toks(x):
    return {t for t in re.findall(r"[a-z0-9]+", str(x or "").lower())
            if len(t) > 2 and t not in STOP}

def jaccard(a,b):
    A,B=toks(a),toks(b)
    if not A or not B: return 0.0
    return len(A&B)/len(A|B)

def conf_score(x):
    x=str(x or "").upper()
    if "HIGH" in x: return 20
    if "MODERATE" in x: return 10
    if "WEAK" in x: return 2
    return 0

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("genbank")
    ap.add_argument("--db", required=True)
    ap.add_argument("--integrated-reaction-evidence", required=True)
    ap.add_argument("--all-candidates", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--context-window-bp", type=int, default=20000)
    ap.add_argument("--winner-margin", type=float, default=25.0)
    a=ap.parse_args()

    db=Path(a.db)
    reactions={r["reaction_id"]:r for r in tsv(db/"reactions.tsv")
               if str(r.get("active","yes")).strip().lower() not in {"no","false","0","inactive"}}
    pathways=tsv(db/"pathways.tsv")
    evrows=tsv(a.integrated_reaction_evidence)
    ev={r["reaction_id"]:r for r in evrows}
    candidates=tsv(a.all_candidates)

    # GenBank locus coordinates and products
    loc={}
    bycontig=defaultdict(list)
    for rec in SeqIO.parse(a.genbank,"genbank"):
        for ft in rec.features:
            if ft.type!="CDS": continue
            q=ft.qualifiers
            tag=(q.get("locus_tag") or q.get("protein_id") or [""])[0]
            if not tag: continue
            row={
                "locus_tag":tag,"contig":rec.id,
                "start":int(ft.location.start)+1,"end":int(ft.location.end),
                "gene":(q.get("gene") or [""])[0],
                "product":(q.get("product") or [""])[0],
            }
            loc[tag]=row
            bycontig[rec.id].append(row)
    for c in bycontig:
        bycontig[c].sort(key=lambda x:x["start"])

    # Candidate metadata by reaction+locus
    cmeta=defaultdict(list)
    for r in candidates:
        cmeta[(r.get("reaction_id",""),r.get("locus_tag",""))].append(r)

    # Pathway membership for context scoring
    rx_to_paths=defaultdict(set)
    path_to_rx=defaultdict(set)
    for p in pathways:
        pid=p.get("pathway_id",""); rid=p.get("reaction_id","")
        if pid and rid:
            rx_to_paths[rid].add(pid)
            path_to_rx[pid].add(rid)

    # Evidence loci by reaction
    rx_loci={}
    for rid,r in ev.items():
        rx_loci[rid]=set(split_tags(r.get("all_evidence_loci","")))
    locus_rx=defaultdict(set)
    for rid,tags in rx_loci.items():
        for tag in tags:
            locus_rx[tag].add(rid)

    # Nearby reaction support map, using existing evidence only.
    locus_supported_rx=defaultdict(set)
    for rid,tags in rx_loci.items():
        for tag in tags:
            if tag in loc:
                locus_supported_rx[tag].add(rid)

    def neighborhood_score(tag,rid):
        if tag not in loc: return 0
        L=loc[tag]
        wanted=set()
        for pid in rx_to_paths.get(rid,set()):
            wanted |= path_to_rx.get(pid,set())
        if not wanted: return 0
        n=0
        for other_tag,ors in locus_supported_rx.items():
            if other_tag==tag or other_tag not in loc: continue
            O=loc[other_tag]
            if O["contig"]!=L["contig"]: continue
            if abs(O["start"]-L["start"]) <= a.context_window_bp:
                if ors & wanted:
                    n += 1
        return min(n,4)*5

    def assignment_score(tag,rid):
        r=ev.get(rid,{})
        score=0.0
        reasons=[]
        ann=set(split_tags(r.get("annotation_loci","")))
        if tag in ann:
            score += 100
            reasons.append("EXACT_ANNOTATION")
        state=r.get("integrated_state","")
        if state=="SUPPORTED_MULTI_SOURCE":
            score += 35; reasons.append("MULTI_SOURCE")
        elif state=="SUPPORTED_RESCUE":
            score += 20; reasons.append("RESCUE_SUPPORTED")
        elif state=="SUPPORTED_ANNOTATION":
            score += 25; reasons.append("ANNOTATION_SUPPORTED")
        elif "CANDIDATE" in state:
            score += 5; reasons.append("CANDIDATE_ONLY")

        cms=cmeta.get((rid,tag),[])
        if cms:
            best=max(conf_score(x.get("diamond_confidence","")) for x in cms)
            score += best
            if best: reasons.append("DIAMOND_"+("HIGH" if best==20 else "MODERATE" if best==10 else "WEAK"))
            # product vs target specificity
            ptxt=loc.get(tag,{}).get("product","")
            ename=reactions.get(rid,{}).get("enzyme_name","")
            sim=jaccard(ptxt,ename)
            if sim >= 0.40:
                score += 20; reasons.append("PRODUCT_SPECIFICITY_HIGH")
            elif sim >= 0.20:
                score += 8; reasons.append("PRODUCT_SPECIFICITY_MODERATE")

        ns=neighborhood_score(tag,rid)
        if ns:
            score += ns; reasons.append(f"PATHWAY_CONTEXT_{ns}")
        return score,reasons

    decisions=[]
    keep=defaultdict(set)
    lose=defaultdict(set)
    ambiguous=defaultdict(set)

    for tag,rids in sorted(locus_rx.items()):
        if len(rids)==1:
            rid=next(iter(rids))
            sc,rs=assignment_score(tag,rid)
            keep[rid].add(tag)
            decisions.append({
                "locus_tag":tag,"reaction_id":rid,"score":f"{sc:.1f}",
                "decision":"UNIQUE_KEEP","winner_reaction":rid,"margin":"",
                "reasons":"|".join(rs),
                "product":loc.get(tag,{}).get("product","")
            })
            continue

        scored=[]
        for rid in sorted(rids):
            sc,rs=assignment_score(tag,rid)
            scored.append((sc,rid,rs))
        scored.sort(reverse=True)
        top=scored[0]
        second=scored[1]
        margin=top[0]-second[0]

        # Exact annotation wins unless two competing reactions both have the same exact locus.
        exact_top="EXACT_ANNOTATION" in top[2]
        exact_second="EXACT_ANNOTATION" in second[2]
        clear=(margin >= a.winner_margin and not (exact_top and exact_second))
        if exact_top and not exact_second:
            clear=True

        if clear:
            winner=top[1]
            for sc,rid,rs in scored:
                if rid==winner:
                    keep[rid].add(tag); dec="WINNER_KEEP"
                else:
                    lose[rid].add(tag); dec="LOSER_REMOVE"
                decisions.append({
                    "locus_tag":tag,"reaction_id":rid,"score":f"{sc:.1f}",
                    "decision":dec,"winner_reaction":winner,
                    "margin":f"{margin:.1f}","reasons":"|".join(rs),
                    "product":loc.get(tag,{}).get("product","")
                })
        else:
            for sc,rid,rs in scored:
                ambiguous[rid].add(tag)
                decisions.append({
                    "locus_tag":tag,"reaction_id":rid,"score":f"{sc:.1f}",
                    "decision":"AMBIGUOUS_REVIEW","winner_reaction":"",
                    "margin":f"{margin:.1f}","reasons":"|".join(rs),
                    "product":loc.get(tag,{}).get("product","")
                })

    # Locus-level evidence provenance. Multi-source is only awarded when
    # annotation and countable rescue evidence converge on the SAME locus.
    rescue_supported=defaultdict(set)
    rescue_review=defaultdict(set)
    for c in candidates:
        rid=c.get("reaction_id",""); tag=c.get("locus_tag","")
        if not rid or not tag: continue
        adj=str(c.get("adjudication","")).upper()
        if adj in {"COUNTABLE","SUPPORTED","KEEP","ACCEPTED"}:
            rescue_supported[rid].add(tag)
        elif "REVIEW" in adj or "AMBIG" in adj:
            rescue_review[rid].add(tag)

    # Reaction-level update.
    out_ev=[]
    for r in evrows:
        rid=r["reaction_id"]
        orig=set(split_tags(r.get("all_evidence_loci","")))
        remaining=(orig - lose[rid] - ambiguous[rid]) | keep[rid]
        removed=(orig & lose[rid])
        amb=(orig & ambiguous[rid])

        nr=dict(r)
        nr["pre_adjudication_state"]=r.get("integrated_state","")

        ann=(set(split_tags(r.get("annotation_loci",""))) & remaining)
        rsc=(rescue_supported[rid] & remaining)
        same_multi=ann & rsc
        supported=(ann | rsc)
        review=((rescue_review[rid] & remaining) | amb) - supported
        candidate_only=remaining - supported - review

        nr["annotation_supported_loci"]="|".join(sorted(ann))
        nr["rescue_supported_loci"]="|".join(sorted(rsc))
        nr["same_locus_multisource_loci"]="|".join(sorted(same_multi))
        nr["supported_loci"]="|".join(sorted(supported))
        nr["review_loci"]="|".join(sorted(review))
        nr["candidate_only_loci"]="|".join(sorted(candidate_only))
        nr["adjudicated_loci"]="|".join(sorted(remaining))
        nr["removed_competing_loci"]="|".join(sorted(removed))
        nr["ambiguous_loci"]="|".join(sorted(amb))

        # Correct reaction-level evidence relationship.
        # Do not call evidence MULTI_SOURCE merely because different loci
        # supplied different source types.
        if same_multi:
            nr["integrated_state"]="SUPPORTED_MULTI_SOURCE"
            nr["evidence_relationship"]="SAME_LOCUS_MULTI_SOURCE"
        elif ann and rsc:
            nr["integrated_state"]="SUPPORTED_MIXED_LOCI"
            nr["evidence_relationship"]="MIXED_LOCI_SOURCES"
            nr["review_flag"]="YES"
        elif ann:
            nr["integrated_state"]="SUPPORTED_ANNOTATION"
            nr["evidence_relationship"]="ANNOTATION_ONLY"
        elif rsc:
            nr["integrated_state"]="SUPPORTED_RESCUE"
            nr["evidence_relationship"]="RESCUE_ONLY"
        else:
            nr["evidence_relationship"]="NO_CLEAN_SUPPORTED_LOCUS"

        if orig and not remaining:
            # Nothing clean remains. Never claim exact support.
            nr["integrated_state"]="CANDIDATE_CONFLICTING"
            nr["candidate_present"]=True
            nr["exact_supported"]=False
            nr["review_flag"]="YES"
            nr["all_evidence_loci"]="|".join(sorted(orig))
            nr["evidence_sources"]=(r.get("evidence_sources","")+"|CROSS_REACTION_CONFLICT").strip("|")
        else:
            nr["all_evidence_loci"]="|".join(sorted(remaining))
            if amb:
                nr["review_flag"]="YES"
                nr["evidence_sources"]=(r.get("evidence_sources","")+"|CROSS_REACTION_REVIEW").strip("|")
        out_ev.append(nr)

    # Alternative-aware pathway scoring.
    aev={r["reaction_id"]:r for r in out_ev}
    pgroup=defaultdict(list)
    pnames={}
    for p in pathways:
        pid=p.get("pathway_id","")
        if not pid: continue
        pnames[pid]=p.get("pathway_name","") or pnames.get(pid,"")
        pgroup[pid].append(p)

    pout=[]
    for pid,rows in pgroup.items():
        # Required units:
        # - ordinary required reaction = one unit
        # - same alternative_group = one OR unit
        units=[]
        seen_alt=set()
        for p in rows:
            required=str(p.get("required","")).lower() in {"yes","true","1","required"}
            if not required: continue
            ag=p.get("alternative_group","").strip()
            if ag:
                if ag in seen_alt: continue
                members=[x.get("reaction_id","") for x in rows
                         if x.get("alternative_group","").strip()==ag]
                units.append(("ALT:"+ag,members))
                seen_alt.add(ag)
            else:
                units.append((p.get("reaction_id",""),[p.get("reaction_id","")]))

        cand_units=0; supp_units=0; review_units=0; missing=[]
        observed_rx=set()
        for uname,members in units:
            cm=False; sm=False
            for rid in members:
                rr=aev.get(rid,{})
                state=rr.get("integrated_state","")
                cand=truth(rr.get("candidate_present")) or state not in {"","NO_CANDIDATE","NOT_DETECTED"}
                # LOGIC A SUPPORT RULE:
                # An already-adjudicated SUPPORTED_* reaction remains countable support.
                # Separate review-level evidence never creates support and does not erase clean support.
                supp=state.startswith("SUPPORTED")
                if cand:
                    cm=True; observed_rx.add(rid)
                if supp:
                    sm=True; observed_rx.add(rid)
            cand_units += int(cm); supp_units += int(sm)
            # A unit is review-level when it has candidate evidence but no countable support
            # and at least one member is explicitly review-required.
            rv=any(aev.get(rid,{}).get("review_flag")=="YES" for rid in members)
            review_units += int(cm and (not sm) and rv)
            if not cm:
                missing.append(uname+"["+ "|".join(members)+"]" if uname.startswith("ALT:") else uname)

        total=len(units)
        cpct=100*cand_units/total if total else 0
        spct=100*supp_units/total if total else 0

        # cluster summary from remaining adjudicated loci
        pts=[]
        for rid in observed_rx:
            for tag in split_tags(aev.get(rid,{}).get("all_evidence_loci","")):
                if tag in loc:
                    L=loc[tag]
                    pts.append((L["contig"],L["start"],L["end"],rid,tag))
        best=None
        for contig in {x[0] for x in pts}:
            cp=sorted([x for x in pts if x[0]==contig],key=lambda x:x[1])
            for i in range(len(cp)):
                rrset=set();mx=cp[i][2]
                for j in range(i,len(cp)):
                    mx=max(mx,cp[j][2]);span=mx-cp[i][1]+1
                    if span>a.context_window_bp: break
                    rrset.add(cp[j][3])
                    tup=(len(rrset),-span,contig,span,sorted(rrset))
                    if best is None or tup[:2]>best[:2]: best=tup

        cluster_n=best[0] if best else 0
        cluster_span=best[3] if best else ""
        cluster_contig=best[2] if best else ""
        cluster_rx="|".join(best[4]) if best else ""

        if total and supp_units==total:
            state="SUPPORTED_COMPLETE"
        elif total and cand_units==total:
            state="COMPLETE_CANDIDATE_REVIEW_REQUIRED"
        elif cluster_n>=3 and supp_units>=2:
            state="PARTIAL_WITH_STRONG_LOCAL_CONTEXT"
        elif cluster_n>=2 and cand_units>=2:
            state="PARTIAL_WITH_LOCAL_CONTEXT"
        elif cand_units:
            state="PARTIAL"
        else:
            state="NOT_DETECTED"

        pout.append({
            "pathway_id":pid,"pathway_name":pnames.get(pid,""),
            "required_units":total,
            "candidate_units":cand_units,"candidate_percent":f"{cpct:.1f}",
            "supported_units":supp_units,"supported_percent":f"{spct:.1f}",
            "review_units":review_units,
            "adjudicated_pathway_state":state,
            "cluster_reaction_count":cluster_n,
            "cluster_contig":cluster_contig,
            "cluster_span_bp":cluster_span,
            "cluster_reactions":cluster_rx,
            "missing_required_units":"|".join(missing),
            "alternative_group_rule":"OR_WITHIN_GROUP",
            "context_rule":"CONTEXT_NEVER_CREATES_REACTIONS"
        })

    efields=list(out_ev[0].keys())
    dfields=["locus_tag","reaction_id","score","decision","winner_reaction","margin","reasons","product"]
    pfields=["pathway_id","pathway_name","required_units","candidate_units","candidate_percent",
             "supported_units","supported_percent","review_units","adjudicated_pathway_state",
             "cluster_reaction_count","cluster_contig","cluster_span_bp","cluster_reactions",
             "missing_required_units","alternative_group_rule","context_rule"]

    wt(a.out_prefix+"_cross_reaction_locus_decisions.tsv",decisions,dfields)
    wt(a.out_prefix+"_adjudicated_reaction_evidence.tsv",out_ev,efields)
    wt(a.out_prefix+"_adjudicated_pathway_summary.tsv",pout,pfields)

    print("Cross-Reaction Adjudicator v0.8.3")
    print("Shared loci examined:",sum(1 for x in locus_rx.values() if len(x)>1))
    print("Winner/loser locus competitions:",
          len({d["locus_tag"] for d in decisions if d["decision"]=="WINNER_KEEP"}))
    print("Ambiguous shared loci:",
          len({d["locus_tag"] for d in decisions if d["decision"]=="AMBIGUOUS_REVIEW"}))
    print("Reaction calls downgraded to conflicting:",
          sum(r.get("integrated_state")=="CANDIDATE_CONFLICTING" and
              r.get("pre_adjudication_state")!="CANDIDATE_CONFLICTING" for r in out_ev))
    print("Alternative pathway groups are scored as OR, not AND.")
    print("Wrote:",a.out_prefix+"_cross_reaction_locus_decisions.tsv")
    print("Wrote:",a.out_prefix+"_adjudicated_reaction_evidence.tsv")
    print("Wrote:",a.out_prefix+"_adjudicated_pathway_summary.tsv")

if __name__=="__main__":
    main()
