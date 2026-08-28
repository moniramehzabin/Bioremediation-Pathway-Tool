#!/usr/bin/env python3
"""
Integrated Evidence v0.6.3

Combines:
  1) specific GenBank annotation evidence
  2) Hybrid Rescue reaction evidence
  3) genomic clustering/context

Safety invariants:
  - context NEVER creates a reaction
  - only UNIQUE exact annotation aliases can create annotation support
  - generic/ambiguous aliases are review-only, not exact support
  - rescue conflicts remain visible
"""
from __future__ import annotations
import argparse,csv,re
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

TRUE={"true","yes","1","present","supported","high"}
FALSE={"false","no","0","inactive"}

def read_tsv(p):
    if not Path(p).exists(): return []
    with open(p,encoding="utf-8") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def write_tsv(p,rows,fields):
    with open(p,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        w.writeheader();w.writerows(rows)

def yes(x): return str(x or "").strip().lower() in TRUE
def active_row(r): return str(r.get("active","yes")).strip().lower() not in FALSE
def norm(x):
    return re.sub(r"\s+"," ",str(x or "").strip().lower())
def split_tags(x):
    if not x:return []
    z=str(x).replace(";","|").replace(",","|")
    return [q.strip() for q in z.split("|") if q.strip()]

ap=argparse.ArgumentParser()
ap.add_argument("genbank")
ap.add_argument("--db",required=True)
ap.add_argument("--reaction-evidence",required=True)
ap.add_argument("--out-prefix",default="integrated_v063")
ap.add_argument("--cluster-window-bp",type=int,default=20000)
a=ap.parse_args()

db=Path(a.db)
reactions={r["reaction_id"]:r for r in read_tsv(db/"reactions.tsv") if active_row(r)}
pathways=read_tsv(db/"pathways.tsv")
aliases=read_tsv(db/"enzyme_aliases.tsv")
rescue={r["reaction_id"]:r for r in read_tsv(a.reaction_evidence)}

# -------- annotation alias policy --------
# Only exact aliases that map uniquely to one active reaction are allowed to
# create exact annotation support. Everything else is review-only.
alias_map=defaultdict(list)
for r in aliases:
    rid=r.get("reaction_id","")
    if rid not in reactions: continue
    mode=norm(r.get("match_mode",""))
    conf=norm(r.get("confidence",""))
    alias=norm(r.get("alias",""))
    atype=norm(r.get("alias_type",""))
    if not alias: continue
    if mode=="exact" and conf in {"high","specific","curated",""}:
        alias_map[(atype,alias)].append(rid)

unique_exact={k:v[0] for k,v in alias_map.items() if len(set(v))==1}
ambiguous_exact={k:sorted(set(v)) for k,v in alias_map.items() if len(set(v))>1}

# -------- parse genome annotations --------
loci={}
annot_hits=defaultdict(list)
for rec in SeqIO.parse(a.genbank,"genbank"):
    for ft in rec.features:
        if ft.type!="CDS":continue
        q=ft.qualifiers
        tag=(q.get("locus_tag") or q.get("protein_id") or [""])[0]
        if not tag:continue
        gene=(q.get("gene") or [""])[0]
        product=(q.get("product") or [""])[0]
        loci[tag]={
          "locus_tag":tag,"contig":rec.id,
          "start":int(ft.location.start)+1,"end":int(ft.location.end),
          "strand":"+" if ft.location.strand==1 else "-" if ft.location.strand==-1 else ".",
          "gene":gene,"product":product,
        }
        candidates=[]
        if gene:
            candidates += [("gene",norm(gene)),("gene_name",norm(gene))]
        if product:
            candidates += [("product",norm(product)),("protein",norm(product)),("enzyme",norm(product))]
        matched=set(); ambiguous=set()
        for key in candidates:
            if key in unique_exact: matched.add(unique_exact[key])
            if key in ambiguous_exact: ambiguous.update(ambiguous_exact[key])
        # If an exact alias uniquely identifies one reaction, count it.
        for rid in matched:
            annot_hits[rid].append({
              **loci[tag],
              "annotation_evidence":"UNIQUE_EXACT_ALIAS",
              "annotation_review":"YES" if rid in ambiguous else "NO",
            })

# -------- integrated reaction states --------
integrated=[]
for rid,rmeta in reactions.items():
    rr=rescue.get(rid,{})
    signal=str(rr.get("reaction_signal","NO_CANDIDATE")).upper()
    rescue_supported=(signal=="SUPPORTED" or int(rr.get("supported_exact_loci") or 0)>0)
    rescue_candidate=(rescue_supported or signal in {"CANDIDATE","CANDIDATE_REVIEW","CANDIDATE_CONFLICTING"}
                      or int(rr.get("candidate_loci") or 0)>0)
    ann=annot_hits.get(rid,[])
    ann_tags=sorted({x["locus_tag"] for x in ann})
    annotation_supported=bool(ann_tags)

    if rescue_supported and annotation_supported:
        state="SUPPORTED_MULTI_SOURCE"
        exact=True;candidate=True;review=False
    elif rescue_supported:
        state="SUPPORTED_RESCUE"
        exact=True;candidate=True;review=False
    elif annotation_supported:
        # Exact unique annotation is enough to establish reaction presence,
        # but provenance remains annotation rather than sequence rescue.
        state="SUPPORTED_ANNOTATION"
        exact=True;candidate=True;review=False
    elif rescue_candidate:
        state=signal if signal!="NO_CANDIDATE" else "CANDIDATE_REVIEW"
        exact=False;candidate=True;review=True
    else:
        state="NO_CANDIDATE"
        exact=False;candidate=False;review=False

    all_tags=sorted(set(ann_tags+split_tags(rr.get("supported_locus_tags"))+split_tags(rr.get("candidate_locus_tags"))))
    integrated.append({
      "reaction_id":rid,
      "enzyme_name":rmeta.get("enzyme_name",""),
      "integrated_state":state,
      "candidate_present":candidate,
      "exact_supported":exact,
      "review_flag":"YES" if review else "NO",
      "annotation_supported":annotation_supported,
      "annotation_loci":"|".join(ann_tags),
      "rescue_signal":signal,
      "rescue_best_sequence_evidence":rr.get("best_sequence_evidence",""),
      "rescue_candidate_loci":rr.get("candidate_loci","0"),
      "rescue_supported_loci":rr.get("supported_exact_loci","0"),
      "all_evidence_loci":"|".join(all_tags),
      "evidence_sources":"|".join(x for x in [
          "GENBANK_EXACT_ALIAS" if annotation_supported else "",
          "HYBRID_RESCUE" if rescue_candidate else ""
      ] if x),
    })

iby={r["reaction_id"]:r for r in integrated}

# -------- pathway integration + context --------
by_path=defaultdict(list);pnames={}
for p in pathways:
    pid=p.get("pathway_id","")
    if not pid:continue
    pnames[pid]=p.get("pathway_name","") or pnames.get(pid,"")
    if str(p.get("required","")).strip().lower() in {"yes","true","1","required"}:
        rid=p.get("reaction_id","")
        if rid and rid in reactions and rid not in by_path[pid]:
            by_path[pid].append(rid)

psum=[]
for pid,rids in by_path.items():
    total=len(rids)
    cand=[x for x in rids if iby.get(x,{}).get("candidate_present")]
    supp=[x for x in rids if iby.get(x,{}).get("exact_supported")]
    missing=[x for x in rids if x not in cand]

    pts=[]
    for rid in cand:
        for tag in split_tags(iby[rid].get("all_evidence_loci","")):
            if tag in loci:
                L=loci[tag];pts.append((L["contig"],L["start"],L["end"],rid,tag))
    best=None
    for contig in {x[0] for x in pts}:
        cp=sorted([x for x in pts if x[0]==contig],key=lambda x:x[1])
        for i in range(len(cp)):
            rrset=set();tags=[];mx=cp[i][2]
            for j in range(i,len(cp)):
                mx=max(mx,cp[j][2]);span=mx-cp[i][1]+1
                if span>a.cluster_window_bp:break
                rrset.add(cp[j][3]);tags.append(cp[j][4])
                tup=(len(rrset),-span,contig,span,sorted(rrset),tags[:])
                if best is None or tup[:2]>best[:2]:best=tup
    cluster_n=best[0] if best else 0
    cluster_span=best[3] if best else ""
    cluster_contig=best[2] if best else ""
    cluster_rx="|".join(best[4]) if best else ""

    cpct=100*len(cand)/total if total else 0
    spct=100*len(supp)/total if total else 0
    if total and len(supp)==total:
        state="SUPPORTED_COMPLETE"
    elif total and len(cand)==total:
        state="COMPLETE_CANDIDATE_REVIEW_REQUIRED"
    elif cluster_n>=3 and len(supp)>=2:
        state="PARTIAL_WITH_STRONG_LOCAL_CONTEXT"
    elif cluster_n>=2 and len(cand)>=2:
        state="PARTIAL_WITH_LOCAL_CONTEXT"
    elif cand:
        state="PARTIAL"
    else:
        state="NOT_DETECTED"

    psum.append({
      "pathway_id":pid,"pathway_name":pnames.get(pid,""),
      "required_reactions":total,
      "candidate_reactions":len(cand),"candidate_percent":f"{cpct:.1f}",
      "supported_reactions":len(supp),"supported_percent":f"{spct:.1f}",
      "integrated_pathway_state":state,
      "cluster_reaction_count":cluster_n,"cluster_contig":cluster_contig,
      "cluster_span_bp":cluster_span,"cluster_reactions":cluster_rx,
      "missing_required_reactions":"|".join(missing),
      "context_rule":"CONTEXT_NEVER_CREATES_REACTIONS",
    })

rfields=["reaction_id","enzyme_name","integrated_state","candidate_present","exact_supported","review_flag",
         "annotation_supported","annotation_loci","rescue_signal","rescue_best_sequence_evidence",
         "rescue_candidate_loci","rescue_supported_loci","all_evidence_loci","evidence_sources"]
pfields=["pathway_id","pathway_name","required_reactions","candidate_reactions","candidate_percent",
         "supported_reactions","supported_percent","integrated_pathway_state","cluster_reaction_count",
         "cluster_contig","cluster_span_bp","cluster_reactions","missing_required_reactions","context_rule"]
write_tsv(a.out_prefix+"_integrated_reaction_evidence.tsv",integrated,rfields)
write_tsv(a.out_prefix+"_integrated_pathway_summary.tsv",psum,pfields)

print("Integrated Evidence v0.6.3")
print("Active reactions:",len(reactions))
print("Unique exact annotation-supported reactions:",sum(1 for r in integrated if r["annotation_supported"]))
for r in psum:
    if r["pathway_id"]=="BPWY_HPA_HOMOPROTOCATECHUATE":
        print("HPA:")
        print(f'  candidate: {r["candidate_reactions"]}/{r["required_reactions"]} ({r["candidate_percent"]}%)')
        print(f'  supported: {r["supported_reactions"]}/{r["required_reactions"]} ({r["supported_percent"]}%)')
        print("  state:",r["integrated_pathway_state"])
        print("  cluster reactions:",r["cluster_reaction_count"])
        print("  cluster:",r["cluster_contig"],"span",r["cluster_span_bp"],"bp")
        print("  missing:",r["missing_required_reactions"])
print("IMPORTANT: context never creates reactions; only unique exact annotation aliases can add annotation support.")
print("Wrote:",a.out_prefix+"_integrated_reaction_evidence.tsv")
print("Wrote:",a.out_prefix+"_integrated_pathway_summary.tsv")
