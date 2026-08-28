#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,re
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

def read_tsv(p):
    with open(p,encoding="utf-8") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def norm(x):
    return re.sub(r"\s+"," ",str(x or "").strip().lower())

def active_row(r):
    return str(r.get("active","yes")).strip().lower() not in {"no","false","0","inactive"}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("genbank")
    ap.add_argument("--db",required=True)
    ap.add_argument("--out",required=True)
    a=ap.parse_args()

    db=Path(a.db)
    reactions={r["reaction_id"]:r for r in read_tsv(db/"reactions.tsv") if active_row(r)}
    aliases=read_tsv(db/"enzyme_aliases.tsv")

    alias_map=defaultdict(list)
    for r in aliases:
        rid=r.get("reaction_id","")
        if rid not in reactions: continue
        mode=norm(r.get("match_mode",""))
        conf=norm(r.get("confidence",""))
        alias=norm(r.get("alias",""))
        atype=norm(r.get("alias_type",""))
        if alias and mode=="exact" and conf in {"high","specific","curated",""}:
            alias_map[(atype,alias)].append(rid)

    unique_exact={k:v[0] for k,v in alias_map.items() if len(set(v))==1}
    ambiguous_exact={k:sorted(set(v)) for k,v in alias_map.items() if len(set(v))>1}

    rows=[]
    for rec in SeqIO.parse(a.genbank,"genbank"):
        for ft in rec.features:
            if ft.type!="CDS": continue
            q=ft.qualifiers
            tag=(q.get("locus_tag") or q.get("protein_id") or [""])[0]
            if not tag: continue
            gene=(q.get("gene") or [""])[0]
            product=(q.get("product") or [""])[0]
            candidates=[]
            if gene:
                candidates += [("gene",norm(gene)),("gene_name",norm(gene))]
            if product:
                candidates += [("product",norm(product)),("protein",norm(product)),("enzyme",norm(product))]

            matched=set()
            ambiguous_rids=set()
            matched_keys=defaultdict(list)
            for key in candidates:
                if key in unique_exact:
                    rid=unique_exact[key]
                    matched.add(rid)
                    matched_keys[rid].append(key)
                if key in ambiguous_exact:
                    ambiguous_rids.update(ambiguous_exact[key])

            for rid in sorted(matched):
                if rid in ambiguous_rids:
                    continue
                keys=matched_keys[rid]
                rows.append({
                    "reaction_id":rid,
                    "locus_tag":tag,
                    "gene":gene,
                    "product":product,
                    "contig":rec.id,
                    "annotation_evidence":"UNIQUE_EXACT_ALIAS",
                    "matched_alias_types":"|".join(sorted({k[0] for k in keys})),
                    "matched_aliases":"|".join(sorted({k[1] for k in keys})),
                    "resolved_for_same_reaction":"YES",
                })

    fields=["reaction_id","locus_tag","gene","product","contig",
            "annotation_evidence","matched_alias_types","matched_aliases",
            "resolved_for_same_reaction"]
    with open(a.out,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows,key=lambda x:(x["reaction_id"],x["locus_tag"])))

    print("Resolved Annotation Scanner v0.8.10")
    print("Unique exact reaction-locus pairs:",len(rows))
    print("Wrote:",a.out)

if __name__=="__main__":
    main()
