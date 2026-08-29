#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, html
from collections import defaultdict, Counter
from pathlib import Path

VERSION = "0.8.9-fixed1"

def tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def esc(x):
    return html.escape(str(x or ""))

def intish(x, default=0):
    try: return int(float(x))
    except Exception: return default

def truth(x):
    return str(x or "").strip().lower() in {"true","yes","1","supported"}

def re_safe(x):
    return ''.join(c if c.isalnum() else '_' for c in str(x))

def rx_state(r):
    if not r:
        return "NOT_DETECTED"
    state = str(r.get("integrated_state","") or "")
    review = str(r.get("review_flag","")).upper() == "YES"
    relationship = str(r.get("evidence_relationship","") or "")
    if state in {"NO_CANDIDATE","NOT_DETECTED",""} and not truth(r.get("candidate_present")):
        return "NOT_DETECTED"
    # A clean adjudicated SUPPORTED_* state remains supported even when separate
    # review-level evidence is also present. Review evidence itself never counts as support.
    if state.startswith("SUPPORTED"):
        return "SUPPORTED"
    if review or state == "CANDIDATE_CONFLICTING":
        return "REVIEW"
    if truth(r.get("candidate_present")) or "CANDIDATE" in state:
        return "CANDIDATE"
    if state.startswith("SUPPORTED"):
        return "SUPPORTED"
    return "NOT_DETECTED"

def rx_label(state):
    return {"SUPPORTED":"Supported","REVIEW":"Review","CANDIDATE":"Candidate only","NOT_DETECTED":"Not detected"}[state]

def rx_class(state):
    return {"SUPPORTED":"s-supported","PARTIAL":"s-partial","REVIEW":"s-review","CANDIDATE":"s-candidate","NOT_DETECTED":"s-missing"}.get(state,"s-missing")

def confidence_text(r):
    if not r: return "No evidence"
    st = rx_state(r)
    if st == "REVIEW": return "Review required"
    if st == "NOT_DETECTED": return "No accepted evidence"
    sources = str(r.get("evidence_sources","") or "")
    relationship = str(r.get("evidence_relationship","") or "")
    if st == "SUPPORTED":
        if "SAME_LOCUS" in relationship: return "Multi-source, same locus"
        if "MIXED" in relationship: return "Supported; mixed-locus provenance"
        if "GENBANK_EXACT_ALIAS" in sources and "HYBRID_RESCUE" in sources: return "Multi-source"
        if "HYBRID_RESCUE" in sources: return "Sequence/domain supported"
        if "GENBANK_EXACT_ALIAS" in sources: return "Annotation supported"
        return "Supported"
    return "Candidate evidence"

def pathway_result(ps, defs, rev):
    req = intish(ps.get("required_units", ps.get("required_reactions", 0)))
    sup = intish(ps.get("supported_units", ps.get("supported_reactions", 0)))
    review = intish(ps.get("review_units", 0))
    cand = intish(ps.get("candidate_units", ps.get("candidate_reactions", 0)))
    extra = max(0, cand - sup)
    if req and sup == req:
        status, label = "SUPPORTED", "Supported"
    elif sup > 0:
        status, label = "PARTIAL", "Partial"
    elif review > 0:
        status, label = "REVIEW", "Review evidence only"
    elif cand > 0:
        status, label = "CANDIDATE", "Candidate only"
    else:
        status, label = "NOT_DETECTED", "Not detected"

    optional_states=[]
    for d in defs:
        required=str(d.get("required","")).strip().lower() in {"yes","true","1","required"}
        if not required:
            optional_states.append(rx_state(rev.get(d.get("reaction_id",""),{})))
    opt_sup=sum(x=="SUPPORTED" for x in optional_states)
    opt_total=len(optional_states)

    if status=="SUPPORTED":
        text=f"All {req} required function(s) are supported."
    elif status=="PARTIAL":
        text=f"{sup} of {req} required function(s) are supported."
    elif status=="REVIEW":
        text=f"No required function is accepted as supported; {review} function(s) remain review-level."
    elif status=="CANDIDATE":
        text="No required function is accepted as supported; candidate evidence exists only."
    else:
        text="No required pathway function is currently supported."
    if extra:
        text += f" {extra} additional required function(s) have candidate/review evidence."
    if opt_total:
        text += f" Optional/supporting functions supported: {opt_sup}/{opt_total}."
    return dict(status=status,label=label,req=req,sup=sup,review=review,extra=extra,opt_sup=opt_sup,opt_total=opt_total,text=text)

def synthetic_capability(name,cid,rid,rev,supported_text,missing_text):
    r=rev.get(rid,{})
    st=rx_state(r)
    if st=="SUPPORTED":
        result,support,interp="Supported","1/1",supported_text
    elif st=="REVIEW":
        result,support,interp="Review evidence only","0/1",missing_text+" Evidence remains review-level."
    elif st=="CANDIDATE":
        result,support,interp="Candidate only","0/1",missing_text+" Candidate evidence is not counted as support."
    else:
        result,support,interp="Not detected","0/1",missing_text
    return dict(id=cid,name=name,supported=support,review="1" if st=="REVIEW" else "0",optional="—",
                result=result,status=st,interpretation=interp,reaction_ids=[rid],synthetic=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--db",required=True)
    ap.add_argument("--pathway-summary",required=True)
    ap.add_argument("--reaction-evidence",required=True)
    ap.add_argument("--cross-reaction-decisions")
    ap.add_argument("--out-html",required=True)
    ap.add_argument("--sample-name",default="")
    a=ap.parse_args()

    db=Path(a.db)
    paths=tsv(db/"pathways.tsv")
    reactions=tsv(db/"reactions.tsv")
    P=tsv(a.pathway_summary)
    R=tsv(a.reaction_evidence)
    D=tsv(a.cross_reaction_decisions) if a.cross_reaction_decisions and Path(a.cross_reaction_decisions).exists() else []

    rdef={r.get("reaction_id",""):r for r in reactions}
    rev={r.get("reaction_id",""):r for r in R}
    psum={r.get("pathway_id",""):r for r in P}
    psteps=defaultdict(list)
    for p in paths: psteps[p.get("pathway_id","")].append(p)
    for pid in psteps:
        psteps[pid].sort(key=lambda x:(intish(x.get("step_order"),9999),x.get("reaction_id","")))

    rows=[]; suppress=set()
    for pid,defs in psteps.items():
        rids={d.get("reaction_id","") for d in defs}
        if "BRXN_ARS_C" in rids and "BRXN_ARS_EFFLUX" in rids:
            suppress.add(pid)
            rows.append(synthetic_capability("Arsenate reduction","CAP_ARSENATE_REDUCTION","BRXN_ARS_C",rev,
                "Arsenate reduction potential is supported.","Arsenate reduction is not currently supported."))
            rows.append(synthetic_capability("Arsenic resistance / efflux","CAP_ARSENIC_EFFLUX","BRXN_ARS_EFFLUX",rev,
                "Arsenic resistance/efflux potential is supported.","Arsenic resistance/efflux is not currently supported."))

    for pid,ps in psum.items():
        if pid in suppress: continue
        defs=psteps.get(pid,[])
        pr=pathway_result(ps,defs,rev)
        rows.append(dict(
            id=pid,
            name=ps.get("pathway_name","") or (defs[0].get("pathway_name","") if defs else pid),
            supported=f"{pr['sup']}/{pr['req']}",
            review=str(pr["review"]),
            optional=f"{pr['opt_sup']}/{pr['opt_total']}" if pr["opt_total"] else "—",
            result=pr["label"],status=pr["status"],interpretation=pr["text"],
            reaction_ids=[d.get("reaction_id","") for d in defs],synthetic=False
        ))

    rank={"SUPPORTED":0,"PARTIAL":1,"REVIEW":2,"CANDIDATE":3,"NOT_DETECTED":4}
    rows.sort(key=lambda x:(rank.get(x["status"],9),x["name"]))
    counts=Counter(r["status"] for r in rows)

    summary_rows=[]
    for r in rows:
        detail_id="detail_"+re_safe(r["id"])
        summary_rows.append(
            f'<tr data-status="{esc(r["status"])}"><td><a href="#{esc(detail_id)}"><b>{esc(r["name"])}</b></a>'
            f'<div class="small">{esc(r["id"])}</div></td><td class="num">{esc(r["supported"])}</td>'
            f'<td class="num">{esc(r["review"])}</td><td class="num">{esc(r["optional"])}</td>'
            f'<td><span class="pill {esc(rx_class(r["status"]))}">{esc(r["result"])}</span></td>'
            f'<td>{esc(r["interpretation"])}</td></tr>'
        )

    reaction_rows=[]
    for r in R:
        rid=r.get("reaction_id",""); st=rx_state(r)
        enzyme=r.get("enzyme_name","") or rdef.get(rid,{}).get("enzyme_name","")
        supported=r.get("supported_loci","")
        review_loci=r.get("review_loci","") or r.get("ambiguous_loci","")
        candidate=r.get("candidate_only_loci","")
        if not supported and st=="SUPPORTED":
            supported=r.get("adjudicated_loci","") or r.get("all_evidence_loci","")
        if not candidate and st=="CANDIDATE":
            candidate=r.get("adjudicated_loci","") or r.get("all_evidence_loci","")
        loci=supported or review_loci or candidate or "—"
        reaction_rows.append(
            f'<tr><td><b>{esc(rid)}</b><div class="small">{esc(enzyme)}</div></td>'
            f'<td><span class="pill {esc(rx_class(st))}">{esc(rx_label(st))}</span></td>'
            f'<td>{esc(loci)}</td><td>{esc(confidence_text(r))}</td>'
            f'<td class="num">{"Yes" if st=="SUPPORTED" else "No"}</td></tr>'
        )

    row_by_id={r["id"]:r for r in rows}
    details=[]
    for r in rows:
        if r["synthetic"]:
            rid=r["reaction_ids"][0]; rr=rev.get(rid,{})
            st=rx_state(rr); enzyme=rr.get("enzyme_name","") or rdef.get(rid,{}).get("enzyme_name","")
            details.append(
                f'<details id="detail_{esc(re_safe(r["id"]))}" class="detail"><summary>{esc(r["name"])} - {esc(r["result"])}</summary>'
                f'<div class="interpret">{esc(r["interpretation"])}</div>'
                f'<div class="mini-step {esc(rx_class(st))}"><b>{esc(rid)}</b><div>{esc(enzyme)}</div>'
                f'<div class="small">{esc(rx_label(st))} · {esc(confidence_text(rr))}</div></div></details>'
            )

    for pid,ps in psum.items():
        if pid in suppress or pid not in row_by_id: continue
        r=row_by_id[pid]; defs=psteps.get(pid,[])
        flow=[]; i=0
        while i<len(defs):
            d=defs[i]; ag=d.get("alternative_group","").strip()
            if ag:
                group=[]; j=i
                while j<len(defs) and defs[j].get("alternative_group","").strip()==ag:
                    group.append(defs[j]); j+=1
                alts=[]
                for g in group:
                    rid=g.get("reaction_id",""); rr=rev.get(rid,{})
                    st=rx_state(rr); enzyme=rr.get("enzyme_name","") or rdef.get(rid,{}).get("enzyme_name","")
                    alts.append(f'<div class="mini-step {esc(rx_class(st))}"><b>{esc(rid)}</b><div>{esc(enzyme)}</div><div class="small">{esc(rx_label(st))}</div></div>')
                flow.append('<div class="alt-box"><div class="alt-title">One of these alternatives</div><div class="alt-grid">'+''.join(alts)+'</div></div>')
                i=j
            else:
                rid=d.get("reaction_id",""); rr=rev.get(rid,{})
                st=rx_state(rr); enzyme=rr.get("enzyme_name","") or rdef.get(rid,{}).get("enzyme_name","")
                required=str(d.get("required","")).strip().lower() in {"yes","true","1","required"}
                flow.append(f'<div class="mini-step {esc(rx_class(st))}"><b>{esc(rid)}</b><div>{esc(enzyme)}</div><div class="small">{"Required" if required else "Optional/supporting"} · {esc(rx_label(st))}</div></div>')
                i+=1
            if i<len(defs): flow.append('<div class="arrow">&rarr;</div>')
        details.append(
            f'<details id="detail_{esc(re_safe(pid))}" class="detail"><summary>{esc(r["name"])} - {esc(r["result"])}</summary>'
            f'<div class="interpret">{esc(r["interpretation"])}</div>'
            f'<div class="flow">{"".join(flow) if flow else "<div class=small>No pathway-flow definition available.</div>"}</div></details>'
        )

    flagged=[d for d in D if d.get("decision") in {"WINNER_KEEP","LOSER_REMOVE","AMBIGUOUS_REVIEW"}]
    review_html=""
    if flagged:
        trs=''.join(
            f'<tr><td>{esc(d.get("locus_tag"))}</td><td>{esc(d.get("reaction_id"))}</td><td>{esc(d.get("decision"))}</td>'
            f'<td>{esc(d.get("winner_reaction") or "—")}</td><td>{esc(d.get("product"))}</td></tr>' for d in flagged
        )
        review_html='<section><h2>Cross-reaction review</h2><p class="muted">Shared proteins are not automatically counted for several competing reactions.</p><div class="table-wrap"><table><thead><tr><th>Locus</th><th>Reaction</th><th>Decision</th><th>Winner</th><th>Product</th></tr></thead><tbody>'+trs+'</tbody></table></div></section>'

    sample=esc(a.sample_name or Path(a.pathway_summary).stem)

    html_doc=f'''<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{sample} - Bioremediation report</title>
<style>
:root{{--bg:#f7f8fa;--card:#fff;--text:#17212b;--muted:#66717d;--border:#d8dee6;--green:#eaf6ee;--greenb:#6fa47f;--blue:#edf4ff;--blueb:#6f92bd;--orange:#fff1eb;--orangeb:#c77f62;--purple:#f4efff;--purpleb:#9276bc;--gray:#f1f3f5;--grayb:#aab2ba}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);font-family:Segoe UI,Arial,sans-serif;color:var(--text);line-height:1.45}}
main{{max-width:1380px;margin:auto;padding:24px 16px 60px}} h1{{margin:3px 0 4px;font-size:28px}} h2{{margin:28px 0 10px;font-size:20px}}
.hero{{background:var(--card);border:1px solid var(--border);border-radius:16px;padding:20px}}
.summary{{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:8px;margin-top:16px}} .summary div{{border:1px solid var(--border);border-radius:10px;padding:10px;background:#fafbfc}} .summary b{{font-size:20px;display:block}} .summary span,.small,.muted{{font-size:12px;color:var(--muted)}}
.controls{{display:flex;gap:8px;flex-wrap:wrap;margin:16px 0}} button{{border:1px solid var(--border);background:white;padding:8px 12px;border-radius:9px;cursor:pointer}}
.table-wrap{{overflow:auto;background:white;border:1px solid var(--border);border-radius:14px}} table{{border-collapse:collapse;width:100%;min-width:980px;font-size:12px}} th,td{{padding:10px;border-bottom:1px solid var(--border);vertical-align:top;text-align:left}} th{{background:#f1f4f7;position:sticky;top:0;z-index:1}} td.num{{text-align:center;white-space:nowrap}}
a{{color:inherit}} .pill{{display:inline-block;padding:4px 8px;border:1px solid;border-radius:999px;font-weight:700;font-size:11px;white-space:nowrap}}
.s-supported{{background:var(--green);border-color:var(--greenb)}} .s-partial{{background:var(--blue);border-color:var(--blueb)}} .s-review{{background:var(--orange);border-color:var(--orangeb)}} .s-candidate{{background:var(--purple);border-color:var(--purpleb)}} .s-missing{{background:var(--gray);border-color:var(--grayb)}}
.detail{{background:white;border:1px solid var(--border);border-radius:12px;padding:12px 14px;margin:8px 0;scroll-margin-top:12px}} summary{{font-weight:800;cursor:pointer}} .interpret{{padding:9px 11px;background:#f7f8fa;border-radius:9px;margin:10px 0;font-size:13px}}
.flow{{display:flex;align-items:center;gap:8px;overflow-x:auto;padding:6px 0}} .mini-step{{min-width:190px;border:2px solid var(--grayb);border-radius:11px;padding:9px;background:var(--gray)}} .mini-step.s-supported{{background:var(--green);border-color:var(--greenb)}} .mini-step.s-review{{background:var(--orange);border-color:var(--orangeb)}} .mini-step.s-candidate{{background:var(--purple);border-color:var(--purpleb)}} .mini-step.s-missing{{background:var(--gray);border-color:var(--grayb)}} .arrow{{font-size:20px;color:var(--muted)}} .alt-box{{border:2px dashed var(--border);border-radius:12px;padding:9px;min-width:410px}} .alt-title{{font-size:11px;font-weight:800;text-transform:uppercase;letter-spacing:.05em;margin-bottom:7px}} .alt-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:7px}}
.hidden{{display:none!important}} @media(max-width:760px){{.summary{{grid-template-columns:repeat(2,minmax(0,1fr))}}}} @media print{{body{{background:white}} .controls{{display:none}}}}
</style></head>
<body><main>
<section class="hero"><div class="muted">Bioremediation Unified Runner - Table-first report v{VERSION}</div><h1>{sample}</h1>
<div class="muted">Pathway completeness and evidence strength are separate. Review/candidate evidence never counts as supported.</div>
<div class="summary"><div><b>{counts["SUPPORTED"]}</b><span>Supported</span></div><div><b>{counts["PARTIAL"]}</b><span>Partial</span></div><div><b>{counts["REVIEW"]}</b><span>Review evidence only</span></div><div><b>{counts["CANDIDATE"]}</b><span>Candidate only</span></div><div><b>{counts["NOT_DETECTED"]}</b><span>Not detected</span></div></div></section>

<section><h2>Results summary</h2><div class="controls"><button data-filter="all">All</button><button data-filter="detected">Hide not detected</button><button id="print">Print / Save PDF</button></div>
<div class="table-wrap"><table id="results"><thead><tr><th>Pathway / capability</th><th>Required supported</th><th>Review units</th><th>Optional supported</th><th>Result</th><th>Biological interpretation</th></tr></thead><tbody>{''.join(summary_rows)}</tbody></table></div></section>

<section><h2>Reaction evidence</h2><p class="muted">Counted = Yes means the reaction is accepted as supported. Review and candidate rows are preserved but do not increase pathway support.</p>
<div class="table-wrap"><table><thead><tr><th>Reaction</th><th>Evidence state</th><th>Locus/loci</th><th>Evidence quality</th><th>Counted?</th></tr></thead><tbody>{''.join(reaction_rows)}</tbody></table></div></section>

<section><h2>Pathway details</h2><p class="muted">Expand only when you want the reaction structure. These diagrams support the table rather than replace it.</p>{''.join(details)}</section>
{review_html}
</main><script>
(function(){{const rows=[...document.querySelectorAll('#results tbody tr')];document.querySelectorAll('[data-filter]').forEach(b=>b.addEventListener('click',()=>{{const mode=b.dataset.filter;rows.forEach(r=>r.classList.toggle('hidden',mode==='detected'&&r.dataset.status==='NOT_DETECTED'));}}));document.getElementById('print').addEventListener('click',()=>window.print());}})();
</script></body></html>'''

    Path(a.out_html).write_text(html_doc,encoding="utf-8")
    print("Table-first Bioremediation Reporter v0.8.9")
    print("Summary rows:",len(rows))
    print("Wrote:",a.out_html)

if __name__=="__main__":
    main()
