#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, html
from collections import defaultdict, Counter
from pathlib import Path

VERSION='0.8.5'

def tsv(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def esc(x): return html.escape(str(x or ''))
def split_pipe(x): return [s.strip() for s in str(x or '').replace(';','|').split('|') if s.strip()]
def truth(x): return str(x or '').strip().lower() in {'true','yes','1','supported'}
def intish(x,d=0):
    try: return int(float(x))
    except Exception: return d

STATE_RANK={'SUPPORTED_COMPLETE':0,'COMPLETE_CANDIDATE_REVIEW_REQUIRED':1,'PARTIAL_WITH_STRONG_LOCAL_CONTEXT':2,'PARTIAL_WITH_LOCAL_CONTEXT':3,'PARTIAL':4,'NOT_DETECTED':5}
STATE_LABEL={'SUPPORTED_COMPLETE':'Supported complete','COMPLETE_CANDIDATE_REVIEW_REQUIRED':'Candidate-complete · review','PARTIAL_WITH_STRONG_LOCAL_CONTEXT':'Partial · strong local context','PARTIAL_WITH_LOCAL_CONTEXT':'Partial · local context','PARTIAL':'Partial','NOT_DETECTED':'Not detected'}

def quality(r):
    state=r.get('integrated_state',''); rel=r.get('evidence_relationship','')
    if r.get('review_flag')=='YES' or 'MIXED' in rel or 'CONFLICT' in state or r.get('ambiguous_loci'):
        return 'review','Review'
    if state=='SUPPORTED_MULTI_SOURCE' and 'SAME_LOCUS' in rel:
        return 'strong','Multi-source · same locus'
    if state in {'SUPPORTED_MULTI_SOURCE','SUPPORTED_RESCUE'}: return 'strong','Strong support'
    if state=='SUPPORTED_ANNOTATION': return 'annotation','Annotation supported'
    if truth(r.get('exact_supported')): return 'supported','Supported'
    if truth(r.get('candidate_present')) or 'CANDIDATE' in state: return 'candidate','Candidate'
    return 'missing','Missing'

def pathway_conf(rows):
    qs=[quality(r)[0] for r in rows if r]
    if not qs: return 'No evidence'
    if 'review' in qs: return 'Review required'
    if 'strong' in qs: return 'Strong evidence'
    if 'supported' in qs: return 'Supported'
    if 'annotation' in qs: return 'Annotation-level evidence'
    if 'candidate' in qs: return 'Candidate-level evidence'
    return 'No evidence'

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--db',required=True)
    ap.add_argument('--pathway-summary',required=True)
    ap.add_argument('--reaction-evidence',required=True)
    ap.add_argument('--cross-reaction-decisions')
    ap.add_argument('--out-html',required=True)
    ap.add_argument('--sample-name',default='')
    a=ap.parse_args()
    db=Path(a.db)
    pdef=tsv(db/'pathways.tsv'); rdefrows=tsv(db/'reactions.tsv')
    P=tsv(a.pathway_summary); R=tsv(a.reaction_evidence)
    D=tsv(a.cross_reaction_decisions) if a.cross_reaction_decisions and Path(a.cross_reaction_decisions).exists() else []
    rdef={r.get('reaction_id',''):r for r in rdefrows}; rev={r.get('reaction_id',''):r for r in R}; ps={p.get('pathway_id',''):p for p in P}
    steps=defaultdict(list)
    for p in pdef: steps[p.get('pathway_id','')].append(p)
    for pid in steps: steps[pid].sort(key=lambda x:(intish(x.get('step_order'),9999),x.get('reaction_id','')))
    pids=sorted(set(ps)|set(steps),key=lambda pid:(STATE_RANK.get(ps.get(pid,{}).get('adjudicated_pathway_state',''),9),ps.get(pid,{}).get('pathway_name',pid)))
    counts=Counter(ps.get(pid,{}).get('adjudicated_pathway_state','UNKNOWN') for pid in pids)
    detected=sum(1 for pid in pids if ps.get(pid,{}).get('adjudicated_pathway_state')!='NOT_DETECTED')
    partial=sum(counts[x] for x in ('PARTIAL','PARTIAL_WITH_LOCAL_CONTEXT','PARTIAL_WITH_STRONG_LOCAL_CONTEXT'))
    review_count=counts['COMPLETE_CANDIDATE_REVIEW_REQUIRED']+sum(1 for r in R if r.get('review_flag')=='YES')

    cards=[]
    for pid in pids:
        s=ps.get(pid,{}); state=s.get('adjudicated_pathway_state','UNKNOWN'); defs=steps.get(pid,[])
        pname=s.get('pathway_name') or (defs[0].get('pathway_name') if defs else pid)
        evidence_rows=[]; flow=[]
        def step_html(d):
            rid=d.get('reaction_id',''); rr=rev.get(rid,{}); q,qlabel=quality(rr); evidence_rows.append(rr)
            enz=rdef.get(rid,{}).get('enzyme_name','') or rr.get('enzyme_name','')
            alt=d.get('alternative_group','').strip(); optional=str(d.get('required','')).lower() not in {'yes','true','1','required'}
            badge='alternative' if alt else ('optional' if optional else 'required')
            loci=rr.get('supported_loci') or rr.get('adjudicated_loci') or ''
            rel=rr.get('evidence_relationship','')
            note=''
            if q=='candidate': note='<div class="plainnote">Candidate detected; not counted as supported.</div>'
            elif q=='review': note='<div class="plainnote">Evidence exists, but review is required before counting.</div>'
            elif q=='missing': note='<div class="plainnote">No countable evidence detected.</div>'
            return f'<div class="step {q}"><div class="rid">{esc(rid)}</div><div class="enz">{esc(enz)}</div><div class="meta"><span>{esc(qlabel)}</span><span>{esc(badge)}</span></div>' + (f'<div class="loci">{esc(loci)}</div>' if loci else '') + note + (f'<div class="rel">{esc(rel)}</div>' if rel else '') + '</div>'
        i=0
        while i < len(defs):
            d=defs[i]; alt=d.get('alternative_group','').strip()
            if alt:
                grp=[]
                while i < len(defs) and defs[i].get('alternative_group','').strip()==alt:
                    grp.append(defs[i]); i+=1
                opts='<div class="orlabel">ONE OF THESE ALTERNATIVES</div><div class="oropts">' + '<div class="orword">OR</div>'.join(step_html(x) for x in grp) + '</div>'
                flow.append('<div class="orgroup">'+opts+'</div>')
            else:
                flow.append(step_html(d)); i+=1
            if i < len(defs): flow.append('<div class="arrow">→</div>')
        if not flow: flow=['<div class="empty">No ordered pathway definition available.</div>']
        trs=[]
        for d in defs:
            rid=d.get('reaction_id',''); rr=rev.get(rid,{})
            enz=rdef.get(rid,{}).get('enzyme_name','') or rr.get('enzyme_name','')
            trs.append('<tr>'+f'<td><b>{esc(rid)}</b><div class="muted">{esc(enz)}</div></td>'+f'<td>{esc(rr.get("integrated_state","NO_CANDIDATE"))}</td>'+f'<td>{esc(rr.get("supported_loci") or "—")}</td>'+f'<td>{esc(rr.get("review_loci") or "—")}</td>'+f'<td>{esc(rr.get("candidate_only_loci") or "—")}</td>'+f'<td>{esc(rr.get("evidence_relationship") or quality(rr)[1])}</td>'+'</tr>')
        req=s.get('required_units',s.get('required_reactions','0')); sup=s.get('supported_units',s.get('supported_reactions','0')); cand=s.get('candidate_units',s.get('candidate_reactions','0'))
        missing=s.get('missing_required_units',s.get('missing_required_reactions',''))
        functional_note=''
        arsenic_text=(pid+' '+pname).lower()
        if 'arsenic' in arsenic_text and ('reduction' in arsenic_text and ('efflux' in arsenic_text or 'detox' in arsenic_text)):
            functional_note='<div class="functionnote"><b>Functional interpretation:</b> arsenate reduction and arsenic resistance/efflux are distinct capabilities. They are shown together here only because this database pathway definition is composite; efflux must not be interpreted as evidence for arsenate reduction.</div>'
        context=''
        if s.get('cluster_reaction_count') not in {'','0',None}:
            context=f'<div class="context"><b>Genomic context:</b> {esc(s.get("cluster_reaction_count"))} reaction(s) on {esc(s.get("cluster_contig"))}' + (f', span {esc(s.get("cluster_span_bp"))} bp' if s.get('cluster_span_bp') else '') + '</div>'
        cards.append(f'''<section class="path-card {'not-detected' if state=='NOT_DETECTED' else 'detected'}">
<div class="head"><div><div class="pid">{esc(pid)}</div><h2>{esc(pname)}</h2></div><div class="badges"><span>{esc(STATE_LABEL.get(state,state))}</span><span>{esc(pathway_conf(evidence_rows))}</span></div></div>{functional_note}
<div class="metrics compactmetrics"><div><b>{esc(sup)}/{esc(req)} ({esc(s.get('supported_percent','0'))}%)</b><span>required functions supported</span></div><div><b>{esc(max(0,intish(cand)-intish(sup)))}</b><span>additional candidate/review function(s)</span></div><div><b>{esc(STATE_LABEL.get(state,state))}</b><span>pathway status</span></div></div>
<div class="flowwrap"><div class="flow">{''.join(flow)}</div></div>{context}{f'<div class="missing"><b>Not supported / still required:</b> {esc(missing)}</div>' if missing else ''}
<details><summary>Reaction & locus evidence</summary><div class="tablewrap"><table><thead><tr><th>Reaction</th><th>State</th><th>Supported loci</th><th>Review loci</th><th>Candidate-only</th><th>Evidence relationship</th></tr></thead><tbody>{''.join(trs) if trs else '<tr><td colspan="6">No reaction rows.</td></tr>'}</tbody></table></div></details></section>''')

    flagged=[d for d in D if d.get('decision') in {'WINNER_KEEP','LOSER_REMOVE','AMBIGUOUS_REVIEW'}]
    review_html=''
    if flagged:
        rows=''.join(f'<tr><td>{esc(d.get("locus_tag"))}</td><td>{esc(d.get("reaction_id"))}</td><td>{esc(d.get("decision"))}</td><td>{esc(d.get("winner_reaction") or "—")}</td><td>{esc(d.get("score"))}</td><td>{esc(d.get("product"))}</td></tr>' for d in flagged)
        review_html=f'<section class="reviewcard"><h2>Cross-reaction review</h2><p class="muted">Shared proteins are not automatically counted for several unrelated reactions.</p><div class="tablewrap"><table><thead><tr><th>Locus</th><th>Reaction</th><th>Decision</th><th>Winner</th><th>Score</th><th>Product</th></tr></thead><tbody>{rows}</tbody></table></div></section>'

    sample=esc(a.sample_name or Path(a.pathway_summary).stem)
    css='''<style>:root{--bg:#f6f7f9;--card:#fff;--text:#18202a;--muted:#657180;--border:#d9dee5;--good:#e9f6ee;--goodb:#73a987;--ann:#fff7e5;--annb:#c9a45a;--review:#fff0eb;--reviewb:#c98265;--cand:#f5f0ff;--candb:#9277bc;--miss:#f2f3f5;--missb:#aab1ba}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text);font-family:Segoe UI,Arial,sans-serif;line-height:1.45}main{max-width:1240px;margin:auto;padding:28px 18px 56px}.hero,.path-card,.reviewcard{background:var(--card);border:1px solid var(--border);border-radius:18px;padding:20px;margin:14px 0}.hero h1{margin:3px 0 5px;font-size:28px}.muted{color:var(--muted);font-size:12px}.summary,.metrics{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px;margin-top:16px}.summary div,.metrics div{border:1px solid var(--border);border-radius:12px;padding:11px;background:#fafbfc}.summary b,.metrics b{display:block;font-size:21px}.summary span,.metrics span{font-size:11px;color:var(--muted)}.controls{display:flex;flex-wrap:wrap;gap:8px;margin:14px 0}.controls button{border:1px solid var(--border);background:#fff;padding:9px 12px;border-radius:10px;cursor:pointer}.controls button.active{outline:2px solid #516779}.head{display:flex;justify-content:space-between;gap:14px}.pid{font-size:11px;color:var(--muted);font-weight:700}.head h2{margin:2px 0;font-size:20px}.badges{display:flex;gap:6px;flex-wrap:wrap;justify-content:flex-end}.badges span{font-size:11px;padding:5px 8px;border:1px solid var(--border);border-radius:999px}.compactmetrics{grid-template-columns:1.25fr 1fr 1.5fr}.flowwrap,.tablewrap{overflow:auto}.flow{display:flex;align-items:stretch;gap:8px;min-width:max-content;padding:8px 0}.step{width:210px;border:2px solid var(--missb);background:var(--miss);border-radius:14px;padding:10px}.step.strong{background:var(--good);border-color:var(--goodb)}.step.supported{background:#eef5ff;border-color:#7196c5}.step.annotation{background:var(--ann);border-color:var(--annb)}.step.review{background:var(--review);border-color:var(--reviewb)}.step.candidate{background:var(--cand);border-color:var(--candb)}.step.missing{opacity:.7}.rid{font-size:11px;font-weight:800}.enz{font-size:12px;margin:5px 0;min-height:34px}.meta{display:flex;justify-content:space-between;gap:4px;font-size:10px;color:var(--muted)}.loci,.rel{font-size:9px;color:var(--muted);margin-top:5px;overflow-wrap:anywhere}.arrow{align-self:center;font-weight:800;color:var(--muted)}.orgroup{border:2px dashed #9aa4af;border-radius:16px;padding:8px;background:#fbfcfd}.orlabel{text-align:center;font-size:10px;font-weight:800;color:var(--muted);margin-bottom:6px}.oropts{display:flex;align-items:stretch;gap:8px}.orword{align-self:center;font-weight:900;color:#5e6874}.plainnote{font-size:9px;margin-top:6px;font-weight:600}.functionnote{font-size:12px;padding:9px 11px;border-radius:10px;background:#eef5ff;margin-top:9px;border-left:4px solid #7196c5}.context,.missing{font-size:12px;padding:8px 10px;border-radius:10px;background:#f7f8fa;margin-top:7px}.missing{background:#fff5f1}details{margin-top:11px}summary{cursor:pointer;font-weight:700;font-size:13px}table{border-collapse:collapse;width:100%;min-width:780px;font-size:11px}th,td{padding:8px;border-bottom:1px solid var(--border);text-align:left;vertical-align:top}th{background:#f4f6f8}.legend{display:flex;flex-wrap:wrap;gap:7px;font-size:11px;margin-top:11px}.legend span{border:1px solid var(--border);border-radius:999px;padding:4px 7px}.hidden{display:none!important}@media(max-width:720px){.summary,.metrics{grid-template-columns:repeat(2,minmax(0,1fr))}.head{display:block}.badges{justify-content:flex-start;margin-top:7px}}@media print{body{background:#fff}.controls{display:none}.path-card,.hero,.reviewcard{break-inside:avoid}details>summary{display:none}details>*{display:block!important}}</style>'''
    script='''<script>(function(){const d=document.getElementById('detected'),a=document.getElementById('all'),p=document.getElementById('print');function mode(m){document.querySelectorAll('.path-card.not-detected').forEach(x=>x.classList.toggle('hidden',m==='detected'));d.classList.toggle('active',m==='detected');a.classList.toggle('active',m==='all')}d.addEventListener('click',()=>mode('detected'));a.addEventListener('click',()=>mode('all'));p.addEventListener('click',()=>window.print());mode('detected')})();</script>'''
    doc='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+sample+' - Bioremediation pathway report</title>'+css+'</head><body><main>'
    doc+=f'<section class="hero"><div class="muted">Bioremediation Unified Runner · Visual report v{VERSION}</div><h1>{sample}</h1><div class="muted">Pathway completeness and evidence confidence are shown separately. Computational evidence requires biological/experimental validation.</div><div class="summary"><div><b>{counts["SUPPORTED_COMPLETE"]}</b><span>supported complete</span></div><div><b>{partial}</b><span>partial pathways/modules</span></div><div><b>{detected}</b><span>non-zero pathway/module calls</span></div><div><b>{review_count}</b><span>review flags / review-complete</span></div></div><div class="legend"><span>Green = strong evidence</span><span>Blue = supported</span><span>Amber = annotation only</span><span>Peach = review/conflict</span><span>Purple = candidate</span><span>Gray = missing</span></div></section>'
    doc+='<div class="controls"><button id="detected" class="active">Detected / partial only</button><button id="all">Show all pathways</button><button id="print">Print / Save PDF</button></div><div id="pathways">'+''.join(cards)+'</div>'+review_html+'</main>'+script+'</body></html>'
    Path(a.out_html).write_text(doc,encoding='utf-8')
    print('Visual Pathway Reporter v0.8.4')
    print('Pathways/modules rendered:',len(pids))
    print('Detected/partial rendered:',detected)
    print('Wrote:',a.out_html)

if __name__=='__main__': main()
