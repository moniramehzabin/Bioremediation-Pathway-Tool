#!/usr/bin/env python3
"""Universal Rescue Engine v0.2.2 — generalized evidence-aware rescue."""
from __future__ import annotations
import argparse, csv, os, subprocess, sys, time
from collections import defaultdict
from pathlib import Path
import requests
from Bio import SeqIO

VERSION = "0.2.2"
CONF_RANK = {"Weak":1,"Moderate":2,"High":3}
IPR_BASE = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
IPR_SUPPORTED="SUPPORTED"; IPR_NOT_SUPPORTED="NOT_SUPPORTED"; IPR_NOT_CURATED="NOT_CURATED"; IPR_NO_RESULT="NO_RESULT"


def read_tsv(path):
    with open(path, encoding="utf-8") as f: return list(csv.DictReader(f, delimiter="\t"))


def run(cmd):
    print("RUN:", " ".join(map(str,cmd)))
    p=subprocess.run(cmd,text=True,capture_output=True)
    if p.stdout.strip(): print(p.stdout.strip())
    if p.returncode!=0:
        if p.stderr.strip(): print(p.stderr,file=sys.stderr)
        raise SystemExit(p.returncode)


def extract_proteome(genbank,out_faa):
    proteins={}; seqs={}; n=0
    with open(out_faa,"w",encoding="utf-8") as out:
        for rec in SeqIO.parse(genbank,"genbank"):
            for feat in rec.features:
                if feat.type!="CDS": continue
                q=feat.qualifiers; aa=q.get("translation",[""])[0]
                if not aa: continue
                locus=q.get("locus_tag",[""])[0] or q.get("protein_id",[""])[0] or f"CDS_{n+1}"
                if locus in proteins:
                    raise ValueError(f"Duplicate protein locus identifier in GenBank: {locus}")
                proteins[locus]={"locus_tag":locus,"gene":q.get("gene",[""])[0],"product":q.get("product",[""])[0],"contig":rec.id,"start":int(feat.location.start)+1,"end":int(feat.location.end),"strand":"+" if feat.location.strand==1 else "-" if feat.location.strand==-1 else "."}
                seqs[locus]=aa
                out.write(f">{locus}\n")
                for i in range(0,len(aa),70): out.write(aa[i:i+70]+"\n")
                n+=1
    return proteins,seqs


def load_reference_metadata(path):
    by_ref={}
    for r in read_tsv(path):
        rid=r.get("reference_id","").strip()
        if not rid: raise ValueError("reference_metadata.tsv requires reference_id")
        if rid in by_ref: raise ValueError(f"Duplicate reference_id: {rid}")
        r.setdefault("expected_interpro_ids",""); r.setdefault("expected_domain_terms","")
        by_ref[rid]=r
    return by_ref


def parse_baseline(path):
    d=defaultdict(list)
    if path:
        for r in read_tsv(path):
            if r.get("reaction_id",""): d[r["reaction_id"]].append(r)
    return d


def parse_systems(path):
    d=defaultdict(list)
    if path:
        for r in read_tsv(path):
            if r.get("reaction_id",""): d[r["reaction_id"]].append(r)
    return d


def coverage(aln_len,full_len):
    try: return min(100.0,max(0.0,100.0*float(aln_len)/float(full_len)))
    except Exception: return 0.0


def parse_interproscan_tsv(path):
    out=defaultdict(list)
    with open(path,encoding="utf-8",errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            c=line.rstrip("\n").split("\t")
            if len(c)<6: continue
            locus=c[0].strip()
            if not locus: continue
            out[locus].append({"analysis":c[3] if len(c)>3 else "","signature_accession":c[4] if len(c)>4 else "","signature_description":c[5] if len(c)>5 else "","interpro_id":c[11] if len(c)>11 and c[11]!="-" else "","interpro_description":c[12] if len(c)>12 and c[12]!="-" else "","go_terms":c[13] if len(c)>13 and c[13]!="-" else "","pathways":c[14] if len(c)>14 and c[14]!="-" else ""})
    return out


def _split_expected(v): return [x.strip() for x in str(v or "").split("|") if x.strip()]


def interpro_evidence_status(locus,meta,ipr):
    ids=set(_split_expected(meta.get("expected_interpro_ids","")))
    terms=[x.lower() for x in _split_expected(meta.get("expected_domain_terms",""))]
    if not ids and not terms: return IPR_NOT_CURATED,"No curated InterPro expectations for this reference",False
    rows=ipr.get(locus,[])
    if not rows: return IPR_NO_RESULT,"No parsed InterPro result for this locus",False
    for r in rows:
        sig=str(r.get("signature_accession","")).strip(); ipid=str(r.get("interpro_id","")).strip()
        if ids and (sig in ids or ipid in ids):
            m=sig if sig in ids else ipid
            return IPR_SUPPORTED,f"Expected accession match: {m}",True
        blob=" ".join([r.get("signature_description",""),r.get("interpro_description",""),r.get("analysis","")]).lower()
        for t in terms:
            if t and t in blob: return IPR_SUPPORTED,f"Expected term match: {t}",True
    return IPR_NOT_SUPPORTED,"InterPro result present but no curated expectation matched",False


def local_interproscan(candidate_faa,executable,out_prefix,applications=None):
    out_tsv=out_prefix+"_interpro.tsv"; cmd=[executable,"-i",candidate_faa,"-f","TSV","-o",out_tsv]
    if applications: cmd += ["-appl",applications]
    run(cmd); return out_tsv


def _request_with_retry(method,url,*,retries=5,base_wait=5,**kwargs):
    last=None
    for attempt in range(1,retries+1):
        try:
            r=requests.request(method,url,**kwargs)
            if r.status_code in {429,500,502,503,504}: raise requests.HTTPError(f"{r.status_code} Server/Rate-limit response for {url}",response=r)
            r.raise_for_status(); return r
        except requests.RequestException as exc:
            last=exc
            if attempt>=retries: break
            wait=min(base_wait*attempt,30); print(f"  transient InterPro web error; retry {attempt}/{retries} in {wait}s: {exc}"); time.sleep(wait)
    raise last


def ebi_submit(seq_fasta,email,title,applications=None):
    p={"email":email,"title":title,"stype":"p","sequence":seq_fasta,"goterms":"false","pathways":"false"}
    if applications: p["appl"]=applications
    return _request_with_retry("POST",IPR_BASE+"/run/",data=p,timeout=60,retries=5,base_wait=5).text.strip()


def ebi_wait(job_id,poll_seconds=5,timeout_seconds=3600):
    start=time.time()
    while True:
        status=_request_with_retry("GET",IPR_BASE+f"/status/{job_id}",timeout=60,retries=5,base_wait=5).text.strip()
        if status=="FINISHED": return
        if status in {"ERROR","FAILURE","NOT_FOUND"}: raise RuntimeError(f"InterProScan web job {job_id} ended with status {status}")
        if time.time()-start>timeout_seconds: raise TimeoutError(f"Timed out waiting for InterProScan job {job_id}")
        time.sleep(poll_seconds)


def ebi_fetch_tsv(job_id):
    return _request_with_retry("GET",IPR_BASE+f"/result/{job_id}/tsv",timeout=120,retries=5,base_wait=5).text


def web_interproscan(candidate_sequences,email,out_prefix,applications=None,max_jobs=200):
    combined=out_prefix+"_interpro_web.tsv"; failed_path=out_prefix+"_interpro_failed.tsv"
    already=set()
    if os.path.exists(combined) and os.path.getsize(combined)>0:
        try: already=set(parse_interproscan_tsv(combined).keys())
        except Exception: already=set()
    pending=[(l,s) for l,s in candidate_sequences.items() if l not in already]
    if len(pending)>max_jobs: raise RuntimeError(f"{len(pending)} new InterPro candidates exceed --interpro-max-jobs={max_jobs}.")
    if already: print(f"InterPro resume: reusing {len(already)} locus/loci from {combined}")
    print(f"New InterPro web jobs required: {len(pending)}")
    failures=[]; mode="a" if os.path.exists(combined) else "w"
    with open(combined,mode,encoding="utf-8") as out:
        for i,(locus,seq) in enumerate(pending,1):
            print(f"InterPro web {i}/{len(pending)}: {locus}"); job=""
            try:
                job=ebi_submit(f">{locus}\n{seq}\n",email,f"UniversalRescue_{locus}",applications); ebi_wait(job); tsv=ebi_fetch_tsv(job); out.write(tsv)
                if tsv and not tsv.endswith("\n"): out.write("\n")
                out.flush()
            except Exception as exc:
                print(f"WARNING: InterPro unavailable for {locus}; retaining DIAMOND candidate. Reason: {exc}",file=sys.stderr)
                failures.append({"locus_tag":locus,"job_id":job,"status":"INTERPRO_UNAVAILABLE_NEEDS_RETRY","error":str(exc)})
    with open(failed_path,"w",newline="",encoding="utf-8") as f:
        fields=["locus_tag","job_id","status","error"]; w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(failures)
    print(f"InterPro completed with {len(failures)} unavailable candidate(s). See: {failed_path}" if failures else "InterPro completed for all submitted candidates.")
    return combined


def choose_interpro_mode(args):
    if args.interpro_mode!="auto": return args.interpro_mode
    if args.interpro_tsv: return "precomputed"
    if args.interproscan: return "local"
    if args.email: return "web"
    return "off"


def context_support(locus_info,reaction_id,systems,baseline,window_bp):
    my=systems.get(reaction_id,[])
    if not my: return False,""
    contig=locus_info.get("contig",""); mid=(int(locus_info.get("start",0))+int(locus_info.get("end",0)))/2
    for sr in my:
        pid=sr.get("pathway_id",""); same=[rid for rid,rows in systems.items() if any(x.get("pathway_id","")==pid for x in rows)]
        same_contig=set(); nearest=None
        for rid in same:
            if rid==reaction_id: continue
            for ev in baseline.get(rid,[]):
                if ev.get("contig","")!=contig: continue
                same_contig.add(rid)
                try:
                    emid=(int(ev.get("start",0))+int(ev.get("end",0)))/2; dist=abs(mid-emid)
                    if nearest is None or dist<nearest[0]: nearest=(dist,rid)
                except Exception: pass
        if len(same_contig)>=2: return True,f"same contig as >=2 reactions in {pid}"
        if nearest and nearest[0]<=window_bp: return True,f"within {int(nearest[0])} bp of {nearest[1]} in {pid}"
    return False,""


def diamond_confidence(pid,qcov,scov,evalue):
    if pid>=40 and qcov>=70 and scov>=70 and evalue<=1e-20: return "High"
    if pid>=25 and qcov>=50 and scov>=50 and evalue<=1e-5: return "Moderate"
    return "Weak"


def classify(pid,qcov,scov,evalue,ip_supported,cx_ok):
    d=diamond_confidence(pid,qcov,scov,evalue)
    if d=="High": return "High"
    if d=="Moderate": return "High" if ip_supported and cx_ok else "Moderate"
    return "Weak"


def main():
    ap=argparse.ArgumentParser(description="Universal Rescue Engine v0.2.2")
    ap.add_argument("genbank"); ap.add_argument("--references",required=True); ap.add_argument("--reference-metadata",required=True); ap.add_argument("--reactions",required=True); ap.add_argument("--pathways"); ap.add_argument("--baseline-evidence"); ap.add_argument("--diamond",default=r".\diamond.exe"); ap.add_argument("--out-prefix",default="rescue")
    ap.add_argument("--interpro-mode",choices=["auto","precomputed","local","web","off"],default="auto"); ap.add_argument("--interpro-tsv"); ap.add_argument("--interproscan"); ap.add_argument("--interpro-applications",default=None); ap.add_argument("--email"); ap.add_argument("--interpro-max-jobs",type=int,default=200); ap.add_argument("--interpro-selection",choices=["moderate","all"],default="moderate")
    ap.add_argument("--evalue",type=float,default=1e-5); ap.add_argument("--min-identity",type=float,default=25.0); ap.add_argument("--min-qcov",type=float,default=50.0); ap.add_argument("--min-scov",type=float,default=50.0); ap.add_argument("--context-window",type=int,default=20000)
    args=ap.parse_args()
    reactions={r["reaction_id"]:r for r in read_tsv(args.reactions)}; ref_meta=load_reference_metadata(args.reference_metadata); baseline=parse_baseline(args.baseline_evidence); systems=parse_systems(args.pathways)
    proteome_faa=args.out_prefix+"_proteome.faa"; proteins,seqs=extract_proteome(args.genbank,proteome_faa); print(f"Extracted translated CDS proteins: {len(proteins)}")
    db=args.out_prefix+"_refdb"; raw=args.out_prefix+"_diamond_raw.tsv"; run([args.diamond,"makedb","--in",args.references,"-d",db])
    fields=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]
    run([args.diamond,"blastp","-d",db,"-q",proteome_faa,"-o",raw,"--evalue",str(args.evalue),"--max-target-seqs","50","--outfmt","6",*fields])
    candidates=[]; all_loci=set(); moderate=set(); high=set()
    with open(raw,encoding="utf-8") as f:
        for line in f:
            if not line.strip(): continue
            h=dict(zip(fields,line.rstrip("\n").split("\t"))); pid=float(h["pident"]); qcov=coverage(h["length"],h["qlen"]); scov=coverage(h["length"],h["slen"]); ev=float(h["evalue"])
            if pid<args.min_identity or qcov<args.min_qcov or scov<args.min_scov: continue
            ref_id=h["sseqid"].split()[0]; meta=ref_meta.get(ref_id)
            if not meta: continue
            rid=meta.get("reaction_id","")
            if rid not in reactions: continue
            dc=diamond_confidence(pid,qcov,scov,ev)
            if dc=="Weak": continue
            row=dict(h); row.update({"reference_id":ref_id,"reaction_id":rid,"identity_pct":pid,"query_coverage_pct":qcov,"reference_coverage_pct":scov,"diamond_confidence":dc}); candidates.append(row); locus=h["qseqid"]; all_loci.add(locus); (high if dc=="High" else moderate).add(locus)
    selected=set(all_loci) if args.interpro_selection=="all" else set(moderate); high_only=high-moderate if args.interpro_selection=="moderate" else set()
    print(f"DIAMOND candidates passing rescue thresholds: {len(candidates)}"); print(f"Unique threshold-passing proteins: {len(all_loci)}"); print(f"Unique proteins with High DIAMOND assignment(s): {len(high)}"); print(f"Unique proteins with Moderate DIAMOND assignment(s): {len(moderate)}"); print(f"High-only proteins skipped from NEW InterPro jobs: {len(high_only)}"); print(f"Unique proteins selected for NEW InterPro jobs: {len(selected)}")
    candidate_faa=args.out_prefix+"_interpro_candidates.faa"
    with open(candidate_faa,"w",encoding="utf-8") as out:
        for locus in sorted(selected):
            out.write(f">{locus}\n"); seq=seqs[locus]
            for i in range(0,len(seq),70): out.write(seq[i:i+70]+"\n")
    mode=choose_interpro_mode(args); print(f"InterPro mode: {mode}"); ipr=defaultdict(list); ipfile=""
    if mode=="precomputed":
        if not args.interpro_tsv: raise SystemExit("--interpro-mode precomputed requires --interpro-tsv")
        ipfile=args.interpro_tsv
    elif mode=="local":
        if not args.interproscan: raise SystemExit("--interpro-mode local requires --interproscan")
        ipfile=local_interproscan(candidate_faa,args.interproscan,args.out_prefix,args.interpro_applications)
    elif mode=="web":
        if not args.email: raise SystemExit("--interpro-mode web requires --email")
        ipfile=web_interproscan({loc:seqs[loc] for loc in sorted(selected)},args.email,args.out_prefix,args.interpro_applications,args.interpro_max_jobs)
    elif mode=="off": print("WARNING: InterPro evidence disabled for this run.")
    if ipfile: ipr=parse_interproscan_tsv(ipfile); print(f"Proteins with parsed InterPro results: {len(ipr)}")
    rows=[]
    for h in candidates:
        locus=h["qseqid"]; ref_id=h["reference_id"]; rid=h["reaction_id"]; meta=ref_meta[ref_id]; pinfo=proteins.get(locus,{})
        ips,ipd,ipok=interpro_evidence_status(locus,meta,ipr); cx,cxd=context_support(pinfo,rid,systems,baseline,args.context_window); pid=float(h["identity_pct"]); qcov=float(h["query_coverage_pct"]); scov=float(h["reference_coverage_pct"]); ev=float(h["evalue"]); conf=classify(pid,qcov,scov,ev,ipok,cx)
        if mode=="precomputed": sub="PRECOMPUTED" if locus in ipr else "NOT_SELECTED"
        elif locus in selected: sub="SUBMITTED_OR_AVAILABLE"
        elif locus in high_only: sub="SKIPPED_HIGH_DIAMOND"
        else: sub="NOT_SELECTED"
        rows.append({"reaction_id":rid,"gene_family":reactions[rid].get("gene_family",""),"enzyme_name":reactions[rid].get("enzyme_name",""),"locus_tag":locus,"gene":pinfo.get("gene",""),"product":pinfo.get("product",""),"contig":pinfo.get("contig",""),"reference_id":ref_id,"reference_name":meta.get("reference_name",""),"reference_accession":meta.get("accession",""),"identity_pct":round(pid,2),"query_coverage_pct":round(qcov,2),"reference_coverage_pct":round(scov,2),"evalue":h["evalue"],"bitscore":h["bitscore"],"diamond_confidence":h["diamond_confidence"],"interpro_support":ipok,"interpro_status":ips,"interpro_submission":sub,"interpro_detail":ipd,"context_support":cx,"context_detail":cxd,"confidence":conf,"evidence_source":"DIAMOND"+("+InterPro" if ipok else "")+("+Context" if cx else "")})
    best={}
    for r in rows:
        k=(r["reaction_id"],r["locus_tag"]); old=best.get(k); newkey=(CONF_RANK.get(r["confidence"],0),1 if r["interpro_status"]==IPR_SUPPORTED else 0,float(r["bitscore"])); oldkey=(-1,-1,-1) if old is None else (CONF_RANK.get(old["confidence"],0),1 if old["interpro_status"]==IPR_SUPPORTED else 0,float(old["bitscore"]));
        if old is None or newkey>oldkey: best[k]=r
    rows=sorted(best.values(),key=lambda r:(r["reaction_id"],-CONF_RANK.get(r["confidence"],0),r["locus_tag"]))
    out_evidence=args.out_prefix+"_rescued_evidence.tsv"; cols=["reaction_id","gene_family","enzyme_name","locus_tag","gene","product","contig","reference_id","reference_name","reference_accession","identity_pct","query_coverage_pct","reference_coverage_pct","evalue","bitscore","diamond_confidence","interpro_support","interpro_status","interpro_submission","interpro_detail","context_support","context_detail","confidence","evidence_source"]
    with open(out_evidence,"w",newline="",encoding="utf-8") as f: w=csv.DictWriter(f,fieldnames=cols,delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary=defaultdict(list)
    for r in rows: summary[r["reaction_id"]].append(r)
    out_summary=args.out_prefix+"_rescue_summary.tsv"; sfields=["reaction_id","gene_family","enzyme_name","rescued","rescued_loci","best_confidence","best_interpro_status","loci"]
    with open(out_summary,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=sfields,delimiter="\t"); w.writeheader()
        for rid in sorted(reactions):
            rr=summary.get(rid,[]); bc=max((x["confidence"] for x in rr),key=lambda c:CONF_RANK.get(c,0)) if rr else ""; st={x["interpro_status"] for x in rr}; bi=IPR_SUPPORTED if IPR_SUPPORTED in st else IPR_NOT_SUPPORTED if IPR_NOT_SUPPORTED in st else IPR_NO_RESULT if IPR_NO_RESULT in st else IPR_NOT_CURATED if IPR_NOT_CURATED in st else ""
            w.writerow({"reaction_id":rid,"gene_family":reactions[rid].get("gene_family",""),"enzyme_name":reactions[rid].get("enzyme_name",""),"rescued":bool(rr),"rescued_loci":len({x["locus_tag"] for x in rr}),"best_confidence":bc,"best_interpro_status":bi,"loci":"|".join(sorted({x["locus_tag"] for x in rr}))})
    print(); print(f"Universal Rescue Engine v{VERSION}"); print(f"Rescued reaction-locus assignments: {len(rows)}"); print(f"Reactions rescued: {len(summary)}"); print(f"Wrote: {out_evidence}"); print(f"Wrote: {out_summary}");
    if ipfile: print(f"InterPro results used: {ipfile}")
    print("Computational assignments require biological/experimental validation.")

if __name__=="__main__": main()
