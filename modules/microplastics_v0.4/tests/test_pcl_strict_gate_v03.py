#!/usr/bin/env python3
from pathlib import Path
import csv,tempfile,importlib.util
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('runner',ROOT/'runner'/'run_microplastics_branch.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
RF=['reaction_id','enzyme_name','integrated_state','candidate_present','exact_supported','review_flag','annotation_supported','annotation_loci','rescue_signal','rescue_best_sequence_evidence','rescue_candidate_loci','rescue_supported_loci','all_evidence_loci','evidence_sources','pre_adjudication_state','annotation_supported_loci','same_locus_multisource_loci','supported_loci','review_loci','candidate_only_loci','adjudicated_loci','removed_competing_loci','ambiguous_loci','evidence_relationship']
PF=['pathway_id','pathway_name','required_units','candidate_units','candidate_percent','supported_units','supported_percent','review_units','adjudicated_pathway_state','cluster_reaction_count','cluster_contig','cluster_span_bp','cluster_reactions','missing_required_units','alternative_group_rule','context_rule']
CF=['reaction_id','locus_tag','identity_pct','query_coverage_pct','reference_coverage_pct','diamond_confidence','policy_interpro_status','count_toward_pathway']
def write(path,fields,rows):
    with path.open('w',encoding='utf-8',newline='') as f: w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def run_case(pid,q,sc,annot=False):
    with tempfile.TemporaryDirectory() as td:
        td=Path(td); r=td/'r.tsv';p=td/'p.tsv';c=td/'c.tsv';a=td/'a.tsv'
        rr={k:'' for k in RF}; rr.update(reaction_id='BRXN_PCL_DEP',enzyme_name='PCL depolymerase',integrated_state='SUPPORTED_RESCUE',candidate_present='True',exact_supported='True',review_flag='NO',annotation_supported='True' if annot else 'False',rescue_signal='SUPPORTED',supported_loci='L1',adjudicated_loci='L1')
        pr={k:'' for k in PF};pr.update(pathway_id='BMOD_PCL_DEPOLYMERIZATION',pathway_name='PCL',required_units='1',candidate_units='1',candidate_percent='100',supported_units='1',supported_percent='100',review_units='0',adjudicated_pathway_state='SUPPORTED')
        cr=dict(reaction_id='BRXN_PCL_DEP',locus_tag='L1',identity_pct=str(pid),query_coverage_pct=str(q),reference_coverage_pct=str(sc),diamond_confidence='High',policy_interpro_status='COMPATIBLE',count_toward_pathway='True')
        write(r,RF,[rr]);write(p,PF,[pr]);write(c,CF,[cr]);m.apply_pcl_strict_gate(r,p,c,a)
        out=m.read_tsv(r)[0]; po=m.read_tsv(p)[0]
        return out,po
x,p=run_case(55,90,90,False); assert x['review_flag']=='YES' and p['supported_units']=='0' and p['review_units']=='1'
x,p=run_case(65,90,90,False); assert x['review_flag']=='NO' and p['supported_units']=='1'
x,p=run_case(45,50,50,True); assert x['review_flag']=='NO' and p['supported_units']=='1'
print('PASS PCL strict gate: below-threshold rescue -> REVIEW; strong rescue -> keep; exact annotation -> keep')
