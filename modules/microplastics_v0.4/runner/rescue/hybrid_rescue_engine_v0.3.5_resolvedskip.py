#!/usr/bin/env python3
"""
Hybrid Rescue Engine v0.3.2
Evidence-aware, context-aware triage for bioremediation reaction rescue.

Design:
- High DIAMOND: usually accepted without InterPro.
- High DIAMOND in ambiguity-prone families: InterPro is still requested when
  curated conflicting-family rules exist.
- Moderate DIAMOND: InterPro is requested only when the reaction family has a
  usable curated policy. If InterPro is not useful/curated, the hit stays
  explicit as sequence-only evidence.
- Weak DIAMOND: disabled by default. If the user deliberately lowers screening
  thresholds AND enables --enable-weak-rescue, weak hits are sent to InterPro
  only when the reaction policy can specifically adjudicate them.

Reaction adjudication:
  SPECIFIC / COMPATIBLE / CONFLICTING / AMBIGUOUS /
  NOT_CURATED / NO_RESULT / NO_MATCH

Pathway-context rule:
- Context can strengthen a Moderate compatible/uncurated hit.
- Context NEVER rescues a CONFLICTING or AMBIGUOUS exact-reaction assignment.
- Missing evidence is never inferred merely because neighboring steps exist.

Outputs:
- *_all_candidates.tsv
- *_adjudicated_rescued_evidence.tsv  (exact-locus-supported evidence only)
- *_reaction_adjudication_summary.tsv
- normal DIAMOND / candidate / InterPro intermediate files

This is a research prototype. Computational assignments require validation.
"""

from __future__ import annotations
import argparse, csv, os, sys
from collections import defaultdict
from pathlib import Path

import rescue_core_v022 as core
import evidence_policy_v054 as ep

VERSION = "0.3.4"
CONF_RANK = {"Weak": 1, "Moderate": 2, "High": 3}

COUNTABLE = "COUNTABLE"
REVIEW = "REVIEW"
REJECT = "REJECT"


def load_policy(path):
    return ep.load_policy(path)


def policy_has_specific(p):
    if not p:
        return False
    return bool(ep.splitv(p.get("specific_ids")) or ep.splitv(p.get("specific_terms")))


def policy_has_conflicts(p):
    if not p:
        return False
    return bool(ep.splitv(p.get("conflicting_ids")) or ep.splitv(p.get("conflicting_terms")))


def policy_is_usable(p):
    return bool(p) and p.get("policy_status", "").strip() != "NOT_CURATED"


def needs_interpro(diamond_conf, reaction_id, policy, enable_weak=False, single_reaction_modules=None):
    """
    Decide whether NEW InterPro work is scientifically useful.

    High:
      only when the reaction has curated conflicting-family rules.
      This catches cases like LadA-like / P-type ATPase overcalls.

    Moderate:
      when there is a usable policy. InterPro can then specifically support,
      broadly support, or contradict the hit.

    Weak:
      only if weak rescue is explicitly enabled AND the policy has specific
      reaction-family rules. Generic-domain confirmation is not enough.
    """
    p = policy.get(reaction_id)
    single_reaction_modules = single_reaction_modules or set()
    # InterPro is a selective referee, not a gate for Moderate evidence.
    # High/Moderate hits are checked only when the family is ambiguity-prone
    # (curated conflicts) or when a Moderate hit represents a stand-alone 1/1
    # functional module with a specific policy (e.g. ChrR-like Cr(VI) reduction).
    if diamond_conf == "High":
        return policy_has_conflicts(p)
    if diamond_conf == "Moderate":
        return policy_has_conflicts(p) or (reaction_id in single_reaction_modules and policy_has_specific(p))
    if diamond_conf == "Weak":
        return bool(enable_weak and policy_has_specific(p))
    return False


def baseline_strong_reactions(baseline):
    strong = set()
    for rid, rows in baseline.items():
        for r in rows:
            if str(r.get("confidence", "")).strip().lower() == "high":
                strong.add(rid)
                break
    return strong


def preliminary_strong_reactions(rows):
    """
    Strong reaction evidence for pathway-context calculation.

    Only High DIAMOND assignments without CONFLICTING/AMBIGUOUS policy evidence
    contribute. This prevents circular reinforcement from uncertain Moderate hits.
    """
    strong = set()
    for r in rows:
        if r["diamond_confidence"] != "High":
            continue
        if r["policy_interpro_status"] in {"CONFLICTING", "AMBIGUOUS"}:
            continue
        strong.add(r["reaction_id"])
    return strong


def pathway_context(reaction_id, pathways_by_reaction, strong_reactions, min_other_steps=2, min_fraction=0.50):
    """
    Evaluate pathway-level support from OTHER reactions only.

    Returns (bool, detail, best_fraction, best_other_count).
    Single-step modules do not create self-support.
    """
    memberships = pathways_by_reaction.get(reaction_id, [])
    best_fraction = 0.0
    best_other = 0
    best_detail = ""

    for membership in memberships:
        pid = membership.get("pathway_id", "")
        if not pid:
            continue
        members = set()
        for rid, rows in pathways_by_reaction.items():
            for rr in rows:
                if rr.get("pathway_id", "") == pid:
                    members.add(rid)
                    break
        others = members - {reaction_id}
        if not others:
            continue
        present = others & strong_reactions
        frac = len(present) / len(others)
        if (len(present), frac) > (best_other, best_fraction):
            best_other = len(present)
            best_fraction = frac
            best_detail = f"{pid}: {len(present)}/{len(others)} other reaction(s) strongly supported"

    ok = best_other >= min_other_steps and best_fraction >= min_fraction
    return ok, best_detail, best_fraction, best_other


def adjudicate(diamond_conf, policy_status, pathway_ok=False, genomic_context=False):
    """
    Evidence-preserving treatment for v0.3.1.

    High and Moderate DIAMOND assignments remain detected candidates.
    Moderate assignments are always visibly REVIEW-flagged, even when InterPro
    supports them. InterPro contradiction is retained as a warning rather than
    silently deleting the candidate.

    Weak assignments are not automatically countable. They require specific
    InterPro support plus either pathway or genomic context. Pathway context is
    optional evidence and is never required for High/Moderate or stand-alone
    modules.
    """
    if diamond_conf == "High":
        if policy_status == "CONFLICTING":
            return REVIEW, "High", False
        if policy_status == "AMBIGUOUS":
            return REVIEW, "High", False
        return COUNTABLE, "High", True

    if diamond_conf == "Moderate":
        # Moderate sequence evidence is always retained as a reaction candidate.
        # Exact-locus support is withheld only when protein-family evidence
        # explicitly conflicts or is ambiguous. The reaction signal itself is
        # never erased; it remains visible in candidate-level pathway evidence.
        if policy_status in {"CONFLICTING", "AMBIGUOUS"}:
            return REVIEW, "Moderate", False
        return REVIEW, "Moderate", True

    if diamond_conf == "Weak":
        if policy_status == "SPECIFIC" and (pathway_ok or genomic_context):
            return REVIEW, "Weak", True
        return REVIEW, "Weak", False

    return REVIEW, diamond_conf or "Weak", False


def main():
    ap = argparse.ArgumentParser(description="Hybrid Rescue Engine v0.3.3 component-aware")
    ap.add_argument("genbank")
    ap.add_argument("--references", required=True)
    ap.add_argument("--reference-metadata", required=True)
    ap.add_argument("--reactions", required=True)
    ap.add_argument("--pathways", required=True)
    ap.add_argument("--reaction-components", help="reaction_components.tsv; defaults to sibling of --reactions")
    ap.add_argument("--baseline-evidence")
    ap.add_argument("--resolved-reaction-loci",
                    help="TSV of reaction_id+locus_tag pairs already resolved by exact annotation; same pair is skipped during rescue")
    ap.add_argument("--policy", default=str(Path(__file__).with_name("reaction_interpro_policy_v054.tsv")))
    ap.add_argument("--diamond", default=r".\diamond.exe")
    ap.add_argument("--out-prefix", default="hybrid_rescue")

    ap.add_argument("--interpro-mode", choices=["auto","precomputed","local","web","off"], default="auto")
    ap.add_argument("--interpro-tsv")
    ap.add_argument("--interproscan")
    ap.add_argument("--interpro-applications", default=None)
    ap.add_argument("--email")
    ap.add_argument("--interpro-max-jobs", type=int, default=200)

    # Standard threshold-passing band.
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--min-identity", type=float, default=25.0)
    ap.add_argument("--min-qcov", type=float, default=50.0)
    ap.add_argument("--min-scov", type=float, default=50.0)

    # Weak rescue is opt-in. User must ALSO deliberately lower one or more
    # screening thresholds if they want hits below the normal Moderate band.
    ap.add_argument("--enable-weak-rescue", action="store_true")
    ap.add_argument("--context-window", type=int, default=20000)
    ap.add_argument("--pathway-context-min-other", type=int, default=2)
    ap.add_argument("--pathway-context-min-fraction", type=float, default=0.50)
    args = ap.parse_args()

    reactions = {r["reaction_id"]: r for r in core.read_tsv(args.reactions) if str(r.get("active","yes")).strip().lower() not in {"no","false","0","inactive"}}
    ref_meta = core.load_reference_metadata(args.reference_metadata)
    # Ignore references that still point to inactive legacy reactions.
    ref_meta = {k:v for k,v in ref_meta.items() if v.get("reaction_id","") in reactions}
    baseline = core.parse_baseline(args.baseline_evidence)
    resolved_pairs=set()
    if args.resolved_reaction_loci:
        for rr in core.read_tsv(args.resolved_reaction_loci):
            rid=str(rr.get("reaction_id","")).strip()
            tag=str(rr.get("locus_tag","")).strip()
            if rid and tag:
                resolved_pairs.add((rid,tag))
    skipped_resolved_pairs=0
    systems = core.parse_systems(args.pathways)
    policy = load_policy(args.policy)
    components_path = args.reaction_components or str(Path(args.reactions).with_name("reaction_components.tsv"))
    component_rows = core.read_tsv(components_path) if Path(components_path).exists() else []
    component_requirements = defaultdict(list)
    for cr in component_rows:
        if str(cr.get("required","")).strip().lower() in {"yes","true","1","required"}:
            component_requirements[cr.get("reaction_id","")].append(cr.get("component_id",""))
    component_requirements = {rid:[c for c in comps if c] for rid,comps in component_requirements.items() if comps}

    # Reactions that occur in at least one one-reaction module. These cannot
    # rely on multi-step pathway context, so Moderate hits can receive targeted
    # InterPro adjudication when a specific policy exists.
    pathway_members = defaultdict(set)
    for rid, memberships in systems.items():
        for m in memberships:
            pid = m.get("pathway_id", "")
            if pid:
                pathway_members[pid].add(rid)
    single_reaction_modules = {
        rid for rid, memberships in systems.items()
        if any(len(pathway_members.get(m.get("pathway_id", ""), set())) == 1 for m in memberships)
    }

    # Extract proteome.
    proteome_faa = args.out_prefix + "_proteome.faa"
    proteins, seqs = core.extract_proteome(args.genbank, proteome_faa)
    print(f"Extracted translated CDS proteins: {len(proteins)}")

    # DIAMOND.
    db = args.out_prefix + "_refdb"
    raw = args.out_prefix + "_diamond_raw.tsv"
    core.run([args.diamond, "makedb", "--in", args.references, "-d", db])
    fields = ["qseqid","sseqid","pident","length","mismatch","gapopen",
              "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]
    core.run([args.diamond, "blastp", "-d", db, "-q", proteome_faa, "-o", raw,
              "--evalue", str(args.evalue), "--max-target-seqs", "50",
              "--outfmt", "6", *fields])

    candidates = []
    selected_loci = set()
    all_loci = set()
    high_n = moderate_n = weak_n = 0

    with open(raw, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            h = dict(zip(fields, line.rstrip("\n").split("\t")))
            pid = float(h["pident"])
            qcov = core.coverage(h["length"], h["qlen"])
            scov = core.coverage(h["length"], h["slen"])
            evalue = float(h["evalue"])

            # CLI screening gate.
            if pid < args.min_identity or qcov < args.min_qcov or scov < args.min_scov:
                continue

            ref_id = h["sseqid"].split()[0]
            meta = ref_meta.get(ref_id)
            if not meta:
                continue
            rid = meta.get("reaction_id", "")
            if rid not in reactions:
                continue

            # Pair-specific skip: exact annotation already resolved THIS reaction at THIS locus.
            # The same locus remains eligible for rescue against DIFFERENT reactions.
            if (rid, h["qseqid"]) in resolved_pairs:
                skipped_resolved_pairs += 1
                continue

            dconf = core.diamond_confidence(pid, qcov, scov, evalue)
            if dconf == "Weak" and not args.enable_weak_rescue:
                continue

            h.update({
                "reference_id": ref_id,
                "reaction_id": rid,
                "identity_pct": pid,
                "query_coverage_pct": qcov,
                "reference_coverage_pct": scov,
                "diamond_confidence": dconf,
            })
            candidates.append(h)
            all_loci.add(h["qseqid"])
            if dconf == "High":
                high_n += 1
            elif dconf == "Moderate":
                moderate_n += 1
            else:
                weak_n += 1

            if needs_interpro(dconf, rid, policy, args.enable_weak_rescue, single_reaction_modules):
                selected_loci.add(h["qseqid"])

    print(f"DIAMOND candidate assignments retained: {len(candidates)}")
    print(f"Unique retained proteins: {len(all_loci)}")
    print(f"Resolved same-reaction annotation hits skipped before rescue: {skipped_resolved_pairs}")
    print(f"Assignment bands: High={high_n}, Moderate={moderate_n}, Weak={weak_n}")
    print(f"Unique proteins selected for InterPro by hybrid policy: {len(selected_loci)}")

    # InterPro candidate FASTA.
    candidate_faa = args.out_prefix + "_interpro_candidates.faa"
    with open(candidate_faa, "w", encoding="utf-8") as out:
        for locus in sorted(selected_loci):
            seq = seqs[locus]
            out.write(f">{locus}\n")
            for i in range(0, len(seq), 70):
                out.write(seq[i:i+70] + "\n")

    mode = core.choose_interpro_mode(args)
    print(f"InterPro mode: {mode}")
    ipr = defaultdict(list)
    ipfile = ""

    if mode == "precomputed":
        if not args.interpro_tsv:
            raise SystemExit("--interpro-mode precomputed requires --interpro-tsv")
        ipfile = args.interpro_tsv
    elif mode == "local":
        if not args.interproscan:
            raise SystemExit("--interpro-mode local requires --interproscan")
        ipfile = core.local_interproscan(candidate_faa, args.interproscan, args.out_prefix, args.interpro_applications)
    elif mode == "web":
        if not args.email:
            raise SystemExit("--interpro-mode web requires --email")
        ipfile = core.web_interproscan(
            {loc: seqs[loc] for loc in sorted(selected_loci)},
            args.email,
            args.out_prefix,
            args.interpro_applications,
            args.interpro_max_jobs,
        )
    elif mode == "off":
        print("WARNING: InterPro disabled; policy states may remain NO_RESULT/NOT_CURATED.")

    if ipfile:
        ipr = core.parse_interproscan_tsv(ipfile)
        print(f"Proteins with parsed InterPro results: {len(ipr)}")

    # First pass: best reference assignment per reaction+locus, with policy evidence.
    assessed = []
    for h in candidates:
        locus = h["qseqid"]
        rid = h["reaction_id"]
        ref_id = h["reference_id"]
        pinfo = proteins.get(locus, {})

        pstatus, pdetail = ep.assess(rid, locus, policy, ipr)
        gctx, gdetail = core.context_support(pinfo, rid, systems, baseline, args.context_window)

        assessed.append({
            "reaction_id": rid,
            "gene_family": reactions[rid].get("gene_family", ""),
            "enzyme_name": reactions[rid].get("enzyme_name", ""),
            "locus_tag": locus,
            "gene": pinfo.get("gene", ""),
            "product": pinfo.get("product", ""),
            "contig": pinfo.get("contig", ""),
            "reference_id": ref_id,
            "reference_name": ref_meta[ref_id].get("reference_name", ""),
            "reference_accession": ref_meta[ref_id].get("accession", ""),
            "component_id": ref_meta[ref_id].get("component_id", ""),
            "component_gene_family": ref_meta[ref_id].get("component_gene_family", ""),
            "identity_pct": round(float(h["identity_pct"]), 2),
            "query_coverage_pct": round(float(h["query_coverage_pct"]), 2),
            "reference_coverage_pct": round(float(h["reference_coverage_pct"]), 2),
            "evalue": h["evalue"],
            "bitscore": h["bitscore"],
            "diamond_confidence": h["diamond_confidence"],
            "policy_interpro_status": pstatus,
            "policy_interpro_detail": pdetail,
            "genomic_context_support": gctx,
            "genomic_context_detail": gdetail,
            "_bitscore_float": float(h["bitscore"]),
        })

    best = {}
    policy_rank = {"SPECIFIC": 5, "COMPATIBLE": 4, "NOT_CURATED": 3, "NO_RESULT": 2,
                   "NO_MATCH": 1, "AMBIGUOUS": 0, "CONFLICTING": -1}
    for r in assessed:
        k = (r["reaction_id"], r.get("component_id",""), r["locus_tag"])
        key = (
            CONF_RANK.get(r["diamond_confidence"], 0),
            policy_rank.get(r["policy_interpro_status"], 0),
            r["_bitscore_float"],
        )
        if k not in best or key > best[k][0]:
            best[k] = (key, r)
    assessed = [x[1] for x in best.values()]

    # Pathway context uses only strong, non-conflicting evidence + direct High baseline.
    strong_reactions = preliminary_strong_reactions(assessed) | baseline_strong_reactions(baseline)

    final_rows = []
    countable_rows = []
    for r in assessed:
        pok, pdetail, pfrac, pother = pathway_context(
            r["reaction_id"], systems, strong_reactions,
            args.pathway_context_min_other,
            args.pathway_context_min_fraction,
        )
        treatment, final_conf, countable = adjudicate(
            r["diamond_confidence"],
            r["policy_interpro_status"],
            pathway_ok=pok,
            genomic_context=bool(r["genomic_context_support"]),
        )
        x = dict(r)
        x.pop("_bitscore_float", None)
        x["pathway_context_support"] = pok
        x["pathway_context_fraction"] = round(pfrac, 3)
        x["pathway_context_detail"] = pdetail
        x["adjudication"] = treatment
        x["count_toward_pathway"] = countable
        x["candidate_detected"] = (r["diamond_confidence"] in {"High", "Moderate"}) or countable
        x["exact_locus_supported"] = countable
        x["confidence"] = final_conf
        warnings = []
        if r["diamond_confidence"] == "Moderate":
            warnings.append("MODERATE_DIAMOND_REVIEW")
        if r["policy_interpro_status"] == "CONFLICTING":
            warnings.append("INTERPRO_CONFLICT")
        elif r["policy_interpro_status"] == "AMBIGUOUS":
            warnings.append("INTERPRO_AMBIGUOUS")
        elif r["policy_interpro_status"] == "SPECIFIC":
            warnings.append("INTERPRO_SPECIFIC_SUPPORT")
        if r["diamond_confidence"] == "Weak":
            warnings.append("WEAK_RESCUE")
        x["evidence_warning"] = "|".join(warnings)
        source = ["DIAMOND"]
        if r["policy_interpro_status"] == "SPECIFIC":
            source.append("InterPro-specific")
        elif r["policy_interpro_status"] == "COMPATIBLE":
            source.append("InterPro-compatible")
        if r["genomic_context_support"]:
            source.append("GenomicContext")
        if pok:
            source.append("PathwayContext")
        x["evidence_source"] = "+".join(source)
        final_rows.append(x)
        if countable:
            countable_rows.append(x)

    # Write complete assessed candidate table.
    all_out = args.out_prefix + "_all_candidates.tsv"
    cols = [
        "reaction_id","gene_family","enzyme_name","locus_tag","gene","product","contig",
        "reference_id","reference_name","reference_accession","component_id","component_gene_family",
        "identity_pct","query_coverage_pct","reference_coverage_pct","evalue","bitscore",
        "diamond_confidence","policy_interpro_status","policy_interpro_detail",
        "genomic_context_support","genomic_context_detail",
        "pathway_context_support","pathway_context_fraction","pathway_context_detail",
        "adjudication","count_toward_pathway","candidate_detected","exact_locus_supported","confidence","evidence_warning","evidence_source",
    ]
    with open(all_out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(final_rows, key=lambda x:(x["reaction_id"],x["locus_tag"])))

    # Safe Pathway Engine input: only countable exact-reaction assignments.
    adj_out = args.out_prefix + "_adjudicated_rescued_evidence.tsv"
    with open(adj_out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(countable_rows, key=lambda x:(x["reaction_id"],x["locus_tag"])))

    # Reaction summary.
    by_rxn = defaultdict(list)
    for r in final_rows:
        by_rxn[r["reaction_id"]].append(r)

    sum_out = args.out_prefix + "_reaction_adjudication_summary.tsv"
    sfields = [
        "reaction_id","gene_family","enzyme_name","candidate_loci","countable_loci",
        "review_loci","rejected_loci","best_confidence","reaction_status",
        "countable_locus_tags","review_locus_tags","rejected_locus_tags",
    ]
    with open(sum_out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=sfields, delimiter="\t")
        w.writeheader()
        for rid in sorted(reactions):
            rr = by_rxn.get(rid, [])
            counted = [x for x in rr if bool(x["count_toward_pathway"])]
            review = [x for x in rr if x["adjudication"] == REVIEW]
            rejected = [x for x in rr if x["adjudication"] == REJECT]
            if counted:
                status = "ACCEPTED"
                best_conf = max((x["confidence"] for x in counted), key=lambda c: CONF_RANK.get(c, 0))
            elif review:
                status = "REVIEW"
                best_conf = max((x["confidence"] for x in review), key=lambda c: CONF_RANK.get(c, 0))
            elif rejected:
                status = "REJECTED"
                best_conf = ""
            else:
                status = "NO_CANDIDATE"
                best_conf = ""
            w.writerow({
                "reaction_id": rid,
                "gene_family": reactions[rid].get("gene_family", ""),
                "enzyme_name": reactions[rid].get("enzyme_name", ""),
                "candidate_loci": len({x["locus_tag"] for x in rr}),
                "countable_loci": len({x["locus_tag"] for x in counted}),
                "review_loci": len({x["locus_tag"] for x in review}),
                "rejected_loci": len({x["locus_tag"] for x in rejected}),
                "best_confidence": best_conf,
                "reaction_status": status,
                "countable_locus_tags": "|".join(sorted({x["locus_tag"] for x in counted})),
                "review_locus_tags": "|".join(sorted({x["locus_tag"] for x in review})),
                "rejected_locus_tags": "|".join(sorted({x["locus_tag"] for x in rejected})),
            })

    # Reaction-level evidence summary: separate biological signal from exact-locus support.
    # This is the key v0.3.2 distinction. A Moderate/conflicting locus may still
    # indicate a reaction candidate, while not being claimed as an exact enzyme ID.
    reaction_evidence = {}
    component_summary_rows = []
    for rid in sorted(reactions):
        rr = by_rxn.get(rid, [])
        req_components = component_requirements.get(rid, [])
        if req_components:
            comp_states={}
            all_candidate_tags=set(); all_supported_tags=set()
            any_conflict=False; any_review=False; seqs=[]
            for comp in req_components:
                cr=[x for x in rr if x.get("component_id","")==comp]
                cc=[x for x in cr if bool(x["candidate_detected"])]
                cs=[x for x in cr if bool(x["exact_locus_supported"])]
                cf=[x for x in cr if x["policy_interpro_status"]=="CONFLICTING"]
                rv=[x for x in cr if x["adjudication"]==REVIEW]
                comp_states[comp]={"candidate":bool(cc),"supported":bool(cs),
                                   "candidate_tags":{x["locus_tag"] for x in cc},
                                   "supported_tags":{x["locus_tag"] for x in cs},
                                   "conflicting":bool(cc) and len(cf)==len(cc)}
                all_candidate_tags |= comp_states[comp]["candidate_tags"]
                all_supported_tags |= comp_states[comp]["supported_tags"]
                any_conflict = any_conflict or comp_states[comp]["conflicting"]
                any_review = any_review or bool(rv)
                seqs += [x["diamond_confidence"] for x in cc]
                component_summary_rows.append({
                    "reaction_id":rid,"component_id":comp,
                    "component_required":"yes","candidate_present":bool(cc),
                    "exact_supported":bool(cs),
                    "candidate_locus_tags":"|".join(sorted(comp_states[comp]["candidate_tags"])),
                    "supported_locus_tags":"|".join(sorted(comp_states[comp]["supported_tags"])),
                })
            parent_candidate=all(comp_states[c]["candidate"] for c in req_components)
            parent_supported=all(comp_states[c]["supported"] for c in req_components)
            if parent_supported:
                state="SUPPORTED"
            elif parent_candidate and any_conflict:
                state="CANDIDATE_CONFLICTING"
            elif parent_candidate:
                state="CANDIDATE_REVIEW"
            else:
                state="NO_CANDIDATE"
            best_seq=max(seqs,key=lambda c:CONF_RANK.get(c,0)) if seqs else ""
            reaction_evidence[rid]={
                "reaction_id":rid,
                "gene_family":reactions[rid].get("gene_family",""),
                "enzyme_name":reactions[rid].get("enzyme_name",""),
                "reaction_signal":state,
                "best_sequence_evidence":best_seq,
                "candidate_loci":len(all_candidate_tags),
                "supported_exact_loci":len(all_supported_tags),
                "review_loci":sum(1 for x in rr if x["adjudication"]==REVIEW),
                "conflicting_loci":sum(1 for x in rr if x["policy_interpro_status"]=="CONFLICTING"),
                "candidate_locus_tags":"|".join(sorted(all_candidate_tags)),
                "supported_locus_tags":"|".join(sorted(all_supported_tags)),
                "component_rule":"ALL_REQUIRED_COMPONENTS",
                "required_components":"|".join(req_components),
                "candidate_components":"|".join(c for c in req_components if comp_states[c]["candidate"]),
                "supported_components":"|".join(c for c in req_components if comp_states[c]["supported"]),
            }
        else:
            cand = [x for x in rr if bool(x["candidate_detected"])]
            supported = [x for x in rr if bool(x["exact_locus_supported"])]
            conflicts = [x for x in rr if x["policy_interpro_status"] == "CONFLICTING"]
            reviews = [x for x in rr if x["adjudication"] == REVIEW]
            if supported:
                state = "SUPPORTED"
            elif cand and conflicts and len(conflicts) == len(cand):
                state = "CANDIDATE_CONFLICTING"
            elif cand:
                state = "CANDIDATE_REVIEW"
            else:
                state = "NO_CANDIDATE"
            best_seq = max((x["diamond_confidence"] for x in cand), key=lambda c: CONF_RANK.get(c,0)) if cand else ""
            reaction_evidence[rid] = {
                "reaction_id": rid,
                "gene_family": reactions[rid].get("gene_family", ""),
                "enzyme_name": reactions[rid].get("enzyme_name", ""),
                "reaction_signal": state,
                "best_sequence_evidence": best_seq,
                "candidate_loci": len({x["locus_tag"] for x in cand}),
                "supported_exact_loci": len({x["locus_tag"] for x in supported}),
                "review_loci": len({x["locus_tag"] for x in reviews}),
                "conflicting_loci": len({x["locus_tag"] for x in conflicts}),
                "candidate_locus_tags": "|".join(sorted({x["locus_tag"] for x in cand})),
                "supported_locus_tags": "|".join(sorted({x["locus_tag"] for x in supported})),
                "component_rule":"",
                "required_components":"",
                "candidate_components":"",
                "supported_components":"",
            }

    # Parent exact-locus evidence for component reactions is safe only when all
    # required components are exact-supported. Remove incomplete component-system
    # rows from the rescued file.
    if component_requirements:
        safe=[]
        for x in countable_rows:
            rid=x["reaction_id"]
            if rid not in component_requirements:
                safe.append(x)
            elif reaction_evidence.get(rid,{}).get("reaction_signal")=="SUPPORTED":
                safe.append(x)
        countable_rows=safe

    comp_out=args.out_prefix+"_component_evidence_summary.tsv"
    if component_summary_rows:
        cfields=["reaction_id","component_id","component_required","candidate_present","exact_supported",
                 "candidate_locus_tags","supported_locus_tags"]
        with open(comp_out,"w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=cfields,delimiter="\t");w.writeheader();w.writerows(component_summary_rows)

    rxev_out = args.out_prefix + "_reaction_evidence_summary.tsv"
    rxev_fields = ["reaction_id","gene_family","enzyme_name","reaction_signal","best_sequence_evidence",
                   "candidate_loci","supported_exact_loci","review_loci","conflicting_loci",
                   "candidate_locus_tags","supported_locus_tags","component_rule","required_components",
                   "candidate_components","supported_components"]
    with open(rxev_out, "w", newline="", encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=rxev_fields, delimiter="\t")
        w.writeheader(); w.writerows(reaction_evidence[r] for r in sorted(reaction_evidence))

    # Pathway-level evidence has TWO completeness values:
    # candidate completeness (High/Moderate sequence signal retained) and
    # supported completeness (at least one exact locus not contradicted/ambiguous).
    pathway_rows = defaultdict(list)
    for rid, memberships in systems.items():
        for m in memberships:
            pid=m.get("pathway_id","")
            if pid:
                pathway_rows[pid].append((rid,m))
    pev_out = args.out_prefix + "_pathway_evidence_summary.tsv"
    pev_fields = ["pathway_id","pathway_name","total_reactions","candidate_reactions","candidate_completeness_pct",
                  "supported_reactions","supported_completeness_pct","review_or_conflicting_reactions",
                  "missing_candidate_reactions","evidence_interpretation"]
    with open(pev_out,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=pev_fields,delimiter="\t"); w.writeheader()
        for pid in sorted(pathway_rows):
            # Score only required/scored reactions. Optional rows are retained in
            # pathways.tsv for biological context but must not inflate the
            # completeness denominator.
            unique={}
            for rid,m in pathway_rows[pid]:
                required=str(m.get("required","")).strip().lower()
                if required in {"yes","true","1","required"}:
                    unique[rid]=m
            # Legacy fallback for modules whose pathway rows do not explicitly
            # define required=yes: score rows that are not explicitly optional.
            if not unique:
                for rid,m in pathway_rows[pid]:
                    required=str(m.get("required","")).strip().lower()
                    if required not in {"no","false","0"}:
                        unique[rid]=m
            rids=list(unique)
            total=len(rids)
            candidate=[rid for rid in rids if reaction_evidence.get(rid,{}).get("reaction_signal") != "NO_CANDIDATE"]
            supported=[rid for rid in rids if reaction_evidence.get(rid,{}).get("reaction_signal") == "SUPPORTED"]
            review=[rid for rid in rids if reaction_evidence.get(rid,{}).get("reaction_signal") in {"CANDIDATE_REVIEW","CANDIDATE_CONFLICTING"}]
            missing=[rid for rid in rids if reaction_evidence.get(rid,{}).get("reaction_signal") == "NO_CANDIDATE"]
            cpc=round(100*len(candidate)/total,1) if total else 0.0
            spc=round(100*len(supported)/total,1) if total else 0.0
            if total and len(candidate)==total and not review:
                interp="COMPLETE_SUPPORTED"
            elif total and len(candidate)==total:
                interp="COMPLETE_CANDIDATE_WITH_REVIEW"
            elif candidate:
                interp="PARTIAL_CANDIDATE"
            else:
                interp="NOT_DETECTED"
            pname=next((m.get("pathway_name","") for _,m in pathway_rows[pid] if m.get("pathway_name","")),"")
            w.writerow({"pathway_id":pid,"pathway_name":pname,"total_reactions":total,
                        "candidate_reactions":len(candidate),"candidate_completeness_pct":cpc,
                        "supported_reactions":len(supported),"supported_completeness_pct":spc,
                        "review_or_conflicting_reactions":"|".join(sorted(review)),
                        "missing_candidate_reactions":"|".join(sorted(missing)),
                        "evidence_interpretation":interp})

    print()
    print(f"Hybrid Rescue Engine v{VERSION}")
    print(f"Assessed reaction-locus candidates: {len(final_rows)}")
    print(f"Countable reaction-locus assignments: {len(countable_rows)}")
    print(f"Rejected exact-reaction assignments: {sum(1 for x in final_rows if x['adjudication']==REJECT)}")
    print(f"Manual-review assignments: {sum(1 for x in final_rows if x['adjudication']==REVIEW)}")
    print(f"Wrote: {all_out}")
    print(f"Wrote: {adj_out}")
    print(f"Wrote: {sum_out}")
    print(f"Wrote: {rxev_out}")
    if component_summary_rows:
        print(f"Wrote: {comp_out}")
    print(f"Wrote: {pev_out}")
    if ipfile:
        print(f"InterPro results used: {ipfile}")
    print("Computational assignments require biological/experimental validation.")


if __name__ == "__main__":
    main()
