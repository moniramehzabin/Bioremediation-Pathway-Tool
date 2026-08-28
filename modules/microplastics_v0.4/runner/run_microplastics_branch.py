#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, subprocess, sys
from pathlib import Path

MICROPLASTIC_IDS = ("BRXN_PET_HYD", "BRXN_MHET_HYD", "BRXN_TPA_DIOX", "BRXN_TPA_DCD_DH", "BRXN_PCL_DEP", "BRXN_PBAT_DEP", "BRXN_PU_DEP")
PET_IDS = MICROPLASTIC_IDS[:4]
PCL_ID = "BRXN_PCL_DEP"
PBAT_ID = "BRXN_PBAT_DEP"
PU_ID = "BRXN_PU_DEP"

def require_file(p: Path, label: str) -> Path:
    if not p.is_file():
        raise SystemExit(f"ERROR: {label} not found: {p}")
    return p

def require_dir(p: Path, label: str) -> Path:
    if not p.is_dir():
        raise SystemExit(f"ERROR: {label} not found: {p}")
    return p

def run(cmd):
    print("\nRUN:", " ".join(map(str, cmd)))
    p = subprocess.run(cmd)
    if p.returncode:
        raise SystemExit(p.returncode)

def read_tsv(p: Path):
    with p.open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def exact_microplastics_summary(reaction_file: Path, out_file: Path):
    rows = {r.get("reaction_id", ""): r for r in read_tsv(reaction_file)}
    with out_file.open("w", encoding="utf-8", newline="") as f:
        fields = ["reaction_id", "enzyme_name", "integrated_state", "review_flag", "supported_loci", "review_loci", "candidate_only_loci"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for rid in MICROPLASTIC_IDS:
            r = rows.get(rid, {})
            w.writerow({k: r.get(k, "") for k in fields} | {"reaction_id": rid})

    print("\nMICROPLASTICS EXACT REACTION SUMMARY")
    print("-" * 72)
    for rid in MICROPLASTIC_IDS:
        r = rows.get(rid, {})
        state = r.get("integrated_state", "MISSING_FROM_OUTPUT")
        print(f"{rid:20s} {state}")
    print("-" * 72)
    print("Exact microplastics summary:", out_file)

def apply_pcl_strict_gate(reaction_file: Path, pathway_file: Path, candidates_file: Path, audit_file: Path):
    """
    Conservative PCL gate for v0.3.
    Exact PCL annotation can remain supported. Sequence-only rescue is accepted as
    supported only when the best countable Q6A0I4-like hit has >=60% identity and
    >=80% query AND reference coverage. Otherwise any detected PCL rescue remains REVIEW.
    This does not alter PET reactions or the protected main database.
    """
    rrows = read_tsv(reaction_file)
    prows = read_tsv(pathway_file)
    crows = [r for r in read_tsv(candidates_file) if r.get("reaction_id") == PCL_ID]
    rr = next((r for r in rrows if r.get("reaction_id") == PCL_ID), None)
    if rr is None:
        raise SystemExit("ERROR: PCL reaction missing from adjudicated output")

    exact_annotation = str(rr.get("annotation_supported", "")).strip().lower() in {"true","yes","1"}
    eligible=[]
    for c in crows:
        try:
            pid=float(c.get("identity_pct") or 0); q=float(c.get("query_coverage_pct") or 0); sc=float(c.get("reference_coverage_pct") or 0)
        except ValueError:
            continue
        if (pid >= 60.0 and q >= 80.0 and sc >= 80.0
            and c.get("diamond_confidence") == "High"
            and str(c.get("policy_interpro_status", "")) not in {"CONFLICTING","AMBIGUOUS"}
            and str(c.get("count_toward_pathway", "")).strip().lower() in {"true","yes","1"}):
            eligible.append((pid,q,sc,c))

    gate_pass = exact_annotation or bool(eligible)
    detected = bool(crows) or str(rr.get("candidate_present", "")).strip().lower() in {"true","yes","1"}
    action="KEEP_SUPPORTED" if gate_pass else ("DOWNGRADE_TO_REVIEW" if detected else "NO_CANDIDATE")

    if not gate_pass and detected:
        loci=sorted({c.get("locus_tag","") for c in crows if c.get("locus_tag")})
        rr["integrated_state"]="CANDIDATE_STRICT_REVIEW"
        rr["review_flag"]="YES"
        rr["exact_supported"]="False"
        rr["supported_loci"]=""
        rr["review_loci"]="|".join(loci)
        rr["candidate_present"]="True"
        rr["evidence_relationship"]="PCL_STRICT_SEQUENCE_GATE"
        rr["rescue_signal"]="REVIEW"

        for pr in prows:
            if pr.get("pathway_id") == "BMOD_PCL_DEPOLYMERIZATION":
                pr["candidate_units"]="1"
                pr["candidate_percent"]="100.0"
                pr["supported_units"]="0"
                pr["supported_percent"]="0.0"
                pr["review_units"]="1"
                pr["adjudicated_pathway_state"]="REVIEW"
                pr["missing_required_units"]=PCL_ID

    # rewrite using original field order
    with reaction_file.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rrows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rrows)
    with pathway_file.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(prows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(prows)

    best=max(eligible, default=None, key=lambda x:(x[0],x[1],x[2]))
    with audit_file.open("w",encoding="utf-8",newline="") as f:
        fields=["reaction_id","exact_annotation_supported","candidate_rows","strict_gate_pass","action","best_identity_pct","best_query_coverage_pct","best_reference_coverage_pct","threshold"]
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerow({
            "reaction_id":PCL_ID,"exact_annotation_supported":exact_annotation,"candidate_rows":len(crows),
            "strict_gate_pass":gate_pass,"action":action,
            "best_identity_pct":best[0] if best else "","best_query_coverage_pct":best[1] if best else "","best_reference_coverage_pct":best[2] if best else "",
            "threshold":"identity>=60;qcov>=80;scov>=80;High;non-conflicting"})
    print("PCL STRICT GATE:", action)
    print("PCL strict-gate audit:", audit_file)

def apply_pbat_strict_gate(reaction_file: Path, pathway_file: Path, candidates_file: Path, audit_file: Path):
    """Conservative PBAT gate: exact PBAT annotation or strong TfCut-like sequence only."""
    rrows = read_tsv(reaction_file); prows = read_tsv(pathway_file)
    crows = [r for r in read_tsv(candidates_file) if r.get("reaction_id") == PBAT_ID]
    rr = next((r for r in rrows if r.get("reaction_id") == PBAT_ID), None)
    if rr is None: raise SystemExit("ERROR: PBAT reaction missing from adjudicated output")
    exact_annotation = str(rr.get("annotation_supported", "")).strip().lower() in {"true","yes","1"}
    eligible=[]
    for c in crows:
        try: pid=float(c.get("identity_pct") or 0); q=float(c.get("query_coverage_pct") or 0); sc=float(c.get("reference_coverage_pct") or 0)
        except ValueError: continue
        if (pid >= 60.0 and q >= 80.0 and sc >= 80.0 and c.get("diamond_confidence") == "High"
            and str(c.get("policy_interpro_status", "")) not in {"CONFLICTING","AMBIGUOUS"}
            and str(c.get("count_toward_pathway", "")).strip().lower() in {"true","yes","1"}): eligible.append((pid,q,sc,c))
    gate_pass = exact_annotation or bool(eligible)
    detected = bool(crows) or str(rr.get("candidate_present", "")).strip().lower() in {"true","yes","1"}
    action="KEEP_SUPPORTED" if gate_pass else ("DOWNGRADE_TO_REVIEW" if detected else "NO_CANDIDATE")
    if not gate_pass and detected:
        loci=sorted({c.get("locus_tag","") for c in crows if c.get("locus_tag")})
        rr["integrated_state"]="CANDIDATE_STRICT_REVIEW"; rr["review_flag"]="YES"; rr["exact_supported"]="False"; rr["supported_loci"]=""
        rr["review_loci"]="|".join(loci); rr["candidate_present"]="True"; rr["evidence_relationship"]="PBAT_STRICT_SEQUENCE_GATE"; rr["rescue_signal"]="REVIEW"
        for pr in prows:
            if pr.get("pathway_id") == "BMOD_PBAT_DEPOLYMERIZATION":
                pr["candidate_units"]="1"; pr["candidate_percent"]="100.0"; pr["supported_units"]="0"; pr["supported_percent"]="0.0"; pr["review_units"]="1"; pr["adjudicated_pathway_state"]="REVIEW"; pr["missing_required_units"]=PBAT_ID
    with reaction_file.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rrows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rrows)
    with pathway_file.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(prows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(prows)
    best=max(eligible, default=None, key=lambda x:(x[0],x[1],x[2]))
    with audit_file.open("w",encoding="utf-8",newline="") as f:
        fields=["reaction_id","exact_annotation_supported","candidate_rows","strict_gate_pass","action","best_identity_pct","best_query_coverage_pct","best_reference_coverage_pct","threshold"]
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerow({"reaction_id":PBAT_ID,"exact_annotation_supported":exact_annotation,"candidate_rows":len(crows),"strict_gate_pass":gate_pass,"action":action,"best_identity_pct":best[0] if best else "","best_query_coverage_pct":best[1] if best else "","best_reference_coverage_pct":best[2] if best else "","threshold":"identity>=60;qcov>=80;scov>=80;High;non-conflicting"})
    print("PBAT STRICT GATE:", action); print("PBAT strict-gate audit:", audit_file)

def main():
    ap = argparse.ArgumentParser(description="Isolated PET+PCL+PBAT+PU microplastics optional-branch runner. Does not modify the main v0.7.6 database.")
    ap.add_argument("genbank", nargs="?", help="Input bacterial GenBank file")
    ap.add_argument("--baseline-db", default="Bioremediation_DB_v0.7.6_dyes", help="Protected main database used only for SHA-256 guard")
    ap.add_argument("--interpro-tsv", help="Precomputed InterPro TSV for the genome")
    ap.add_argument("--diamond", default="diamond.exe")
    ap.add_argument("--out-prefix", default="PET_branch_run")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--self-test", action="store_true", help="Validate package + baseline only; do not analyze genome")
    a = ap.parse_args()

    runner = Path(__file__).resolve().parent
    branch = runner.parent
    db = require_dir(branch / "database", "microplastics branch database")
    guard = require_file(branch / "guard" / "verify_baseline_unchanged.py", "baseline guard")
    baseline = require_dir(Path(a.baseline_db).resolve(), "protected v0.7.6 baseline database")

    # Always protect the main database before doing anything else.
    run([a.python, str(guard), str(baseline)])

    expected = {
        "genbank validator": runner / "genbank_input_validator_v082.py",
        "resolved annotation scanner": runner / "resolved_annotation_scanner_v0810.py",
        "hybrid rescue": runner / "rescue" / "hybrid_rescue_engine_v0.3.5_resolvedskip.py",
        "integrated evidence": runner / "integrated_evidence_v063.py",
        "cross-reaction adjudicator": runner / "cross_reaction_adjudicator_v088.py",
        "table reporter": runner / "table_first_reporter_v089.py",
        "microplastics InterPro policy": branch / "reaction_interpro_policy_microplastics_v04.tsv",
    }
    for label, p in expected.items(): require_file(p, label)
    for fn in ["reference_proteins.faa", "reference_metadata.tsv", "reactions.tsv", "pathways.tsv", "reaction_components.tsv", "enzyme_aliases.tsv"]:
        require_file(db / fn, f"microplastics database {fn}")

    # Hard isolation checks.
    rids = {r.get("reaction_id", "") for r in read_tsv(db / "reactions.tsv")}
    if rids != set(MICROPLASTIC_IDS):
        raise SystemExit(f"ERROR: microplastics branch reaction set is not exact: {sorted(rids)}")
    if "BRXN_COP_ATPASE" in rids:
        raise SystemExit("ERROR: copper ATPase contamination detected in PET branch")

    print("\nMICROPLASTICS BRANCH SELF-CHECK: PASS")
    print("Main v0.7.6 baseline: hash-locked and unchanged")
    print("Reaction IDs:", ", ".join(MICROPLASTIC_IDS))
    print("PCL strict sequence gate: >=60% identity and >=80% query/reference coverage")
    print("PBAT: curated TfCut seed; generic cutinase wording alone is insufficient")
    print("PU: exact pueA/pueB/polyurethanase annotation only in v0.4; sequence rescue held out pending checksum verification")
    print("Copper ATPase contamination: NO")
    print("Main database write access: NONE")

    if a.self_test:
        return
    if not a.genbank or not a.interpro_tsv:
        raise SystemExit("ERROR: genbank and --interpro-tsv are required unless --self-test is used")

    gb = require_file(Path(a.genbank).resolve(), "GenBank")
    ipr = require_file(Path(a.interpro_tsv).resolve(), "InterPro TSV")
    diamond = require_file(Path(a.diamond).resolve(), "DIAMOND executable")
    prefix = Path(a.out_prefix)

    norm = Path(str(prefix) + "_normalized_input.gbk")
    audit = Path(str(prefix) + "_GENBANK_INPUT_AUDIT.txt")
    resolved = Path(str(prefix) + "_resolved_annotation_pairs.tsv")
    rpre = str(prefix) + "_rescue"
    ipre = str(prefix) + "_integrated"
    apre = str(prefix) + "_adjudicated"

    run([a.python, str(expected["genbank validator"]), str(gb), "--out", str(norm), "--audit", str(audit)])
    used = norm if norm.exists() and norm.stat().st_size else gb

    run([a.python, str(expected["resolved annotation scanner"]), str(used), "--db", str(db), "--out", str(resolved)])

    run([a.python, str(expected["hybrid rescue"]), str(used),
         "--references", str(db / "reference_proteins.faa"),
         "--reference-metadata", str(db / "reference_metadata.tsv"),
         "--reactions", str(db / "reactions.tsv"),
         "--pathways", str(db / "pathways.tsv"),
         "--reaction-components", str(db / "reaction_components.tsv"),
         "--resolved-reaction-loci", str(resolved),
         "--policy", str(expected["microplastics InterPro policy"]),
         "--diamond", str(diamond),
         "--interpro-mode", "precomputed", "--interpro-tsv", str(ipr),
         "--out-prefix", rpre])

    rsum = Path(rpre + "_reaction_evidence_summary.tsv")
    rall = Path(rpre + "_all_candidates.tsv")
    require_file(rsum, "microplastics rescue reaction summary")
    require_file(rall, "microplastics rescue candidate table")

    run([a.python, str(expected["integrated evidence"]), str(used), "--db", str(db),
         "--reaction-evidence", str(rsum), "--out-prefix", ipre])

    ire = Path(ipre + "_integrated_reaction_evidence.tsv")
    require_file(ire, "microplastics integrated reaction evidence")

    run([a.python, str(expected["cross-reaction adjudicator"]), str(used), "--db", str(db),
         "--integrated-reaction-evidence", str(ire), "--all-candidates", str(rall), "--out-prefix", apre])

    are = Path(apre + "_adjudicated_reaction_evidence.tsv")
    aps = Path(apre + "_adjudicated_pathway_summary.tsv")
    dec = Path(apre + "_cross_reaction_locus_decisions.tsv")
    for p, label in [(are,"microplastics adjudicated reaction evidence"),(aps,"microplastics adjudicated pathway summary"),(dec,"microplastics cross-reaction decisions")]: require_file(p,label)

    pcl_audit = Path(str(prefix) + "_PCL_STRICT_GATE.tsv")
    apply_pcl_strict_gate(are, aps, rall, pcl_audit)
    pbat_audit = Path(str(prefix) + "_PBAT_STRICT_GATE.tsv")
    apply_pbat_strict_gate(are, aps, rall, pbat_audit)

    html = Path(str(prefix) + "_MICROPLASTICS_TABLE_REPORT.html")
    run([a.python, str(expected["table reporter"]), "--db", str(db), "--pathway-summary", str(aps),
         "--reaction-evidence", str(are), "--cross-reaction-decisions", str(dec),
         "--out-html", str(html), "--sample-name", gb.stem])

    exact = Path(str(prefix) + "_MICROPLASTICS_EXACT_REACTIONS.tsv")
    exact_microplastics_summary(are, exact)

    manifest = {
        "branch": "Microplastics Optional Branch v0.4 PET+PCL+PBAT+PU",
        "main_database_modified": False,
        "baseline_database": str(baseline),
        "microplastics_database": str(db),
        "genbank": str(gb),
        "outputs": {"pathway_summary":str(aps),"reaction_evidence":str(are),"decisions":str(dec),"table_report":str(html),"exact_microplastics_reactions":str(exact),"pcl_strict_gate":str(pcl_audit),"pbat_strict_gate":str(pbat_audit)}
    }
    Path(str(prefix) + "_MICROPLASTICS_RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("\nMICROPLASTICS OPTIONAL-BRANCH RUN COMPLETE")
    print("Main v0.7.6 database was not modified.")
    print("Table report:", html)

if __name__ == "__main__":
    main()
