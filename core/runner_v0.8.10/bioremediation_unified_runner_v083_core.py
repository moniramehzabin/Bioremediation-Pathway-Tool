#!/usr/bin/env python3
"""
Bioremediation Unified Runner core pipeline.

GenBank
 -> exact-annotation resolved-pair scan
 -> Hybrid Rescue v0.3.5_resolvedskip
 -> Integrated Evidence v0.6.3
 -> Cross-Reaction Locus Adjudication v0.8.8
 -> alternative-aware pathway scoring
 -> reports

This runner only wires validated components together. It does not change
biological scoring thresholds or evidence interpretation.
"""
from __future__ import annotations
import argparse, csv, json, subprocess, sys
from pathlib import Path
from datetime import datetime

VERSION = "0.8.10"

def check(p, label, dir=False):
    p = Path(p)
    ok = p.is_dir() if dir else p.is_file()
    if not ok:
        raise SystemExit(f"ERROR: {label} not found: {p}")
    return p

def run(cmd, log):
    print("\nRUN:", " ".join(map(str, cmd)))
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    print(p.stdout, end="")
    with open(log, "a", encoding="utf-8") as f:
        f.write("\nRUN: " + " ".join(map(str, cmd)) + "\n" + p.stdout)
    if p.returncode:
        raise SystemExit(f"ERROR: command failed ({p.returncode}); see {log}")

def tsv(p):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def report(prefix, pathway, reaction, decisions):
    P = tsv(pathway)
    R = tsv(reaction)
    D = tsv(decisions)
    rank = {
        "SUPPORTED_COMPLETE": 0,
        "COMPLETE_CANDIDATE_REVIEW_REQUIRED": 1,
        "PARTIAL_WITH_STRONG_LOCAL_CONTEXT": 2,
        "PARTIAL_WITH_LOCAL_CONTEXT": 3,
        "PARTIAL": 4,
        "NOT_DETECTED": 5,
    }
    P = sorted(
        P,
        key=lambda r: (
            rank.get(r.get("adjudicated_pathway_state", ""), 9),
            r.get("pathway_name", ""),
        ),
    )
    fp = Path(prefix + "_REPORT.txt")
    with open(fp, "w", encoding="utf-8") as f:
        f.write("BIOREMEDIATION PATHWAY REPORT — v0.8.10\n")
        f.write("=" * 76 + "\n\n")
        f.write("Cross-reaction specificity competition is enabled.\n")
        f.write("One shared protein is not automatically counted for several unrelated reactions.\n")
        f.write("Alternative reaction groups are scored as OR.\n")
        f.write("Genomic context never creates a missing reaction.\n\n")
        for r in P:
            f.write(f"{r.get('pathway_name')} [{r.get('pathway_id')}]\n")
            f.write(f"  State: {r.get('adjudicated_pathway_state')}\n")
            f.write(
                f"  Candidate units: {r.get('candidate_units')}/{r.get('required_units')} "
                f"({r.get('candidate_percent')}%)\n"
            )
            f.write(
                f"  Supported units: {r.get('supported_units')}/{r.get('required_units')} "
                f"({r.get('supported_percent')}%)\n"
            )
            if r.get("cluster_reaction_count") not in {"", "0"}:
                f.write(
                    f"  Local cluster: {r.get('cluster_reaction_count')} reactions; "
                    f"{r.get('cluster_contig')}; span {r.get('cluster_span_bp')} bp\n"
                )
            if r.get("missing_required_units"):
                f.write(f"  Missing: {r.get('missing_required_units')}\n")
            f.write("\n")

        f.write("\nINFORMATIVE REACTION CALLS\n" + "=" * 76 + "\n")
        for r in sorted(R, key=lambda x: x.get("reaction_id", "")):
            if r.get("integrated_state") == "NO_CANDIDATE":
                continue
            f.write(f"{r.get('reaction_id')} | {r.get('enzyme_name','')}\n")
            f.write(f"  State: {r.get('integrated_state')}\n")
            f.write(f"  Supported loci: {r.get('supported_loci') or '-'}\n")
            f.write(f"  Evidence relationship: {r.get('evidence_relationship') or '-'}\n")
            if r.get("review_loci"):
                f.write(f"  Review loci: {r.get('review_loci')}\n")
            if r.get("candidate_only_loci"):
                f.write(f"  Candidate-only loci: {r.get('candidate_only_loci')}\n")
            if r.get("removed_competing_loci"):
                f.write(f"  Competing loci removed: {r.get('removed_competing_loci')}\n")
            if r.get("ambiguous_loci"):
                f.write(f"  Ambiguous loci: {r.get('ambiguous_loci')}\n")
            if r.get("review_flag") == "YES":
                f.write("  REVIEW REQUIRED\n")
            f.write("\n")

        f.write("\nCROSS-REACTION COMPETITIONS\n" + "=" * 76 + "\n")
        for d in D:
            if d.get("decision") in {"WINNER_KEEP", "LOSER_REMOVE", "AMBIGUOUS_REVIEW"}:
                f.write(
                    f"{d.get('locus_tag')} | {d.get('reaction_id')} | "
                    f"{d.get('decision')} | score={d.get('score')}"
                )
                if d.get("winner_reaction"):
                    f.write(f" | winner={d.get('winner_reaction')}")
                f.write("\n")
    return fp

def add_interpro_cli(ap):
    ap.add_argument(
        "--interpro-mode",
        choices=["auto", "precomputed", "local", "web", "off"],
        default="auto",
    )
    ap.add_argument("--interpro-tsv")
    ap.add_argument("--interproscan")
    ap.add_argument("--interpro-applications", default=None)
    ap.add_argument("--email")
    ap.add_argument("--interpro-max-jobs", type=int, default=200)
    ap.add_argument("--enable-weak-rescue", action="store_true")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--min-identity", type=float, default=25.0)
    ap.add_argument("--min-qcov", type=float, default=50.0)
    ap.add_argument("--min-scov", type=float, default=50.0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genbank")
    ap.add_argument("--db", default=None)
    add_interpro_cli(ap)
    ap.add_argument("--diamond", required=True)
    ap.add_argument(
        "--rescue-engine",
        default=str(
            Path(__file__).resolve().parent
            / "Bioremediation_Hybrid_Rescue_v0.3.5_resolvedskip"
            / "hybrid_rescue_engine_v0.3.5_resolvedskip.py"
        ),
    )
    ap.add_argument("--integrated-engine", default=None)
    ap.add_argument(
        "--adjudicator",
        default=str(Path(__file__).resolve().parent / "cross_reaction_adjudicator_v088.py"),
    )
    ap.add_argument(
        "--policy",
        default=str(
            Path(__file__).resolve().parent
            / "Bioremediation_Hybrid_Rescue_v0.3.5_resolvedskip"
            / "reaction_interpro_policy_v054.tsv"
        ),
    )
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()

    here = Path(__file__).resolve().parent
    repo_root = here.parent.parent
    if a.db is None:
        a.db = str(repo_root / "core" / "database_v0.7.6")
    if a.integrated_engine is None:
        a.integrated_engine = str(
            repo_root / "core" / "integrated_evidence_v0.6.3" / "integrated_evidence_v063.py"
        )

    if a.interpro_mode == "precomputed" and not a.interpro_tsv:
        raise SystemExit("ERROR: --interpro-mode precomputed requires --interpro-tsv.")
    if a.interpro_mode == "local" and not a.interproscan:
        raise SystemExit("ERROR: --interpro-mode local requires --interproscan.")
    if a.interpro_mode == "web" and not a.email:
        raise SystemExit("ERROR: --interpro-mode web requires --email.")

    gb = check(a.genbank, "GenBank")
    db = check(a.db, "database", True)
    if a.interpro_tsv:
        check(a.interpro_tsv, "InterPro TSV")
    rescue = check(a.rescue_engine, "rescue engine")
    integ = check(a.integrated_engine, "integrated engine")
    adjud = check(a.adjudicator, "adjudicator")
    policy = check(a.policy, "policy")
    for fn in [
        "reference_proteins.faa",
        "reference_metadata.tsv",
        "reactions.tsv",
        "pathways.tsv",
        "reaction_components.tsv",
        "enzyme_aliases.tsv",
    ]:
        check(db / fn, f"database {fn}")

    prefix = a.out_prefix
    rpre = prefix + "_rescue"
    ipre = prefix + "_integrated"
    apre = prefix + "_adjudicated"
    log = prefix + "_RUN.log"

    resolved = Path(prefix + "_resolved_exact_annotation.tsv")
    resolver = here / "resolved_annotation_scanner_v0810.py"
    run(
        [
            a.python,
            str(resolver),
            str(gb),
            "--db",
            str(db),
            "--out",
            str(resolved),
        ],
        log,
    )

    rescue_cmd = [
        a.python,
        str(rescue),
        str(gb),
        "--references",
        str(db / "reference_proteins.faa"),
        "--reference-metadata",
        str(db / "reference_metadata.tsv"),
        "--reactions",
        str(db / "reactions.tsv"),
        "--pathways",
        str(db / "pathways.tsv"),
        "--reaction-components",
        str(db / "reaction_components.tsv"),
        "--policy",
        str(policy),
        "--diamond",
        a.diamond,
        "--resolved-reaction-loci",
        str(resolved),
        "--interpro-mode",
        a.interpro_mode,
        "--interpro-max-jobs",
        str(a.interpro_max_jobs),
        "--evalue",
        str(a.evalue),
        "--min-identity",
        str(a.min_identity),
        "--min-qcov",
        str(a.min_qcov),
        "--min-scov",
        str(a.min_scov),
        "--out-prefix",
        rpre,
    ]
    if a.interpro_tsv:
        rescue_cmd += ["--interpro-tsv", a.interpro_tsv]
    if a.interproscan:
        rescue_cmd += ["--interproscan", a.interproscan]
    if a.interpro_applications:
        rescue_cmd += ["--interpro-applications", a.interpro_applications]
    if a.email:
        rescue_cmd += ["--email", a.email]
    if a.enable_weak_rescue:
        rescue_cmd.append("--enable-weak-rescue")
    run(rescue_cmd, log)

    rsum = Path(rpre + "_reaction_evidence_summary.tsv")
    rall = Path(rpre + "_all_candidates.tsv")
    run(
        [
            a.python,
            str(integ),
            str(gb),
            "--db",
            str(db),
            "--reaction-evidence",
            str(rsum),
            "--out-prefix",
            ipre,
        ],
        log,
    )

    ire = Path(ipre + "_integrated_reaction_evidence.tsv")
    run(
        [
            a.python,
            str(adjud),
            str(gb),
            "--db",
            str(db),
            "--integrated-reaction-evidence",
            str(ire),
            "--all-candidates",
            str(rall),
            "--out-prefix",
            apre,
        ],
        log,
    )

    are = Path(apre + "_adjudicated_reaction_evidence.tsv")
    aps = Path(apre + "_adjudicated_pathway_summary.tsv")
    dec = Path(apre + "_cross_reaction_locus_decisions.tsv")
    for p in [are, aps, dec]:
        check(p, "final output")

    rep = report(prefix, aps, are, dec)

    visual = Path(prefix + "_VISUAL_REPORT.html")
    visual_reporter = here / "visual_pathway_reporter_v086.py"
    run(
        [
            a.python,
            str(visual_reporter),
            "--db",
            str(db),
            "--pathway-summary",
            str(aps),
            "--reaction-evidence",
            str(are),
            "--cross-reaction-decisions",
            str(dec),
            "--out-html",
            str(visual),
            "--sample-name",
            Path(a.genbank).stem,
        ],
        log,
    )

    manifest = {
        "runner_version": VERSION,
        "finished": datetime.now().isoformat(timespec="seconds"),
        "genbank": str(gb.resolve()),
        "database": str(db.resolve()),
        "interpro": {
            "mode": a.interpro_mode,
            "precomputed_tsv": a.interpro_tsv or "",
            "local_interproscan": a.interproscan or "",
            "web_email_supplied": bool(a.email),
            "max_jobs": a.interpro_max_jobs,
        },
        "weak_rescue_enabled": bool(a.enable_weak_rescue),
        "outputs": {
            "final_pathway_summary": str(aps),
            "final_reaction_evidence": str(are),
            "cross_reaction_decisions": str(dec),
            "report": str(rep),
            "visual_report": str(visual),
            "log": log,
        },
    }
    with open(prefix + "_RUN_MANIFEST.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print("\n" + "=" * 76)
    print("UNIFIED RUN v0.8.10 COMPLETE")
    print("=" * 76)
    print("Final pathway summary:", aps)
    print("Final reaction evidence:", are)
    print("Cross-reaction decisions:", dec)
    print("Human-readable report:", rep)
    print("Visual pathway report:", visual)

if __name__ == "__main__":
    main()
