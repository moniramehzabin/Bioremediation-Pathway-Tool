#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys, shutil
from pathlib import Path

def run(cmd):
    print("RUN:", " ".join(map(str, cmd)))
    p = subprocess.run(cmd)
    if p.returncode:
        raise SystemExit(p.returncode)

def _append_interpro_args(cmd, a):
    cmd += ["--interpro-mode", a.interpro_mode]
    if a.interpro_tsv:
        cmd += ["--interpro-tsv", a.interpro_tsv]
    if a.interproscan:
        cmd += ["--interproscan", a.interproscan]
    if a.interpro_applications:
        cmd += ["--interpro-applications", a.interpro_applications]
    if a.email:
        cmd += ["--email", a.email]
    cmd += ["--interpro-max-jobs", str(a.interpro_max_jobs)]

    if a.enable_weak_rescue:
        cmd.append("--enable-weak-rescue")
    cmd += [
        "--evalue", str(a.evalue),
        "--min-identity", str(a.min_identity),
        "--min-qcov", str(a.min_qcov),
        "--min-scov", str(a.min_scov),
    ]
    return cmd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genbank")
    ap.add_argument("--db", default=None)
    ap.add_argument(
        "--interpro-mode",
        choices=["auto", "precomputed", "local", "web", "off"],
        default="auto",
        help=(
            "InterPro routing mode. auto uses a supplied --interpro-tsv first, "
            "then --interproscan, then --email for EMBL-EBI web jobs; if none is supplied, the user is prompted for an email for selective web InterProScan."
        ),
    )
    ap.add_argument("--interpro-tsv", help="Precomputed InterProScan TSV.")
    ap.add_argument("--interproscan", help="Path/name of local InterProScan executable/script.")
    ap.add_argument("--interpro-applications", default=None)
    ap.add_argument("--email", help="Email required for EMBL-EBI InterProScan web mode.")
    ap.add_argument("--interpro-max-jobs", type=int, default=200)
    ap.add_argument(
        "--enable-weak-rescue",
        action="store_true",
        help="Allow weak DIAMOND hits to reach InterPro only under the curated weak-rescue policy.",
    )
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--min-identity", type=float, default=25.0)
    ap.add_argument("--min-qcov", type=float, default=50.0)
    ap.add_argument("--min-scov", type=float, default=50.0)
    ap.add_argument(
        "--diamond",
        default=None,
        help="Path/name of DIAMOND executable. If omitted, search PATH then repository root.",
    )
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()

    here = Path(__file__).resolve().parent
    repo_root = here.parent.parent
    if a.db is None:
        a.db = str(repo_root / "core" / "database_v0.7.6")

    integrated = repo_root / "core" / "integrated_evidence_v0.6.3" / "integrated_evidence_v063.py"

    if a.interpro_mode == "precomputed" and not a.interpro_tsv:
        raise SystemExit("ERROR: --interpro-mode precomputed requires --interpro-tsv.")
    if a.interpro_mode == "local" and not a.interproscan:
        raise SystemExit("ERROR: --interpro-mode local requires --interproscan.")
    if a.interpro_mode == "web" and not a.email:
        a.email = input("Enter email for EMBL-EBI InterProScan web jobs: ").strip()
        if not a.email:
            raise SystemExit("ERROR: an email is required for EMBL-EBI InterProScan web mode.")
    if a.interpro_mode == "auto" and not a.interpro_tsv and not a.interproscan and not a.email:
        a.email = input("Enter email for selective EMBL-EBI InterProScan web jobs: ").strip()
        if not a.email:
            raise SystemExit("ERROR: InterPro is required for a normal run. Supply an email, --interpro-tsv, --interproscan, or explicitly choose --interpro-mode off.")

    diamond = None
    if a.diamond:
        explicit = Path(a.diamond).expanduser()
        if explicit.is_file():
            diamond = str(explicit.resolve())
        else:
            diamond = shutil.which(a.diamond)
        if not diamond:
            raise SystemExit(f"ERROR: DIAMOND executable not found: {a.diamond}")
    else:
        diamond = shutil.which("diamond") or shutil.which("diamond.exe")
        if not diamond:
            for candidate in (repo_root / "diamond.exe", repo_root / "diamond"):
                if candidate.is_file():
                    diamond = str(candidate.resolve())
                    break
        if not diamond:
            raise SystemExit(
                "ERROR: DIAMOND was not found. "
                "Install DIAMOND and add it to PATH, "
                "place diamond.exe in the repository root, "
                "or use --diamond <path>."
            )

    print("DIAMOND:", diamond)

    norm = Path(a.out_prefix + "_normalized_input.gbk")
    audit = Path(a.out_prefix + "_GENBANK_INPUT_AUDIT.txt")

    run([
        a.python,
        str(here / "genbank_input_validator_v082.py"),
        a.genbank,
        "--out", str(norm),
        "--audit", str(audit),
    ])

    used = norm if norm.exists() and norm.stat().st_size else Path(a.genbank)

    cmd = [
        a.python,
        str(here / "bioremediation_unified_runner_v083_core.py"),
        str(used),
        "--db", a.db,
        "--diamond", diamond,
        "--integrated-engine", str(integrated),
        "--out-prefix", a.out_prefix,
    ]
    run(_append_interpro_args(cmd, a))

    pathway_summary = Path(a.out_prefix + "_adjudicated_adjudicated_pathway_summary.tsv")
    reaction_evidence = Path(a.out_prefix + "_adjudicated_adjudicated_reaction_evidence.tsv")
    cross_decisions = Path(a.out_prefix + "_adjudicated_cross_reaction_locus_decisions.tsv")
    visual_report = Path(a.out_prefix + "_VISUAL_REPORT.html")

    run([
        a.python,
        str(here / "table_first_reporter_v089.py"),
        "--db", a.db,
        "--pathway-summary", str(pathway_summary),
        "--reaction-evidence", str(reaction_evidence),
        "--cross-reaction-decisions", str(cross_decisions),
        "--out-html", str(visual_report),
        "--sample-name", Path(a.genbank).stem,
    ])

if __name__ == "__main__":
    main()
