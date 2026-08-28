#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys, shutil
from pathlib import Path

def run(cmd):
    print("RUN:", " ".join(map(str, cmd)))
    p = subprocess.run(cmd)
    if p.returncode:
        raise SystemExit(p.returncode)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genbank")
    ap.add_argument("--db", default="Bioremediation_DB_v0.7.6_dyes")
    ap.add_argument("--interpro-tsv", required=True)
    ap.add_argument("--diamond", default=None, help="Path/name of DIAMOND executable. If omitted, search PATH then repository root.")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()

    here = Path(__file__).resolve().parent
    repo_root = here.parent.parent
    integrated = repo_root / "core" / "integrated_evidence_v0.6.3" / "integrated_evidence_v063.py"
    diamond = None
    if a.diamond:
        explicit = Path(a.diamond).expanduser()
        if explicit.is_file():
            diamond = str(explicit.resolve())
        else:
            diamond = shutil.which(a.diamond)
        if not diamond:
            raise SystemExit(
                f"ERROR: DIAMOND executable not found: {a.diamond}"
            )
    else:
        diamond = shutil.which("diamond") or shutil.which("diamond.exe")

        if not diamond:
            for candidate in (
                repo_root / "diamond.exe",
                repo_root / "diamond",
            ):
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

    run([
        a.python,
        str(here / "bioremediation_unified_runner_v083_core.py"),
        str(used),
        "--db", a.db,
        "--interpro-tsv", a.interpro_tsv,
        "--diamond", diamond,
        "--integrated-engine", str(integrated),
        "--out-prefix", a.out_prefix,
    ])

    # Final user-facing report:
    # Re-render the already-adjudicated TSV outputs with the validated
    # table-first reporter. This does not rerun or alter biological scoring.
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
