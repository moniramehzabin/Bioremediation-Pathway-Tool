# Pathway Tool local app

This interface runs the existing Bioremediation Pathway Tool analysis engine on
your computer. It opens in your browser; it is not a publicly hosted service.
The v0.2 app preview adds a separate DyP annotation-seeded family-validation
route. The reference database and DIAMOND thresholds are unchanged. This route
nominates specifically annotated proteins even without a DIAMOND match, then
requires conventional, checksum-matched family evidence and both ordered DyP
terminal domains before supporting a family-level call. It cannot establish
activity against a particular dye. AI-only or incomplete results remain review.

## Start the app

On Windows, double-click **Open_Pathway_Tool.bat**. It first checks native Python
and DIAMOND. If these are unavailable or a native Python module is blocked, it
tries Ubuntu-24.04 and the existing Conda environment named `bpt`.

For a configured Linux or macOS environment, run `python app.py`, or use
`bash Open_Pathway_Tool.sh` from that environment. macOS and the Windows launcher
have not been execution-tested in the development container.

Keep the launcher window open while using the app. Closing the browser tab does
not cancel an analysis. Use **Cancel run** to stop it, or press Ctrl+C in the
launcher to stop the app and its active analysis.

## Run an analysis

1. Choose an annotated GenBank genome with translated CDS features.
2. Choose online InterPro validation and enter your email, or upload an existing
   InterProScan TSV with matching protein identifiers.
3. Click **Run analysis**. Progress and error messages appear in the log panel.
4. When the process finishes, open the report, download the Excel workbook or
   download the complete result archive.

Excel output contains separate sheets for pathways, reactions, candidates,
family validation, cross-reaction review and failed InterPro checks. Tables have
filters and frozen headers. The original TSV outputs remain in the archive.
No extra Python package is required for Excel export.

Skipping InterPro is available explicitly for preliminary exploration only.
Completed engine execution is not a guarantee that every InterPro check
succeeded. The app surfaces the failed-job count when present; review the report,
logs and evidence tables before interpreting results. Moderate or conflicting
assignments must not be treated as experimentally verified functions.

The stable core runs in this first interface. The optional microplastics branch,
advanced threshold controls, arbitrary custom database selection and local
InterProScan installation are not exposed in this first version.

## Setup on a new computer

The interface itself uses Python's standard library. The analysis still requires
Biopython, requests and a compatible DIAMOND executable. It is not yet a bundled
installer: launching the interface does not silently install system software or
change security settings.

For an existing Conda installation on Linux/WSL:

```bash
conda create -n bpt -c conda-forge -c bioconda python=3.12 biopython requests diamond -y
conda activate bpt
python app.py
```

For a native Python installation, install `requirements.txt` and make DIAMOND
available on PATH or place the platform-appropriate executable in the app folder.
Windows Smart App Control can block native compiled packages; the launcher can
use the configured WSL environment without disabling that protection. WSL and
Conda must already be installed for that fallback.

## Files and privacy

Results are saved under `BPT_results` in the home directory of the environment
running the app. For WSL, this is the Linux home directory. Download the ZIP to
save a Windows-accessible copy. Each run uses a new directory and does not
overwrite previous results. Archives include input files, logs and settings;
review them before sharing because the InterPro email can appear in logs.

The interface listens only on `127.0.0.1` with a per-launch private URL. Do not
publish this server or use it as a shared multiuser service. Online InterPro mode
sends policy-selected protein sequences and the supplied email to EMBL-EBI.
Genomes are not uploaded to a separate app-hosting server.

Results remain on disk after restarting, but the interface currently lists only
runs from its active session. Open older downloaded HTML reports directly or
inspect the saved run directory.
