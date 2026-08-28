# Bioremediation Pathway Tool v0.1.0

A bacterial-genome evidence-mining workflow for curated bioremediation reactions/pathways, with an isolated optional microplastics branch.

This repository is **Bioremediation Pathway Tool**, a separate project from **Bioremediation Gene Miner**. The two projects should not be merged or treated as versions of one another.

## v0.1.0 scope

- Stable core database v0.7.6: **34 pathway/module IDs, 81 reaction IDs, 379 reference proteins**.
- Core runner v0.8.10, with the missing reporter reference repaired in this release candidate (`visual_pathway_reporter_v086.py`).
- Optional microplastics v0.4: PET, PCL, PBAT and polyester-PU evidence models (**7 exact reaction IDs**).
- P450 is deliberately **not** merged into the stable core; it remains a separate development track and is not included in v0.1.0.

## Evidence interpretation

A computational hit is evidence for a curated reaction/function, not proof that an organism degrades a polymer or pollutant in vivo. Experimental validation remains required. Generic esterase/lipase/cutinase annotations are not automatically accepted as PET/PCL/PBAT/PU degradation evidence.

`NO_CANDIDATE` means no candidate met the current branch's detection rules; it does not prove biological incapacity.

## Why a report may show 6 summary rows but 7 microplastics reactions

The table reporter summarizes **pathway/module IDs**, not unique reaction IDs. Microplastics v0.4 contains 6 pathway/module IDs but 7 unique reaction IDs. PET reactions are reused across multiple PET pathway/module definitions. The exact-reaction TSV is the authoritative seven-reaction view.

## Requirements

- Python 3.10+
- Biopython
- requests
- DIAMOND executable available locally for sequence-rescue runs
- Precomputed InterPro TSV for normal integrated runs

Install Python dependencies:

```bash
python -m pip install -r requirements.txt
```

## Release audit

```bash
python tests/release_audit.py
```

The audit checks the stable-core counts, microplastics reaction set, PET regression lock, PCL strict gate, synthetic positive annotation control, ATPase contamination, and personal absolute paths.

## Microplastics self-test

Windows CMD, from the repository root:

```cmd
python "modules\microplastics_v0.4\runner\run_microplastics_branch.py" --baseline-db "core\database_v0.7.6" --self-test
```

## Microplastics genome run

```cmd
python "modules\microplastics_v0.4\runner\run_microplastics_branch.py" "sample.gbk" --baseline-db "core\database_v0.7.6" --interpro-tsv "sample_interpro.tsv" --diamond "diamond.exe" --out-prefix "sample_microplastics"
```

## Stable-core run

```cmd
python "core\runner_v0.8.10\bioremediation_unified_runner_v0810.py" "sample.gbk" --db "core\database_v0.7.6" --interpro-tsv "sample_interpro.tsv" --diamond "diamond.exe" --out-prefix "sample_core"
```

## Synthetic positive control

`examples/microplastics_annotation_positive_control.gbk` is deliberately synthetic and tests **annotation resolution only**. It is not biological sequence evidence and is not a substitute for end-to-end positive genomes.

## Citation`r`n`r`nBioremediation Pathway Tool v0.1.0 is permanently archived on Zenodo.`r`n`r`n**DOI:** 10.5281/zenodo.22148201`r`n`r`n## Release information

Bioremediation Pathway Tool v0.1.0 is released under the MIT License. Software authors are Monira Mehzabin and Khandoker Md Rezwan; citation metadata is provided in `CITATION.cff`. End-to-end positive biological controls can be expanded as the validation set grows.
