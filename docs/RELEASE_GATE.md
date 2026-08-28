# Release gate

## Passed in this candidate

- Stable core database: 26 protected files present.
- Stable core counts: 34 pathway/module IDs, 81 reaction IDs, 379 reference proteins.
- Microplastics reaction set: exactly 7 reaction IDs.
- PET v0.2 regression lock passes.
- PCL strict-gate regression passes.
- Synthetic annotation-positive control resolves all 7 reaction IDs and all 3 TPADO components.
- Microplastics ATPase contamination check passes.
- Personal absolute-path scan passes.
- The 6-vs-7 report discrepancy is explained: 6 pathway/module IDs vs 7 unique reactions.

## Still required before calling the software fully validated

- Owner chooses an open-source license.
- Owner fills author/ORCID in CITATION.cff.
- Add end-to-end biological positive controls for each optional polymer module as verified datasets are selected.
- CI should be added after the public repository structure is finalized.
