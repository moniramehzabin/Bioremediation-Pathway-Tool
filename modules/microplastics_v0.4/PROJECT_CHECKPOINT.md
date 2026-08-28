# Project checkpoint — Microplastics optional branch v0.3

- Protected main database: Bioremediation_DB_v0.7.6_dyes (34 pathways/modules, 81 reactions, 379 reference proteins; 26 files hash-locked).
- Main S180 validated baseline remains: 9 Supported / 7 Partial / 2 Review-only / 0 Candidate-only / 17 Not detected.
- P450 remains optional and is NOT merged into the main database.
- PET v0.2 remains the rollback checkpoint; S180 PET v0.2 result: all four PET/TPA reactions NO_CANDIDATE.
- v0.3 adds PCL only to the optional microplastics branch.
- PCL direct seed: Q6A0I4, Thermobifida fusca Cut2; experimentally reported PCL depolymerization.
- PCL generic cutinase/lipase/esterase annotation is not accepted as direct support.
- PCL sequence-only support must pass >=60% identity, >=80% query coverage, >=80% reference coverage, High DIAMOND, and non-conflicting policy evidence.
- WP_004373894.1 and WP_003239806.1 are documented as direct PCLases but held out until exact sequences are independently checksum-verified.
- Main v0.7.6 is never modified by this branch.
