# Changelog

## Release candidate — 2026-08-28

- Preserved stable core database v0.7.6 (34 pathway/module IDs, 81 reaction IDs, 379 references, 26 protected files).
- Repaired core runner reference from missing visual reporter v0.8.4 to bundled v0.8.6.
- Added isolated microplastics v0.4 branch: PET, PCL, PBAT and polyester-PU.
- Added strict PCL/PBAT evidence gating and PET regression lock.
- Added synthetic annotation positive control for all seven microplastics reaction IDs.
- Added release audit and personal-path scan.
- Kept P450 outside the stable core by design.

## Release-candidate component propagation fix
- Fixed a code-only metadata propagation defect for existing multi-component reactions.
- Older curated references that already carry an exact component target/gene identity are now mapped to the matching `reaction_components.tsv` component when and only when the mapping is unique.
- No reaction definitions, thresholds, reference sequences, or protected v0.7.6 database files were changed.
- Added regression test covering PcaI -> A and PcaJ -> B.

## v0.1.0 - 2026-08-28
- Finalized the first public Bioremediation Pathway Tool package.
- Preserved the validated core database and biological scoring rules.
- Retained the component-aware metadata propagation fix for multi-component reactions.
- Corrected table-first reporter handling so adjudicated supported calls are not demoted solely because of mixed-locus provenance.
- Integrated the validated table-first v0.8.9 report into the unified runner.
- Completed end-to-end M. terreus regression testing and final release audit.
