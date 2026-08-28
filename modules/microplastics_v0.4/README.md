# Microplastics Optional Branch v0.4 — PET + PCL + PBAT + PU

This is an optional branch. It does not modify the protected main `Bioremediation_DB_v0.7.6_dyes`.

## Modules
- PET/TPA: inherited unchanged from validated PET v0.2
- PCL: inherited from validated v0.3 with strict sequence gate
- PBAT: experimentally grounded TfCut seed; generic cutinase/lipase evidence alone is not PBAT support
- Polyester PU: conservative exact-annotation module (pueA/pueB/polyurethanase). Sequence rescue is intentionally held out until exact reference FASTAs are checksum-verified.

## Safety
Run `runner/run_microplastics_branch.py --baseline-db Bioremediation_DB_v0.7.6_dyes --self-test` first. The baseline guard must report 26 protected files unchanged.
