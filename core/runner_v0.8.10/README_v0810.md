# Unified Runner v0.8.10 - skip already-resolved annotation pairs

The runner now scans exact GenBank annotation evidence before Hybrid Rescue.

If a unique exact alias already resolves reaction X at locus A:
- X/A is kept as annotation evidence;
- X/A is skipped during rescue for reaction X;
- locus A can still be evaluated for another reaction Y.

This avoids wasting DIAMOND/InterPro effort rediscovering the same reaction at
the same already-resolved locus, while preserving cross-reaction checking.

New audit:
`<prefix>_resolved_exact_annotation.tsv`

New rescue console line:
`Resolved same-reaction annotation hits skipped before rescue: N`

For B. diminuta, `GBFBEFJJ_01512 / BRXN_ALK_B` should be present in the
resolved annotation TSV and absent from rescue candidates for BRXN_ALK_B.
