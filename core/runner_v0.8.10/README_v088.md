# Unified Runner v0.8.8 — strict pathway support scoring

This revision fixes a real scoring inconsistency found in B. diminuta:

- Green / accepted support counts toward `supported_units`.
- Orange / `review_flag=YES` does **not** count toward `supported_units`.
- Purple / candidate does **not** count toward `supported_units`.
- `SUPPORTED_MIXED_LOCI` remains review-level and does not count as supported.
- Review-level units are reported separately as `review_units`.

This makes colors and pathway percentages obey the same rule everywhere.

## Correct an existing B. diminuta pathway summary without rerunning DIAMOND or InterPro

```cmd
py Bioremediation_Unified_Runner_v0.8.8\rescore_existing_pathways_v088.py --db Bioremediation_DB_v0.7.6_dyes --reaction-evidence B_diminuta_unified_v082_adjudicated_adjudicated_reaction_evidence.tsv --old-pathway-summary B_diminuta_unified_v082_adjudicated_adjudicated_pathway_summary.tsv --out B_diminuta_unified_v082_pathway_summary_STRICT_v088.tsv
```

Then compare Styrene and Naphthalene:

```cmd
powershell -NoProfile -Command "Import-Csv 'B_diminuta_unified_v082_pathway_summary_STRICT_v088.tsv' -Delimiter \"`t\" | Where-Object {$_.pathway_id -in @('BPWY_STYRENE_SIDECHAIN','BPWY_NAPHTHALENE_CANONICAL')} | Format-List *"
```
