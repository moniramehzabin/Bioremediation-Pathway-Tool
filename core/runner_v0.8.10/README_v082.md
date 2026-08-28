# Unified Runner v0.8.2 — GenBank input compatibility

The runner now validates GenBank input before analysis.

Design:
- standards-compliant GenBank files (including ordinary NCBI/PGAP output) are used unchanged;
- known non-standard Prokka/SPAdes LOCUS spacing is repaired only in an internal copy;
- the user's original GBK is never edited;
- the repaired copy is parsed again and audited before the scientific pipeline runs.

Run with the same precomputed InterPro workflow:

```cmd
py Bioremediation_Unified_Runner_v0.8.2\bioremediation_unified_runner_v082.py SAMPLE.gbk --db Bioremediation_DB_v0.7.6_dyes --interpro-tsv SAMPLE_interpro.tsv --diamond .\diamond.exe --out-prefix SAMPLE_v082
```

The run writes `<prefix>_GENBANK_INPUT_AUDIT.txt`.

## v0.8.3 evidence-provenance patch
`bioremediation_unified_runner_v083.py` preserves v0.8.2 GenBank normalization and adds locus-aware evidence provenance. `SUPPORTED_MULTI_SOURCE` now requires annotation and countable rescue support on the same locus. If annotation and rescue support occur on different loci, the reaction is labeled `SUPPORTED_MIXED_LOCI` and flagged for review. Final reaction evidence also separates `supported_loci`, `review_loci`, `candidate_only_loci`, and removed/ambiguous loci.
