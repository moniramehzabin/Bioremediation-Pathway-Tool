# Integrated Evidence v0.6.3

Combines three layers:

1. **Specific GenBank annotation**
2. **Hybrid Rescue evidence**
3. **Genomic context**

Safety rules:

- Context never creates a reaction.
- Only an exact annotation alias that maps uniquely to one active reaction may
  create `SUPPORTED_ANNOTATION`.
- Generic or ambiguous annotation aliases do not become exact support.
- Rescue evidence and annotation provenance remain visible separately.
- Inactive legacy reaction IDs are ignored.

For the cleaned HPA model this should allow an exact `hpaD` annotation to support
canonical `BRXN_HPA_D`, while still leaving genuinely missing HPA_BC/F/G absent.

Run:

```cmd
py Bioremediation_Integrated_Evidence_v0.6.3\integrated_evidence_v063.py M_terreus_scaffolds_above500.bgpipe.output_3693362.gb --db Bioremediation_DB_v0.6.1_hpa_cleanup --reaction-evidence M_terreus_hybrid_v034_hpa_clean_reaction_evidence_summary.tsv --out-prefix M_terreus_integrated_v063
```
