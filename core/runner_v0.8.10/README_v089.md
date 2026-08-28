# Bioremediation Unified Runner v0.8.9 - table-first report

Reporting-layer redesign built on strict v0.8.8 scoring.

Order:
1. Results summary table
2. Reaction evidence table
3. Expandable pathway details
4. Cross-reaction review

Rules:
- Supported = counted.
- Review = not counted.
- Candidate only = not counted.
- NO_CANDIDATE = Not detected.
- Required-core support is separate from optional/supporting functions.
- Arsenate reduction and arsenic efflux are shown as separate capabilities.

Render B. diminuta using the strict v0.8.8 pathway summary:

```cmd
py Bioremediation_Unified_Runner_v0.8.9\render_existing_run_v089.py --db Bioremediation_DB_v0.7.6_dyes --run-prefix B_diminuta_unified_v082 --pathway-summary B_diminuta_unified_v082_pathway_summary_STRICT_v088.tsv --sample-name B_diminuta_S180
```

Then:

```cmd
start B_diminuta_unified_v082_TABLE_REPORT_v089.html
```
