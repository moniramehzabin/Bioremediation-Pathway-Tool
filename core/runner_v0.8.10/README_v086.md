# Visual Pathway Report v0.8.6

This revision fixes two presentation problems found during B. diminuta review.

- Arsenate reduction (`BRXN_ARS_C`) and arsenic resistance/efflux
  (`BRXN_ARS_EFFLUX`) are displayed as separate functional cards and scored
  from their own reaction evidence.
- Efflux is explicitly not counted as reduction evidence.
- OR pathway alternatives remain grouped rather than looking sequential.
- Candidate-only reactions explicitly say they are not counted as supported.
- Pathway arrows use HTML entities to avoid terminal/browser encoding artifacts.

Render the existing B. diminuta result without rerunning DIAMOND or InterPro:

```cmd
py Bioremediation_Unified_Runner_v0.8.6\render_existing_run_v086.py --db Bioremediation_DB_v0.7.6_dyes --run-prefix B_diminuta_unified_v082 --sample-name B_diminuta_S180
```

Then:

```cmd
start B_diminuta_unified_v082_VISUAL_REPORT.html
```
