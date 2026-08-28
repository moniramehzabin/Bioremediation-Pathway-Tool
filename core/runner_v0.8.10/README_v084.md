# Bioremediation Unified Runner v0.8.4

v0.8.4 adds a presentation-only visual pathway report. Scientific calls are unchanged from v0.8.3.

New output:

`<out-prefix>_VISUAL_REPORT.html`

The report shows pathway/module results first, ordered reaction flows from `step_order`, completeness separately from evidence confidence, genomic context, missing reactions, locus-level provenance, and cross-reaction review. It also has Show all pathways and Print / Save PDF controls.

To visualize an already completed run without rerunning DIAMOND/InterPro:

```cmd
py Bioremediation_Unified_Runner_v0.8.4\render_existing_run_v084.py --db Bioremediation_DB_v0.7.6_dyes --run-prefix B_diminuta_unified_v083 --sample-name B_diminuta_S180
```

Then open:

`B_diminuta_unified_v083_VISUAL_REPORT.html`

For a new analysis, use `bioremediation_unified_runner_v084.py` with the same arguments as v0.8.3. The visual HTML will be generated automatically at the end.
