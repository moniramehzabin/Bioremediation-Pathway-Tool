# Bioremediation Unified Runner v0.8.1

v0.8.1 fixes the main scientific issue exposed by the first unified M. terreus run:

**one broadly conserved locus is no longer automatically allowed to support several competing pollutant-specific reactions.**

New layer:
- scores all reaction assignments sharing the same locus;
- exact GenBank annotation has highest specificity weight;
- multi-source and DIAMOND evidence contribute;
- coherent pathway neighborhood contributes;
- strong winner keeps the locus;
- losing assignments lose that locus;
- close competitions remain `AMBIGUOUS_REVIEW`;
- a reaction is only downgraded when no clean supporting locus remains.

It also fixes pathway alternatives:
`alternative_group` members are now scored as **OR**, not AND.

## Run

From the project root:

```cmd
py Bioremediation_Unified_Runner_v0.8.1\bioremediation_unified_runner_v081.py M_terreus_scaffolds_above500.bgpipe.output_3693362.gb --db Bioremediation_DB_v0.7.6_dyes --interpro-tsv M_terreus_rescue_v053_chromium_interpro_web.tsv --diamond .\diamond.exe --out-prefix M_terreus_unified_v081
```

Then inspect:

```cmd
type M_terreus_unified_v081_REPORT.txt
```

and, for the competition table:

```cmd
powershell -NoProfile -Command "Import-Csv 'M_terreus_unified_v081_adjudicated_cross_reaction_locus_decisions.tsv' -Delimiter \"`t\" | Where-Object {$_.decision -ne 'UNIQUE_KEEP'} | Format-Table -Wrap -AutoSize"
```
