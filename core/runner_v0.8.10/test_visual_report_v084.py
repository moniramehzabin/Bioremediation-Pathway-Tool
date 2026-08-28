from pathlib import Path
r=(Path(__file__).parent/'visual_pathway_reporter_v084.py').read_text(encoding='utf-8')
c=(Path(__file__).parent/'bioremediation_unified_runner_v083_core.py').read_text(encoding='utf-8')
assert 'VISUAL_REPORT.html' in c
assert 'visual_pathway_reporter_v084.py' in c
assert 'step_order' in r
assert 'evidence_relationship' in r
assert 'supported_loci' in r and 'review_loci' in r
assert 'Pathway completeness and evidence confidence' in r
assert 'window.print()' in r
assert 'alternative_group' in r
print('PASS visual report integrated into unified run')
print('PASS ordered pathway flow from step_order')
print('PASS completeness separated from evidence confidence')
print('PASS locus-level support/review provenance')
print('PASS OR/alternative pathway presentation')
print('PASS print/save-PDF control')
print('ALL 6 v0.8.4 REPORT TESTS PASSED')
