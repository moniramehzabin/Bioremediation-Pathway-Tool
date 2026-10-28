import hashlib
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / 'core/runner_v0.8.10/Bioremediation_Hybrid_Rescue_v0.3.5_resolvedskip/family_seed_evidence.py'
spec = importlib.util.spec_from_file_location('family_seed', MODULE)
family = importlib.util.module_from_spec(spec); spec.loader.exec_module(family)


class FamilyValidationTests(unittest.TestCase):
    def evaluate(self, rows, sequence='A'*468):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'domains.tsv'
            path.write_text('\n'.join('\t'.join(row) for row in rows))
            return family.evaluate_dyp(path, {'test': sequence}, ['test'])[0]

    def rows(self):
        digest = hashlib.md5(('A'*468).encode()).hexdigest()
        return [['test',digest,'468',db,sig,'domain',str(start),str(end),score,db,
                 'date',ipr,'description'] for db,sig,start,end,score,ipr in [
            ('NCBIFAM','TIGR01413',117,456,'5.5E-64','IPR006314'),
            ('Pfam','PF04261',117,259,'2.1E-24','IPR048327'),
            ('Pfam','PF20628',271,455,'1.2E-28','IPR048328')]]

    def test_conventional_family_and_ordered_domains_support(self):
        self.assertTrue(self.evaluate(self.rows())['countable_family'])

    def test_single_domain_or_family_only_stays_review(self):
        self.assertFalse(self.evaluate(self.rows()[:2])['countable_family'])
        self.assertFalse(self.evaluate(self.rows()[1:])['countable_family'])
        self.assertFalse(self.evaluate([])['countable_family'])

    def test_ai_predictions_cannot_create_support(self):
        rows=self.rows()
        for row in rows:row[9]='InterPro-N'
        self.assertFalse(self.evaluate(rows)['countable_family'])

    def test_wrong_sequence_results_cannot_create_support(self):
        self.assertFalse(self.evaluate(self.rows(), sequence='G'*468)['countable_family'])

    def test_weak_or_disordered_domains_stay_review(self):
        rows=self.rows();rows[1][8]='0.1'
        self.assertFalse(self.evaluate(rows)['countable_family'])
        rows=self.rows();rows[1][6:8]=['300','440']
        self.assertFalse(self.evaluate(rows)['countable_family'])

    def test_general_peroxidase_annotation_not_nominated(self):
        proteins={'a':{'product':'peroxidase'},'b':{'product':'Dyp-type peroxidase'},
                  'c':{'product':'multicopper oxidase'}}
        self.assertEqual(family.nominated_dyp(proteins, {'BRXN_DYP':{}}), ['b'])
        self.assertEqual(family.nominated_dyp(proteins, {'BRXN_DYP':{}}, {('BRXN_DYP','b')}), [])

    def test_provenance_does_not_fabricate_diamond_scores(self):
        result=self.evaluate(self.rows())
        row=family.candidate_row(result, {'test':{}}, {'BRXN_DYP':{}})
        self.assertEqual(row['identity_pct'], '')
        self.assertEqual(row['reference_id'], '')
        self.assertNotIn('DIAMOND',row['evidence_source'])


if __name__=='__main__':unittest.main()
