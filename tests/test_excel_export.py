import importlib.util
from pathlib import Path
import tempfile
import unittest
import zipfile
from xml.etree import ElementTree as ET

spec=importlib.util.spec_from_file_location('excel_export',Path(__file__).resolve().parents[1]/'excel_export.py')
export=importlib.util.module_from_spec(spec);spec.loader.exec_module(export)


class ExcelTests(unittest.TestCase):
    def test_numeric_values_identifiers_and_formula_safety(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'analysis_rescue_all_candidates.tsv'
            p.write_text('locus_tag\tevalue\tidentity_pct\tproduct\nAC6666_06575\t5.5E-64\t32.4\t=1+1\n')
            workbook=export.export_results(folder)
            with zipfile.ZipFile(workbook) as z:
                for name in z.namelist():
                    if name.endswith('.xml') or name.endswith('.rels'):ET.fromstring(z.read(name))
                text=z.read('xl/worksheets/sheet2.xml').decode()
                self.assertIn('<v>5.5E-64</v>',text)
                self.assertIn('<v>32.4</v>',text)
                self.assertNotIn('<f>',text)
                self.assertIn('=1+1',text)
                self.assertIn('autoFilter',text)
                self.assertIn('state="frozen"',text)


if __name__=='__main__':unittest.main()
