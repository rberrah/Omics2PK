"""Fail-closed original Jarmusch source/permissions, never disclose donor-level rows."""
import io
import json
import unittest
from zipfile import ZipFile
from scripts.audit_jarmusch14_original_omics_targeted_pk import inspect_meta,inspect_source


class JarmuschOriginalSourceTests(unittest.TestCase):
    def test_correct_original_article_and_oa_supplement_only(self):
        meta={"pmcid":"PMC7485946","license_code":"CC BY",
              "media_urls":[
                  "s3://pmc-oa-opendata/PMC7485946.1/media/cts12785-sup-0002-TableS2.xlsx",
                  "https://example.org/fake_exfiltration.xlsx"]}
        result,docs=inspect_meta(json.dumps(meta).encode())
        self.assertEqual(result["n_official_original_office_supplements"],1)
        self.assertEqual(len(docs),1)

    def test_different_pmc_is_never_a_valid_pk_join(self):
        with self.assertRaisesRegex(ValueError,"Incorrect"):
            inspect_meta(b'{"pmcid":"PMC999","media_urls":[]}')

    def test_wrong_file_must_not_qualify(self):
        with self.assertRaises(ValueError):
            inspect_source(b"<html>fake</html>","cts12785-sup-0002-TableS2.xlsx")

    def test_metadata_source_lists_are_not_a_patient_join(self):
        with io.BytesIO() as buffer:
            with ZipFile(buffer,"w") as z:
                z.writestr("xl/workbook.xml",'<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheets><sheet name="Sheet1" sheetId="1"/></sheets></workbook>')
                z.writestr("xl/worksheets/sheet1.xml",'<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>ID</t></is></c><c r="B1" t="inlineStr"><is><t>Summary</t></is></c></row><row r="2"><c r="A2" t="inlineStr"><is><t>PERSON_PRIVATE_ABC</t></is></c><c r="B2"><v>7.1</v></c></row></sheetData></worksheet>')
            result=inspect_source(buffer.getvalue(),"table.xlsx")
        self.assertEqual(result["structural_evidence"]["n_sheets"],1)
        self.assertFalse(result["structural_evidence"]["subject_join_is_proven"])
        self.assertNotIn("PERSON_PRIVATE_ABC",repr(result))


if __name__=="__main__":
    unittest.main()
