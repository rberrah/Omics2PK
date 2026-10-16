"""Synthetic-only source archive/container QA; real studies never included in CI artifacts."""
import io
import unittest
from unittest.mock import patch
from zipfile import ZipFile,ZIP_DEFLATED

from scripts import audit_public_predose_metabolomics_tac_pk_sources as t


def fake_zip(contents):
    buf=io.BytesIO()
    with ZipFile(buf,"w",ZIP_DEFLATED) as z:
        for name,blob in contents.items():z.writestr(name,blob)
    return buf.getvalue()


class PrimaryOmicsPKStudiesTests(unittest.TestCase):
    def test_original_zip_aggregate_shape_no_raw_patient_values(self):
        worksheet='''<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
         <sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>PERSON_SECRET</t></is></c>
         <c r="B1"><v>5</v></c></row>
         <row r="2"><c r="A2"><v>1</v></c><c r="B2"><v>3</v></c></row></sheetData></worksheet>'''
        xlsx=fake_zip({"xl/worksheets/sheet1.xml":worksheet})
        zipfile=fake_zip({"SUPPLEMENT_original.xlsx":xlsx,"method.txt":"MODEL"})
        result=t._inspect_archive(zipfile)
        self.assertEqual(result["n_contained_spreadsheets"],1)
        self.assertEqual(result["spreadsheet_matrix_shapes_without_original_values"][0]["worksheet_shapes_no_source_values"][0]["n_original_source_rows"],2)
        self.assertNotIn("PERSON_SECRET",repr(result))
        self.assertFalse(result["original_patient_identifiers_or_values_exported"])

    def test_reject_path_traversal(self):
        result=t._inspect_archive(fake_zip({"../clinical_secret.csv":b"STUFF"}))
        self.assertEqual(result["status"],"ORIGINAL_ZIP_UNSAFE_PATH_ABSTAIN")

    def test_reject_authority_based_fabricated_Phapale_pk_crosswalk(self):
        with patch.object(t,"_fetch",return_value=(None,{"status":"SOURCE_HTTP_ACCESS_DENIED_OR_UNAVAILABLE","http_status":403})):
            report=t.run()
        src=report["sources"]["Phapale_2010_high_priority"]
        self.assertFalse(src["individual_omics_x_PK_source_matrix_retrieved"])
        self.assertEqual(src["original_serial_plasma_TAC_horizon_h"],72)
        self.assertEqual(src["original_n_human_volunteers"],29)
        for k in ("Burghelea_2022_Biomedicines","He_2022_Molecules"):
            self.assertFalse(report["sources"][k]["individual_NCA_AUC_Cmax_Tmax_dataset_qualified"])


if __name__=="__main__":unittest.main()
