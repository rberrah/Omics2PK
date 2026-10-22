"""Synthetic tests for the Zhang same-100-HLM companion source audit."""
import io
import unittest
import zipfile

from scripts.audit_zhang100_companion_pmc4986620 import audit


def fake_xlsx():
    wb=(
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="Data" sheetId="1" r:id="rId1"/></sheets></workbook>'
    )
    rel=(
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Target="worksheets/sheet1.xml" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"/>'
        '</Relationships>'
    )
    rows=[]
    for i in range(1,102):
        rows.append(
            f'<row r="{i}"><c r="A{i}"><v>{i}</v></c>'
            f'<c r="B{i}"><v>{i+1}</v></c><c r="C{i}"><v>{i+2}</v></c></row>'
        )
    sheet=(
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<sheetData>'+''.join(rows)+'</sheetData></worksheet>'
    )
    mem=io.BytesIO()
    with zipfile.ZipFile(mem,"w") as z:
        z.writestr("xl/workbook.xml",wb)
        z.writestr("xl/_rels/workbook.xml.rels",rel)
        z.writestr("xl/worksheets/sheet1.xml",sheet)
    return mem.getvalue()


class ZhangCompanionAuditTests(unittest.TestCase):
    def test_numeric_workbook_is_only_structural_candidate(self):
        blob=fake_xlsx()
        out=audit(fetch=lambda url:(blob,"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"))
        self.assertEqual(out["access"],"OFFICIAL_COMPANION_SUPPLEMENT_INSPECTED_STRUCTURE_ONLY")
        self.assertTrue(out["structural_candidate_for_donor_level_source"])
        self.assertFalse(out["source_rows_acquired"])
        self.assertFalse(out["rna_protein_activity_crosswalk_qualified"])
        self.assertFalse(out["model_fitted"])
        self.assertTrue(out["same_100_HLM_panel_as_primary_Zhang_publication_level"])

    def test_nonarchive_abstains_from_crosswalk(self):
        out=audit(fetch=lambda url:(b"<xml>no supplement archive</xml>","application/xml"))
        self.assertEqual(out["access"],"OFFICIAL_ENDPOINT_NOT_NUMERICAL_ARCHIVE")
        self.assertFalse(out["structural_candidate_for_donor_level_source"])
        self.assertFalse(out["source_rows_acquired"])


if __name__=="__main__":
    unittest.main()
