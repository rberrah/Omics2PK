"""Synthetic-only schema audit tests; no public donors or assay values."""
import io
import unittest
import zipfile

from scripts.audit_wegler_source_structure import workbook_structure_only


def fake_xlsx() -> bytes:
    xml = (
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="Hepatocytes" sheetId="1" r:id="rId1"/></sheets>'
        '</workbook>'
    )
    rel = (
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Target="worksheets/sheet1.xml" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"/>'
        '</Relationships>'
    )
    sheet = (
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>DO_NOT_LEAK</t></is></c>'
        '<c r="B1" t="inlineStr"><is><t>DONOR_ID</t></is></c></row>'
        '<row r="2"><c r="A2"><v>234.5</v></c><c r="B2"><v>456.7</v></c></row>'
        '</sheetData></worksheet>'
    )
    memory = io.BytesIO()
    with zipfile.ZipFile(memory, "w") as z:
        z.writestr("xl/workbook.xml", xml)
        z.writestr("xl/_rels/workbook.xml.rels", rel)
        z.writestr("xl/worksheets/sheet1.xml", sheet)
    return memory.getvalue()


class WeglerStructuralTests(unittest.TestCase):
    def test_sheet_counts_do_not_publish_values_or_donor_labels(self):
        data = workbook_structure_only(fake_xlsx(), ("Hepatocytes",))
        self.assertEqual(data["Hepatocytes"]["xml_rows"], 2)
        self.assertEqual(data["Hepatocytes"]["first_four_row_cell_counts"], [2, 2])
        self.assertNotIn("DO_NOT_LEAK", repr(data))
        self.assertNotIn("DONOR_ID", repr(data))
        self.assertNotIn("234.5", repr(data))

    def test_wrong_sheet_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "No required compartments"):
            workbook_structure_only(fake_xlsx(), ("Microsomes",))


if __name__ == "__main__":
    unittest.main()
