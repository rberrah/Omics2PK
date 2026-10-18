"""Synthetic-only tests for exhaustive Rosuvastatin original publisher sheets."""
import unittest
from io import BytesIO
from zipfile import ZipFile, ZIP_DEFLATED

from scripts.audit_rosuvastatin40_all_excel_source_cells import (
    inspect_original_xlsx, inspect_supplement_package,
)


def workbook(n=40, *, hidden_heading=False, omit_cmax=False):
    headers = "AUC" if omit_cmax else "AUC Cmax"
    def c(ref, s):
        return f'<c r="{ref}" t="inlineStr"><is><t>{s}</t></is></c>'
    def ncell(ref, value):
        return f'<c r="{ref}"><v>{value}</v></c>'
    rows = []
    if hidden_heading:
        rows.append('<row r="1">' + c("A1","Dataset note") + "</row>")
    header_row = 2 if hidden_heading else 1
    rows.append(f'<row r="{header_row}">'+c(f"A{header_row}","AUC")+(
        "" if omit_cmax else c(f"B{header_row}","Cmax"))+c(f"C{header_row}","Subject")+"</row>")
    for i in range(n):
        r = i + header_row + 1
        rows.append(f'<row r="{r}">'+ncell(f"A{r}",100+i)+(
            "" if omit_cmax else ncell(f"B{r}",10+i)) + c(f"C{r}","FAKE")+"</row>")
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    xml = f'<worksheet xmlns="{ns}"><sheetData>{"".join(rows)}</sheetData></worksheet>'
    mem=BytesIO()
    with ZipFile(mem,"w",ZIP_DEFLATED) as z:
        z.writestr("xl/worksheets/sheet1.xml",xml)
    return mem.getvalue()


class OriginalWorkbookSourceCellTests(unittest.TestCase):
    def test_40_person_table_with_heading_in_nonfirst_row_detected(self):
        d=inspect_original_xlsx(workbook(hidden_heading=True))
        self.assertTrue(d["worksheet_aggregate_structure"][0]["candidate_40_original_subject_PK_outcome_table"])
        self.assertEqual(d["worksheet_aggregate_structure"][0]["max_both_PK_outcome_numeric_rows_below_same_header"],40)
        self.assertNotIn("FAKE",str(d))

    def test_aggregate_table_does_not_masquerade_as_original_pk_matrix(self):
        d=inspect_original_xlsx(workbook(n=3))
        self.assertFalse(d["worksheet_aggregate_structure"][0]["candidate_40_original_subject_PK_outcome_table"])

    def test_missing_cmax_is_not_validated(self):
        d=inspect_original_xlsx(workbook(omit_cmax=True))
        self.assertFalse(d["worksheet_aggregate_structure"][0]["candidate_40_original_subject_PK_outcome_table"])

    def test_supplement_wrapper_reports_only_hash_and_metadata(self):
        mem=BytesIO()
        with ZipFile(mem,"w") as z:
            z.writestr("Table1.xlsx", workbook())
            z.writestr("Figure.pdf", b"SYNTHETIC")
        d=inspect_supplement_package(mem.getvalue())
        self.assertEqual(d["n_original_excel_files"],1)
        self.assertTrue(d["no_patient_table_or_omics_PK_link_claim"])
        self.assertNotIn("FAKE",str(d))


if __name__=="__main__":
    unittest.main()
