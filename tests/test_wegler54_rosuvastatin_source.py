"""Synthetic-only source classification: no clinical observations bundled."""
from __future__ import annotations
from io import BytesIO
from zipfile import ZipFile, ZIP_DEFLATED
import unittest

from scripts.audit_wegler54_rosuvastatin_proteomics_pk_source import (
    inspect_xlsx_schema,inspect_supplement_package
)


def workbook(n=54,*,with_pk=True):
    headers=["Donor","OATP1B1","OATP1B3","OATP2B1","NTCP"]
    if with_pk:headers+=["AUC","Cmax"]
    rows=[]
    def cell(ref,val):
        return f'<c r="{ref}" t="inlineStr"><is><t>{val}</t></is></c>'
    rows.append('<row r="1">'+''.join(cell(chr(65+j)+"1",h) for j,h in enumerate(headers))+'</row>')
    for k in range(n):
        ri=k+2
        values=[f"FAKE_PERSON_{k}",3+k/100,2+k/100,1+k/100,0.5+k/100]
        if with_pk:values+=[90+k,10+k/2]
        rows.append(f'<row r="{ri}">'+''.join(cell(chr(65+j)+str(ri),str(v)) for j,v in enumerate(values))+'</row>')
    content=f'<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>{"".join(rows)}</sheetData></worksheet>'
    o=BytesIO()
    with ZipFile(o,"w",ZIP_DEFLATED) as z:z.writestr("xl/worksheets/sheet1.xml",content)
    return o.getvalue()


class WeglerPublicSourceContractTests(unittest.TestCase):
    def test_matching_true_pk_protein_matrix_detected_only_when_paired(self):
        result=inspect_xlsx_schema(workbook())
        schema=result["source_sheets_schema_only"][0]
        self.assertTrue(schema["measured_54_person_auc_cmax_and_proteome_matrix_verified"])
        self.assertEqual(schema["max_numeric_auc_cmax_subject_records_under_joint_headers"],54)
        self.assertTrue(schema["source_cell_values_or_subject_names_exported"] is False)
        self.assertNotIn("FAKE_PERSON",str(result))

    def test_missing_pk_cannot_pass(self):
        schema=inspect_xlsx_schema(workbook(with_pk=False))["source_sheets_schema_only"][0]
        self.assertFalse(schema["measured_54_person_auc_cmax_and_proteome_matrix_verified"])
        self.assertFalse(schema["role_keyword_present_anywhere"]["auc"])

    def test_only_53_person_records_cannot_pass(self):
        schema=inspect_xlsx_schema(workbook(n=53))["source_sheets_schema_only"][0]
        self.assertFalse(schema["measured_54_person_auc_cmax_and_proteome_matrix_verified"])

    def test_publisher_supplement_bundle_schema_only(self):
        b=BytesIO()
        with ZipFile(b,"w") as z:
            z.writestr("SuppinfoS1.xlsx",workbook())
            z.writestr("Supplement.pdf",b"FICTIONAL")
        result=inspect_supplement_package(b.getvalue())
        self.assertEqual(result["public_supplement_xlsx_count"],1)
        self.assertTrue(result["any_original_54_person_auc_cmax_proteomics_table_qualified"])
        self.assertNotIn("FAKE_PERSON",str(result))


if __name__=="__main__":unittest.main()
