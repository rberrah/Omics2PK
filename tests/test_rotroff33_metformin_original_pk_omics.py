"""Synthetic-only original Rotroff n33 public supplement table structural gates."""
import unittest
from io import BytesIO
from zipfile import ZipFile

from scripts.audit_rotroff33_metformin_original_pk_omics import (
    analyze_rows, audit_package, source_xlsx_schema
)


def rows(n=33, *, pk=True, baseline=True, subject=True):
    h={}
    if subject:h[0]="Subject ID"
    if baseline:h[1]="Baseline"
    for j in range(25):h[j+2]=f"metabolite_feature_{j}"
    if pk:
        h[30]="AUC ng h/mL"
        h[31]="Cmax ng/mL"
    table=[h]
    for i in range(n):
        r={}
        if subject:r[0]=f"SYNTHETIC_{i}"
        if baseline:r[1]="predose"
        for j in range(25):r[j+2]=str(1+(i+j)%7)
        if pk:
            r[30]=str(350+i)
            r[31]=str(15+i)
        table.append(r)
    return table


def make_xlsx(source_rows):
    def col_ref(i):
        s=""
        i+=1
        while i:
            i,rem=divmod(i-1,26)
            s=chr(65+rem)+s
        return s
    lines=[]
    for ix,row in enumerate(source_rows,1):
        cells="".join(f'<c r="{col_ref(j)}{ix}" t="inlineStr"><is><t>{val}</t></is></c>'
                      for j,val in row.items())
        lines.append(f'<row r="{ix}">{cells}</row>')
    xml=('<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
         '<sheetData>'+"".join(lines)+'</sheetData></worksheet>')
    buf=BytesIO()
    with ZipFile(buf,"w") as z:z.writestr("xl/worksheets/sheet1.xml",xml)
    return buf.getvalue()


class PublicPKPairingTests(unittest.TestCase):
    def test_true_selfcontained_33_original_individual_matrix(self):
        d=analyze_rows(rows())
        self.assertEqual(d["max_subjects_with_numeric_auc_and_cmax"],33)
        self.assertEqual(d["max_subjects_with_20_plus_features_and_baseline_keyword"],33)
        self.assertTrue(d["self_contained_33_person_baseline_omics_and_pk_table"])
        self.assertNotIn("SYNTHETIC_",str(d))

    def test_no_pk_stops_clinical_claim(self):
        d=analyze_rows(rows(pk=False))
        self.assertFalse(d["self_contained_33_person_baseline_omics_and_pk_table"])

    def test_no_proven_predose_stops_clinical_claim(self):
        d=analyze_rows(rows(baseline=False))
        self.assertFalse(d["self_contained_33_person_baseline_omics_and_pk_table"])

    def test_only_32_participants_stops_clinical_claim(self):
        d=analyze_rows(rows(n=32))
        self.assertFalse(d["self_contained_33_person_baseline_omics_and_pk_table"])

    def test_absent_original_subject_id_stops_clinical_claim(self):
        d=analyze_rows(rows(subject=False))
        self.assertFalse(d["self_contained_33_person_baseline_omics_and_pk_table"])

    def test_original_supplement_wrapper_never_exports_person_info(self):
        buf=BytesIO()
        with ZipFile(buf,"w") as z:
            z.writestr("dataset.xlsx",make_xlsx(rows()))
            z.writestr("Supplement.docx",b"FAKE_DOCX")
        r=audit_package(buf.getvalue())
        self.assertEqual(r["n_original_workbooks"],1)
        self.assertEqual(r["n_source_workbooks_containing_selfcontained_33_patient_E3_table"],1)
        self.assertFalse(r["actual_baseline_omics_pk_patient_crosswalk_qualified"])
        self.assertNotIn("SYNTHETIC_",str(r))


if __name__=="__main__":unittest.main()
