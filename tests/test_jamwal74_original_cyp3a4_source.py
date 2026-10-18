"""Synthetic-only original n74 mRNA/protein/activity release-gate tests."""
from __future__ import annotations
from io import BytesIO
from zipfile import ZipFile
import unittest

from scripts.audit_jamwal74_original_cyp3a4_rna_protein_activity import (
    _row_audit, inspect_original_docx,audit_official_original_supplement
)

def fake_docx(n=74,three=True,has_ids=True):
    cols=["Donor ID","CYP3A4 mRNA","CYP3A4 protein","midazolam hydroxylation activity"] if three else [
        "Donor ID","CYP3A4 mRNA","CYP3A4 protein"]
    if not has_ids:cols[0]="Group"
    def row(cells):
        return '<w:tr>'+''.join('<w:tc><w:p><w:r><w:t>'+str(v)+'</w:t></w:r></w:p></w:tc>' for v in cells)+'</w:tr>'
    rs=[row(cols)]
    for i in range(n):
        values=[f"FICTIONAL_{i:03d}",str(4+i/10),str(5+i/10)]
        if three:values.append(str(1+i/11))
        rs.append(row(values))
    doc='<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body><w:tbl>'+''.join(rs)+'</w:tbl></w:body></w:document>'
    f=BytesIO()
    with ZipFile(f,"w") as z:z.writestr("word/document.xml",doc)
    return f.getvalue()

class OriginalSourceContractTests(unittest.TestCase):
    def test_positive_original_format_structure_only(self):
        report=inspect_original_docx(fake_docx())
        self.assertEqual(report["n_original_docx_tables"],1)
        self.assertEqual(report["n_strict_source_74person_numeric_RNA_protein_activity_tables"],1)
        self.assertNotIn("FICTIONAL",str(report))

    def test_missing_mdz_denies_source_matrix(self):
        report=inspect_original_docx(fake_docx(three=False))
        self.assertEqual(report["n_strict_source_74person_numeric_RNA_protein_activity_tables"],0)

    def test_wrong_donor_count_denies_source_matrix(self):
        report=inspect_original_docx(fake_docx(n=24))
        self.assertEqual(report["n_strict_source_74person_numeric_RNA_protein_activity_tables"],0)

    def test_no_source_subject_key_denies(self):
        report=inspect_original_docx(fake_docx(has_ids=False))
        self.assertEqual(report["n_strict_source_74person_numeric_RNA_protein_activity_tables"],0)

    def test_original_package_never_qualifies_source_without_actual_row_semantics(self):
        b=BytesIO()
        with ZipFile(b,"w") as z:
            z.writestr("supplement.docx",fake_docx())
            z.writestr("figure.png",b"FICTIONAL")
        report=audit_official_original_supplement(b.getvalue())
        self.assertEqual(report["n_candidate_original_numeric_74person_RNA_protein_activity_tables"],1)
        self.assertFalse(report["actual_original_74donor_RNA_vs_multimodal_clearance_comparison_justified"])
        self.assertFalse(report["verified_individual_74donor_RNA_protein_activity_values_available"])
        self.assertEqual(report["original_supplement_member_format_counts"][".docx"],1)

if __name__=="__main__":unittest.main()
