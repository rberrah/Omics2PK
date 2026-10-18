"""Public hepatic source gates: all source observations synthetic; no DOI raw data."""
from __future__ import annotations
from io import BytesIO
from zipfile import ZipFile
from unittest import TestCase

from scripts.audit_all_hepatic_rna_protein_original_public_sources import (
    STUDIES, _row_metadata,inspect_package
)


def fabricated_docx_matrix():
    def tc(x):
        return "<w:tc><w:p><w:r><w:t>"+x+"</w:t></w:r></w:p></w:tc>"
    names=["Donor ID","CYP3A4 mRNA","CYP3A4 protein abundance",
           "midazolam 1 hydroxylation activity"]
    doc="<w:document xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\"><w:body><w:tbl>"
    doc+="<w:tr>"+"".join(map(tc,names))+"</w:tr>"
    for i in range(100):
        doc+="<w:tr>"+"".join(map(tc,[
            "SYNTHETIC_PERSON_"+str(i),
            str(10+i/30),str(12+i/45),str(25+i/55),
        ]))+"</w:tr>"
    doc+="</w:tbl></w:body></w:document>"
    output=BytesIO()
    with ZipFile(output,"w") as z:z.writestr("word/document.xml",doc)
    return output.getvalue()


class PublicPrimaryLiverSourceTests(TestCase):
    def test_all_six_studies_are_distinct_true_published_sources(self):
        self.assertEqual(len(STUDIES),6)
        doi=[row["doi"] for row in STUDIES.values()]
        self.assertEqual(len(set(doi)),len(doi))
        self.assertEqual(STUDIES["Bhatt_2018_UGT2B17_n455"]["published_n_activity"],346)
        self.assertEqual(STUDIES["Bhatt_2018_UGT2B17_n455"]["published_n_RNA"],230)
        self.assertFalse(STUDIES["Wenzel_2024_jejunum_liver_n11"]["drug_turnover_assayed"])

    def test_hundred_genuine_rows_look_eligible_but_not_source_certified(self):
        matrix=[
            ["Donor ID","CYP3A4 mRNA","CYP3A4 protein abundance","midazolam hydroxylation activity"]
        ]+[[f"FICTIONAL-{i}",str(i+2.),str(i+3.),str(i+4.)] for i in range(100)]
        out=_row_metadata(matrix)
        self.assertEqual(out["n_rows"],101)
        self.assertTrue(all(out["source_header_role_flags"].values()))
        self.assertEqual(out["n_source_data_rows_with_three_numeric_measurements"],100)
        self.assertNotIn("FICTIONAL",str(out))

    def test_real_original_study_requires_public_data_not_synthetic_table(self):
        b=BytesIO()
        with ZipFile(b,"w") as outer:
            outer.writestr("fake_original_study.docx",fabricated_docx_matrix())
        summary=inspect_package(b.getvalue())
        self.assertTrue(summary["archive_usable"])
        self.assertEqual(summary["documents"][0]["max_source_table_rows"],101)
        self.assertEqual(summary["documents"][0]["n_rough_original_paired_table_candidates"],1)
        self.assertNotIn("SYNTHETIC_PERSON",str(summary))
        # A potential donor/feature numerical row is not a qualified original
        # DOI-and-same-human joint matrix. NO model should be fitted here.

    def test_xml_error_never_recast_as_empty_patient_matrix(self):
        out=inspect_package(b"<?xml version=\"1.0\"?><error>No files</error>")
        self.assertFalse(out["archive_usable"])

    def test_member_path_traversal_refused(self):
        b=BytesIO()
        with ZipFile(b,"w") as z:z.writestr("../malicious.csv",b"bad")
        self.assertEqual(inspect_package(b.getvalue())["reason"],"UNSAFE_SOURCE_MEMBER")


if __name__=="__main__":
    import unittest
    unittest.main()
