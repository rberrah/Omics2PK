"""Never promote a published same-human design into numeric E3 without a join."""
import json
import unittest
from scripts.audit_metformin_17_pharmacometabolomic_public import (
    PMCID, source_metadata, supplement_structure,
)


class MetforminClinicalProofTests(unittest.TestCase):
    def test_exact_original_article_and_no_numeric_omic_claim(self):
        source={"pmcid":PMCID,"version":1,"license_code":"CC BY",
                "media_urls":["s3://pmc-oa-opendata/PMC9231303.1/media/supplement.docx"]}
        result=source_metadata(json.dumps(source).encode())
        self.assertEqual(result["pmcid"],PMCID)
        self.assertEqual(result["n_declared_supplementary_media"],1)

    def test_incorrect_article_denied(self):
        with self.assertRaisesRegex(ValueError,"Incorrect original"):
            source_metadata(b'{"pmcid":"PMC999","media_urls":[]}')

    def test_structured_html_or_graph_does_not_count_as_individual_omics(self):
        result=supplement_structure(b"feature-level plots", "supplement.pdf")
        self.assertFalse(result["has_person_level_omics_matrix_verified"])
        self.assertFalse(result["has_original_numeric_patient_PK_verified"])

    def test_bogus_external_source_excluded(self):
        source={"pmcid":PMCID,"media_urls":[
            "http://example.org/unsafe.xlsx",
            "s3://pmc-oa-opendata/PMC9231303.1/media/tableS5.xlsx"]}
        result=source_metadata(json.dumps(source).encode())
        self.assertEqual(result["n_declared_supplementary_media"],1)


if __name__=="__main__":
    unittest.main()
