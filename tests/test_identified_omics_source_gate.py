"""Qualify discovery metadata without ever inventing donor-level crosswalks."""
from __future__ import annotations
import unittest
from unittest.mock import patch

from scripts.audit_identified_omics_new_public_sources_20261106 import (
    _pride_file_list,_workbench_study,SOURCES,run,
)


class SourceMetadataSafetyContract(unittest.TestCase):
    def test_exact_original_repositories_and_no_fake_pk_result(self):
        self.assertEqual(SOURCES["imatinib_genotype_target_protein"],("PRIDE","PXD066482"))
        self.assertEqual(SOURCES["tucatinib_CYP2C8_CYP3A4_CYP3A5"],("PRIDE","PXD057282"))
        self.assertEqual(SOURCES["pharmmet_liver_s9_drug_metabolite_library"],("WORKBENCH","ST003991"))

    def test_mock_real_shape_pride_index_is_metadata_only(self):
        sample={"_embedded":{"files":[
            {"fileName":"CYP3A4_targeted_protein_quant.txt","fileCategory":"RESULT"},
            {"fileName":"HLM_proteomics.wiff","fileCategory":"RAW"},
        ]}}
        with patch("scripts.audit_identified_omics_new_public_sources_20261106._fetch_json",return_value=sample):
            result=_pride_file_list("PXD057282")
        self.assertEqual(result["status"],"REAL_PUBLIC_PRIDE_FILE_INDEX")
        self.assertEqual(result["file_extension_counts"],{"txt":1,"wiff":1})
        self.assertFalse(result["vetted_source_numeric_donor_PGx_omics_function_join"])
        self.assertTrue(result["cannot_claim_patient_PK_from_raw_mass_spectra"])
        self.assertNotIn("HLM_proteomics.wiff",str(result))

    def test_empty_original_source_api_never_invented(self):
        with patch("scripts.audit_identified_omics_new_public_sources_20261106._fetch_json",return_value={}):
            result=_pride_file_list("PXD066482")
        self.assertEqual(result["status"],"ABSTAIN_PRIDE_NO_QUALIFIED_PUBLIC_FILE_INDEX")

    def test_workbench_does_not_fabricate_subject_pk(self):
        with patch("scripts.audit_identified_omics_new_public_sources_20261106._fetch_json",
                   return_value={"study_title":"PharmMet","project_id":"PR002497"}):
            result=_workbench_study("ST003991")
        self.assertEqual(result["status"],"PUBLIC_ORIGINAL_WORKBENCH_STUDY_SUMMARY")
        self.assertTrue(result["not_clinical_PK"])
        self.assertTrue(result["source_donor_level_response_or_PK_not_verified"])


if __name__=="__main__":
    unittest.main()
