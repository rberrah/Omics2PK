"""Synthetic individual Omics→PK network/signature export tests."""
import unittest
from omics2pk.omic_signature_contract import OmicSignatureError, qualify_individual_signatures

def good():
    return [{
      "subject_id":"SYN-1","sample_id":"syn_sample1","assay_id":"RNA-A1",
      "omic":"transcriptomics","signature_id":"CYP3A_pathway_score",
      "score":0.88,"score_scale":"z_training","fit_scope":"training_frozen_test_projection",
      "reference_model_sha256":"1"*64,"feature_definition_sha256":"2"*64,
      "source_exchange_analysis_id":"SYN-ANALYSIS"
    }]

class ScoresTests(unittest.TestCase):
    def test_one_valid_individual_molecular_score_is_not_PopPK_claim(self):
        q=qualify_individual_signatures(good())
        self.assertEqual(q["n_named_signatures"],1)
        self.assertTrue(q["no_patient_or_score_values_exported"])
        self.assertNotIn("SYN-1",repr(q))
        self.assertEqual(q["status"],"QUALIFIED_INDIVIDUAL_SIGNATURE_SCHEMA_NOT_PREDICTIVE_EFFECT")

    def test_group_only_and_transductive_scores_prohibited(self):
        for scope in ("group_contrast_pvalue","transductive_cohort_network","in_sample_projection"):
            bad=good();bad[0]["fit_scope"]=scope
            with self.assertRaisesRegex(OmicSignatureError,"In-sample"):
                qualify_individual_signatures(bad)

    def test_missing_fit_hash_and_nonfinite_rejected(self):
        for key,value in (("reference_model_sha256","bad"),("score","NaN")):
            bad=good();bad[0][key]=value
            with self.assertRaises(OmicSignatureError):qualify_individual_signatures(bad)

    def test_duplicate_assay_signature_rejected(self):
        with self.assertRaisesRegex(OmicSignatureError,"Duplicate"):
            qualify_individual_signatures(good()+good())

if __name__=="__main__":
    unittest.main()
