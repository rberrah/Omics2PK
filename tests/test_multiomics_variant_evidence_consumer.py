"""Pure synthetic optional sidecar interoperability, no human plasma PK."""
import unittest
from omics2pk.multiomics_variant_evidence_consumer import (
    qualify_variant_evidence_sidecar,OmicsVariantEvidenceHandoffError,
)
def example():
    return {
        "schema_version":"variant_omic_evidence_sidecar_v1",
        "status":"QUALIFIED_ORIGINAL_ASSAY_DESCRIPTIVE_HYPOTHESES_NOT_CAUSAL_PROOF",
        "gene":"CYP3A5","genomic_consequence":"missense",
        "genotype_coverage_state":"original_variant_called",
        "original_per_assay_signals_measured_not_inferred":{
            "total_RNA":"reference_like","protein_full_length":"decreased"},
        "assays_share_source_biological_unit_tissue_and_time":True,
        "hypotheses_to_test":[
            "HYPOTHESIS_POSTTRANSCRIPTIONAL_TRANSLATION_OR_PROTEIN_INSTABILITY"],
        "causal_variant_mechanism_proven_by_this_report":False,
        "protein_decay_observed_directly":False,
        "in_vivo_drug_PK_or_PD_improvement_estimated":False,
        "can_export_as_optional_separate_annotation_not_exchange_V1_schema_change":True,
    }
class MultiOmicsHypothesisHandoffTest(unittest.TestCase):
    def test_accept_optional_sidecar_as_only_prospective_PopPK_hypothesis(self):
        v=qualify_variant_evidence_sidecar(example(),drug="tacrolimus",endpoint="AUC")
        self.assertEqual(v["status"],
                         "QUALIFIED_HYPOTHESIS_PRIOR_ONLY_PATIENT_PK_EFFECT_NOT_ESTABLISHED")
        self.assertEqual(v["mechanism_gene"],"CYP3A5")
        self.assertEqual(len(v["observed_omic_layers_not_fitted_as_PK_predictors"]),2)
        self.assertFalse(v["eligible_as_patient_level_PK_predictor"])
        self.assertFalse(v["any_superiority_vs_PGx_only_proven"])
        self.assertFalse(v["original_PK_PGx_omic_linkage_separately_verified"])
        self.assertTrue(v["original_MultiOmics_DOWNSTREAM_EXCHANGE_V1_stays_unchanged"])

    def test_genuine_pk_original_join_cannot_turn_sidecar_into_PK_predictor(self):
        v=qualify_variant_evidence_sidecar(
            example(),drug="tacrolimus",endpoint="AUC",
            orig_PGx_pk_join_verified=True)
        self.assertTrue(v["original_PK_PGx_omic_linkage_separately_verified"])
        self.assertFalse(v["eligible_as_patient_level_PK_predictor"])
        self.assertFalse(v["any_superiority_vs_PGx_only_proven"])
        self.assertTrue(v["need_original_per_subject_timed_omic_matrices"])

    def test_forbid_causal_variant_or_measured_decay_fabrication(self):
        for field in ("causal_variant_mechanism_proven_by_this_report",
                      "protein_decay_observed_directly",
                      "in_vivo_drug_PK_or_PD_improvement_estimated"):
            with self.subTest(field=field),self.assertRaises(OmicsVariantEvidenceHandoffError):
                qualify_variant_evidence_sidecar(
                    {**example(),field:True},drug="tacrolimus",endpoint="AUC")

    def test_invalid_PGx_coding_organ_hypothesis_or_readout_abstains(self):
        for change in (
            {"schema_version":"DOWNSTREAM_EXCHANGE_V1"},
            {"gene":"SLC22A2"},
            {"hypotheses_to_test":["PROVEN_CAUSAL_CYP3A5_CLEARANCE_LOSS"]},
            {"original_per_assay_signals_measured_not_inferred":{"nonexistent_RNA":"reference_like"}},
            {"genotype_coverage_state":"fully_reference_normal_function"},
            {"assays_share_source_biological_unit_tissue_and_time":False},
            {"status":"ABSTAIN_UNMATCHED_SUBJECT_TISSUE_OR_SAMPLE_LANDMARK"},
        ):
            with self.subTest(change=change),self.assertRaises(OmicsVariantEvidenceHandoffError):
                qualify_variant_evidence_sidecar({**example(),**change},
                                                 drug="tacrolimus",endpoint="AUC")
        with self.assertRaises(OmicsVariantEvidenceHandoffError):
            qualify_variant_evidence_sidecar(example(),drug="metformin",endpoint="CLr")

if __name__=="__main__":
    unittest.main()
