"""Synthetic original-source-like model comparisons. No human data or PK claims."""
import unittest
from omics2pk.omic_incremental_endpoint_audit import (
    audit_omic_incremental_PK_prediction,OmicIncrementalEvidenceError,
)


def synthetic(n=50):
    observations=[]
    for i in range(n):
        actual=4.0+(i%11)*0.5
        context=("known_reduced_function_variant" if i%2 else
                 "no_reported_variant_in_assayed_sites")
        observations.append({
            "subject_id":f"SYNTHETIC_SOURCE_{i:03d}",
            "fold_id":str(i%5),
            "genotype_stratum":context,
            "observed":actual,
            "predictions":{
                "clinical_plus_genotype":actual*1.5,
                "add_predose_RNA":actual*1.11,
                "add_predose_protein":actual*0.97,
                "add_predose_activity":actual*1.01,
            },
        })
    return observations


def assess(rows=None,**kwargs):
    defaults=dict(
        records=synthetic() if rows is None else rows,
        endpoint="human_plasma_AUC",
        original_units="ng h/mL",
        genotype_qualifier="source_complete_diplotype",
        evaluation_regime="true_out_of_fold_per_biological_subject",
        predeclared_primary_metric="mean_absolute_log_error",
        source_omic_sample_prelandmark_verified=True,
        source_prediction_training_is_subject_disjoint=True,
        original_source_omics_and_PK_same_subject_verified=True,
    )
    defaults.update(kwargs)
    return audit_omic_incremental_PK_prediction(**defaults)


class OmicAddedValueResearchTest(unittest.TestCase):
    def test_four_predictor_arms_single_PK_endpoint_many_metrics(self):
        q=assess()
        self.assertEqual(q["n_original_unique_biological_subjects"],50)
        self.assertEqual(q["n_original_heldout_groups"],5)
        self.assertEqual(q["endpoint_class"],"CLINICAL_PK")
        self.assertEqual(q["predeclared_primary_loss"],"mean_absolute_log_error")
        a=q["heldout_descriptives_by_predeclared_model_arm"]
        self.assertEqual(len(a),4)
        self.assertLess(a["add_predose_protein"]["mean_absolute_log_error"],
                        a["clinical_plus_genotype"]["mean_absolute_log_error"])
        self.assertGreater(a["add_predose_protein"]["P20"],
                           a["clinical_plus_genotype"]["P20"])
        self.assertTrue(q["baseline_clinical_plus_genotype"])
        self.assertFalse(q["selected_best_omic_claimed"])
        self.assertFalse(q["nested_hyperparameter_or_modality_selection_independently_verified"])
        self.assertFalse(q["clinical_dose_recommendation_or_improved_patient_PK_proven"])
        self.assertNotIn("SYNTHETIC_SOURCE_000",str(q))

    def test_bootstraps_strictly_paired_and_deterministic(self):
        q=assess()
        self.assertEqual(q,assess())
        contrasts=q["added_omic_minus_genotype_primary_loss_paired_bootstrap"]
        self.assertEqual(set(contrasts),{
            "add_predose_RNA","add_predose_protein","add_predose_activity"})
        x=contrasts["add_predose_protein"]
        self.assertLess(x["delta_primary_loss_negative_favors_added_omics"],0)
        self.assertFalse(x["interval_includes_zero"])
        self.assertFalse(x["bootstrap_refitted_any_model"])
        self.assertFalse(x["bootstrap_accounts_for_omic_selection"])

    def test_known_variant_strata_require_minimum_independent_subjects(self):
        q=assess()
        states=q["within_variant_function_groups"]
        self.assertEqual(len(states),2)
        self.assertEqual(states["known_reduced_function_variant"]["n_observed_source_subjects"],25)
        self.assertEqual(states["known_reduced_function_variant"]["status"],
                         "DESCRIPTIVE_WITHIN_KNOWN_VARIANT_STRATUM_NOT_CAUSAL")
        short=assess(synthetic(15))["within_variant_function_groups"]
        self.assertEqual(short["known_reduced_function_variant"]["status"],
                         "ABSTAIN_TOO_FEW_INDEPENDENT_VARIANT_STRATUM_SUBJECTS")
        self.assertTrue(short["known_reduced_function_variant"]["no_omic_winner_inferred"])

    def test_mae_log_rmse_p20_primary_metric_separately_declared(self):
        for metric in ("MAE","RMSE","mean_absolute_log_error","P20"):
            with self.subTest(metric=metric):
                x=assess(predeclared_primary_metric=metric)
                self.assertEqual(x["predeclared_primary_loss"],metric)
                self.assertEqual(
                    x["added_omic_minus_genotype_primary_loss_paired_bootstrap"]
                    ["add_predose_protein"]["n_bootstrap"],499)
                self.assertLess(
                    x["added_omic_minus_genotype_primary_loss_paired_bootstrap"]
                    ["add_predose_protein"]["delta_primary_loss_negative_favors_added_omics"],0)

    def test_no_original_pgx_forces_abstention(self):
        with self.assertRaisesRegex(OmicIncrementalEvidenceError,"no original genotype"):
            assess(genotype_qualifier="not_assessed")
        with self.assertRaisesRegex(OmicIncrementalEvidenceError,"Endpoint"):
            assess(endpoint="unknown_tacrolimus_PD")
        with self.assertRaisesRegex(OmicIncrementalEvidenceError,"metric"):
            assess(predeclared_primary_metric="MAE_wins_by_definition")

    def test_no_postdose_omics_as_predose_or_late_join(self):
        for key in (
            "source_omic_sample_prelandmark_verified",
            "source_prediction_training_is_subject_disjoint",
            "original_source_omics_and_PK_same_subject_verified",
        ):
            with self.subTest(key=key),self.assertRaisesRegex(
                OmicIncrementalEvidenceError,"ABSTAIN"):
                assess(**{key:False})

    def test_all_original_modalities_must_have_same_biological_subject_fold_and_outcome(self):
        rows=synthetic()
        rows[1]["subject_id"]=rows[0]["subject_id"]
        with self.assertRaisesRegex(OmicIncrementalEvidenceError,"Duplicate"):
            assess(rows)
        rows=synthetic()
        rows[1]["predictions"].pop("add_predose_RNA")
        with self.assertRaisesRegex(OmicIncrementalEvidenceError,"exactly the same subject"):
            assess(rows)
        rows=synthetic()
        rows[1]["observed"]=0.0
        with self.assertRaisesRegex(OmicIncrementalEvidenceError,"positive and finite"):
            assess(rows)
        rows=synthetic()
        rows[1]["predictions"]["add_predose_RNA"]=-10
        with self.assertRaisesRegex(OmicIncrementalEvidenceError,"positive and finite"):
            assess(rows)

    def test_n13_PHh_is_invitro_and_never_clinical_PopPK(self):
        q=assess(synthetic(13),
                 endpoint="PHH_imatinib_invitro_clearance",
                 genotype_qualifier="one_measured_SNP",
                 original_units="microlitre/million cells/min")
        self.assertEqual(q["endpoint_class"],"IN_VITRO_DRUG_METABOLISM_NOT_CLINICAL_PK")
        self.assertFalse(q["geno_reference_was_complete_diplotype"])
        self.assertFalse(q["selected_best_omic_claimed"])
        self.assertTrue(all(
            value["status"]=="ABSTAIN_TOO_FEW_INDEPENDENT_VARIANT_STRATUM_SUBJECTS"
            for value in q["within_variant_function_groups"].values()))


if __name__=="__main__":
    unittest.main()
