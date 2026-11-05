"""Synthetic heldout predictive accuracy, bias, ranks and PK endpoint guards."""
import unittest
try:
    import numpy as np
    from omics2pk.predictive_diagnostics import (
        DiagnosticsError, summarize_heldout, paired_heldout_contrast,
        multiplicative_pk_diagnostics,
    )
except ImportError:
    np=None


@unittest.skipUnless(np is not None, "NumPy optional scientific extra not installed")
class HeldoutDiagnosticsTests(unittest.TestCase):
    def test_mae_and_rmse_have_different_loss_preferences(self):
        y=np.linspace(0.1,3.0,30)
        regular=np.vstack([y+.5,y+.5,y+.5])
        outlier=np.vstack([y+.1,y+.1,y+.1])
        outlier[:,-1]+=5.9
        a=summarize_heldout(y,regular)
        b=summarize_heldout(y,outlier)
        self.assertLess(b["MAE_original_source_scale"],a["MAE_original_source_scale"])
        self.assertGreater(b["RMSE_original_source_scale"],a["RMSE_original_source_scale"])
        self.assertEqual(a["n_independent_source_cell_lines"],30)
        self.assertEqual(a["n_repeated_predictions_not_independent"],90)
        self.assertAlmostEqual(a["mean_signed_error_prediction_minus_observed"],.5)
        self.assertAlmostEqual(a["Spearman_per_unique_cell_line_repeat_mean_prediction"],1.0)
        self.assertTrue(a["source_has_not_proven_physical_GI50_unit_or_log_transform"])

    def test_baseline_overshoot_and_exploratory_paired_resampling(self):
        y=np.linspace(1,10,30)
        good=np.stack([y+.2]*3)
        bad=np.stack([y+.7]*3)
        contrast=paired_heldout_contrast(y,good,bad)
        self.assertAlmostEqual(contrast["delta_MAE_candidate_minus_reference_negative_is_better"],-.5)
        self.assertEqual(contrast["n_independent_line_units_resampled"],30)
        self.assertEqual(contrast["n_bootstrap_draws"],1999)
        self.assertTrue(contrast["repeat_errors_averaged_within_cell_line_before_resampling"])
        self.assertEqual(contrast,paired_heldout_contrast(y,good,bad))

    def test_no_fake_PK_metrics_on_negative_GI50_or_unverified_endpoint(self):
        y=np.linspace(1,10,30)
        p=np.stack([y*1.1,y*1.1])
        with self.assertRaises(DiagnosticsError):
            multiplicative_pk_diagnostics(y,p,positive_endpoint="AUC")
        with self.assertRaises(DiagnosticsError):
            multiplicative_pk_diagnostics(y,p,positive_endpoint="GI50",original_units_verified=True)
        with self.assertRaises(DiagnosticsError):
            multiplicative_pk_diagnostics(-y,-p,positive_endpoint="AUC",original_units_verified=True)
        with self.assertRaises(DiagnosticsError):
            summarize_heldout(y,p,endpoint_kind="clinical_human_AUC")
        valid=multiplicative_pk_diagnostics(y,p,positive_endpoint="AUC",
                                            original_units_verified=True,
                                            clinically_relevant_ratio=1.2)
        self.assertAlmostEqual(valid["geometric_mean_fold_error_pred_over_observed"],1.1)
        self.assertAlmostEqual(valid["geometric_mean_absolute_fold_error"],1.1)
        self.assertEqual(valid["within_prespecified_fold_ratio_fraction"],1.0)

    def test_constant_prediction_calibration_slope_abstains(self):
        y=np.linspace(-2,2,30)
        p=np.zeros((2,30))
        res=summarize_heldout(y,p)
        self.assertIsNone(res["pooled_calibration_slope_observed_on_predicted_descriptive"])
        self.assertIsNone(res["Spearman_per_unique_cell_line_repeat_mean_prediction"])
        self.assertLess(res["out_of_fold_R2_original_source_scale"],1.0)

    def test_same_biological_unit_and_finite_predictions_required(self):
        y=np.arange(30,dtype=float)
        with self.assertRaises(DiagnosticsError):
            summarize_heldout(y,np.ones((2,31)))
        with self.assertRaises(DiagnosticsError):
            summarize_heldout(y,np.full((2,30),np.nan))
        with self.assertRaises(DiagnosticsError):
            paired_heldout_contrast(y,np.ones((3,30)),np.ones((2,30)))


if __name__=="__main__":
    unittest.main()
