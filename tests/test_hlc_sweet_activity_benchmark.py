"""Synthetic qualification for the HLC site-disjoint SWEET activity benchmark."""
import unittest

try:
    import numpy as np
    from omics2pk.hlc_sweet_activity_benchmark import (
        HLCSweetBenchmarkError,
        evaluate_hlc_sweet_activity,
    )
except ModuleNotFoundError as exc:
    if exc.name != "numpy":
        raise
    np = None


@unittest.skipIf(np is None, "NumPy optional extra not installed")
class HLCSweetActivityBenchmarkTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261022)
        self.genes = tuple(["CYP3A4"] + [f"G{i:02}" for i in range(39)])
        self.ids = tuple(f"D{i:03}" for i in range(90))
        self.sites = tuple(["SITE_A"] * 30 + ["SITE_B"] * 30 + ["SITE_C"] * 30)
        self.age = rng.uniform(18, 85, 90)
        self.sex = np.asarray([i % 2 for i in range(90)], dtype=float)
        self.x = rng.normal(size=(40, 90))
        # Deliberately introduce site shifts so leave-site-out QC is meaningful.
        self.x[5:10, 30:60] += 0.5
        self.x[10:15, 60:] -= 0.4
        latent = (
            0.3 * self.x[0]
            + 0.15 * self.x[5] * self.x[6]
            + 0.005 * self.age
            + 0.1 * self.sex
        )
        self.activity = np.exp(latent - latent.min() + 0.3)

    def run_benchmark(self, x=None, activity=None):
        return evaluate_hlc_sweet_activity(
            self.x if x is None else x,
            genes=self.genes,
            subject_ids=self.ids,
            site_ids=self.sites,
            ages=self.age,
            sexes=self.sex,
            activities=self.activity if activity is None else activity,
            direct_cyp_gene="CYP3A4",
            activity_label="synthetic_midazolam_turnover",
            activity_unit="synthetic_units",
            edge_top_k=12,
            internal_crossfit_folds=5,
        )

    def test_site_disjoint_benchmark_runs_and_emits_only_aggregate_metrics(self):
        out = self.run_benchmark()
        self.assertEqual(out["status"], "HLC_SITE_DISJOINT_E2_BENCHMARK_EXECUTED")
        self.assertEqual(out["n_complete_same_donor_subjects"], 90)
        self.assertEqual(out["n_original_sites"], 3)
        self.assertEqual(out["validation"], "leave_one_original_site_out")
        self.assertFalse(out["site_used_as_predictor"])
        self.assertTrue(out["site_used_only_for_outer_holdout"])
        self.assertFalse(out["pk_pd_target_used_to_build_network"])
        self.assertFalse(out["heldout_expression_used_to_fit_preprocessing_or_gene_selection"])
        self.assertFalse(out["individual_predictions_exported"])
        self.assertTrue(out["not_clinical_PK_or_PopPK"])
        self.assertFalse(out["equivalent_to_official_sweet"])
        self.assertEqual(set(out["site_results"]), {"SITE_A", "SITE_B", "SITE_C"})
        self.assertNotIn("D000", str(out))
        for site in out["site_results"].values():
            self.assertEqual(site["n_test"], 30)
            self.assertEqual(site["n_train"], 60)
            self.assertEqual(site["n_edge_genes"], 12)
            self.assertGreaterEqual(site["heldout_weight_clamp_fraction"], 0)
            self.assertLessEqual(site["heldout_weight_clamp_fraction"], 1)

    def test_heldout_site_expression_cannot_change_its_training_selection_or_reference(self):
        a = self.run_benchmark()
        changed = self.x.copy()
        changed[:, :30] += np.linspace(-50, 50, 40)[:, None]
        b = self.run_benchmark(x=changed)
        self.assertEqual(
            a["site_results"]["SITE_A"]["edge_selection_fingerprint"],
            b["site_results"]["SITE_A"]["edge_selection_fingerprint"],
        )
        self.assertEqual(
            a["site_results"]["SITE_A"]["outer_reference_fingerprint"],
            b["site_results"]["SITE_A"]["outer_reference_fingerprint"],
        )

    def test_heldout_outcomes_cannot_change_network_fit_fingerprints(self):
        a = self.run_benchmark()
        activity = self.activity.copy()
        activity[:30] *= 50
        b = self.run_benchmark(activity=activity)
        self.assertEqual(
            a["site_results"]["SITE_A"]["edge_selection_fingerprint"],
            b["site_results"]["SITE_A"]["edge_selection_fingerprint"],
        )
        self.assertEqual(
            a["site_results"]["SITE_A"]["outer_reference_fingerprint"],
            b["site_results"]["SITE_A"]["outer_reference_fingerprint"],
        )

    def test_requires_three_sites_and_source_annotated_direct_gene(self):
        with self.assertRaisesRegex(HLCSweetBenchmarkError, "at least three"):
            evaluate_hlc_sweet_activity(
                self.x,
                genes=self.genes,
                subject_ids=self.ids,
                site_ids=("A",) * 45 + ("B",) * 45,
                ages=self.age,
                sexes=self.sex,
                activities=self.activity,
                direct_cyp_gene="CYP3A4",
                activity_label="x",
                activity_unit="u",
                edge_top_k=12,
            )
        with self.assertRaisesRegex(HLCSweetBenchmarkError, "source-annotated"):
            evaluate_hlc_sweet_activity(
                self.x,
                genes=self.genes,
                subject_ids=self.ids,
                site_ids=self.sites,
                ages=self.age,
                sexes=self.sex,
                activities=self.activity,
                direct_cyp_gene="CYP_NOT_PRESENT",
                activity_label="x",
                activity_unit="u",
                edge_top_k=12,
            )

    def test_activity_and_demographic_contracts_fail_closed(self):
        bad = self.activity.copy()
        bad[0] = -1
        with self.assertRaisesRegex(HLCSweetBenchmarkError, "cannot be negative"):
            self.run_benchmark(activity=bad)
        sex = self.sex.copy()
        sex[0] = 2
        with self.assertRaisesRegex(HLCSweetBenchmarkError, "0/1"):
            evaluate_hlc_sweet_activity(
                self.x,
                genes=self.genes,
                subject_ids=self.ids,
                site_ids=self.sites,
                ages=self.age,
                sexes=sex,
                activities=self.activity,
                direct_cyp_gene="CYP3A4",
                activity_label="x",
                activity_unit="u",
                edge_top_k=12,
            )


if __name__ == "__main__":
    unittest.main()
