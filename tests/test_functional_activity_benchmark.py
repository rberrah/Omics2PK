"""Synthetic and adversarial tests; NO original enzyme or PK results."""
from datetime import datetime, timedelta, timezone
from math import exp
from random import Random
import unittest

from omics2pk.functional_activity_benchmark import (
    FunctionalBenchmarkError, benchmark_cyp3a5_activity,
)


def fixture(n=60, include_metabolome=True):
    rng = Random(13)
    ids = [f"SYN-{i:03d}" for i in range(n)]
    allele = ["*3/*3", "*1/*3", "*1/*1"]
    start = datetime(2026, 1, 1, tzinfo=timezone.utc)
    times = [(start + timedelta(days=i)).isoformat() for i in range(n)]
    source = {
        "donor_ids": ids,
        "diplotypes": [allele[i % 3] for i in range(n)],
        "genotype_call_statuses": ["complete"] * n,
        "measured_activity": [],
        "activity_assay_at": times,
        "activity_unit": "pmol/min/mg",
        "activity_assay_validated": True,
        "source_kind": "synthetic",
        "source_sha256": "a" * 64,
        "modalities": {},
    }
    gen = []
    rn = []
    protein = []
    ms = []
    for i in range(n):
        g = i % 3
        functional_noise = rng.gauss(0, 1)
        metab_signal = functional_noise + 0.12 * rng.gauss(0, 1)
        canonical = .5 * g + functional_noise + .15 * rng.gauss(0, 1)
        full_protein = .65 * g + functional_noise + .2 * rng.gauss(0, 1)
        rn.append([rng.gauss(0, 2), canonical])
        protein.append([rng.gauss(0, 2), full_protein])
        ms.append([rng.gauss(0, 2), metab_signal])
        gen.append(exp(-.7 + .5 * g + .68 * functional_noise + rng.gauss(0, .08)))
    source["measured_activity"] = gen
    base = dict(tissue="liver", preprocessing_scope="per_sample",
                assay_relation_to_outcome="independent_assay")
    source["modalities"]["transcriptomics"] = dict(
        **base, donor_ids=list(ids), collected_at=list(times),
        feature_ids=["rna_total", "rna_canonical"], values=rn,
        annotations={
            "rna_total": {"target": "total_gene", "assay_validated": True},
            "rna_canonical": {"target": "canonical_splice_junction",
                              "assay_validated": True, "external_evidence": "CURATED_TEST_ONLY"},
        })
    source["modalities"]["proteomics"] = dict(
        **base, donor_ids=list(ids), collected_at=list(times),
        feature_ids=["pan_cyp", "full_length"], values=protein,
        annotations={
            "pan_cyp": {"target": "pan_cyp3a", "assay_validated": True},
            "full_length": {"target": "cyp3a5_full_length", "assay_validated": True,
                            "external_evidence": "CURATED_TEST_ONLY"},
        })
    if include_metabolome:
        source["modalities"]["metabolomics"] = dict(
            **base, donor_ids=list(ids), collected_at=list(times),
            feature_ids=["shared_metabolite", "specific_metabolite"],
            values=ms, annotations={
                "shared_metabolite": {"target": "metabolite_abundance",
                                      "assay_validated": True,
                                      "external_evidence": "CURATED_TEST_ONLY",
                                      "candidate_enzymes": ["CYP3A4", "CYP3A5"],
                                      "validated_as_cyp3a5_specific": False},
                "specific_metabolite": {"target": "metabolite_abundance",
                                        "assay_validated": True,
                                        "external_evidence": "CURATED_TEST_ONLY",
                                        "candidate_enzymes": ["CYP3A5"],
                                        "validated_as_cyp3a5_specific": True},
            })
    return source


class FunctionalActivityBenchmarkTests(unittest.TestCase):
    def test_five_frozen_arms_are_compared_on_same_donors(self):
        result = benchmark_cyp3a5_activity(**fixture())
        self.assertEqual(result["n_independent_donors"], 60)
        self.assertEqual(result["n_outer_donor_heldout_folds"], 5)
        self.assertEqual(result["n_knowledge_curated_features"],
                         {"metabolomics": 1, "proteomics": 1, "transcriptomics": 1})
        self.assertEqual(len(result["MAE_by_arm_observed_activity_units"]), 6)
        self.assertIn("genotype_only", result["MAE_by_arm_observed_activity_units"])
        self.assertFalse(result["in_vivo_PK_prediction_demonstrated"])
        self.assertFalse(result["uncertainty_calibrated"])
        self.assertIn("SYNTHETIC", result["status"])
        # Only synthetic ground truth deliberately constructed to contain signal.
        self.assertLess(result["MAE_by_arm_observed_activity_units"]["genotype_plus_curated_omics"],
                        result["MAE_by_arm_observed_activity_units"]["genotype_only"])

    def test_permutation_of_source_rows_invariant_when_aligned(self):
        d = fixture()
        r = benchmark_cyp3a5_activity(**d)
        reverse = list(reversed(range(len(d["donor_ids"]))))
        for key in ("donor_ids", "diplotypes", "genotype_call_statuses",
                    "measured_activity", "activity_assay_at"):
            d[key] = [d[key][i] for i in reverse]
        for mod in d["modalities"].values():
            for key in ("donor_ids", "values", "collected_at"):
                mod[key] = [mod[key][i] for i in reverse]
        again = benchmark_cyp3a5_activity(**d)
        self.assertEqual(r["MAE_by_arm_observed_activity_units"],
                         again["MAE_by_arm_observed_activity_units"])

    def test_mismatched_donor_id_order_fails(self):
        d = fixture()
        d["modalities"]["metabolomics"]["donor_ids"][0] = "INVALID"
        with self.assertRaisesRegex(FunctionalBenchmarkError, "donor ID order"):
            benchmark_cyp3a5_activity(**d)

    def test_post_outcome_metabolomics_fails(self):
        d = fixture()
        d["modalities"]["metabolomics"]["collected_at"][3] = "2026-12-31T01:00:00+00:00"
        with self.assertRaisesRegex(FunctionalBenchmarkError, "Post-outcome"):
            benchmark_cyp3a5_activity(**d)

    def test_missing_genotype_or_uncovered_allele_fails(self):
        for dip, status in [("*3/*3", "partial"), ("*4/*8", "complete")]:
            d = fixture()
            d["diplotypes"][0] = dip
            d["genotype_call_statuses"][0] = status
            with self.subTest(dip=dip), self.assertRaisesRegex(
                    FunctionalBenchmarkError, "Incomplete CYP3A5"):
                benchmark_cyp3a5_activity(**d)

    def test_unvalidated_activity_fails(self):
        d = fixture()
        d["activity_assay_validated"] = False
        with self.assertRaisesRegex(FunctionalBenchmarkError, "validated CYP3A5"):
            benchmark_cyp3a5_activity(**d)

    def test_missing_curated_evidence_fails(self):
        d = fixture()
        del d["modalities"]["proteomics"]["annotations"]["full_length"]["external_evidence"]
        with self.assertRaisesRegex(FunctionalBenchmarkError, "source evidence"):
            benchmark_cyp3a5_activity(**d)

    def test_curated_metabolome_is_not_automatically_shared_cyp3a(self):
        d = fixture()
        a = d["modalities"]["metabolomics"]["annotations"]["specific_metabolite"]
        a["candidate_enzymes"] = ["CYP3A4", "CYP3A5"]
        r = benchmark_cyp3a5_activity(**d)
        self.assertEqual(r["n_knowledge_curated_features"]["metabolomics"], 0)
        self.assertIn("genotype_plus_metabolomics", r["MAE_by_arm_observed_activity_units"])

    def test_external_feature_transform_not_assumed_when_preprocessing_unknown(self):
        d = fixture()
        d["modalities"]["transcriptomics"]["preprocessing_scope"] = "all_donor_batch_correction"
        with self.assertRaisesRegex(FunctionalBenchmarkError, "externally frozen"):
            benchmark_cyp3a5_activity(**d)

    def test_omic_missing_values_rejected(self):
        d = fixture()
        d["modalities"]["metabolomics"]["values"][12][1] = float("nan")
        with self.assertRaisesRegex(FunctionalBenchmarkError, "Missing or infinite"):
            benchmark_cyp3a5_activity(**d)

    def test_no_patient_data_exported(self):
        d = fixture()
        r = benchmark_cyp3a5_activity(**d)
        self.assertNotIn("SYN-000", str(r))
        self.assertNotIn("predictions_by_donor", r)

    def test_requested_intervals_use_independent_calibration_and_never_export_rows(self):
        d = fixture()
        out = benchmark_cyp3a5_activity(**d, interval_alpha=0.1)
        audit = out["prediction_interval_audit"]
        self.assertEqual(audit["status"], "EXPLORATORY_SPLIT_CONFORMAL_EXCHANGEABILITY_UNVERIFIED")
        self.assertEqual(audit["requested_nominal_coverage"], 0.9)
        self.assertTrue(audit["interval_centers_fit_on_proper_training_only"])
        self.assertFalse(audit["independent_external_calibration_tested"])
        self.assertEqual(set(audit["model_interval_metrics"]),
                         set(out["MAE_by_arm_observed_activity_units"]))
        for result in audit["model_interval_metrics"].values():
            self.assertGreaterEqual(result["min_calibration_donors_per_fold"], 9)
            self.assertGreater(result["mean_prediction_interval_width_in_endpoint_units"], 0)
            self.assertTrue(0 <= result["empirical_heldout_coverage"] <= 1)
        self.assertNotIn("SYN-", str(audit))

    def test_intervals_do_not_change_prespecified_full_outer_fit_maes(self):
        d = fixture()
        baseline = benchmark_cyp3a5_activity(**d)
        with_intervals = benchmark_cyp3a5_activity(**d, interval_alpha=.1)
        self.assertEqual(baseline["MAE_by_arm_observed_activity_units"],
                         with_intervals["MAE_by_arm_observed_activity_units"])
        self.assertEqual(baseline["prediction_interval_audit"]["status"], "NOT_REQUESTED")

    def test_intervals_order_invariant_and_fixed_seed(self):
        d = fixture()
        r = benchmark_cyp3a5_activity(**d, interval_alpha=.1)["prediction_interval_audit"]
        reverse = list(reversed(range(len(d["donor_ids"]))))
        for key in ("donor_ids", "diplotypes", "genotype_call_statuses",
                    "measured_activity", "activity_assay_at"):
            d[key] = [d[key][i] for i in reverse]
        for mod in d["modalities"].values():
            for key in ("donor_ids", "values", "collected_at"):
                mod[key] = [mod[key][i] for i in reverse]
        again = benchmark_cyp3a5_activity(**d, interval_alpha=.1)["prediction_interval_audit"]
        self.assertEqual(r, again)

    def test_small_n_high_coverage_abstains_without_fake_interval(self):
        d = fixture(n=30)
        r = benchmark_cyp3a5_activity(**d, interval_alpha=.01)
        audit = r["prediction_interval_audit"]
        self.assertEqual(audit["status"], "ABSTAIN_INSUFFICIENT_PROPER_TRAIN_OR_CALIBRATION_DONORS")
        self.assertIsNone(audit["model_interval_metrics"])
        self.assertIn("genotype_only", r["MAE_by_arm_observed_activity_units"])

    def test_invalid_interval_level_rejected(self):
        for bad in (-1, 0, 1, True, float("nan")):
            with self.subTest(bad=bad), self.assertRaisesRegex(FunctionalBenchmarkError, "interval_alpha"):
                benchmark_cyp3a5_activity(**fixture(), interval_alpha=bad)

    def test_activity_derived_metabolite_blocked_even_if_measured_at_landmark(self):
        d = fixture()
        d["modalities"]["metabolomics"]["assay_relation_to_outcome"] = "same_activity_assay_derived"
        with self.assertRaisesRegex(FunctionalBenchmarkError, "independent assay"):
            benchmark_cyp3a5_activity(**d)

    def test_unknown_assay_independence_blocked(self):
        d = fixture()
        d["modalities"]["proteomics"]["assay_relation_to_outcome"] = "unknown"
        with self.assertRaisesRegex(FunctionalBenchmarkError, "independent assay"):
            benchmark_cyp3a5_activity(**d)

    def test_works_without_metabolome_cannot_claim_metabolite_gain(self):
        d = fixture(include_metabolome=False)
        r = benchmark_cyp3a5_activity(**d)
        self.assertNotIn("genotype_plus_metabolomics", r["MAE_by_arm_observed_activity_units"])
        self.assertEqual(len(r["MAE_by_arm_observed_activity_units"]), 5)


if __name__ == "__main__":
    unittest.main()
