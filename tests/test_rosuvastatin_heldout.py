"""Entirely synthetic tests; never publish original patient observations."""
from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from omics2pk.rosuvastatin_heldout import (
    ATTESTATIONS, BenchmarkGateError, evaluate_original_combined,
)


class PrivateE3FixedBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "PRIVATE_DO_NOT_COMMIT.csv"
        self.headers = ["patient", "AUC", "Cmax", "age", "BMI", "metab_a", "metab_b", "metab_c"]
        self.rows = [
            {
                "patient": f"PRIVATE_DONOR_{i:03d}",
                "AUC": str(100 + i + (i % 3) * 3),
                "Cmax": str(20 + i / 3 + i % 4),
                "age": str(20 + i % 50),
                "BMI": str(19 + i % 12),
                "metab_a": str(1000 + i * 17),
                "metab_b": str(400 + i * 3 + i % 7),
                "metab_c": "" if i % 7 == 0 else str(60 + i % 9),
            }
            for i in range(40)
        ]
        self.attested = {name: True for name in ATTESTATIONS}
        self.write()

    def write(self):
        with self.path.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=self.headers)
            w.writeheader()
            w.writerows(self.rows)
        self.sha = hashlib.sha256(self.path.read_bytes()).hexdigest()

    def evaluate(self, **kwargs):
        opts = dict(
            person_id_col="patient", auc_col="AUC", cmax_col="Cmax",
            omic_columns=["metab_a", "metab_b", "metab_c"],
            clinical_columns=["age", "BMI"],
            expected_sha256=self.sha,
            attestations=self.attested,
            bootstrap_replicates=100,
        )
        opts.update(kwargs)
        return evaluate_original_combined(self.path, **opts)

    def test_all_forty_real_humans_holdout_without_source_rows_export(self):
        got = self.evaluate()
        self.assertEqual(got["n_heldout_subject_folds"], 40)
        self.assertEqual(got["status"], "PRIVATE_SOURCE_ATTESTED_EXPLORATORY_E3_LODO_NOT_EXTERNAL_VALIDATION")
        self.assertEqual(set(got["outcomes"]), {"AUC", "Cmax"})
        self.assertEqual(got["n_selected_predose_omic_features"], 3)
        self.assertTrue(got["reference_is_clinical"])
        self.assertFalse(got["hyperparameter_search_performed"])
        self.assertFalse(got["input_features_prespecification_independently_certified"])
        serialized = json.dumps(got)
        self.assertNotIn("PRIVATE_DONOR_", serialized)
        self.assertNotIn("metab_a", serialized)
        self.assertNotIn('"patient"', serialized)
        self.assertTrue(all(
            key in got["outcomes"][outcome]["clinical_plus_predose_omics"]
            for outcome in ("AUC", "Cmax") for key in ("mae", "rmse", "bias")
        ))

    def test_shuffle_file_rows_does_not_change_heldout_predictions(self):
        old = self.evaluate()
        self.rows.reverse()
        self.write()
        new = self.evaluate()
        self.assertEqual(old["outcomes"], new["outcomes"])

    def test_source_hash_mismatch_abstains(self):
        with self.assertRaisesRegex(BenchmarkGateError, "SHA256 differs"):
            self.evaluate(expected_sha256="0" * 64)

    def test_manual_attestations_are_never_silently_assumed(self):
        a = dict(self.attested)
        a["all_selected_omic_signals_predose_verified"] = False
        with self.assertRaisesRegex(BenchmarkGateError, "four documented"):
            self.evaluate(attestations=a)

    def test_source_pk_outcome_cannot_be_a_molecular_predictor(self):
        with self.assertRaisesRegex(BenchmarkGateError, "PK outcome appears"):
            self.evaluate(omic_columns=["AUC", "metab_a"])

    def test_non_declared_postdose_concentration_column_blocks_all_fits(self):
        self.headers.append("rosuvastatin_concentration_postdose")
        for row in self.rows:
            row["rosuvastatin_concentration_postdose"] = "5"
        self.write()
        with self.assertRaisesRegex(BenchmarkGateError, "extra potential postdose PK"):
            self.evaluate()

    def test_duplicated_human_does_not_produce_forty_fold_fake_sample_size(self):
        self.rows[39]["patient"] = self.rows[0]["patient"]
        self.write()
        with self.assertRaisesRegex(BenchmarkGateError, "Duplicate"):
            self.evaluate()

    def test_feature_missingness_is_train_only_and_abstains_when_too_high(self):
        for row in self.rows[:15]:
            row["metab_c"] = "NA"
        self.write()
        with self.assertRaisesRegex(BenchmarkGateError, "80%"):
            self.evaluate()

    def test_no_source_clinical_comparator_does_not_claim_clinical_gain(self):
        result = self.evaluate(clinical_columns=[])
        self.assertFalse(result["reference_is_clinical"])
        self.assertEqual(result["n_selected_clinical_covariates"], 0)
        self.assertIn("clinical_or_constant_baseline", result["outcomes"]["AUC"])

    def test_no_predose_intensities_abstains(self):
        for row in self.rows:
            row["metab_a"] = "annotation"
            row["metab_b"] = "annotation"
        self.write()
        with self.assertRaises(BenchmarkGateError):
            self.evaluate()


if __name__ == "__main__":
    unittest.main()
