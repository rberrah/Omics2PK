"""Synthetic only: NEVER include or serialize source patient records."""
import json
import unittest

try:
    import numpy as np
    import sklearn
    from omics2pk.predose_representation_benchmark import (
        compare_predose_rna_pk, RepresentationBenchmarkError)
except ModuleNotFoundError as exc:
    if exc.name not in ("numpy", "sklearn"):
        raise
    np = None


@unittest.skipIf(np is None, "Requires optional NumPy/sklearn benchmark dependencies")
class PredoseRepresentationBenchmarkTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261019)
        n = 12
        self.inputs = dict(
            subject_ids=[f"FAKE_PERSON_{i:02d}" for i in range(n)],
            gene_ids=[f"GENE_{j}" for j in range(7)],
            expression=rng.normal(5, 1, size=(7, n)),
            clinical=rng.normal(size=(n, 2)),
            measured_pk=np.exp(rng.normal(2, 0.2, size=n)),
            genotype_dosages=rng.integers(0, 3, size=(n, 2)).astype(float),
            collected_at=["2025-01-02T09:00:00+00:00"] * n,
            first_dose_at=["2025-01-03T09:00:00+00:00"] * n,
            tissue="liver", target_name="synthetic_cl", target_unit="fake_L_per_h",
            source_kind="synthetic", source_sha256="a" * 64,
            prespecified_gene_weights={"GENE_0": 0.7, "GENE_1": -0.3},
        )

    def test_complete_competing_models_share_same_persons(self):
        report = compare_predose_rna_pk(**self.inputs)
        self.assertEqual(report["n_independent_subjects"], 12)
        self.assertEqual(report["status"], "SYNTHETIC_SOFTWARE_CHECK_NOT_PK_EVIDENCE")
        self.assertEqual(report["primary_reference"], "clinical_pgx")
        self.assertEqual(len(report["models"]), 7)
        self.assertTrue(all(m["n_biological_test_people"] == 12
                            for m in report["models"].values()))
        self.assertTrue(all(np.isfinite(m["mae"]) for m in report["models"].values()))
        text = json.dumps(report)
        self.assertNotIn("FAKE_PERSON_", text)
        self.assertNotIn("GENE_0", text)
        self.assertNotIn(str(self.inputs["measured_pk"][0]), text)
        self.assertIn("not official SWEET", report["graph_qualifier"])

    def test_optional_genotype_and_biology_weights(self):
        self.inputs["genotype_dosages"] = None
        self.inputs["prespecified_gene_weights"] = None
        r = compare_predose_rna_pk(**self.inputs)
        self.assertEqual(r["primary_reference"], "clinical")
        self.assertEqual(len(r["models"]), 5)

    def test_duplicate_person_rejected_before_modeling(self):
        self.inputs["subject_ids"][1] = self.inputs["subject_ids"][0]
        with self.assertRaisesRegex(RepresentationBenchmarkError, "unique"):
            compare_predose_rna_pk(**self.inputs)

    def test_after_dose_omics_rejected(self):
        self.inputs["collected_at"][0] = "2025-01-04T08:00:00+00:00"
        with self.assertRaisesRegex(RepresentationBenchmarkError, "Postdose"):
            compare_predose_rna_pk(**self.inputs)

    def test_no_source_date_time_zone_rejected(self):
        self.inputs["collected_at"][0] = "2025-01-02T09:00:00"
        with self.assertRaisesRegex(RepresentationBenchmarkError, "timezones"):
            compare_predose_rna_pk(**self.inputs)

    def test_real_genotype_must_be_valid_dosage(self):
        self.inputs["genotype_dosages"][0, 0] = 1.5
        with self.assertRaisesRegex(RepresentationBenchmarkError, "0/1/2"):
            compare_predose_rna_pk(**self.inputs)

    def test_bad_predefined_score_weights_rejected(self):
        self.inputs["prespecified_gene_weights"] = {"NOT_ASSAYED": 1.0}
        with self.assertRaisesRegex(RepresentationBenchmarkError, "prespecified"):
            compare_predose_rna_pk(**self.inputs)

    def test_permuting_outcome_has_no_effect_on_graph_source_contract(self):
        # Null controls should be explicitly run in real prospective protocols.
        other = dict(self.inputs)
        other["measured_pk"] = self.inputs["measured_pk"][::-1].copy()
        r1 = compare_predose_rna_pk(**self.inputs)
        r2 = compare_predose_rna_pk(**other)
        self.assertNotEqual(r1["models"]["clinical_pgx"]["mae"],
                            r2["models"]["clinical_pgx"]["mae"])


if __name__ == "__main__":
    unittest.main()
