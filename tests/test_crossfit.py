"""Synthetic-only split/leakage controls for experimental SWEET-inspired features."""
import unittest
from unittest.mock import patch

try:
    import numpy as np
    from omics2pk.inductive_sweet import (
        CONSTRUCTION_V0,
        CONSTRUCTION_V1,
        InductiveNetworkError,
        build_outer_split_features,
        crossfit_training_features,
        fit_frozen_reference,
        gene_strength_features,
        predict_network,
    )
except ModuleNotFoundError as exc:
    if exc.name != "numpy":
        raise
    np = None


@unittest.skipIf(np is None, "NumPy is optional: qualified in inductive-sweet workflow")
class CrossfitTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261016)
        self.genes = tuple(f"G{i}" for i in range(16))
        self.ids = tuple(f"SUB{i:02}" for i in range(18))
        self.x = rng.normal(size=(len(self.genes), len(self.ids)))
        self.train = self.ids[:13]
        self.test = self.ids[13:]

    def split(self, x=None, train=None, test=None, folds=0):
        return build_outer_split_features(
            self.x if x is None else x, genes=self.genes,
            subject_ids=self.ids,
            train_ids=self.train if train is None else train,
            test_ids=self.test if test is None else test,
            tissue="liver", folds=folds, construction=CONSTRUCTION_V1)

    def test_leave_one_out_excludes_the_actual_subject(self):
        r = self.split()
        self.assertEqual(len(r["training"]), 13)
        self.assertEqual(len(r["test"]), 5)
        self.assertTrue(all(x["reference_n"] == 12 for x in r["training"]))
        self.assertTrue(all(x["reference_n"] == 13 for x in r["test"]))
        self.assertTrue(all(x["construction"] == CONSTRUCTION_V1 for x in r["training"]+r["test"]))
        self.assertEqual({x["effective_cohort_n_for_weight"] for x in r["training"]},{13})
        self.assertEqual({x["effective_cohort_n_for_weight"] for x in r["test"]},{14})
        record = r["training"][0]
        loo = fit_frozen_reference(
            self.x[:, 1:13], genes=self.genes,
            subject_ids=self.train[1:], tissue="liver")
        self.assertEqual(loo.fingerprint, record["reference_fingerprint"])
        graph, _ = predict_network(
            loo, self.x[:, 0], subject_id=self.train[0],
            genes=self.genes, tissue="liver", construction=CONSTRUCTION_V1)
        expected = gene_strength_features(graph, genes=self.genes)
        for gene in self.genes:
            self.assertAlmostEqual(record["features"][gene], expected[gene], places=12)
        self.assertFalse(r["trained_with_pk_pd_target"])
        self.assertFalse(r["equivalent_to_official_sweet"])


    def test_build_outer_default_remains_legacy_v0(self):
        legacy = build_outer_split_features(
            self.x, genes=self.genes, subject_ids=self.ids,
            train_ids=self.train, test_ids=self.test, tissue="liver")
        self.assertEqual(legacy["construction"],CONSTRUCTION_V0)
        self.assertTrue(all(r["construction"] == CONSTRUCTION_V0
                            for r in legacy["training"]+legacy["test"]))

    def test_test_patient_and_other_test_patients_do_not_change_reference(self):
        first = self.split()
        shifted = self.x.copy()
        shifted[:, -1] += np.linspace(0, 6, len(self.genes))
        second = self.split(x=shifted)
        self.assertEqual(first["outer_reference_fingerprint"],
                         second["outer_reference_fingerprint"])
        self.assertEqual(first["training"], second["training"])
        self.assertEqual(first["test"][:-1], second["test"][:-1])

    def test_outer_training_covariates_do_not_see_any_outer_test_subject(self):
        first = self.split()
        shifted = self.x.copy()
        shifted[:, 13:] += np.random.default_rng(42).normal(size=(16, 5)) * 2
        second = self.split(x=shifted)
        self.assertEqual(first["training"], second["training"])
        self.assertEqual(first["outer_reference_fingerprint"],
                         second["outer_reference_fingerprint"])

    def test_fold_crossfit_fits_one_frozen_reference_per_fold(self):
        with patch(
            "omics2pk.inductive_sweet.fit_frozen_reference",
            wraps=fit_frozen_reference,
        ) as mocked_fit:
            records = crossfit_training_features(
                self.x[:, :13],
                genes=self.genes,
                subject_ids=self.train,
                tissue="liver",
                folds=3,
                construction=CONSTRUCTION_V1,
            )
        self.assertEqual(len(records), 13)
        self.assertEqual(mocked_fit.call_count, 3)
        groups = {}
        for record in records:
            groups.setdefault(record["excluded_fold"], set()).add(
                record["reference_fingerprint"]
            )
        self.assertEqual(len(groups), 3)
        self.assertTrue(all(len(fingerprints) == 1 for fingerprints in groups.values()))

    def test_crossfit_references_exclude_complete_training_fold(self):
        r = self.split(folds=3)
        groups = {}
        for record in r["training"]:
            groups.setdefault(record["excluded_fold"], []).append(record)
        self.assertEqual(len(groups), 3)
        for group in groups.values():
            self.assertEqual(len({e["reference_fingerprint"] for e in group}), 1)
            self.assertEqual(group[0]["reference_n"], 13 - len(group))
        self.assertEqual({item["reference_n"] for item in r["test"]}, {13})

    def test_fold_allocation_does_not_depend_on_input_column_order(self):
        original = self.split(folds=3)
        order = list(reversed(range(18)))
        reordered = build_outer_split_features(
            self.x[:, order], genes=self.genes,
            subject_ids=tuple(self.ids[i] for i in order),
            train_ids=self.train, test_ids=self.test, tissue="liver", folds=3,
            construction=CONSTRUCTION_V1)
        a = {r["subject_id"]: r["excluded_fold"] for r in original["training"]}
        b = {r["subject_id"]: r["excluded_fold"] for r in reordered["training"]}
        self.assertEqual(a, b)

    def test_duplicate_and_overlapping_subject_ids_fail_closed(self):
        with self.assertRaisesRegex(InductiveNetworkError, "partition"):
            self.split(test=(self.train[0],) + self.test[1:])
        with self.assertRaisesRegex(InductiveNetworkError, "partition"):
            self.split(test=self.test[:-1])
        with self.assertRaisesRegex(InductiveNetworkError, "partition"):
            self.split(train=self.train + (self.train[0],))

    def test_loo_requires_sufficient_training_subjects(self):
        with self.assertRaisesRegex(InductiveNetworkError, "six"):
            crossfit_training_features(
                self.x[:, :5], genes=self.genes, subject_ids=self.ids[:5],
                tissue="liver")
        with self.assertRaisesRegex(InductiveNetworkError, "Too few"):
            crossfit_training_features(
                self.x[:, :6], genes=self.genes, subject_ids=self.ids[:6],
                tissue="liver", folds=2)

    def test_folds_boolean_and_invalid_count_rejected(self):
        for folds in (True, -1, 1, 100):
            with self.subTest(folds=folds):
                with self.assertRaises(InductiveNetworkError):
                    self.split(folds=folds)


    def test_full_weight_gene_universe_is_separate_and_query_is_subject_specific(self):
        rng = np.random.default_rng(20261017)
        all_genes = self.genes + tuple(f"EXTRA_{j}" for j in range(24))
        all_expr = np.vstack([self.x, rng.normal(size=(24, len(self.ids)))])
        r = build_outer_split_features(
            self.x, genes=self.genes, subject_ids=self.ids,
            train_ids=self.train, test_ids=self.test, tissue="liver",
            weight_expression=all_expr, weight_genes=all_genes,
            construction=CONSTRUCTION_V1)
        self.assertEqual(r["n_weight_genes"], 40)
        self.assertEqual(len(r["training"]), 13)
        self.assertEqual(len(r["test"]), 5)
        self.assertTrue(all(rec["reference_n"] == 12 for rec in r["training"]))
        self.assertTrue(all(rec["reference_n"] == 13 for rec in r["test"]))
        # A perturbation in another held-out subject's full-gene profile cannot
        # alter any training feature or another held-out prediction.
        changed = all_expr.copy()
        changed[19:, -1] += np.linspace(-2.0, 3.0, 21)
        r2 = build_outer_split_features(
            self.x, genes=self.genes, subject_ids=self.ids,
            train_ids=self.train, test_ids=self.test, tissue="liver",
            weight_expression=changed, weight_genes=all_genes,
            construction=CONSTRUCTION_V1)
        self.assertEqual(r["training"], r2["training"])
        self.assertEqual(r["test"][:-1], r2["test"][:-1])
        self.assertEqual(r["outer_reference_fingerprint"], r2["outer_reference_fingerprint"])

    def test_full_gene_alignment_and_missing_query_are_rejected(self):
        rng = np.random.default_rng(123)
        wgenes = self.genes + ("EXTRA",)
        wx = np.vstack([self.x[:, :13], rng.normal(size=(1, 13))])
        ref = fit_frozen_reference(
            self.x[:, :13], genes=self.genes, subject_ids=self.train,
            tissue="liver", weight_expression=wx, weight_genes=wgenes)
        with self.assertRaisesRegex(InductiveNetworkError, "full aligned"):
            predict_network(ref, self.x[:, 13], subject_id=self.test[0],
                            genes=self.genes, tissue="liver")
        with self.assertRaisesRegex(InductiveNetworkError, "full aligned"):
            predict_network(ref, self.x[:, 13], subject_id=self.test[0],
                            genes=self.genes, tissue="liver",
                            weight_query=np.r_[self.x[:, 13], 0.7],
                            weight_genes=tuple(reversed(wgenes)))
        with self.assertRaisesRegex(InductiveNetworkError, "must be supplied together"):
            fit_frozen_reference(
                self.x[:, :13], genes=self.genes, subject_ids=self.train,
                tissue="liver", weight_expression=wx)
        with self.assertRaisesRegex(InductiveNetworkError, "include edge genes"):
            fit_frozen_reference(
                self.x[:, :13], genes=self.genes, subject_ids=self.train,
                tissue="liver", weight_expression=wx,
                weight_genes=self.genes[:-1] + ("EXTRA", "OTHER"))

    def test_weighted_network_v1_equation_matches_independent_numpy_oracle(self):
        rng = np.random.default_rng(881)
        full_genes = self.genes + tuple(f"OTHER_{i}" for i in range(36))
        weight_values = np.vstack([self.x, rng.normal(size=(36, len(self.ids)))])
        ref = fit_frozen_reference(
            self.x[:, :13], genes=self.genes, subject_ids=self.train,
            tissue="liver", weight_expression=weight_values[:, :13],
            weight_genes=full_genes)
        predicted, info = predict_network(
            ref, self.x[:, 13], subject_id=self.test[0],
            genes=self.genes, tissue="liver",
            weight_query=weight_values[:, 13], weight_genes=full_genes,
            construction=CONSTRUCTION_V1)
        # Independent translation of the query-anchored v1 formula.
        anchor = np.corrcoef(
            np.column_stack([self.x[:, :13], self.x[:, 13]]))
        contrasted = np.corrcoef(
            np.column_stack([self.x[:, :13], self.x[:, 13], self.x[:, 13]]))
        corr = np.corrcoef(weight_values[:, :13].T)
        ref_means = (corr.sum(axis=1) - 1) / 12
        augmented_w = np.corrcoef(
            np.column_stack([weight_values[:, :13], weight_values[:, 13]]).T)
        query_similarity = augmented_w[-1, :-1].mean()
        clipped_position = np.clip(
            (query_similarity - ref_means.min() + 0.01) /
            (ref_means.max() - ref_means.min() + 0.01), 0, 1)
        manual_weight = 0.1 * 14 * clipped_position
        expected = anchor + manual_weight * (contrasted - anchor)
        expected = (expected + expected.T) * 0.5
        np.fill_diagonal(expected, 0.0)
        np.testing.assert_allclose(predicted, expected, atol=1e-12, rtol=1e-12)
        self.assertAlmostEqual(info["weight"], float(manual_weight), places=12)
        self.assertEqual(info["effective_cohort_n_for_weight"],14)
        self.assertTrue(info["uses_full_weight_expression"])

    def test_weighted_network_v0_equation_remains_reproducible(self):
        rng = np.random.default_rng(882)
        full_genes = self.genes + tuple(f"OTHER_{i}" for i in range(12))
        weight_values = np.vstack([self.x, rng.normal(size=(12, len(self.ids)))])
        ref = fit_frozen_reference(
            self.x[:, :13], genes=self.genes, subject_ids=self.train,
            tissue="liver", weight_expression=weight_values[:, :13],
            weight_genes=full_genes)
        predicted, info = predict_network(
            ref, self.x[:, 13], subject_id=self.test[0],
            genes=self.genes, tissue="liver",
            weight_query=weight_values[:, 13], weight_genes=full_genes,
            construction=CONSTRUCTION_V0)
        base = np.corrcoef(self.x[:, :13])
        contrasted = np.corrcoef(
            np.column_stack([self.x[:, :13], self.x[:, 13]]))
        corr = np.corrcoef(weight_values[:, :13].T)
        ref_means = (corr.sum(axis=1) - 1) / 12
        augmented_w = np.corrcoef(
            np.column_stack([weight_values[:, :13], weight_values[:, 13]]).T)
        query_similarity = augmented_w[-1, :-1].mean()
        clipped_position = np.clip(
            (query_similarity - ref_means.min() + 0.01) /
            (ref_means.max() - ref_means.min() + 0.01), 0, 1)
        manual_weight = 0.1 * 13 * clipped_position
        expected = base + manual_weight * (contrasted - base)
        expected = (expected + expected.T) * 0.5
        np.fill_diagonal(expected, 0.0)
        np.testing.assert_allclose(predicted, expected, atol=1e-12, rtol=1e-12)
        self.assertAlmostEqual(info["weight"], float(manual_weight), places=12)

if __name__ == "__main__":
    unittest.main()
