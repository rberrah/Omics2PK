"""Synthetic-only tests for a proposed, NOT validated, inductive SWEET variant."""
import unittest

try:
    import numpy as np
    from omics2pk.inductive_sweet import (
        CONSTRUCTION_V0,
        CONSTRUCTION_V1,
        InductiveNetworkError,
        fit_frozen_reference,
        gene_strength_features,
        predict_network,
    )
except ModuleNotFoundError as exc:
    if exc.name != "numpy":
        raise
    np = None


@unittest.skipIf(np is None, "NumPy optional extra not installed; tested in inductive workflow")
class InductiveSWEETTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261016)
        self.genes = tuple(f"GENE_{i}" for i in range(9))
        self.subjects = tuple(f"TRAIN_{i}" for i in range(12))
        self.matrix = rng.normal(size=(len(self.genes), len(self.subjects)))
        self.query = rng.normal(size=len(self.genes))
        self.reference = fit_frozen_reference(
            self.matrix, genes=self.genes, subject_ids=self.subjects, tissue="blood")

    def predict(self, x, patient="NEW_1"):
        return predict_network(self.reference, x, subject_id=patient,
                               genes=self.genes, tissue="blood")

    def test_network_and_strength_are_finite_and_symmetric(self):
        graph, info = self.predict(self.query)
        self.assertEqual(graph.shape, (9, 9))
        self.assertTrue(np.isfinite(graph).all())
        np.testing.assert_allclose(graph, graph.T, atol=1e-12)
        np.testing.assert_array_equal(np.diag(graph), np.zeros(9))
        self.assertEqual(len(gene_strength_features(graph, genes=self.genes)), 9)
        self.assertFalse(info["equivalent_to_official_sweet"])
        self.assertFalse(info["validated_for_pk"])
        self.assertEqual(info["construction"],CONSTRUCTION_V1)
        self.assertTrue(info["query_specific_anchor_includes_query_once"])
        self.assertTrue(info["perturbation_duplicates_query_once"])
        self.assertEqual(info["effective_cohort_n_for_weight"],len(self.subjects)+1)


    def test_v1_matches_query_anchor_then_duplicate_formula(self):
        graph, info = self.predict(self.query)
        anchor = np.corrcoef(np.column_stack((self.reference.expression, self.query)))
        perturbed = np.corrcoef(np.column_stack(
            (self.reference.expression, self.query, self.query)))
        expected = anchor + info["weight"] * (perturbed - anchor)
        expected = (expected + expected.T) * 0.5
        np.fill_diagonal(expected,0.0)
        np.testing.assert_allclose(graph,expected,rtol=0,atol=1e-12)

    def test_v0_is_retained_for_reproducibility_but_not_default(self):
        v1, i1 = self.predict(self.query)
        v0, i0 = predict_network(
            self.reference,self.query,subject_id="NEW_1",genes=self.genes,
            tissue="blood",construction=CONSTRUCTION_V0)
        self.assertEqual(i0["construction"],CONSTRUCTION_V0)
        self.assertEqual(i0["effective_cohort_n_for_weight"],len(self.subjects))
        self.assertFalse(i0["query_specific_anchor_includes_query_once"])
        self.assertFalse(np.array_equal(v0,v1))

    def test_unknown_construction_rejected(self):
        with self.assertRaisesRegex(InductiveNetworkError,"Unknown"):
            predict_network(
                self.reference,self.query,subject_id="NEW_1",genes=self.genes,
                tissue="blood",construction="invented")

    def test_other_test_subjects_cannot_change_frozen_prediction(self):
        first, a = self.predict(self.query)
        self.predict(self.query + 0.8, patient="NEW_2")
        second, b = self.predict(self.query)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(a["reference_fingerprint"], b["reference_fingerprint"])
        self.assertFalse(self.reference.expression.flags.writeable)

    def test_query_subject_inside_reference_rejected(self):
        with self.assertRaisesRegex(InductiveNetworkError, "new"):
            self.predict(self.query, patient="TRAIN_1")

    def test_gene_order_or_tissue_mismatch_rejected(self):
        with self.assertRaisesRegex(InductiveNetworkError, "Gene order"):
            predict_network(self.reference, self.query, subject_id="NEW",
                            genes=tuple(reversed(self.genes)), tissue="blood")
        with self.assertRaisesRegex(InductiveNetworkError, "tissue"):
            predict_network(self.reference, self.query, subject_id="NEW",
                            genes=self.genes, tissue="liver")

    def test_constant_reference_gene_rejected(self):
        a = self.matrix.copy()
        a[0, :] = 2.0
        with self.assertRaisesRegex(InductiveNetworkError, "Constant"):
            fit_frozen_reference(a, genes=self.genes, subject_ids=self.subjects, tissue="blood")

    def test_constant_query_and_nan_rejected(self):
        with self.assertRaisesRegex(InductiveNetworkError, "constant"):
            self.predict(np.full(len(self.genes), 1.0))
        with self.assertRaisesRegex(InductiveNetworkError, "Invalid query"):
            self.predict(np.full(len(self.genes), np.nan))

    def test_degenerate_small_reference_rejected(self):
        with self.assertRaisesRegex(InductiveNetworkError, "five"):
            fit_frozen_reference(self.matrix[:, :4], genes=self.genes,
                                 subject_ids=self.subjects[:4], tissue="blood")

    def test_fingerprint_tracks_training_values(self):
        alter = self.matrix.copy()
        alter[0, 0] += 0.00001
        fresh = fit_frozen_reference(alter, genes=self.genes,
                                     subject_ids=self.subjects, tissue="blood")
        self.assertNotEqual(self.reference.fingerprint, fresh.fingerprint)


if __name__ == "__main__":
    unittest.main()
