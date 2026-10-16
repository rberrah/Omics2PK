"""Synthetic-only tests for a proposed, NOT validated, inductive SWEET variant."""
import unittest

import numpy as np

from omics2pk.inductive_sweet import (
    FrozenReference,
    InductiveNetworkError,
    fit_frozen_reference,
    gene_strength_features,
    predict_network,
)


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
