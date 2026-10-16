"""Synthetic pipeline integrity; real outcome evaluation in separate public CI."""
import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
    from omics2pk.inductive_sweet import InductiveNetworkError
    from scripts.benchmark_public_multidrug import model_predict, network_split
except ModuleNotFoundError as exc:
    if exc.name not in {"numpy", "sklearn", "xgboost"}:
        raise
    np = None


@unittest.skipIf(np is None, "Optional measured-response benchmark dependencies")
class MeasuredBenchmarkTests(unittest.TestCase):
    def test_identical_outer_splits_and_frozen_network(self):
        rng = np.random.default_rng(42)
        x = rng.normal(size=(22, 8))
        ids = tuple(f"CELL{i}" for i in range(22))
        genes = tuple(f"ABCG{i}" for i in range(8))
        a, b, nt, nv, split = network_split(
            x, ids, genes, np.arange(15), np.arange(15,22))
        self.assertEqual(a.shape, (15,8))
        self.assertEqual(nt.shape, (15,8))
        self.assertEqual(nv.shape, (7,8))
        self.assertEqual(split["n_train"], 15)
        xx = x.copy()
        xx[-1, :] += rng.normal(size=8)
        _, _, nt2, nv2, split2 = network_split(
            xx, ids, genes, np.arange(15), np.arange(15,22))
        np.testing.assert_allclose(nt, nt2)
        np.testing.assert_allclose(nv[:-1], nv2[:-1])
        self.assertEqual(split["outer_reference_fingerprint"],
                         split2["outer_reference_fingerprint"])

    def test_ridge_and_xgboost_fixed_configurations_are_finite(self):
        rng = np.random.default_rng(43)
        x = rng.normal(size=(20, 8))
        y = rng.normal(size=20)
        for kind in ("ridge", "pca_ridge", "xgb"):
            with self.subTest(kind=kind):
                p = model_predict(x[:15], y[:15], x[15:], kind)
                self.assertEqual(p.shape, (5,))
                self.assertTrue(np.isfinite(p).all())

    def test_missing_expression_imputation_uses_train_only(self):
        rng = np.random.default_rng(44)
        x = rng.normal(size=(22, 8))
        x[0,0] = np.nan
        ids = tuple(f"CELL{i}" for i in range(22))
        genes = tuple(f"ABC{i}" for i in range(8))
        raw_train, raw_test, *_ = network_split(
            x, ids, genes, np.arange(15), np.arange(15,22))
        expected = np.nanmedian(x[:15,0])
        self.assertAlmostEqual(raw_train[0,0], expected)
        self.assertFalse(np.isnan(raw_test).any())


if __name__ == "__main__":
    unittest.main()
