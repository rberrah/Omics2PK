"""Contract tests for the H2a selection comparison.

These test the machinery, not the result: fold disjointness, that the shared
decomposition gives the same numbers as fitting each penalty separately, that
every per-fold choice is made on training donors only, and that the decoy and
destroyed-signal constructions are what the protocol says they are. None of
them needs the source archive, so they run wherever the repository runs.
"""
from __future__ import annotations

import importlib.util
import types
import unittest
from pathlib import Path

try:
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"pytest is required for this module: {exc}") from exc

try:
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"numpy and scikit-learn are required: {exc}") from exc

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "run_h2a_hlc_selection_comparison.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("h2a_runner", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H2A = load_runner()


class TestFrozenParameters(unittest.TestCase):
    def test_the_parameters_match_the_protocol(self):
        self.assertEqual(H2A.SEED, 20261108)
        self.assertEqual(H2A.OUTER_FOLDS, 5)
        self.assertEqual(H2A.OUTER_REPEATS, 2)
        self.assertEqual(H2A.INNER_FOLDS, 3)
        self.assertEqual(H2A.PANEL_SIZES, (1, 2, 3, 5, 10, 20, 50, 100))
        self.assertGreaterEqual(H2A.RANDOM_REPLICATIONS, 100)
        self.assertEqual(H2A.N_DECOYS, 200)
        self.assertAlmostEqual(H2A.NON_INFERIORITY_MARGIN, 0.05)
        self.assertAlmostEqual(H2A.TRIM_FRACTION, 0.10)
        self.assertEqual(H2A.N_BOOTSTRAP, 2000)

    def test_ten_endpoints_are_declared(self):
        self.assertEqual(len(H2A.ENDPOINTS), 10)
        self.assertEqual(len(set(H2A.ENDPOINTS)), 10)

    def test_the_wrong_mechanism_rotation_never_maps_an_endpoint_to_itself(self):
        for endpoint, wrong in H2A.WRONG_SEED_ROTATION.items():
            self.assertNotEqual(endpoint, wrong)
        # A rotation, so every endpoint is used exactly once as a wrong seed.
        self.assertEqual(sorted(H2A.WRONG_SEED_ROTATION.values()),
                         sorted(H2A.ENDPOINTS))

    def test_only_the_two_documented_substrates_are_declared(self):
        self.assertEqual(set(H2A.DOCUMENTED_SUBSTRATE),
                         {"CYP3A4M_activity", "CYP3A4T_activity"})


class TestFolds(unittest.TestCase):
    def test_every_donor_is_tested_once_per_repeat(self):
        folds = H2A.repeated_folds(97, n_splits=5, n_repeats=2, seed=1)
        self.assertEqual(len(folds), 10)
        for repeat in range(2):
            tested = np.concatenate([folds[repeat * 5 + i][1] for i in range(5)])
            np.testing.assert_array_equal(np.sort(tested), np.arange(97))

    def test_train_and_test_never_share_a_donor(self):
        for train, test in H2A.repeated_folds(97, n_splits=5, n_repeats=2, seed=1):
            self.assertEqual(np.intersect1d(train, test).size, 0)
            self.assertEqual(train.size + test.size, 97)

    def test_folds_are_reproducible_and_seed_dependent(self):
        first = H2A.repeated_folds(50, n_splits=5, n_repeats=1, seed=3)
        again = H2A.repeated_folds(50, n_splits=5, n_repeats=1, seed=3)
        other = H2A.repeated_folds(50, n_splits=5, n_repeats=1, seed=4)
        for (a, b), (c, d) in zip(first, again):
            np.testing.assert_array_equal(a, c)
            np.testing.assert_array_equal(b, d)
        self.assertFalse(all(np.array_equal(a[1], b[1])
                             for a, b in zip(first, other)))

    def test_the_two_repeats_are_not_the_same_partition(self):
        folds = H2A.repeated_folds(97, n_splits=5, n_repeats=2, seed=1)
        self.assertFalse(np.array_equal(np.sort(folds[0][1]),
                                        np.sort(folds[5][1])))

    def test_leaving_one_centre_out_gives_one_fold_per_centre(self):
        centre = np.array([1, 1, 1, 2, 2, 3])
        folds = H2A.centre_folds(centre)
        self.assertEqual(len(folds), 3)
        for (train, test), value in zip(folds, (1, 2, 3)):
            np.testing.assert_array_equal(np.unique(centre[test]), [value])
            self.assertNotIn(value, centre[train])


class TestRidgePath(unittest.TestCase):
    def reference(self, x_train, y_train, x_test, alpha):
        scaler = StandardScaler().fit(x_train)
        model = Ridge(alpha=alpha).fit(scaler.transform(x_train), y_train)
        return model.predict(scaler.transform(x_test))

    def test_agrees_with_a_separate_fit_when_donors_outnumber_variables(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=(120, 20))
        y = x[:, 0] * 2.0 - x[:, 1] + rng.normal(size=120)
        x_test = rng.normal(size=(15, 20))
        alphas = (0.01, 1.0, 100.0)
        for alpha, predicted in zip(
                alphas, H2A.ridge_path_predictions(x, y, x_test, alphas)):
            np.testing.assert_allclose(
                predicted, self.reference(x, y, x_test, alpha), atol=1e-9)

    def test_agrees_with_a_separate_fit_when_variables_outnumber_donors(self):
        rng = np.random.default_rng(1)
        x = rng.normal(size=(40, 500))
        y = x[:, 3] - x[:, 7] + rng.normal(size=40)
        x_test = rng.normal(size=(9, 500))
        alphas = (0.1, 10.0, 1000.0)
        for alpha, predicted in zip(
                alphas, H2A.ridge_path_predictions(x, y, x_test, alphas)):
            np.testing.assert_allclose(
                predicted, self.reference(x, y, x_test, alpha), atol=1e-8)

    def test_a_constant_variable_does_not_divide_by_zero(self):
        rng = np.random.default_rng(2)
        x = rng.normal(size=(40, 5))
        x[:, 2] = 7.0
        y = rng.normal(size=40)
        predicted = H2A.ridge_path_predictions(x, y, x[:5], (1.0,))[0]
        self.assertTrue(np.all(np.isfinite(predicted)))

    def test_a_larger_penalty_shrinks_towards_the_training_mean(self):
        rng = np.random.default_rng(3)
        x = rng.normal(size=(60, 4))
        y = x[:, 0] * 5.0 + rng.normal(size=60)
        x_test = rng.normal(size=(10, 4))
        light, heavy = H2A.ridge_path_predictions(x, y, x_test, (0.01, 1e8))
        self.assertLess(float(np.std(heavy)), float(np.std(light)))
        np.testing.assert_allclose(heavy, np.full(10, y.mean()), atol=1e-3)

    def test_the_penalty_is_chosen_on_the_training_fold_only(self):
        rng = np.random.default_rng(4)
        x = rng.normal(size=(80, 6))
        y = x[:, 0] + rng.normal(size=80) * 0.1
        # A test block of absurd values must not change the chosen penalty.
        quiet = np.zeros((10, 6))
        wild = np.full((10, 6), 1e6)
        _, alpha_quiet = H2A.ridge_prediction(x, y, quiet, inner_folds=3, seed=1)
        _, alpha_wild = H2A.ridge_prediction(x, y, wild, inner_folds=3, seed=1)
        self.assertEqual(alpha_quiet, alpha_wild)

    def test_the_chosen_penalty_comes_from_the_frozen_grid(self):
        rng = np.random.default_rng(5)
        x = rng.normal(size=(60, 4))
        y = rng.normal(size=60)
        _, alpha = H2A.ridge_prediction(x, y, x[:6], inner_folds=3, seed=1)
        self.assertIn(alpha, H2A.RIDGE_ALPHAS)


class TestLosses(unittest.TestCase):
    def test_trimming_removes_the_worst_errors(self):
        values = np.array([1.0, 2.0, 3.0, 100.0])
        self.assertAlmostEqual(H2A.trimmed_mean(values, 0.0), 26.5)
        self.assertAlmostEqual(H2A.trimmed_mean(values, 0.25), 2.0)

    def test_trimming_is_one_sided(self):
        values = np.array([0.0, 1.0, 2.0, 3.0, 400.0])
        self.assertAlmostEqual(H2A.trimmed_mean(values, 0.20), 1.5)

    def test_the_log_loss_skips_non_positive_pairs_and_counts_them(self):
        y_true = np.array([1.0, 2.0, -1.0, 4.0])
        y_pred = np.array([1.0, 2.0, 1.0, 0.0])
        errors = np.abs(y_pred - y_true)
        summary = H2A.loss_summary(errors, y_true, y_pred)
        self.assertEqual(summary["log_mae_pairs"], 2)
        self.assertAlmostEqual(summary["log_mae"], 0.0)
        self.assertEqual(summary["n"], 4)
        self.assertAlmostEqual(summary["mae"], float(errors.mean()))


class TestPairedContrast(unittest.TestCase):
    def test_a_negative_difference_means_the_first_arm_is_better(self):
        rng = np.random.default_rng(6)
        better = rng.normal(10.0, 1.0, size=200)
        worse = better + 5.0
        result = H2A.paired_contrast(better, worse, seed=1, label="x")
        self.assertLess(result["delta_mae"], 0)
        self.assertTrue(result["interval_excludes_zero"])
        self.assertIn("Negative means the first arm", result["how_to_read"])

    def test_no_difference_is_not_reported_as_equivalence(self):
        rng = np.random.default_rng(7)
        errors = rng.normal(10.0, 3.0, size=120)
        result = H2A.paired_contrast(errors, errors.copy(), seed=1, label="x")
        self.assertAlmostEqual(result["delta_mae"], 0.0)
        self.assertFalse(result["interval_excludes_zero"])
        self.assertIn("never that there is none", result["how_to_read"])

    def test_non_inferiority_needs_the_whole_interval_inside_the_margin(self):
        rng = np.random.default_rng(8)
        reference = rng.normal(10.0, 0.2, size=400)
        close = reference + 0.001
        far = reference * 1.5
        self.assertTrue(H2A.paired_contrast(close, reference, seed=1,
                                            label="x")["non_inferior_within_margin"])
        self.assertFalse(H2A.paired_contrast(far, reference, seed=1,
                                             label="x")["non_inferior_within_margin"])

    def test_the_relative_difference_uses_the_reference_arm(self):
        a = np.full(50, 9.0)
        b = np.full(50, 10.0)
        result = H2A.paired_contrast(a, b, seed=1, label="x")
        self.assertAlmostEqual(result["reference_mae"], 10.0)
        self.assertAlmostEqual(result["relative"], -0.1)


class TestSpearmanStrength(unittest.TestCase):
    def test_matches_a_reference_rank_correlation(self):
        from scipy.stats import spearmanr
        rng = np.random.default_rng(9)
        matrix = rng.normal(size=(60, 8))
        y = matrix[:, 0] * 2.0 + rng.normal(size=60)
        strength = H2A._spearman_strength(matrix, y)
        for column in range(matrix.shape[1]):
            expected = abs(spearmanr(matrix[:, column], y).statistic)
            self.assertAlmostEqual(strength[column], expected, places=9)

    def test_a_constant_column_scores_zero_rather_than_failing(self):
        matrix = np.column_stack([np.arange(20.0), np.full(20, 3.0)])
        strength = H2A._spearman_strength(matrix, np.arange(20.0))
        self.assertAlmostEqual(strength[0], 1.0)
        self.assertEqual(strength[1], 0.0)


def fake_data(values, symbols, *, age=None, sex=None):
    n = values.shape[1]
    return {
        "values": values,
        "symbols": np.asarray(symbols, dtype=str),
        "donors": np.array([f"{1 + i % 3}11{i:04d}" for i in range(n)]),
        "age": np.full(n, 50.0) if age is None else np.asarray(age, dtype=float),
        "sex": np.zeros(n) if sex is None else np.asarray(sex, dtype=float),
        "centre": np.array([1 + i % 3 for i in range(n)]),
        "phenotype": None,
    }


class TestFoldMatrix(unittest.TestCase):
    def universe(self, data):
        by_symbol = {}
        for position, symbol in enumerate(data["symbols"]):
            if symbol:
                by_symbol.setdefault(str(symbol), []).append(position)
        return {"by_symbol": by_symbol, "genes_on_network": [],
                "decoy_symbols": [], "decoy_values": np.zeros((0, data["values"].shape[1]))}

    def test_the_probe_with_the_highest_training_mean_represents_its_gene(self):
        values = np.array([
            [1.0, 1.0, 1.0, 1.0],      # gene A, low probe
            [5.0, 5.0, 5.0, 5.0],      # gene A, high probe
            [2.0, 2.0, 2.0, 2.0],      # gene B
        ])
        data = fake_data(values, ["AAA1", "AAA1", "BBB2"])
        matrix, names, is_decoy = H2A.fold_matrix(
            data, self.universe(data), np.array([0, 1, 2]), np.arange(4))
        self.assertEqual(sorted(names), ["AAA1", "BBB2"])
        column = matrix[:, list(names).index("AAA1")]
        np.testing.assert_allclose(column, 5.0)
        self.assertFalse(is_decoy.any())

    def test_the_choice_of_probe_uses_training_donors_only(self):
        # On the training donors the first probe is higher; over all donors the
        # second would be. The first must win.
        values = np.array([
            [9.0, 9.0, 0.0, 0.0],
            [1.0, 1.0, 99.0, 99.0],
        ])
        data = fake_data(values, ["AAA1", "AAA1"])
        matrix, names, _ = H2A.fold_matrix(
            data, self.universe(data), np.array([0, 1]), np.arange(4))
        np.testing.assert_allclose(matrix[:2, 0], [9.0, 9.0])

    def test_an_unannotated_probe_is_excluded(self):
        values = np.ones((2, 4))
        data = fake_data(values, ["AAA1", ""])
        _, names, _ = H2A.fold_matrix(
            data, self.universe(data), np.array([0, 1, 2]), np.arange(4))
        self.assertEqual(list(names), ["AAA1"])

    def test_a_missing_cell_takes_the_training_mean_of_its_probe(self):
        values = np.array([[2.0, 4.0, np.nan, 10.0]])
        data = fake_data(values, ["AAA1"])
        matrix, _, _ = H2A.fold_matrix(
            data, self.universe(data), np.array([0, 1]), np.arange(4))
        self.assertAlmostEqual(float(matrix[2, 0]), 3.0)
        self.assertAlmostEqual(float(matrix[3, 0]), 10.0)

    def test_a_probe_missing_across_most_of_a_training_fold_is_dropped(self):
        values = np.array([
            [np.nan, np.nan, np.nan, 1.0],
            [1.0, 2.0, 3.0, 4.0],
        ])
        data = fake_data(values, ["AAA1", "BBB2"])
        _, names, _ = H2A.fold_matrix(
            data, self.universe(data), np.array([0, 1, 2]), np.arange(4))
        self.assertEqual(list(names), ["BBB2"])

    def test_decoys_are_appended_and_flagged(self):
        values = np.array([[1.0, 2.0, 3.0, 4.0]])
        data = fake_data(values, ["AAA1"])
        universe = self.universe(data)
        universe["decoy_symbols"] = ["ZZZ9"]
        universe["decoy_values"] = np.array([[4.0, 3.0, 2.0, 1.0]])
        matrix, names, is_decoy = H2A.fold_matrix(
            data, universe, np.array([0, 1, 2]), np.arange(4))
        self.assertEqual(list(names), ["AAA1", "ZZZ9"])
        np.testing.assert_array_equal(is_decoy, [False, True])
        self.assertEqual(int(is_decoy.sum()), 1)


class TestBaselineMatrix(unittest.TestCase):
    def test_a_missing_age_takes_the_training_median(self):
        data = fake_data(np.ones((1, 5)), ["AAA1"],
                         age=[10.0, 20.0, 60.0, np.nan, 99.0])
        baseline = H2A.baseline_matrix(data, np.arange(5), np.array([0, 1, 2]))
        self.assertAlmostEqual(float(baseline[3, 0]), 20.0)

    def test_a_missing_sex_takes_the_training_majority(self):
        data = fake_data(np.ones((1, 5)), ["AAA1"],
                         sex=[1.0, 1.0, 1.0, np.nan, 0.0])
        baseline = H2A.baseline_matrix(data, np.arange(5), np.array([0, 1, 2]))
        self.assertAlmostEqual(float(baseline[3, 1]), 1.0)

    def test_a_covariate_absent_throughout_a_training_fold_is_an_error(self):
        data = fake_data(np.ones((1, 4)), ["AAA1"],
                         age=[np.nan, np.nan, 40.0, 50.0])
        with self.assertRaisesRegex(H2A.H2aError, "absent throughout"):
            H2A.baseline_matrix(data, np.arange(4), np.array([0, 1]))


class TestUniverse(unittest.TestCase):
    def graph(self, symbols):
        return types.SimpleNamespace(index={s: i for i, s in enumerate(symbols)})

    def test_decoys_carry_symbols_the_array_does_not_measure(self):
        H2A_decoys = H2A.N_DECOYS
        try:
            H2A.N_DECOYS = 3
            values = np.arange(40.0).reshape(4, 10)
            data = fake_data(values, ["AAA1", "BBB2", "CCC3", "DDD4"])
            graph = self.graph(["AAA1", "BBB2", "OFF1", "OFF2", "OFF3", "OFF4"])
            universe = H2A.build_universe(data, graph, np.random.default_rng(0))
            self.assertEqual(len(universe["decoy_symbols"]), 3)
            for symbol in universe["decoy_symbols"]:
                self.assertNotIn(symbol, universe["by_symbol"])
                self.assertIn(symbol, graph.index)
            self.assertEqual(universe["genes_on_network"], ["AAA1", "BBB2"])
        finally:
            H2A.N_DECOYS = H2A_decoys

    def test_a_decoy_keeps_a_real_probe_s_values_in_a_shuffled_order(self):
        saved = H2A.N_DECOYS
        try:
            H2A.N_DECOYS = 1
            values = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]])
            data = fake_data(values, ["AAA1"])
            graph = self.graph(["AAA1", "OFF1", "OFF2"])
            universe = H2A.build_universe(data, graph, np.random.default_rng(1))
            decoy = universe["decoy_values"][0]
            np.testing.assert_array_equal(np.sort(decoy), np.sort(values[0]))
        finally:
            H2A.N_DECOYS = saved

    def test_too_few_off_array_genes_is_refused_not_worked_around(self):
        saved = H2A.N_DECOYS
        try:
            H2A.N_DECOYS = 5
            data = fake_data(np.ones((2, 4)), ["AAA1", "BBB2"])
            graph = self.graph(["AAA1", "BBB2", "OFF1"])
            with self.assertRaisesRegex(H2A.H2aError, "too few network genes"):
                H2A.build_universe(data, graph, np.random.default_rng(0))
        finally:
            H2A.N_DECOYS = saved


class TestFullOrder(unittest.TestCase):
    def test_seeds_and_candidates_are_ordered_by_proximity_not_by_name(self):
        ranking = types.SimpleNamespace(
            seeds=("ZZZ9",), ranked=("AAA1", "MMM5"),
            scores={"ZZZ9": 0.5, "AAA1": 0.1, "MMM5": 0.3})
        self.assertEqual(H2A._full_order(ranking), ["ZZZ9", "MMM5", "AAA1"])

    def test_a_tie_falls_back_on_the_name_for_reproducibility(self):
        ranking = types.SimpleNamespace(
            seeds=("BBB2", "AAA1"), ranked=(),
            scores={"BBB2": 0.5, "AAA1": 0.5})
        self.assertEqual(H2A._full_order(ranking), ["AAA1", "BBB2"])

    def test_a_many_seed_panel_is_not_ordered_alphabetically(self):
        # The defect this function exists to prevent: a pathway seed set can
        # hold a hundred genes, and the alphabet must not decide its top-k.
        scores = {f"G{i:03d}": float(i) for i in range(100)}
        ranking = types.SimpleNamespace(
            seeds=tuple(sorted(scores)), ranked=(), scores=scores)
        order = H2A._full_order(ranking)
        self.assertEqual(order[0], "G099")
        self.assertNotEqual(order, sorted(scores))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
