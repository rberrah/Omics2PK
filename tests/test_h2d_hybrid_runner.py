"""Contract tests for the H2d runner and its erratum.

These guard the three corrections the external review asked for, so that none of
them can be lost: agreement must be reported three ways with the collapse
verdict on exact set equality, the hybrid's smaller panels must be compared at
matched dimension as well as at budget, and the nine remaining endpoints must
not be described as nine independent replications.

They also guard the thing that is easiest to lose quietly -- that the rule was
designed after seeing one endpoint, so that endpoint cannot confirm it.
"""
from __future__ import annotations

import importlib.util
import inspect
import subprocess
import sys
import unittest
from pathlib import Path

try:
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"pytest is required for this module: {exc}") from exc

try:
    import numpy  # noqa: F401  (import guard only)
    import sklearn  # noqa: F401  (import guard only)
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"numpy and scikit-learn are required: {exc}") from exc

REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = REPO_ROOT / "scripts" / "run_h2d_hybrid_selection.py"
PRESPEC = (REPO_ROOT / "docs"
           / "OMICS2PK_H2D_HYBRID_SELECTION_PRESPECIFICATION_20261108.md")
ERRATUM = REPO_ROOT / "docs" / "OMICS2PK_H2D_ERRATUM_20261108.md"


def load_runner():
    spec = importlib.util.spec_from_file_location("h2d_runner", RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H2D = load_runner()


class TestTheGeneratingEndpointCannotConfirm(unittest.TestCase):
    def test_it_is_excluded_from_the_confirmatory_set(self):
        self.assertNotIn(H2D.GENERATING_ENDPOINT, H2D.CONFIRMATORY_ENDPOINTS)

    def test_the_confirmatory_set_is_the_other_nine(self):
        self.assertEqual(len(H2D.CONFIRMATORY_ENDPOINTS), 9)
        self.assertEqual(set(H2D.CONFIRMATORY_ENDPOINTS) | {H2D.GENERATING_ENDPOINT},
                         set(H2D.H2A.ENDPOINTS))

    def test_the_generating_endpoint_is_the_one_that_produced_the_hypothesis(self):
        self.assertEqual(H2D.GENERATING_ENDPOINT, "CYP2D6_activity")

    def test_the_runner_caps_the_generating_endpoint_at_exploratory(self):
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn('"exploratory" if exploratory else', source)
        self.assertIn("why_exploratory", source)
        self.assertIn("cannot be confirmed here whatever it", source)


class TestAgreementIsReportedThreeWays(unittest.TestCase):
    """E1: containment alone was reported as convergence, and it is not."""

    def source(self):
        return RUNNER.read_text(encoding="utf-8")

    def test_the_collapse_verdict_uses_exact_set_equality(self):
        self.assertEqual(H2D.COLLAPSE_THRESHOLD, 0.90)
        self.assertIn("exact set equality", self.source())

    def test_the_verdict_excludes_panel_sizes_where_equality_is_expected(self):
        self.assertEqual(H2D.COLLAPSE_VERDICT_MIN_K, 5)
        self.assertIn("why_small_k_is_excluded", self.source())

    def test_all_three_measures_are_recorded(self):
        source = self.source()
        for field in ("mean_containment", "mean_jaccard",
                      "identical_set_fraction"):
            self.assertIn(field, source)

    def test_the_threshold_value_was_not_weakened(self):
        # The erratum changes the definition, never the number.
        self.assertIn("0.90", f"{H2D.COLLAPSE_THRESHOLD:.2f}")
        self.assertIn("threshold is unchanged in value",
                      ERRATUM.read_text(encoding="utf-8"))


class TestMatchedDimensionIsReported(unittest.TestCase):
    """E2: the hybrid may retain fewer than k, so budget alone is not enough."""

    def source(self):
        return RUNNER.read_text(encoding="utf-8")

    def test_losses_are_reported_at_budget_and_at_matched_dimension(self):
        source = self.source()
        self.assertIn("losses_at_budget", source)
        self.assertIn("losses_at_matched_dimension", source)

    def test_the_comparator_is_cut_to_what_the_hybrid_actually_used(self):
        source = self.source()
        self.assertIn("statistical[arm][: len(kept)]", source)
        self.assertIn("neither side is credited for being smaller", source)

    def test_the_number_actually_retained_is_reported_with_abstentions(self):
        source = self.source()
        for field in ("median_retained", "min_retained", "max_retained",
                      "abstentions"):
            self.assertIn(field, source)

    def test_non_inferiority_is_judged_at_both_dimensions(self):
        source = self.source()
        self.assertIn('"at_budget"', source)
        self.assertIn('"at_matched_dimension"', source)

    def test_the_margin_must_hold_against_every_statistical_arm(self):
        source = self.source()
        self.assertIn("non_inferiority_against_every_statistical_arm", source)
        self.assertEqual(H2D.STATISTICAL,
                         ("univariate_spearman", "lasso_entry_order",
                          "xgboost_gain"))


class TestTheValidationSetIsNotCalledIndependent(unittest.TestCase):
    """E4: the nine endpoints share a cohort and the enzymes are co-regulated."""

    def source(self):
        return RUNNER.read_text(encoding="utf-8")

    def test_the_runner_records_the_measured_dependence(self):
        caveat = H2D.VALIDATION_SET_INDEPENDENCE
        self.assertAlmostEqual(caveat["pairwise_spearman_median"], 0.446)
        self.assertAlmostEqual(caveat["first_component_share_of_variance"], 0.496)
        self.assertEqual(caveat["of_endpoints"], 10)
        # The key itself says what the number is not, so a reader of the output
        # cannot take it for a count of replications.
        self.assertAlmostEqual(
            caveat["dependence_index_not_a_replication_count"], 3.55)
        self.assertNotIn("effective_independent_endpoints", caveat)
        self.assertIn("validation_set_independence", self.source())

    def test_it_refuses_the_phrase_that_would_overstate_the_result(self):
        # Asserted on the value, not on the source text, so a line break inside
        # a string literal cannot make the guard pass vacuously.
        caveat = H2D.VALIDATION_SET_INDEPENDENCE
        self.assertIn("NOT independent replication", caveat["claim"])
        self.assertIn("No text may say 'confirmed on nine endpoints'",
                      caveat["counting_rule"])
        self.assertIn("a dependence index of about three and a half",
                      caveat["counting_rule"])
        self.assertIn("NOT a count of",
                      caveat["dependence_index_meaning"])
        self.assertIn("never be described as one",
                      caveat["dependence_index_meaning"])
        self.assertIn("would need another cohort", caveat["why_not_replication"])
        self.assertTrue(caveat["same_cohort"])

    def test_the_report_carries_the_caveat_by_reference(self):
        self.assertIn('"validation_set_independence": VALIDATION_SET_INDEPENDENCE',
                      self.source())

    def test_the_two_substrates_of_one_enzyme_count_once(self):
        self.assertEqual(H2D.SAME_ENZYME_PAIR,
                         ("CYP3A4M_activity", "CYP3A4T_activity"))
        source = self.source()
        self.assertIn("0.841", source)
        self.assertIn("count as one endpoint", source)


class TestNothingDivergesFromH2a(unittest.TestCase):
    """Only the retention step may differ, so the shared parameters are shared."""

    def test_the_runner_imports_the_h2a_module_rather_than_copying_it(self):
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn("def load_h2a()", source)
        self.assertIn("so every shared choice is literally the same code", source)

    def test_the_inherited_parameters_are_the_h2a_ones(self):
        self.assertEqual(H2D.H2A.SEED, 20261108)
        self.assertEqual(H2D.H2A.OUTER_FOLDS, 5)
        self.assertEqual(H2D.H2A.OUTER_REPEATS, 2)
        self.assertEqual(H2D.H2A.INNER_FOLDS, 3)
        self.assertEqual(H2D.H2A.PANEL_SIZES, (1, 2, 3, 5, 10, 20, 50, 100))
        self.assertAlmostEqual(H2D.H2A.NON_INFERIORITY_MARGIN, 0.05)

    def test_the_runner_reuses_the_h2a_fold_and_matrix_builders(self):
        source = RUNNER.read_text(encoding="utf-8")
        for call in ("H2A.repeat_assignment", "H2A.fold_matrix",
                     "H2A.baseline_matrix", "H2A.ridge_prediction",
                     "H2A.paired_contrast", "H2A.repeated_folds",
                     "H2A.statistical_orders", "H2A.build_universe"):
            self.assertIn(call, source)

    def test_the_destroyed_signal_check_is_carried_over(self):
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn("destroyed_signal_repair", source)
        self.assertIn("hybrid_kept_the_dead_gene", source)


class TestTheQualifiedPoolIsSharedWithH2a(unittest.TestCase):
    """H2a and H2d are comparable only on one universe.

    The probe behind a selection is chosen per training fold, so restricting the
    universe is not enough on its own: the feature each arm actually used has to
    be checked against the pool, fold by fold, rather than assumed admissible
    because its gene symbol was.
    """

    def source(self):
        return RUNNER.read_text(encoding="utf-8")

    def test_the_runner_offers_the_same_pool_flags_as_h2a(self):
        done = subprocess.run([sys.executable, str(RUNNER), "--help"],
                              capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)
        for flag in ("--probe-pool", "--pool-name", "--restrict-to-network"):
            self.assertIn(flag, done.stdout, flag)

    def test_the_default_pool_name_matches_the_h2a_default(self):
        done = subprocess.run([sys.executable, str(RUNNER), "--help"],
                              capture_output=True, text=True)
        self.assertIn("A_qualified_human", done.stdout)

    def test_the_endpoint_runner_receives_the_admissible_feature_set(self):
        parameters = inspect.signature(H2D.run_endpoint).parameters
        self.assertIn("pool_features", parameters)
        self.assertIsNone(parameters["pool_features"].default)

    def test_the_restriction_reuses_the_h2a_code_rather_than_a_copy(self):
        source = self.source()
        self.assertIn("H2A.qualified_pool_features", source)
        self.assertIn("H2A.restrict_universe_to_qualified_pool", source)
        self.assertIn("H2A.restrict_universe_to_network", source)

    def test_the_feature_used_for_each_selection_is_traced(self):
        source = self.source()
        self.assertIn("chosen=chosen_feature", source)
        self.assertIn("probe_traceability", source)
        self.assertIn("features_used_by_the_hybrid_and_the_target", source)
        self.assertIn("not traceable from the gene symbol alone", source)

    def test_a_selection_outside_the_pool_stops_the_run(self):
        source = self.source()
        self.assertIn("outside the qualified pool", source)
        self.assertIn("raise H2dError", source)

    def test_the_report_records_the_pool_so_a_comparison_can_be_checked(self):
        source = self.source()
        self.assertIn('"qualified_pool": universe.get("qualified_pool")', source)
        self.assertIn("same_pool_as_h2a", source)
        self.assertIn("an unrestricted H2d", source)

    def test_a_decoy_selection_is_counted_and_not_called_inadmissible(self):
        # A decoy has no feature at all; counting it as a pool violation would
        # confuse a selection error with an annotation problem.
        source = self.source()
        self.assertIn("H2A.DECOY_FEATURE_ID", source)
        self.assertIn("decoy_selections", source)


class TestTheErratumExistsAndWeakensNothing(unittest.TestCase):
    def text(self):
        return ERRATUM.read_text(encoding="utf-8")

    def test_the_erratum_exists_and_predates_execution(self):
        self.assertTrue(ERRATUM.is_file())
        self.assertIn("before the hybrid rule has been run on any cohort endpoint",
                      self.text())

    def test_it_covers_all_four_items(self):
        text = self.text()
        for item in ("E1", "E2", "E3", "E4"):
            self.assertIn(f"## {item}", text)

    def test_it_states_that_no_threshold_became_easier(self):
        text = self.text().lower()
        self.assertIn("no success threshold is changed", text)
        self.assertIn("makes a criterion easier to satisfy", text)

    def test_the_calibration_table_covers_several_regimes(self):
        text = self.text()
        for regime in ("pure noise", "correlated nulls",
                       "driven entirely by the baseline", "weak real signal"):
            self.assertIn(regime, text)

    def test_it_declines_to_claim_general_false_discovery_control(self):
        text = self.text().lower()
        self.assertIn("no claim of general false-discovery control is made", text)

    def test_it_does_not_claim_the_rule_is_calibrated(self):
        text = self.text()
        self.assertIn("do **not** demonstrate exact calibration", text)
        self.assertIn("claimed more than 120 draws can support", text)
        self.assertNotIn("The rule is calibrated on these designs", text)

    def test_it_calls_the_index_an_index_and_not_a_replication_count(self):
        flat = " ".join(self.text().split())
        self.assertIn("spectral index of dependence between endpoints", flat)
        self.assertIn("not a *count* of independent replications", flat)
        self.assertIn("no interpretation as", flat)

    def test_the_prespecification_it_amends_exists(self):
        self.assertTrue(PRESPEC.is_file())
        self.assertIn("1.0.0", PRESPEC.read_text(encoding="utf-8"))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()


class TestTheThreeRunnerDefectsStayFixed(unittest.TestCase):
    """Regression guards for three defects found by review before execution.

    All three are independent of the selection method, so none of them changes a
    hypothesis, a margin or a decision rule.
    """

    def source(self):
        return RUNNER.read_text(encoding="utf-8")

    # 1. The log loss needs predictions, not absolute errors.
    def test_out_of_fold_predictions_are_kept_for_every_arm(self):
        source = self.source()
        self.assertIn('predictions = {name: blank() for name in', source)
        self.assertIn('("hybrid_at_budget", "blind_mechanistic") + STATISTICAL',
                      source)
        self.assertIn('predictions["statistical_at_matched_dimension"]', source)

    def test_the_loss_summary_receives_each_arm_s_own_predictions(self):
        source = self.source()
        self.assertIn("per_donor(predictions[name][k])", source)
        self.assertNotIn(
            'per_donor(predictions["hybrid_at_budget"][k])\n'
            '                if name == "hybrid_at_budget" else', source)

    def test_the_matched_dimension_losses_are_full_summaries_not_a_bare_mae(self):
        source = self.source()
        self.assertIn(
            'per_donor(predictions["statistical_at_matched_dimension"][arm][k])',
            source)
        self.assertNotIn('{"mae": float(np.nanmean(per_donor(', source)

    def test_the_reason_is_recorded_beside_the_fix(self):
        # Comment markers survive whitespace flattening, so the phrases checked
        # here are ones that sit on a single source line.
        source = self.source()
        self.assertIn(
            "absolute errors in place of predictions leaves the plain MAE",
            source)
        self.assertIn("robustness co-criterion rather than a decoration", source)

    # 2. The destroyed-signal permutation must respect block separation.
    def test_the_permutation_happens_within_each_block(self):
        source = self.source()
        self.assertIn("for block_rows in (train, test):", source)
        self.assertIn("block_rows[shuffler.permutation(block_rows.size)]", source)
        self.assertIn('"within each block separately, per-fold seed"', source)

    def test_a_single_permutation_over_all_donors_is_gone(self):
        source = self.source()
        self.assertNotIn(").permutation(y.size)\n            wrecked[:,", source)
        self.assertNotIn("wrecked[:, position_of[target_gene]] = matrix[\n", source)

    def test_the_permutation_seed_is_per_fold_and_reproducible(self):
        source = self.source()
        self.assertIn("shuffler = np.random.default_rng(H2A.SEED + 31337 + fold_number)",
                      source)

    def test_within_block_permutation_keeps_each_block_s_own_values(self):
        # Demonstrated on a small array: the fix must be a permutation of each
        # block, so every block keeps exactly the multiset it started with.
        matrix = numpy.arange(20.0).reshape(10, 2)
        train = numpy.array([0, 1, 2, 3, 4, 5])
        test = numpy.array([6, 7, 8, 9])
        wrecked = matrix.copy()
        shuffler = numpy.random.default_rng(1)
        for block in (train, test):
            wrecked[block, 0] = matrix[block[shuffler.permutation(block.size)], 0]
        numpy.testing.assert_array_equal(numpy.sort(wrecked[train, 0]),
                                         numpy.sort(matrix[train, 0]))
        numpy.testing.assert_array_equal(numpy.sort(wrecked[test, 0]),
                                         numpy.sort(matrix[test, 0]))

    # 3. The repair criterion needs the univariate arm on the same matrix.
    def test_the_univariate_arm_is_computed_on_the_wrecked_matrix(self):
        source = self.source()
        self.assertIn("H2A._spearman_strength(wrecked[train], y[train])", source)
        self.assertIn('destroyed["univariate"]', source)
        self.assertIn("univariate_kept_the_dead_gene", source)

    def test_both_halves_of_the_repair_criterion_are_reported(self):
        source = self.source()
        self.assertIn("contrast_against_blind", source)
        self.assertIn("contrast_against_univariate", source)
        self.assertIn("repair_criterion_met", source)

    def test_the_criterion_requires_both_halves_and_says_so(self):
        source = self.source()
        self.assertIn('against_univariate["non_inferior_within_margin"]', source)
        self.assertIn('and against_blind["interval_excludes_zero"]', source)
        self.assertIn('and against_blind["delta_mae"] < 0.0', source)
        flat = " ".join(source.split())
        self.assertIn("Either one alone is not the repair", flat)
