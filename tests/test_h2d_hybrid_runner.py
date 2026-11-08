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
        source = self.source()
        self.assertIn("validation_set_independence", source)
        self.assertIn("effective_independent_endpoints", source)
        self.assertIn("3.55", source)
        self.assertIn("0.496", source)

    def test_it_refuses_the_phrase_that_would_overstate_the_result(self):
        # Asserted on the value, not on the source text, so a line break inside
        # a string literal cannot make the guard pass vacuously.
        caveat = H2D.VALIDATION_SET_INDEPENDENCE
        self.assertIn("NOT independent replication", caveat["claim"])
        self.assertIn("No text may say 'confirmed on nine endpoints'",
                      caveat["counting_rule"])
        self.assertIn("effective independent dimensionality about three and a half",
                      caveat["counting_rule"])
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
        self.assertIn("no claim of general fdr control is made", text)

    def test_the_prespecification_it_amends_exists(self):
        self.assertTrue(PRESPEC.is_file())
        self.assertIn("1.0.0", PRESPEC.read_text(encoding="utf-8"))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
