"""Contract tests for the hybrid selection rule.

The rule's whole point is that knowledge proposes and training data disposes, so
the tests that matter are: a shortlisted variable with no information must be
demoted, a shortlisted variable with information must be kept, nothing may be
decided from donors outside the training block, and an empty answer must be a
declared abstention rather than a silent fallback.

The prespecification document's own provenance is also checked, because the rule
was designed after seeing a result and any later text that forgets to say so
would be a misrepresentation.
"""
from __future__ import annotations

import re
import unittest
from pathlib import Path

try:
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"pytest is required for this module: {exc}") from exc

try:
    import numpy as np
    from scipy.stats import spearmanr  # noqa: F401  (import guard only)
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"numpy and scipy are required: {exc}") from exc

from omics2pk.hybrid_mechanistic_selection import (
    FROZEN,
    HybridSelectionError,
    agreement_with,
    partial_rank_strength,
    select_hybrid,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
PRESPEC = (REPO_ROOT / "docs"
           / "OMICS2PK_H2D_HYBRID_SELECTION_PRESPECIFICATION_20261108.md")


def build(n=240, seed=0):
    """A block where one shortlisted gene informs and one is pure noise."""
    rng = np.random.default_rng(seed)
    baseline = np.column_stack([rng.normal(50, 10, size=n),
                                rng.integers(0, 2, size=n).astype(float)])
    informative = rng.normal(size=n)
    endpoint = (3.0 * informative + 0.4 * baseline[:, 0]
                + 2.0 * baseline[:, 1] + rng.normal(size=n))
    dead = rng.normal(size=n)
    # A variable that only repeats the baseline carries nothing extra.
    echo = baseline[:, 0] + rng.normal(scale=0.01, size=n)
    matrix = np.column_stack([informative, dead, echo, rng.normal(size=(n, 6))])
    names = ["INFORMS", "DEAD", "ECHO"] + [f"NOISE{i}" for i in range(6)]
    available = {name: index for index, name in enumerate(names)}
    return matrix, endpoint, baseline, names, available


class TestFrozenParameters(unittest.TestCase):
    def test_the_frozen_values_match_the_prespecification(self):
        self.assertEqual(FROZEN["oversample"], 3)
        self.assertEqual(FROZEN["null_permutations"], 200)
        self.assertAlmostEqual(FROZEN["target_fdr"], 0.05)
        self.assertEqual(FROZEN["permutation_seed"], 20261108)

    def test_the_oversampling_factor_is_small_on_purpose(self):
        # A wide shortlist would let the data step do all the work, which is the
        # declared failure mode.
        self.assertLessEqual(FROZEN["oversample"], 5)


class TestPartialRankStrength(unittest.TestCase):
    def test_a_variable_that_only_repeats_the_baseline_scores_near_zero(self):
        matrix, endpoint, baseline, names, _ = build()
        strengths, _ = partial_rank_strength(matrix, endpoint, baseline)
        self.assertLess(strengths[names.index("ECHO")], 0.2)

    def test_an_informative_variable_scores_high(self):
        matrix, endpoint, baseline, names, _ = build()
        strengths, _ = partial_rank_strength(matrix, endpoint, baseline)
        self.assertGreater(strengths[names.index("INFORMS")], 0.6)

    def test_noise_scores_below_signal(self):
        matrix, endpoint, baseline, names, _ = build()
        strengths, _ = partial_rank_strength(matrix, endpoint, baseline)
        self.assertGreater(strengths[names.index("INFORMS")],
                           strengths[names.index("DEAD")])

    def test_the_residualised_endpoint_is_returned_for_the_null(self):
        from scipy.stats import rankdata

        matrix, endpoint, baseline, _, _ = build()
        _, residual = partial_rank_strength(matrix, endpoint, baseline)
        self.assertEqual(residual.shape, endpoint.shape)
        # The adjustment is made in rank space, so the residual is orthogonal to
        # the ranked baseline. Checking against the raw baseline instead would
        # be checking the weaker statistic this rule deliberately does not use.
        design = np.column_stack([np.ones(baseline.shape[0]),
                                  rankdata(baseline, axis=0)])
        np.testing.assert_allclose(design.T @ residual,
                                   np.zeros(design.shape[1]), atol=1e-8)

    def test_an_empty_block_is_handled(self):
        strengths, residual = partial_rank_strength(
            np.zeros((20, 0)), np.arange(20.0), np.ones((20, 1)))
        self.assertEqual(strengths.size, 0)
        self.assertEqual(residual.size, 20)

    def test_misaligned_input_is_refused(self):
        with self.assertRaisesRegex(HybridSelectionError, "must align"):
            partial_rank_strength(np.zeros((10, 2)), np.zeros(9), np.zeros((10, 1)))
        with self.assertRaisesRegex(HybridSelectionError, "two-dimensional"):
            partial_rank_strength(np.zeros(10), np.zeros(10), np.zeros((10, 1)))


class TestTheRuleDemotesWhatCarriesNothing(unittest.TestCase):
    def test_an_informative_gene_at_the_top_of_the_ordering_is_kept(self):
        matrix, endpoint, baseline, _, available = build()
        selection = select_hybrid(["INFORMS", "DEAD", "NOISE0"], available,
                                  matrix, endpoint, baseline, budget=1)
        self.assertEqual(selection.retained, ("INFORMS",))
        self.assertFalse(selection.abstained)

    def test_a_dead_gene_at_the_top_of_the_ordering_is_demoted(self):
        # This is the case the whole rule exists for: knowledge puts a gene
        # first, the data say it carries nothing, and it must not be selected.
        matrix, endpoint, baseline, _, available = build()
        selection = select_hybrid(["DEAD", "NOISE0", "NOISE1"], available,
                                  matrix, endpoint, baseline, budget=1)
        self.assertNotIn("DEAD", selection.retained)
        self.assertIn("DEAD", selection.demoted)

    def test_a_shortlist_of_pure_noise_mostly_produces_an_abstention(self):
        # Abstention on a null shortlist is a stochastic event, not a certainty:
        # the rule controls a false-discovery rate, so it must fire sometimes.
        # What matters is that the rate is near its nominal level rather than
        # that one particular draw abstains, so this measures the rate.
        noise = ["DEAD", "NOISE0", "NOISE1", "NOISE2", "NOISE3"]
        fired = 0
        trials = 40
        for seed in range(trials):
            matrix, endpoint, baseline, _, available = build(seed=seed)
            selection = select_hybrid(noise, available, matrix, endpoint,
                                      baseline, budget=2)
            fired += int(bool(selection.retained))
            self.assertLessEqual(len(selection.retained), 1,
                                 msg=f"two null variables retained at seed {seed}")
        rate = fired / trials
        self.assertLessEqual(rate, 0.25, msg=(
            f"a null shortlist produced a selection in {fired}/{trials} draws, "
            "which is too often for the frozen rate"))

    def test_an_abstention_demotes_the_whole_shortlist_rather_than_hiding_it(self):
        # Constant columns carry nothing with certainty, so abstention here is
        # not a matter of chance and the bookkeeping can be checked exactly.
        rng = np.random.default_rng(3)
        n = 200
        baseline = np.column_stack([rng.normal(50, 10, size=n),
                                    rng.integers(0, 2, size=n).astype(float)])
        endpoint = rng.normal(size=n)
        matrix = np.column_stack([np.full(n, 4.0), np.full(n, -1.0)])
        available = {"FLAT1": 0, "FLAT2": 1}
        selection = select_hybrid(["FLAT1", "FLAT2"], available, matrix,
                                  endpoint, baseline, budget=2)
        self.assertTrue(selection.abstained)
        self.assertEqual(selection.retained, ())
        self.assertEqual(set(selection.demoted), set(selection.shortlist))
        self.assertEqual(set(selection.statistics), {"FLAT1", "FLAT2"})
        self.assertIn("allowed answer", selection.how_to_read)

    def test_the_informative_gene_is_found_even_when_ranked_below_noise(self):
        matrix, endpoint, baseline, _, available = build()
        selection = select_hybrid(["NOISE0", "NOISE1", "INFORMS"], available,
                                  matrix, endpoint, baseline, budget=1)
        self.assertEqual(selection.retained, ("INFORMS",))


class TestShortlistAndBudget(unittest.TestCase):
    def test_the_shortlist_is_the_budget_times_the_oversampling_factor(self):
        matrix, endpoint, baseline, names, available = build()
        selection = select_hybrid(names, available, matrix, endpoint, baseline,
                                  budget=2, oversample=3)
        self.assertEqual(len(selection.shortlist), 6)

    def test_genes_absent_from_the_platform_are_skipped_without_comment(self):
        matrix, endpoint, baseline, _, available = build()
        selection = select_hybrid(["NOWHERE1", "INFORMS", "NOWHERE2"], available,
                                  matrix, endpoint, baseline, budget=1,
                                  oversample=3)
        self.assertEqual(selection.shortlist, ("INFORMS",))

    def test_the_retained_set_never_exceeds_the_budget(self):
        matrix, endpoint, baseline, names, available = build()
        selection = select_hybrid(names, available, matrix, endpoint, baseline,
                                  budget=1, oversample=5)
        self.assertLessEqual(len(selection.retained), 1)

    def test_a_zero_budget_abstains(self):
        matrix, endpoint, baseline, names, available = build()
        selection = select_hybrid(names, available, matrix, endpoint, baseline,
                                  budget=0)
        self.assertTrue(selection.abstained)
        self.assertEqual(selection.retained, ())

    def test_an_empty_ordering_abstains_rather_than_failing(self):
        matrix, endpoint, baseline, _, available = build()
        selection = select_hybrid([], available, matrix, endpoint, baseline,
                                  budget=3)
        self.assertTrue(selection.abstained)

    def test_bad_parameters_are_refused(self):
        matrix, endpoint, baseline, names, available = build()
        with self.assertRaisesRegex(HybridSelectionError, "budget cannot be negative"):
            select_hybrid(names, available, matrix, endpoint, baseline, budget=-1)
        with self.assertRaisesRegex(HybridSelectionError, "oversampling factor"):
            select_hybrid(names, available, matrix, endpoint, baseline,
                          budget=1, oversample=0)
        with self.assertRaisesRegex(HybridSelectionError, "permutation"):
            select_hybrid(names, available, matrix, endpoint, baseline,
                          budget=1, null_permutations=0)


class TestNothingIsDecidedOutsideTheTrainingBlock(unittest.TestCase):
    def test_the_rule_sees_only_the_rows_it_is_given(self):
        matrix, endpoint, baseline, names, available = build(n=240)
        train = np.arange(160)
        on_train = select_hybrid(names, available, matrix[train], endpoint[train],
                                 baseline[train], budget=2)
        # Corrupting the rows outside the training block cannot change anything,
        # because the rule never receives them.
        wrecked = matrix.copy()
        wrecked[160:] = 1e6
        again = select_hybrid(names, available, wrecked[train], endpoint[train],
                              baseline[train], budget=2)
        self.assertEqual(on_train.retained, again.retained)

    def test_the_decision_is_reproducible(self):
        matrix, endpoint, baseline, names, available = build()
        first = select_hybrid(names, available, matrix, endpoint, baseline,
                              budget=2, random_seed=11)
        second = select_hybrid(names, available, matrix, endpoint, baseline,
                               budget=2, random_seed=11)
        self.assertEqual(first.retained, second.retained)

    def test_the_parameters_that_produced_the_decision_travel_with_it(self):
        matrix, endpoint, baseline, names, available = build()
        selection = select_hybrid(names, available, matrix, endpoint, baseline,
                                  budget=2, random_seed=11)
        self.assertEqual(selection.parameters["budget"], 2)
        self.assertEqual(selection.parameters["permutation_seed"], 11)
        self.assertEqual(selection.parameters["shortlist_size"],
                         len(selection.shortlist))


class TestAgreementMeasure(unittest.TestCase):
    """Containment alone was reported as convergence, and it is not.

    A hybrid retaining one gene that a statistical arm's top five also contains
    scored a perfect 1.0 although the panels were plainly different sets. The
    declared failure mode is about the panels being the same, so exact set
    equality has to be measured separately.
    """

    def hybrid_of_one(self):
        matrix, endpoint, baseline, _, available = build()
        return select_hybrid(["INFORMS"], available, matrix, endpoint,
                             baseline, budget=1)

    def test_a_subset_is_contained_but_not_identical(self):
        agreement = agreement_with(self.hybrid_of_one(),
                                   ["INFORMS", "NOISE2", "NOISE3", "DEAD"])
        self.assertEqual(agreement["containment"], 1.0)
        self.assertFalse(agreement["identical"])
        self.assertAlmostEqual(agreement["jaccard"], 0.25)
        self.assertEqual(agreement["hybrid_size"], 1)
        self.assertEqual(agreement["comparator_size"], 4)

    def test_equal_panels_are_identical_on_every_measure(self):
        agreement = agreement_with(self.hybrid_of_one(), ["INFORMS"])
        self.assertEqual(agreement["containment"], 1.0)
        self.assertEqual(agreement["jaccard"], 1.0)
        self.assertTrue(agreement["identical"])

    def test_disjoint_panels_score_zero(self):
        agreement = agreement_with(self.hybrid_of_one(), ["NOISE0", "NOISE1"])
        self.assertEqual(agreement["containment"], 0.0)
        self.assertEqual(agreement["jaccard"], 0.0)
        self.assertFalse(agreement["identical"])

    def test_an_abstention_has_nothing_to_compare(self):
        rng = np.random.default_rng(3)
        n = 200
        baseline = np.column_stack([rng.normal(50, 10, size=n),
                                    rng.integers(0, 2, size=n).astype(float)])
        selection = select_hybrid(["FLAT1"], {"FLAT1": 0},
                                  np.full((n, 1), 2.0), rng.normal(size=n),
                                  baseline, budget=1)
        agreement = agreement_with(selection, ["FLAT1"])
        self.assertTrue(selection.abstained)
        self.assertTrue(agreement["abstained"])
        self.assertTrue(np.isnan(agreement["containment"]))
        self.assertTrue(np.isnan(agreement["jaccard"]))
        # Not False: an abstention is not a disagreement, and scoring it as one
        # would flatter the method.
        self.assertIsNone(agreement["identical"])

    def test_containment_and_equality_can_disagree_sharply(self):
        # The case that motivated the fix, asserted so it cannot return.
        agreement = agreement_with(
            self.hybrid_of_one(),
            ["INFORMS", "NOISE0", "NOISE1", "NOISE2", "NOISE3"])
        self.assertEqual(agreement["containment"], 1.0)
        self.assertFalse(agreement["identical"])
        self.assertLess(agreement["jaccard"], 0.25)


class TestCalibrationAcrossRegimes(unittest.TestCase):
    """The single null design was too thin to call the rule calibrated.

    Three regimes, all synthetic and none touching the cohort: pure noise,
    correlated candidates where a false positive is easiest to produce, and a
    weak signal where the question is sensitivity rather than specificity.
    """

    PERMUTATIONS = 200
    TRIALS = 120

    def run_regime(self, make, budget=2):
        fired = 0
        sizes = []
        for seed in range(self.TRIALS):
            matrix, endpoint, baseline, names, available = make(seed)
            selection = select_hybrid(names, available, matrix, endpoint,
                                      baseline, budget=budget,
                                      null_permutations=self.PERMUTATIONS)
            fired += int(bool(selection.retained))
            sizes.append(len(selection.retained))
        return fired / self.TRIALS, float(np.mean(sizes))

    def test_pure_noise_fires_near_the_nominal_rate(self):
        def make(seed):
            rng = np.random.default_rng(1000 + seed)
            n = 240
            baseline = np.column_stack([rng.normal(50, 10, size=n),
                                        rng.integers(0, 2, size=n).astype(float)])
            endpoint = rng.normal(size=n)
            matrix = rng.normal(size=(n, 6))
            names = [f"NULL{i}" for i in range(6)]
            return matrix, endpoint, baseline, names, {
                name: index for index, name in enumerate(names)}

        rate, mean_size = self.run_regime(make)
        self.assertLessEqual(rate, 0.20, msg=f"pure noise fired at {rate:.3f}")
        self.assertLess(mean_size, 0.5)

    def test_correlated_candidates_do_not_inflate_the_rate(self):
        # Correlated nulls share their fluctuations, which is where a
        # permutation null that ignored the dependence would break down.
        def make(seed):
            rng = np.random.default_rng(2000 + seed)
            n = 240
            baseline = np.column_stack([rng.normal(50, 10, size=n),
                                        rng.integers(0, 2, size=n).astype(float)])
            endpoint = rng.normal(size=n)
            shared = rng.normal(size=n)
            matrix = np.column_stack(
                [0.9 * shared + 0.436 * rng.normal(size=n) for _ in range(6)])
            names = [f"CORR{i}" for i in range(6)]
            return matrix, endpoint, baseline, names, {
                name: index for index, name in enumerate(names)}

        rate, _ = self.run_regime(make)
        self.assertLessEqual(rate, 0.20, msg=f"correlated nulls fired at {rate:.3f}")

    def test_an_endpoint_driven_entirely_by_the_baseline_is_not_mistaken(self):
        # Nothing in the candidate block adds information, so the rule must
        # abstain even though the endpoint is highly predictable.
        def make(seed):
            rng = np.random.default_rng(3000 + seed)
            n = 240
            age = rng.normal(50, 10, size=n)
            sex = rng.integers(0, 2, size=n).astype(float)
            baseline = np.column_stack([age, sex])
            endpoint = 5.0 * age + 20.0 * sex + rng.normal(size=n)
            matrix = np.column_stack([age + rng.normal(scale=0.01, size=n),
                                      sex, rng.normal(size=(n, 4))])
            names = ["ECHO_AGE", "ECHO_SEX"] + [f"NULL{i}" for i in range(4)]
            return matrix, endpoint, baseline, names, {
                name: index for index, name in enumerate(names)}

        rate, _ = self.run_regime(make)
        self.assertLessEqual(rate, 0.25, msg=(
            f"baseline echoes were selected in {rate:.3f} of draws"))

    def test_a_weak_signal_is_detected_more_often_than_noise(self):
        def make(seed):
            rng = np.random.default_rng(4000 + seed)
            n = 240
            baseline = np.column_stack([rng.normal(50, 10, size=n),
                                        rng.integers(0, 2, size=n).astype(float)])
            weak = rng.normal(size=n)
            endpoint = 0.45 * weak + rng.normal(size=n)
            matrix = np.column_stack([weak, rng.normal(size=(n, 5))])
            names = ["WEAK"] + [f"NULL{i}" for i in range(5)]
            return matrix, endpoint, baseline, names, {
                name: index for index, name in enumerate(names)}

        rate, _ = self.run_regime(make)
        self.assertGreater(rate, 0.60, msg=(
            f"a weak but real signal was only found in {rate:.3f} of draws"))

    def test_a_strong_signal_is_always_detected(self):
        def make(seed):
            rng = np.random.default_rng(5000 + seed)
            n = 200
            baseline = np.column_stack([rng.normal(50, 10, size=n),
                                        rng.integers(0, 2, size=n).astype(float)])
            strong = rng.normal(size=n)
            endpoint = 4.0 * strong + rng.normal(size=n)
            matrix = np.column_stack([strong, rng.normal(size=(n, 5))])
            names = ["STRONG"] + [f"NULL{i}" for i in range(5)]
            return matrix, endpoint, baseline, names, {
                name: index for index, name in enumerate(names)}

        self.TRIALS = 40
        rate, _ = self.run_regime(make, budget=1)
        self.assertEqual(rate, 1.0)


class TestThePrespecificationCarriesItsOwnProvenance(unittest.TestCase):
    """The rule was designed after seeing a result. The document must say so."""

    def text(self):
        return PRESPEC.read_text(encoding="utf-8")

    def test_the_document_exists(self):
        self.assertTrue(PRESPEC.is_file())

    def test_it_names_the_generating_endpoint_and_caps_it(self):
        text = self.text()
        self.assertIn("CYP2D6_activity", text)
        self.assertIn("generating endpoint", text)
        self.assertRegex(text, r"capped at \*\*`?exploratory`?\*\*|capped at `exploratory`")

    def test_it_states_that_the_rule_followed_the_observation(self):
        text = self.text().lower()
        self.assertIn("designed after seeing", text)
        self.assertIn("hypothesis-generating", text)

    def test_the_confirmatory_set_excludes_the_generating_endpoint(self):
        text = self.text().lower()
        self.assertIn("confirmatory set is the other nine endpoints", text)

    def test_it_declares_the_collapse_failure_mode_in_advance(self):
        text = self.text().lower()
        self.assertIn("collapses onto the univariate selector", text)
        self.assertIn("90", text)

    def test_it_resolves_the_label_collision_explicitly(self):
        text = self.text()
        self.assertIn("H2d", text)
        self.assertRegex(text, r'called it "H2c"|review\'s "H2c"')

    def test_it_refuses_to_make_mae_the_criterion(self):
        text = self.text().lower()
        self.assertIn("what is not a criterion", text)

    def test_the_frozen_values_in_the_document_match_the_code(self):
        text = self.text()
        for key, value in (("oversample", 3), ("null_permutations", 200),
                           ("target_fdr", 0.05),
                           ("permutation_seed", 20261108)):
            self.assertRegex(text, rf"`{key} = {re.escape(str(value))}`",
                             msg=f"{key} not stated as {value}")
            self.assertEqual(FROZEN[key], value)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
