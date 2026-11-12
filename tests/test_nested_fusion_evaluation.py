"""Contracts for the nested fusion evaluation and the two leaks it closes.

The commissioning review of commit 5eda7c4 named four defects. Three are
testable and are tested here:

1. the context probe chose its compared arm pair from the whole cohort before
   splitting -> `TestArmPairChosenInsideEachFold`;
2. the stacking weights were fitted on a cached prediction layer whose base
   models had trained on the rows the weights were then scored on ->
   `TestNothingLeaksAcrossTheOuterSplit`, which makes the disjointness an
   executable guarantee rather than a comment;
3. an adjusted interval excluding zero was named `CONTEXT_CARRIES_SIGNAL` even
   at rho 0.17, below the declared 0.30 margin -> `TestProbeVerdictLadder` and
   `TestContrastVerdictLadder`.

The fourth — that the three endpoints were selected on these data — cannot be
repaired by code, only labelled, and `TestRunnerDeclarations` checks the label
is present and that the control endpoint names actually exist.
"""
from __future__ import annotations

import importlib.util
import pathlib

import pytest

np = pytest.importorskip("numpy")

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _load(name):
    path = ROOT / "src" / "omics2pk" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"omics2pk.{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def nfe():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import nested_fusion_evaluation
    return nested_fusion_evaluation


@pytest.fixture(scope="module")
def fhd():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import fusion_headroom_diagnostics
    return fusion_headroom_diagnostics


@pytest.fixture(scope="module")
def gateable(nfe):
    arm_a, arm_b, truth, context = nfe.gateable_cohort()
    arms = {"arm_a": arm_a, "arm_b": arm_b}

    def fit_predict(train_index, test_index):
        del train_index
        return {name: values[test_index] for name, values in arms.items()}

    return fit_predict, truth, context, list(arms)


# --------------------------------------------------------------------------
# 1. the probe's compared pair is a training-set decision
# --------------------------------------------------------------------------
class TestArmPairChosenInsideEachFold:
    def test_the_pair_is_flagged_as_chosen_inside_each_fold(self, fhd):
        predictions, truth, context = \
            fhd.counterexample_constant_weights_cannot_see_a_gate(200)
        probe = fhd.context_informativeness(predictions, truth, context)
        assert probe["arm_pair_chosen_inside_each_fold"] is True

    def test_one_pair_is_recorded_per_fold(self, fhd):
        predictions, truth, context = \
            fhd.counterexample_constant_weights_cannot_see_a_gate(200)
        probe = fhd.context_informativeness(predictions, truth, context,
                                            folds=5)
        assert len(probe["compared_arms_per_fold"]) == 5

    def test_nearly_tied_arms_make_the_pair_vary_between_folds(self, fhd):
        """A whole-cohort choice CANNOT produce two different pairs."""
        predictions, truth, context = \
            fhd.counterexample_constant_weights_cannot_see_a_gate(200)
        probe = fhd.context_informativeness(predictions, truth, context)
        distinct = {tuple(pair) for pair in probe["compared_arms_per_fold"]}
        assert len(distinct) > 1, (
            "with two exactly tied arms the fold-local choice must differ "
            "somewhere; an identical pair in every fold means the choice is "
            "still being made on the whole cohort")
        assert probe["arm_pair_stability"] < 1.0

    def test_a_dominant_arm_gives_a_perfectly_stable_pair(self, fhd):
        rng = np.random.default_rng(3)
        n = 300
        truth = rng.normal(size=n)
        predictions = np.column_stack([
            truth + rng.normal(scale=0.2, size=n),
            truth + rng.normal(scale=0.5, size=n),
            truth + rng.normal(scale=3.0, size=n)])
        probe = fhd.context_informativeness(predictions, truth,
                                            rng.normal(size=(n, 2)))
        assert probe["arm_pair_stability"] == 1.0
        assert probe["compared_arms"] == [0, 1]

    def test_an_explicit_pair_is_flagged_as_not_chosen_in_fold(self, fhd):
        """The old behaviour stays reachable, and is labelled as the old one."""
        predictions, truth, context = \
            fhd.counterexample_constant_weights_cannot_see_a_gate(200)
        probe = fhd.context_informativeness(predictions, truth, context,
                                            arm_pair=(0, 1))
        assert probe["arm_pair_chosen_inside_each_fold"] is False
        assert probe["arm_pair_stability"] == 1.0

    def test_the_leak_has_no_assumed_sign(self, fhd):
        """On null cohorts the whole-cohort pair biases rho DOWNWARD.

        This is why the module refuses to describe the leak as optimistic in
        general: the direction has to be measured on the data at hand.
        """
        gaps = []
        for seed in range(12):
            rng = np.random.default_rng(seed)
            n = 120
            truth = rng.normal(size=n)
            predictions = np.column_stack(
                [truth + rng.normal(scale=0.7, size=n) for _ in range(8)])
            context = rng.normal(size=(n, 4))
            order = np.abs(predictions - truth[:, None]).mean(axis=0).argsort()
            in_fold = fhd.context_informativeness(
                predictions, truth, context)["spearman_out_of_sample"]
            whole = fhd.context_informativeness(
                predictions, truth, context,
                arm_pair=(int(order[0]), int(order[1]))
            )["spearman_out_of_sample"]
            gaps.append(whole - in_fold)
        assert float(np.mean(gaps)) < 0.0


# --------------------------------------------------------------------------
# 2. the probe verdict ladder separates "above zero" from "big enough"
# --------------------------------------------------------------------------
class TestProbeVerdictLadder:
    def test_reaching_the_margin_needs_the_lower_bound_at_the_floor(self, fhd):
        assert fhd.probe_verdict([0.30, 0.70], [0.32, 0.68]) == \
            "SIGNAL_REACHES_MARGIN"

    def test_positive_but_entirely_below_the_floor_is_its_own_outcome(self, fhd):
        assert fhd.probe_verdict([0.08, 0.26], [0.05, 0.28]) == \
            "SIGNAL_POSITIVE_BUT_BELOW_MARGIN"

    def test_positive_and_straddling_the_floor_is_undecided(self, fhd):
        assert fhd.probe_verdict([0.10, 0.55], [0.05, 0.60]) == \
            "SIGNAL_POSITIVE_MARGIN_UNDECIDED"

    def test_nominal_only_keeps_its_own_two_names(self, fhd):
        assert fhd.probe_verdict([0.05, 0.20], [-0.02, 0.25]) == \
            "SIGNAL_POSITIVE_UNADJUSTED_ONLY_BELOW_MARGIN"
        assert fhd.probe_verdict([0.05, 0.60], [-0.02, 0.70]) == \
            "SIGNAL_POSITIVE_UNADJUSTED_ONLY_MARGIN_UNDECIDED"

    def test_refutation_still_needs_the_upper_bound_below_the_floor(self, fhd):
        assert fhd.probe_verdict([-0.10, 0.20], [-0.15, 0.25]) == \
            "REFUTED_CONTEXT_UNINFORMATIVE"

    def test_a_wide_null_interval_is_underpowered_not_refuted(self, fhd):
        assert fhd.probe_verdict([-0.20, 0.80], [-0.30, 0.90]) == \
            "NO_DETECTABLE_SIGNAL_UNDERPOWERED"

    @pytest.mark.parametrize("low,high", [(0.087, 0.366), (0.046, 0.428),
                                          (0.030, 0.427), (0.063, 0.431)])
    def test_the_values_actually_observed_do_not_reach_the_margin(self, fhd,
                                                                  low, high):
        """The reported associations on HLC, by their own numbers.

        None of them licenses a claim that the declared 0.30 association was
        met, which is exactly the overstatement the review caught.
        """
        assert fhd.probe_verdict([low, high], [low, high]) == \
            "SIGNAL_POSITIVE_MARGIN_UNDECIDED"

    def test_the_overstating_name_is_gone(self, fhd):
        source = (ROOT / "src" / "omics2pk"
                  / "fusion_headroom_diagnostics.py").read_text(
                      encoding="utf-8")
        assert '"CONTEXT_CARRIES_SIGNAL"' not in source
        assert '"CONTEXT_CARRIES_SIGNAL_UNADJUSTED_ONLY"' not in source

    def test_the_probe_reports_whether_it_reached_the_margin(self, fhd):
        predictions, truth, context = \
            fhd.counterexample_constant_weights_cannot_see_a_gate(200)
        probe = fhd.context_informativeness(predictions, truth, context)
        assert probe["reaches_minimum_interesting_association"] is True


# --------------------------------------------------------------------------
# 3. a non-linear basis exists, because a linear null refutes only linearity
# --------------------------------------------------------------------------
class TestNonLinearBasis:
    def test_the_quadratic_basis_adds_squares_and_products(self, fhd):
        block = np.array([[2.0, 3.0]])
        expanded = fhd.expand_basis(block, "quadratic")
        assert expanded.shape == (1, 5)
        assert list(expanded[0]) == [2.0, 3.0, 4.0, 6.0, 9.0]

    def test_the_linear_basis_is_the_identity(self, fhd):
        block = np.arange(6.0).reshape(3, 2)
        assert np.array_equal(fhd.expand_basis(block, "linear"), block)

    def test_an_unknown_basis_is_refused(self, fhd):
        with pytest.raises(fhd.FusionDiagnosticsError, match="unknown basis"):
            fhd.expand_basis(np.zeros((3, 2)), "cubic")

    def test_a_purely_quadratic_signal_is_missed_by_the_linear_basis(self, fhd):
        """The review's point, made executable.

        The advantage depends on the SQUARE of the context, so a linear probe
        sees nothing and a quadratic one sees it. A linear null therefore
        refutes linear gating, not contextual fusion.
        """
        rng = np.random.default_rng(5)
        n = 400
        context = rng.normal(size=(n, 1))
        side = (context[:, 0] ** 2 > 1.0)
        truth = rng.normal(size=n)
        arm_a = np.where(side, truth, truth + 1.0)
        arm_b = np.where(side, truth + 1.0, truth)
        predictions = np.column_stack([arm_a, arm_b])
        linear = fhd.context_informativeness(predictions, truth, context,
                                             basis="linear")
        quadratic = fhd.context_informativeness(predictions, truth, context,
                                                basis="quadratic")
        assert abs(linear["spearman_out_of_sample"]) < 0.2
        assert quadratic["spearman_out_of_sample"] > 0.5


# --------------------------------------------------------------------------
# 4. NOTHING crosses the outer split, and that is enforced, not asserted
# --------------------------------------------------------------------------
class TestNothingLeaksAcrossTheOuterSplit:
    def test_no_fit_ever_sees_a_row_it_will_predict(self, nfe, gateable):
        _, truth, context, arm_names = gateable
        arm_a, arm_b, _, _ = nfe.gateable_cohort()
        arms = {"arm_a": arm_a, "arm_b": arm_b}
        seen = []

        def spy(train_index, test_index):
            overlap = np.intersect1d(train_index, test_index)
            assert overlap.size == 0, (
                f"{overlap.size} rows are in both the training and the "
                "predicted set of a single fit")
            seen.append((set(train_index.tolist()), set(test_index.tolist())))
            return {name: values[test_index] for name, values in arms.items()}

        nfe.nested_fusion_evaluation(spy, truth, context, arm_names=arm_names)
        assert len(seen) == (nfe.FROZEN["outer_folds"]
                             * (nfe.FROZEN["inner_folds"] + 1))

    def test_an_arm_that_memorises_its_training_rows_gains_nothing(self, nfe):
        """The leak, encoded as a cohort that would expose it.

        `memoriser` returns the exact truth for any row it was trained on and
        a uselessly large constant otherwise. Under a correct nested design it
        is never trained on a row it predicts, so it is useless on every
        scored row and must attract essentially no weight. Under the previous
        cached-layer design the stacking weights were fitted on rows where
        such an arm looks perfect.

        The constant has to be far from the outcome's mean: a first version
        used 0.0, which on a standard-normal outcome is the mean predictor, so
        an eleven per cent weight on it was legitimate shrinkage rather than
        the leak the test was looking for.
        """
        rng = np.random.default_rng(13)
        n = 300
        truth = rng.normal(size=n)
        honest = truth + rng.normal(scale=0.4, size=n)
        context = rng.normal(size=(n, 2))

        def fit_predict(train_index, test_index):
            trained_on = set(train_index.tolist())
            memorised = np.array(
                [truth[i] if i in trained_on else 1.0e3 for i in test_index])
            return {"honest": honest[test_index], "memoriser": memorised}

        result = nfe.nested_fusion_evaluation(
            fit_predict, truth, context, arm_names=["honest", "memoriser"])
        assert result["mean_constant_weights"]["memoriser"] < 0.01
        assert result["mean_gate_weights"]["memoriser"] < 0.25
        assert result["mean_loss"]["K"] < 1.5 * result["mean_loss"]["S"]

    def test_the_outer_folds_are_shared_by_every_comparator(self, nfe,
                                                            gateable):
        fit_predict, truth, context, arm_names = gateable
        result = nfe.nested_fusion_evaluation(fit_predict, truth, context,
                                              arm_names=arm_names)
        sizes = {name: result["losses"][name].size
                 for name in nfe.COMPARATORS}
        assert set(sizes.values()) == {truth.size}, (
            "a paired contrast requires every comparator scored on every "
            "individual, on the same partition")

    def test_it_refuses_a_cohort_too_small_to_nest(self, nfe):
        rng = np.random.default_rng(1)
        truth = rng.normal(size=12)
        predictions = {"a": truth, "b": truth}
        with pytest.raises(nfe.NestedEvaluationError, match="too small"):
            nfe.nested_fusion_evaluation(
                lambda tr, te: {k: v[te] for k, v in predictions.items()},
                truth, rng.normal(size=(12, 2)), arm_names=["a", "b"])

    def test_it_refuses_a_misaligned_context(self, nfe, gateable):
        fit_predict, truth, _, arm_names = gateable
        with pytest.raises(nfe.NestedEvaluationError, match="disagree on n"):
            nfe.nested_fusion_evaluation(fit_predict, truth,
                                         np.zeros((truth.size - 1, 2)),
                                         arm_names=arm_names)

    def test_it_refuses_a_single_arm(self, nfe, gateable):
        fit_predict, truth, context, _ = gateable
        with pytest.raises(nfe.NestedEvaluationError, match="at least two"):
            nfe.nested_fusion_evaluation(fit_predict, truth, context,
                                         arm_names=["arm_a"])

    def test_a_missing_prediction_is_refused(self, nfe, gateable):
        _, truth, context, arm_names = gateable
        arm_a, arm_b, _, _ = nfe.gateable_cohort()

        def broken(train_index, test_index):
            del train_index
            values = arm_a[test_index].copy()
            values[0] = np.nan
            return {"arm_a": values, "arm_b": arm_b[test_index]}

        with pytest.raises(nfe.NestedEvaluationError, match="missing"):
            nfe.nested_fusion_evaluation(broken, truth, context,
                                         arm_names=arm_names)


# --------------------------------------------------------------------------
# 5. the harness can see a gate that exists — the positive control
# --------------------------------------------------------------------------
class TestPositiveControl:
    def test_the_self_test_passes(self, nfe):
        report = nfe.self_test()
        assert report["with_real_context"]["X_vs_K"] == "GAIN_REACHES_MARGIN"
        assert not report["with_permuted_context"]["X_vs_K"].startswith("GAIN_")

    def test_the_gate_nearly_eliminates_the_loss_when_context_is_perfect(
            self, nfe, gateable):
        fit_predict, truth, context, arm_names = gateable
        result = nfe.nested_fusion_evaluation(fit_predict, truth, context,
                                              arm_names=arm_names)
        assert result["mean_loss"]["X"] < 0.2 * result["mean_loss"]["K"]

    def test_stacking_cannot_help_on_a_cohort_only_a_gate_can_solve(
            self, nfe, gateable):
        fit_predict, truth, context, arm_names = gateable
        result = nfe.nested_fusion_evaluation(fit_predict, truth, context,
                                              arm_names=arm_names)
        assert result["mean_loss"]["K"] > 0.8 * result["mean_loss"]["S"], (
            "this cohort is built so constant weights cannot help; if K wins "
            "big here the cohort no longer isolates contextual gating")

    @pytest.mark.parametrize("mode", ["clinical", "predictions", "both"])
    def test_every_context_mode_finds_the_gate(self, nfe, gateable, mode):
        fit_predict, truth, context, arm_names = gateable
        result = nfe.nested_fusion_evaluation(
            fit_predict, truth, context, arm_names=arm_names,
            context_mode=mode)
        assert nfe.contrast(result, "X", "K")["verdict"] == "GAIN_REACHES_MARGIN"

    @pytest.mark.parametrize("mode", ["clinical", "predictions", "both"])
    def test_permuting_the_context_destroys_the_gain_in_every_mode(
            self, nfe, gateable, mode):
        fit_predict, truth, context, arm_names = gateable
        result = nfe.nested_fusion_evaluation(
            fit_predict, truth, context, arm_names=arm_names,
            context_mode=mode, permute_context=True)
        assert not nfe.contrast(result, "X", "K")["verdict"].startswith("GAIN_")

    def test_an_unknown_context_mode_is_refused(self, nfe, gateable):
        fit_predict, truth, context, arm_names = gateable
        with pytest.raises(nfe.NestedEvaluationError,
                           match="unknown context mode"):
            nfe.nested_fusion_evaluation(fit_predict, truth, context,
                                         arm_names=arm_names,
                                         context_mode="telepathy")


# --------------------------------------------------------------------------
# 6. the gate stays inside the family the convex ceiling bounds
# --------------------------------------------------------------------------
class TestTheGateStaysConvex:
    def test_the_gate_never_beats_the_per_individual_convex_ceiling(
            self, nfe, fhd, gateable):
        fit_predict, truth, context, arm_names = gateable
        result = nfe.nested_fusion_evaluation(fit_predict, truth, context,
                                              arm_names=arm_names)
        ceiling = fhd.individual_convex_oracle(
            result["outer_predictions"], truth)["mean_loss"]
        for name in nfe.COMPARATORS:
            assert result["mean_loss"][name] >= ceiling - 1e-9, (
                f"{name} beat the ceiling that bounds the whole family, which "
                "means either the ceiling or the comparator is wrong")

    def test_the_gate_weights_are_a_probability_vector(self, nfe):
        scores = np.array([[1.0, 2.0, 3.0], [0.0, 0.0, 0.0]])
        weights = nfe._softmax(scores)
        assert np.allclose(weights.sum(axis=1), 1.0)
        assert (weights >= 0.0).all()

    def test_a_huge_temperature_approaches_uniform_weights(self, nfe):
        losses = np.array([[0.1, 0.9], [0.8, 0.2], [0.5, 0.5]])
        design = np.array([[0.0], [1.0], [0.5]])
        weights, _ = nfe._gate_weights(design, losses, design, 1e6,
                                       penalty=1.0)
        assert np.allclose(weights, 0.5, atol=1e-3)


# --------------------------------------------------------------------------
# 7. the contrast ladder, which must not repeat the probe's mistake
# --------------------------------------------------------------------------
class TestContrastVerdictLadder:
    @staticmethod
    def _result(better, worse):
        return {"losses": {"X": np.asarray(better, dtype=float),
                           "K": np.asarray(worse, dtype=float)}}

    def test_a_real_but_tiny_gain_is_not_called_a_gain(self, nfe):
        rng = np.random.default_rng(21)
        worse = 100.0 + rng.normal(scale=0.01, size=600)
        better = worse - 0.5
        cell = nfe.contrast(self._result(better, worse), "X", "K")
        assert cell["verdict"] == "GAIN_PRESENT_BUT_BELOW_MARGIN", (
            "a 0.5 % gain against a declared 5 % margin is detectable and "
            "uninteresting, and must not be reported under the same name as "
            "an interesting one")

    def test_a_large_gain_is_called_a_gain(self, nfe):
        rng = np.random.default_rng(22)
        worse = 100.0 + rng.normal(scale=1.0, size=600)
        better = worse - 30.0
        cell = nfe.contrast(self._result(better, worse), "X", "K")
        assert cell["verdict"] == "GAIN_REACHES_MARGIN"

    def test_being_worse_is_a_named_result(self, nfe):
        rng = np.random.default_rng(23)
        worse = 100.0 + rng.normal(scale=1.0, size=600)
        better = worse + 20.0
        cell = nfe.contrast(self._result(better, worse), "X", "K")
        assert cell["verdict"] == "DETECTABLY_WORSE"

    def test_a_tight_null_is_a_refutation(self, nfe):
        rng = np.random.default_rng(24)
        worse = 100.0 + rng.normal(scale=0.05, size=800)
        better = worse + rng.normal(scale=0.05, size=800)
        cell = nfe.contrast(self._result(better, worse), "X", "K")
        assert cell["verdict"] == "REFUTED_GAIN_BELOW_MARGIN"

    def test_a_wide_null_is_not_a_refutation(self, nfe):
        rng = np.random.default_rng(25)
        worse = 100.0 + rng.normal(scale=60.0, size=25)
        better = worse + rng.normal(scale=60.0, size=25)
        cell = nfe.contrast(self._result(better, worse), "X", "K")
        assert cell["verdict"] == "NO_DETECTABLE_GAIN_UNDERPOWERED"
        assert cell["an_interval_containing_zero_is_not_a_refutation"] is True

    def test_an_unknown_comparator_is_refused(self, nfe):
        with pytest.raises(nfe.NestedEvaluationError, match="unknown"):
            nfe.contrast(self._result([1.0], [1.0]), "F", "K")

    def test_a_positive_gain_straddling_the_margin_is_undecided(self, nfe):
        rng = np.random.default_rng(26)
        # Reference mean 100 -> a 5 % margin of 5.0. The paired difference is
        # centred on 5.0 with a standard error near 1.0, so the interval
        # straddles the margin while excluding zero.
        worse = 100.0 + rng.normal(scale=1.0, size=400)
        better = worse - 5.0 + rng.normal(scale=20.0, size=400)
        cell = nfe.contrast(self._result(better, worse), "X", "K")
        assert cell["verdict"] == "GAIN_PRESENT_MARGIN_UNDECIDED"

    def test_the_contrast_ladder_mirrors_the_probe_ladder(self, nfe, fhd):
        """Both must separate below-margin, undecided and reaching-margin."""
        source = (ROOT / "src" / "omics2pk"
                  / "nested_fusion_evaluation.py").read_text(encoding="utf-8")
        for rung in ("GAIN_REACHES_MARGIN", "GAIN_PRESENT_BUT_BELOW_MARGIN",
                     "GAIN_PRESENT_MARGIN_UNDECIDED"):
            assert f'"{rung}"' in source
        assert '"GAIN_PRESENT"' not in source, (
            "a bare GAIN_PRESENT invites exactly the overstatement the "
            "review caught in the probe")

    def test_the_margin_matches_the_headroom_diagnostics(self, nfe, fhd):
        assert (nfe.FROZEN["minimum_interesting_relative_gain"]
                == fhd.FROZEN["minimum_interesting_relative_gain"])


# --------------------------------------------------------------------------
# 8. what the runner must declare out loud
# --------------------------------------------------------------------------
class TestRunnerDeclarations:
    @pytest.fixture(scope="class")
    def source(self):
        return (ROOT / "scripts"
                / "run_nested_fusion_evaluation.py").read_text(
                    encoding="utf-8")

    def test_the_endpoint_names_are_the_cohort_s_own(self, source):
        for name in ("CYP3A4T_activity", "CYP2D6_activity",
                     "CYP2C9_activity", "CYP2E1_activity",
                     "CYP3A4M_activity"):
            assert name in source

    def test_a_mistyped_endpoint_name_cannot_pass_silently(self, source):
        assert "declared endpoints absent from the cohort" in source

    def test_the_selection_caveat_is_stated(self, source):
        assert "prespecification_caveat" in source
        assert "selected_on_these_data" in source

    def test_both_controls_are_run(self, source):
        assert "negative_control_all_arms_at_baseline" in source
        assert "substrate_specificity_control" in source

    def test_the_decoy_and_permutation_controls_are_run(self, source):
        assert "decoy_arm_control" in source
        assert "permuted_context_control" in source

    def test_the_grid_sensitivity_exists(self, source):
        assert "SENSITIVITY_TEMPERATURE_GRID" in source
        assert "grid_sensitivity" in source

    def test_q2_is_declared_untested(self, source):
        assert "Q2_mechanistic_prior" in source

    def test_the_artefact_carries_no_individual_level_data(self, source):
        assert "no_individual_level_data_in_this_artefact" in source
        assert 'k not in ("losses", "outer_predictions")' in source


# --------------------------------------------------------------------------
# 9. the results document must keep saying what the run measured
# --------------------------------------------------------------------------
RESULTS = ROOT / "docs" / "OMICS2PK_NESTED_FUSION_Q1_RESULTS_20261111.md"
HEADROOM = ROOT / "docs" / "OMICS2PK_FUSION_HEADROOM_RESULTS_20261111.md"


class TestResultsDocument:
    @pytest.fixture(scope="class")
    def text(self):
        return RESULTS.read_text(encoding="utf-8")

    def test_it_exists_and_is_versioned(self, text):
        assert "**Version 1.0.0 \u2014 2026-11-11.**" in text

    def test_it_states_the_primary_answer_as_a_count(self, text):
        assert "0** (25 refuted, 5 underpowered)" in text

    def test_it_reports_the_null_surviving_the_wider_grid(self, text):
        assert "widened temperature grid" in text
        assert "not a grid artefact" in text

    def test_it_separates_the_four_defects_by_consequence(self, text):
        assert "negligible here" in text
        assert "this is the one that mattered" in text

    def test_it_quantifies_the_cached_layer_leak(self, text):
        assert "+167.5" in text and "+58.2" in text

    def test_it_records_my_own_contradicted_claim(self, text):
        assert "struck, not deleted" in text
        assert "no reliable sign" in text

    def test_it_admits_reproducing_the_same_error_in_new_code(self, text):
        assert "The same error, found in my own new code" in text
        assert "+0.62 %" in text

    def test_it_names_the_gate_s_own_defect_as_a_limitation(self, text):
        assert "limitation of the gate I built, not a measurement of" in text

    def test_it_reports_the_positive_control(self, text):
        assert "91.7 %" in text
        assert "a null" in text and "is a statement about HLC" in text

    def test_it_reports_that_the_gate_chose_differently_and_was_not_better(
            self, text):
        assert "34 % to 76 %" in text
        assert "they are not better choices" in text

    def test_it_keeps_the_unrepairable_defect_visible(self, text):
        assert "not repairable, only labelled" in text
        assert "would not have been admissible as confirmation" in text

    def test_it_refuses_to_open_q2(self, text):
        assert "Q2 should not be opened on this cohort" in text
        # the phrase wraps in the document, so match the unwrapped fragment
        assert "prior's flexibility, not its knowledge" in text

    def test_it_lists_what_the_null_does_not_cover(self, text):
        for item in ("context this cohort does not publish", "Q2",
                     "clinical PK", "larger n"):
            assert item in text

    def test_it_states_the_absence_of_gain_is_a_valid_result(self, text):
        assert "absence of gain is a valid scientific result" in text

    def test_the_headroom_document_records_its_withdrawal(self):
        text = HEADROOM.read_text(encoding="utf-8")
        assert "**1.1.0**" in text
        assert "WITHDRAWN" in text
        assert "2.9" in text
        assert "SUPERSEDED, see \u00a70" in text

    def test_the_headroom_document_keeps_what_still_stands(self):
        text = HEADROOM.read_text(encoding="utf-8")
        assert "headroom present on 13 of 13" in text.lower()
        assert "**stands**" in text
