"""Tests for the cohort-size simulation.

A sizing module is dangerous in a way an analysis module is not: its output is a
single number that reads like a fact. These tests check that the number cannot
be produced without its assumptions, that the forbidden route from a rank
correlation to a variance fraction is absent, and that the pipeline being sized
really is the pipeline the project runs -- selection and penalty inside the
training folds.
"""
from __future__ import annotations

import json
import math
import pathlib
import re
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PRESPEC = ROOT / "docs" / "OMICS2PK_COHORT_SIZING_PRESPECIFICATION_20261110.md"
RESULTS = ROOT / "docs" / "OMICS2PK_COHORT_SIZING_RESULTS_20261110.md"
ARTEFACT = ROOT / "reports" / "data" / "cohort_size_simulation.json"


@pytest.fixture(scope="module")
def mod():
    from omics2pk import cohort_size_simulation
    return cohort_size_simulation


@pytest.fixture(scope="module")
def payload():
    if not ARTEFACT.is_file():
        pytest.skip("the simulation has not been run in this checkout")
    return json.loads(ARTEFACT.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------
# 1. the one permitted route from a published rank correlation
# --------------------------------------------------------------------------
class TestSpearmanToPearson:
    def test_zero_maps_to_zero(self, mod):
        assert mod.spearman_to_pearson(0.0) == pytest.approx(0.0)

    def test_one_maps_to_one(self, mod):
        assert mod.spearman_to_pearson(1.0) == pytest.approx(1.0)

    def test_it_is_the_gaussian_copula_formula(self, mod):
        for rho in (0.2, 0.53, 0.9):
            assert mod.spearman_to_pearson(rho) == pytest.approx(
                2.0 * math.sin(math.pi * rho / 6.0))

    def test_it_is_monotone(self, mod):
        values = [mod.spearman_to_pearson(r) for r in (0.1, 0.3, 0.6, 0.9)]
        assert values == sorted(values)

    def test_it_is_above_the_rank_value_but_not_by_much(self, mod):
        # The copula mapping lifts a rank correlation slightly; a reader who
        # mistook one for the other would be wrong by a knowable amount.
        for rho in (0.2, 0.5, 0.8):
            lifted = mod.spearman_to_pearson(rho)
            assert rho < lifted < rho + 0.05

    def test_an_impossible_correlation_is_refused(self, mod):
        with pytest.raises(mod.CohortSizeSimulationError, match=r"\[-1, 1\]"):
            mod.spearman_to_pearson(1.4)

    def test_the_docstring_records_that_it_is_an_assumption(self, mod):
        text = mod.spearman_to_pearson.__doc__
        assert "ASSUMPTION" in text
        assert "not permitted" in text

    def test_the_module_never_squares_a_rank_correlation(self):
        source = (ROOT / "src" / "omics2pk"
                  / "cohort_size_simulation.py").read_text(encoding="utf-8")
        for pattern in (r"rho_spearman\s*\*\*\s*2", r"rho_s\s*\*\*\s*2",
                        r"spearman\w*\s*\*\*\s*2"):
            assert not re.search(pattern, source), (
                "a squared rank correlation appears in the sizing module")


# --------------------------------------------------------------------------
# 2. the pipeline being sized is the project's pipeline
# --------------------------------------------------------------------------
class TestFoldSafety:
    def test_selection_uses_only_training_rows(self, mod):
        rng = np.random.default_rng(1)
        # A block whose only informative column is informative ONLY for the
        # held-out rows cannot help, if selection is fold-safe.
        n = 60
        target = rng.normal(size=n)
        block = rng.normal(size=(n, 20))
        block[n // 2:, 0] = target[n // 2:] * 3.0
        errors = mod.cross_validated_absolute_errors(
            block, target, q=3, rng=rng)
        baseline = mod.cross_validated_absolute_errors(
            np.zeros((n, 0)), target, q=0, rng=rng)
        # It may help a little by luck, but must not dominate.
        assert errors.mean() > 0.4 * baseline.mean()

    def test_a_real_signal_is_recovered(self, mod):
        rng = np.random.default_rng(2)
        n = 120
        signal = rng.normal(size=n)
        target = signal + 0.3 * rng.normal(size=n)
        block = np.column_stack([signal + 0.1 * rng.normal(size=n),
                                 rng.normal(size=(n, 19))])
        errors = mod.cross_validated_absolute_errors(
            block, target, q=3, rng=rng)
        baseline = mod.cross_validated_absolute_errors(
            np.zeros((n, 0)), target, q=0, rng=rng)
        assert errors.mean() < 0.6 * baseline.mean()

    def test_an_empty_block_is_the_mean_predictor(self, mod):
        rng = np.random.default_rng(3)
        target = rng.normal(size=40) + 5.0
        errors = mod.cross_validated_absolute_errors(
            np.zeros((40, 0)), target, q=0, rng=rng)
        assert errors.mean() == pytest.approx(
            np.abs(target - target.mean()).mean(), rel=0.25)


# --------------------------------------------------------------------------
# 3. detection requires the right sign
# --------------------------------------------------------------------------
class TestPairedDetection:
    def test_a_clear_improvement_is_detected(self, mod):
        baseline = np.full(60, 1.0)
        arm = baseline - 0.3
        assert mod.paired_detection(arm, baseline, alpha=0.05,
                                    family_size=1) is True

    def test_a_clear_degradation_is_not_detected(self, mod):
        baseline = np.full(60, 1.0)
        arm = baseline + 0.3
        assert mod.paired_detection(arm, baseline, alpha=0.05,
                                    family_size=1) is False

    def test_no_difference_is_not_detected(self, mod):
        rng = np.random.default_rng(4)
        baseline = rng.normal(size=80) + 3.0
        assert mod.paired_detection(baseline, baseline.copy(), alpha=0.05,
                                    family_size=1) is False

    def test_a_larger_family_makes_detection_harder(self, mod):
        rng = np.random.default_rng(5)
        baseline = rng.normal(size=40) + 2.0
        arm = baseline - 0.12
        lenient = mod.paired_detection(arm, baseline, alpha=0.05,
                                       family_size=1)
        strict = mod.paired_detection(arm, baseline, alpha=0.05,
                                      family_size=50)
        assert lenient >= strict

    def test_the_normal_quantile_is_right(self, mod):
        assert mod._two_sided_normal_quantile(0.05) == pytest.approx(
            1.959964, abs=1e-4)
        assert mod._two_sided_normal_quantile(0.01) == pytest.approx(
            2.575829, abs=1e-4)


# --------------------------------------------------------------------------
# 4. power and the search
# --------------------------------------------------------------------------
class TestPowerAndSearch:
    def test_power_rises_with_n_for_a_real_effect(self, mod):
        scenario = {"rho_pearson_log": 0.7, "analytical_cv_log": 0.1,
                    "missing_rate": 0.0}
        small = mod.power_at(15, mod.simulate_enzyme_activity, scenario,
                             replicates=60)
        large = mod.power_at(120, mod.simulate_enzyme_activity, scenario,
                             replicates=60)
        assert large > small

    def test_power_is_near_the_nominal_level_under_the_null(self, mod):
        scenario = {"rho_pearson_log": 0.0, "analytical_cv_log": 0.1,
                    "missing_rate": 0.0}
        observed = mod.power_at(80, mod.simulate_enzyme_activity, scenario,
                                replicates=200)
        # A one-sided detection at alpha=0.05 two-sided is nominally 0.025.
        assert observed < 0.12, (
            "the detection rule rejects far too often with no effect, so any "
            "required n it produces would be meaningless")

    def test_a_stronger_correlation_needs_fewer_individuals(self, mod):
        weak = mod.required_sample_size(
            mod.simulate_enzyme_activity,
            {"rho_pearson_log": 0.3, "analytical_cv_log": 0.1,
             "missing_rate": 0.0}, replicates=40)
        strong = mod.required_sample_size(
            mod.simulate_enzyme_activity,
            {"rho_pearson_log": 0.8, "analytical_cv_log": 0.1,
             "missing_rate": 0.0}, replicates=40)
        assert strong["required_n"] is not None
        assert weak["required_n"] is None or (
            strong["required_n"] < weak["required_n"])

    def test_an_unreachable_scenario_returns_none_and_says_so(self, mod):
        out = mod.required_sample_size(
            mod.simulate_enzyme_activity,
            {"rho_pearson_log": 0.0, "analytical_cv_log": 0.1,
             "missing_rate": 0.0}, replicates=20, high=40)
        assert out["required_n"] is None
        assert "not answerable by enlarging the cohort" in out["verdict"]

    def test_the_search_records_every_n_it_evaluated(self, mod):
        out = mod.required_sample_size(
            mod.simulate_enzyme_activity,
            {"rho_pearson_log": 0.8, "analytical_cv_log": 0.1,
             "missing_rate": 0.0}, replicates=30)
        assert len(out["searched"]) >= 2
        assert all(0.0 <= v <= 1.0 for v in out["searched"].values())


# --------------------------------------------------------------------------
# 5. redundancy and noise must be able to destroy a gain
# --------------------------------------------------------------------------
class TestIncrementalLayer:
    def test_full_redundancy_leaves_nothing_to_add(self, mod):
        rng = np.random.default_rng(6)
        detections = sum(
            mod.simulate_incremental_layer(
                150, baseline_r2=0.3, target_mae_reduction=0.3,
                n_features=20, n_informative=5, redundancy=0.999,
                analytical_cv=0.1, missing_rate=0.0, selected=5,
                rng=rng, alpha=0.05, family_size=1)
            for _ in range(60))
        assert detections / 60 < 0.5, (
            "a layer that duplicates the baseline is being credited with a "
            "gain, so the redundancy parameter is not doing its job")

    def test_a_zero_reduction_is_detected_at_about_the_nominal_rate(self, mod):
        rng = np.random.default_rng(8)
        detections = sum(
            mod.simulate_incremental_layer(
                120, baseline_r2=0.3, target_mae_reduction=0.0,
                n_features=50, n_informative=5, redundancy=0.0,
                analytical_cv=0.2, missing_rate=0.0, selected=5,
                rng=rng, alpha=0.05, family_size=1)
            for _ in range(120))
        assert detections / 120 < 0.15

    def test_a_large_clean_effect_is_detected(self, mod):
        rng = np.random.default_rng(9)
        detections = sum(
            mod.simulate_incremental_layer(
                200, baseline_r2=0.2, target_mae_reduction=0.35,
                n_features=20, n_informative=5, redundancy=0.0,
                analytical_cv=0.05, missing_rate=0.0, selected=5,
                rng=rng, alpha=0.05, family_size=1)
            for _ in range(60))
        assert detections / 60 > 0.5

    def test_bad_parameters_are_refused(self, mod):
        rng = np.random.default_rng(10)
        common = dict(n_features=20, n_informative=5, redundancy=0.0,
                      analytical_cv=0.1, missing_rate=0.0, selected=5,
                      rng=rng, alpha=0.05, family_size=1)
        with pytest.raises(mod.CohortSizeSimulationError, match="baseline_r2"):
            mod.simulate_incremental_layer(50, baseline_r2=1.0,
                                           target_mae_reduction=0.1, **common)
        with pytest.raises(mod.CohortSizeSimulationError, match="fraction"):
            mod.simulate_incremental_layer(50, baseline_r2=0.3,
                                           target_mae_reduction=1.0, **common)

    def test_more_informative_features_than_features_is_refused(self, mod):
        rng = np.random.default_rng(11)
        with pytest.raises(mod.CohortSizeSimulationError,
                           match="non-empty subset"):
            mod.simulate_incremental_layer(
                50, baseline_r2=0.3, target_mae_reduction=0.1,
                n_features=5, n_informative=50, redundancy=0.0,
                analytical_cv=0.1, missing_rate=0.0, selected=5,
                rng=rng, alpha=0.05, family_size=1)

    def test_the_optimistic_simplifications_are_declared(self, mod):
        text = mod.simulate_incremental_layer.__doc__
        assert "declared scenario, not" in text
        source = (ROOT / "src" / "omics2pk"
                  / "cohort_size_simulation.py").read_text(encoding="utf-8")
        assert "OPTIMISTIC for the layer" in source


# --------------------------------------------------------------------------
# 6. confirmation is not discovery
# --------------------------------------------------------------------------
class TestValidationCohortSize:
    def test_a_bigger_effect_needs_fewer_individuals(self, mod):
        small = mod.validation_cohort_size(0.03, 0.2)["required_n"]
        large = mod.validation_cohort_size(0.10, 0.2)["required_n"]
        assert large < small

    def test_it_scales_as_the_inverse_square_of_the_effect(self, mod):
        one = mod.validation_cohort_size(0.05, 0.2)["required_n"]
        half = mod.validation_cohort_size(0.025, 0.2)["required_n"]
        assert half / one == pytest.approx(4.0, rel=0.05)

    def test_a_multiplicity_family_increases_it(self, mod):
        alone = mod.validation_cohort_size(0.05, 0.2, family_size=1)[
            "required_n"]
        among = mod.validation_cohort_size(0.05, 0.2, family_size=10)[
            "required_n"]
        assert among > alone

    def test_it_says_it_is_not_a_discovery_size(self, mod):
        out = mod.validation_cohort_size(0.05, 0.2)
        assert "not a discovery size" in out["what_this_is_not"]

    def test_it_says_the_dispersion_from_a_small_sample_is_biased_low(self, mod):
        out = mod.validation_cohort_size(0.05, 0.2)
        assert "biased low" in out["what_this_assumes"]
        assert "optimistic" in out["what_this_assumes"]

    def test_a_zero_effect_cannot_be_sized(self, mod):
        with pytest.raises(mod.CohortSizeSimulationError,
                           match="zero effect cannot be sized"):
            mod.validation_cohort_size(0.0, 0.2)

    def test_a_non_positive_dispersion_is_refused(self, mod):
        with pytest.raises(mod.CohortSizeSimulationError,
                           match="dispersion must be positive"):
            mod.validation_cohort_size(0.05, 0.0)


# --------------------------------------------------------------------------
# 7. sensitivity is mandatory, not decorative
# --------------------------------------------------------------------------
class TestSensitivity:
    def test_it_ranks_the_variations(self, mod):
        out = mod.one_at_a_time_sensitivity(
            mod.simulate_enzyme_activity,
            {"rho_pearson_log": 0.7, "analytical_cv_log": 0.1,
             "missing_rate": 0.0},
            {"rho_pearson_log": [0.3, 0.9]},
            replicates=30)
        assert len(out["ordered_by_required_n_descending"]) == 2
        assert "reference" in out and "variations" in out

    def test_it_explains_how_to_read_the_ranking(self, mod):
        out = mod.one_at_a_time_sensitivity(
            mod.simulate_enzyme_activity,
            {"rho_pearson_log": 0.8, "analytical_cv_log": 0.1,
             "missing_rate": 0.0},
            {"analytical_cv_log": [0.4]}, replicates=20)
        assert "widest range of required" in out["how_to_read"]


# --------------------------------------------------------------------------
# 8. the module refuses to let a required n travel without its assumptions
# --------------------------------------------------------------------------
class TestFraming:
    def test_the_reading_note_names_the_three_assumptions(self, mod):
        text = mod.how_to_read_any_required_n()
        for word in ("effect size", "redundancy", "noise", "not about biology"):
            assert word in text

    def test_the_validation_scheme_difference_is_declared(self, mod):
        assert "leave-one-donor-out" in mod.VALIDATION_NOTE
        assert "artefact of the scheme" in mod.VALIDATION_NOTE

    def test_the_two_evidence_levels_are_frozen(self, mod):
        assert mod.FROZEN["evidence_levels"] == (
            "transferable_literature", "declared_scenario")


# --------------------------------------------------------------------------
# 9. the executed artefact and its document
# --------------------------------------------------------------------------
class TestArtefact:
    def test_the_prespecification_exists(self):
        assert PRESPEC.is_file()

    def test_question_a_is_level_one_and_b_and_c_are_level_two(self, payload):
        assert payload["question_a_enzyme_activity"]["evidence_level"] == \
            "transferable_literature"
        assert payload["question_b_multiomic_complementarity"][
            "evidence_level"] == "declared_scenario"
        assert payload["question_c_clinical_pk"]["evidence_level"] == \
            "declared_scenario"

    def test_no_in_vitro_effect_is_transferred_to_the_clinical_question(
            self, payload):
        text = payload["question_c_clinical_pk"]["refusal"]
        assert "no effect measured on in-vitro" in text
        assert "is transferred to a clinical" in text
        assert "in-vitro question only" in text

    def test_the_published_inputs_are_named_with_their_source(self, payload):
        block = payload["question_a_enzyme_activity"]["published_inputs"]
        assert "PMC12735975" in block["source"]
        assert "never read as performance" in block["why_these_are_transferable"]
        assert "cannot validate anything" in block["why_they_are_not_a_result"]

    def test_the_three_clinical_endpoints_are_sized_separately(self, payload):
        block = payload["question_c_clinical_pk"]["by_endpoint"]
        assert set(block) == {"single_trough_concentration",
                              "auc_from_dense_sampling",
                              "clearance_from_a_population_model"}
        for endpoint, cell in block.items():
            assert "why_this_endpoint_differs" in cell

    def test_it_holds_no_individual_level_data(self, payload):
        assert payload["no_individual_level_data_in_this_artefact"] is True

    def test_the_confirmatory_sizing_is_separate_from_the_simulation(
            self, payload):
        block = payload["question_a_enzyme_activity"][
            "measured_confirmatory_sizing"]
        for label, cell in block.items():
            assert cell["n_donors_observed"] == 15
            assert "not a discovery size" in cell["confirmatory_sizing"][
                "what_this_is_not"]
