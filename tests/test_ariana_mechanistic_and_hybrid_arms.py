"""Contract tests for the ARIANA mechanistic (M4/M5), PopPK (P1–P3) and HY arms.

Offline: constructed data with known truth, pinned provenance constants, and
the governance rules the analysis must not break. The source workbook is never
touched and no network call is made.

The rules these tests exist above all to defend:

* ``lnMR OH-metoprolol/metoprolol`` shares its denominator with the target
  concentration and must appear in **no** design, in any arm;
* the key secondary must run on the **same folds** as the primary contrast,
  which is only true while the primary runner's frozen constants are unchanged;
* every mechanistic coefficient must be transferred from a named publication,
  never estimated on ARIANA, and every arm must declare how many parameters it
  does fit on the cohort;
* an activity score absent from the published Rüdesheim grid must be reported
  as interpolated, not silently aligned.
"""
from __future__ import annotations

import importlib.util
import json
import re
import unittest
from pathlib import Path

try:  # The `offline` CI job installs no dependencies, so this must skip
    # rather than error during unittest discovery.
    import numpy as np
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(
        f"numpy and pytest are required for this module: {exc}"
    ) from exc

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "src" / "omics2pk" / "ariana_mechanistic_arms.py"
RUNNER = ROOT / "scripts" / "run_h1_ariana_mechanistic_and_hybrid_arms.py"
REPORT = ROOT / "reports" / "data" / "h1_ariana_mechanistic_hybrid.json"
PROHIBITED = "lnMR OH-metoprolol/metoprolol"


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mech():
    return _load(MODULE, "ariana_mechanistic_arms")


@pytest.fixture(scope="module")
def runner():
    try:
        return _load(RUNNER, "run_h1_mechanistic")
    except ImportError as exc:  # pragma: no cover - needs sklearn/xgboost
        pytest.skip(f"the runner needs the modelling stack: {exc}")


@pytest.fixture(scope="module")
def report():
    if not REPORT.is_file():
        pytest.skip("the run artefact is not present in this checkout")
    return json.loads(REPORT.read_text(encoding="utf-8"))


def _flat(path):
    """Source text with runs of whitespace collapsed.

    A phrase wrapped across two source lines is still the same phrase, so a
    test that greps the source must not depend on where the line breaks fell.
    """
    return re.sub(r"\s+", " ", path.read_text(encoding="utf-8"))


# --- the prohibited predictor --------------------------------------------
class TestProhibitedPredictor:
    def test_absent_from_every_hybrid_residual_design(self, runner):
        for name, features in runner.HY_RESIDUAL_DESIGNS.items():
            assert PROHIBITED not in features, name

    def test_absent_from_every_omic_panel(self, mech):
        for name, columns in mech.OMIC_PANELS.items():
            assert PROHIBITED not in columns, name

    def test_preflight_passes_as_shipped(self, runner):
        runner.assert_no_prohibited_predictor()

    def test_preflight_catches_an_injected_design(self, runner):
        original = dict(runner.HY_RESIDUAL_DESIGNS)
        runner.HY_RESIDUAL_DESIGNS["injected"] = [PROHIBITED]
        try:
            with pytest.raises(runner.H1.AnalysisError):
                runner.assert_no_prohibited_predictor()
        finally:
            runner.HY_RESIDUAL_DESIGNS.clear()
            runner.HY_RESIDUAL_DESIGNS.update(original)

    def test_preflight_catches_an_injected_panel(self, runner):
        # Patch the module instance the RUNNER holds: the `mech` fixture loads
        # the file a second time, so patching that object would prove nothing.
        original = dict(runner.MECH.OMIC_PANELS)
        patched = dict(original)
        patched["injected"] = (PROHIBITED,)
        runner.MECH.OMIC_PANELS = patched
        try:
            with pytest.raises(runner.H1.AnalysisError):
                runner.assert_no_prohibited_predictor()
        finally:
            runner.MECH.OMIC_PANELS = original

    def test_the_runner_says_why_it_is_prohibited(self):
        text = _flat(RUNNER)
        assert "shares its denominator" in text


# --- fold identity with the primary contrast -----------------------------
class TestFoldIdentity:
    def test_folds_come_from_the_primary_runner(self, runner):
        folds = runner.H1.make_folds(47, runner.H1.OUTER_FOLDS,
                                     runner.H1.OUTER_REPEATS, runner.H1.SEED)
        assert len(folds) == runner.H1.OUTER_FOLDS * runner.H1.OUTER_REPEATS
        covered = set()
        for _, test in folds[: runner.H1.OUTER_FOLDS]:
            covered |= set(int(i) for i in test)
        assert covered == set(range(47))

    def test_train_and_test_are_disjoint(self, runner):
        for train, test in runner.H1.make_folds(47, 5, 2, runner.H1.SEED):
            assert not (set(int(i) for i in train) & set(int(i) for i in test))

    def test_machinery_check_passes_as_shipped(self, runner):
        runner.check_primary_machinery_unchanged()

    @pytest.mark.parametrize(
        "name,value",
        [("SEED", 1), ("OUTER_FOLDS", 4), ("OUTER_REPEATS", 1),
         ("INNER_FOLDS", 5), ("TRIM_FRACTION", 0.2),
         ("RELATIVE_IMPROVEMENT_THRESHOLD", 0.05)],
    )
    def test_machinery_check_catches_a_changed_constant(self, runner, name, value):
        original = getattr(runner.H1, name)
        setattr(runner.H1, name, value)
        try:
            with pytest.raises(runner.H1.AnalysisError):
                runner.check_primary_machinery_unchanged()
        finally:
            setattr(runner.H1, name, original)

    def test_fingerprint_is_stable_and_sensitive(self, runner):
        folds = runner.H1.make_folds(47, 5, 2, runner.H1.SEED)
        assert runner.fold_fingerprint(folds) == runner.fold_fingerprint(folds)
        other = runner.H1.make_folds(47, 5, 2, runner.H1.SEED + 1)
        assert runner.fold_fingerprint(folds) != runner.fold_fingerprint(other)


# --- pinned provenance ----------------------------------------------------
class TestTransferredConstants:
    def test_rudesheim_grid_is_the_published_one(self, mech):
        assert mech.RUDESHEIM_AS_GRID == (0.0, 0.5, 1.25, 1.5, 2.0, 3.0)
        assert mech.RUDESHEIM_RELATIVE_ACTIVITY == (0.0, 0.19, 0.64, 0.72,
                                                    1.00, 2.13)

    def test_fraction_metabolised(self, mech):
        assert mech.FM_CYP2D6 == 0.80

    def test_phenoconversion_decrement(self, mech):
        assert mech.DELTA_AS_PER_INHIBITOR == 0.54

    def test_sex_effect(self, mech):
        assert mech.SEX_LOG_EFFECT == 0.572

    def test_allometric_exponents(self, mech):
        assert mech.ALLOMETRIC_EXPONENT_CL == 0.75
        assert mech.ALLOMETRIC_EXPONENT_V == 1.0

    def test_omic_exponent_and_weight(self, mech):
        assert mech.OMIC_GAMMA == 1.0
        assert mech.BLEND_WEIGHT == 0.5

    def test_every_constant_names_a_source_and_is_not_fitted(self, mech):
        assert len(mech.LITERATURE_CONSTANTS) >= 6
        for module, coefficient, source in mech.LITERATURE_CONSTANTS:
            assert module and coefficient and source
            assert re.search(r"doi|PMC|FDA|Anderson", source), module

    def test_each_backbone_names_its_source(self, mech):
        for key, backbone in mech.BACKBONES.items():
            assert backbone["source"], key
            assert backbone["structure"] in {"1C", "2C"}, key

    def test_the_weight_is_deliberately_not_the_optimum(self, mech):
        text = _flat(MODULE)
        assert "optimum near 0.75" in text or "NOT the optimum" in text


class TestPrespecifiedChoices:
    def test_primary_panel_is_the_one_with_an_authentic_standard(self, mech):
        assert mech.PRIMARY_OMIC_PANEL == "oh_only"
        assert mech.OMIC_PANELS["oh_only"] == (
            "lnMR OH-solanidine/solanidine",
        )

    def test_the_result_guided_panel_is_demoted_to_a_sensitivity(self, mech):
        text = _flat(MODULE)
        assert "on the observed loss" in text
        assert "oh_and_ssda" in mech.OMIC_PANELS

    def test_primary_backbone_is_the_extended_release_scenario(self, mech):
        assert mech.PRIMARY_BACKBONE == "FDA2006_ER"
        assert mech.BACKBONES["FDA2006_ER"]["role"] == "primary scenario"

    def test_the_other_backbones_are_sensitivities(self, mech):
        for key in ("Chai2026_IR", "Eugene2016_IR"):
            assert mech.BACKBONES[key]["role"] == "structural sensitivity"

    def test_key_secondary_is_m5_against_m4(self, mech):
        assert mech.KEY_SECONDARY_CONTRAST == ("M5_mechanistic_pgx_omic",
                                               "M4_mechanistic_pgx")

    def test_free_parameters_are_declared_for_every_arm(self, mech):
        for arm in mech.MECHANISTIC_ARMS:
            assert arm in mech.FREE_PARAMETERS
        assert mech.FREE_PARAMETERS["M4_mechanistic_pgx"] == 1
        assert mech.FREE_PARAMETERS["M5_mechanistic_pgx_omic"] == 1
        assert mech.FREE_PARAMETERS["P1_poppk_published_coefficients"] == 1
        assert mech.FREE_PARAMETERS["P2_poppk_estimated_coefficients"] == 3
        assert mech.FREE_PARAMETERS["P3_poppk_estimated_plus_omic"] == 4

    def test_hybrid_primary_residual_design_is_declared(self, runner):
        assert runner.HY_PRIMARY_RESIDUAL_DESIGN == "clinical_pgx_omic"
        assert runner.HY_PRIMARY_RESIDUAL_DESIGN in runner.HY_RESIDUAL_DESIGNS


# --- the genotype module -------------------------------------------------
class TestGenotypeModule:
    def test_grid_points_return_the_published_values(self, mech):
        got = mech.relative_activity_from_activity_score(mech.RUDESHEIM_AS_GRID)
        assert np.allclose(got, mech.RUDESHEIM_RELATIVE_ACTIVITY)

    def test_interpolation_is_monotone_between_grid_points(self, mech):
        scores = np.linspace(0, 3, 61)
        activity = mech.relative_activity_from_activity_score(scores)
        assert np.all(np.diff(activity) >= -1e-12)

    def test_beyond_the_grid_is_capped_not_extrapolated(self, mech):
        assert mech.relative_activity_from_activity_score([5.0])[0] == pytest.approx(2.13)

    def test_a_negative_score_is_refused(self, mech):
        with pytest.raises(mech.MechanisticError):
            mech.relative_activity_from_activity_score([-0.5])

    def test_a_non_finite_score_is_refused(self, mech):
        with pytest.raises(mech.MechanisticError):
            mech.relative_activity_from_activity_score([np.nan])

    def test_poor_metaboliser_multiplier_is_a_fivefold_exposure(self, mech):
        multiplier = mech.clearance_multiplier([0.0])[0]
        assert multiplier == pytest.approx(0.20)
        assert 1.0 / multiplier == pytest.approx(5.0)

    def test_normal_metaboliser_multiplier_is_one(self, mech):
        assert mech.clearance_multiplier([1.0])[0] == pytest.approx(1.0)

    def test_fraction_metabolised_out_of_domain_is_refused(self, mech):
        for bad in (0.0, 1.0, -0.1, 1.5):
            with pytest.raises(mech.MechanisticError):
                mech.clearance_multiplier([1.0], fm=bad)

    def test_negative_activity_is_refused(self, mech):
        with pytest.raises(mech.MechanisticError):
            mech.clearance_multiplier([-1.0])

    def test_phenoconversion_subtracts_only_for_inhibitor_users(self, mech):
        got = mech.effective_activity_score([2.0, 2.0], [0, 1])
        assert got[0] == pytest.approx(2.0)
        assert got[1] == pytest.approx(2.0 - 0.54)

    def test_phenoconversion_floors_at_zero(self, mech):
        assert mech.effective_activity_score([0.25], [1])[0] == 0.0


class TestActivityScoreAudit:
    def test_an_on_grid_cohort_reports_no_interpolation(self, mech):
        audit = mech.activity_score_grid_audit([0.0, 0.5, 2.0, 3.0])
        assert audit["n_patients_interpolated"] == 0
        assert audit["scores_absent_from_the_published_grid"] == []

    def test_off_grid_scores_are_named_and_counted(self, mech):
        audit = mech.activity_score_grid_audit([0.25, 0.25, 1.0, 2.0])
        assert audit["scores_absent_from_the_published_grid"] == [0.25, 1.0]
        assert audit["n_patients_interpolated"] == 3
        assert audit["n_patients_on_grid"] == 1
        assert audit["fraction_interpolated"] == pytest.approx(0.75)

    def test_the_audit_says_an_interpolated_value_is_not_published(self, mech):
        audit = mech.activity_score_grid_audit([1.0])
        assert "not a value Rüdesheim 2020 reports" in (
            audit["this_is_a_decision_not_a_published_value"]
        )

    def test_the_audit_carries_the_published_table(self, mech):
        audit = mech.activity_score_grid_audit([2.0])
        assert audit["published_grid"] == list(mech.RUDESHEIM_AS_GRID)
        assert audit["published_relative_activity"] == list(
            mech.RUDESHEIM_RELATIVE_ACTIVITY
        )


# --- the omic module -----------------------------------------------------
def _ratio_column(values):
    return np.asarray(values, dtype=float).reshape(-1, 1)


class TestOmicModule:
    def test_reference_uses_only_inhibitor_free_normal_metabolisers(self, mech):
        scores = np.array([2.0, 2.0, 1.25, 2.0, 0.0])
        inhibitor = np.array([0.0, 0.0, 0.0, 1.0, 0.0])
        ratios = _ratio_column([1.0, 2.0, 3.0, 99.0, 99.0])
        got = mech.omic_reference(scores, inhibitor, ratios)
        assert got["n_reference_patients"] == 3
        assert got["reference"][0] == pytest.approx(2.0)
        assert got["inhibitor_free_requirement_relaxed"] is False

    def test_reference_relaxes_the_inhibitor_rule_when_it_must(self, mech):
        scores = np.array([2.0, 2.0, 2.0, 0.0])
        inhibitor = np.array([1.0, 1.0, 0.0, 0.0])
        ratios = _ratio_column([1.0, 2.0, 3.0, 50.0])
        got = mech.omic_reference(scores, inhibitor, ratios)
        assert got["inhibitor_free_requirement_relaxed"] is True
        assert got["n_reference_patients"] == 3
        assert got["reference"][0] == pytest.approx(2.0)

    def test_reference_refuses_fewer_than_three_normal_metabolisers(self, mech):
        with pytest.raises(mech.MechanisticError):
            mech.omic_reference(np.array([0.0, 0.25, 1.0]),
                                np.array([0.0, 0.0, 0.0]),
                                _ratio_column([1.0, 2.0, 3.0]))

    def test_reference_never_takes_the_concentration_as_an_argument(self, mech):
        import inspect

        parameters = inspect.signature(mech.omic_reference).parameters
        assert set(parameters) == {"train_activity_scores", "train_inhibitor",
                                   "train_ratios"}

    def test_a_patient_at_the_reference_has_activity_one(self, mech):
        got = mech.relative_activity_from_omic(_ratio_column([1.4]), [1.4])
        assert got[0] == pytest.approx(1.0)

    def test_the_ratio_is_proportional_to_activity_at_gamma_one(self, mech):
        got = mech.relative_activity_from_omic(
            _ratio_column([1.4 + np.log(2.0)]), [1.4]
        )
        assert got[0] == pytest.approx(2.0)

    def test_two_ratios_combine_as_a_geometric_mean(self, mech):
        ratios = np.array([[np.log(4.0), np.log(1.0)]])
        got = mech.relative_activity_from_omic(ratios, [0.0, 0.0])
        assert got[0] == pytest.approx(2.0)

    def test_no_weight_is_fitted_when_combining_ratios(self, mech):
        text = _flat(MODULE)
        assert "No weight is fitted on the data" in text

    def test_activity_is_capped_at_the_published_ultrarapid_value(self, mech):
        got = mech.relative_activity_from_omic(_ratio_column([50.0]), [0.0])
        assert got[0] == pytest.approx(2.13)

    def test_a_smaller_gamma_attenuates_the_module(self, mech):
        high = mech.relative_activity_from_omic(_ratio_column([1.0]), [0.0],
                                                gamma=1.0)[0]
        low = mech.relative_activity_from_omic(_ratio_column([1.0]), [0.0],
                                               gamma=0.5)[0]
        assert 1.0 < low < high

    def test_a_non_finite_ratio_is_refused(self, mech):
        with pytest.raises(mech.MechanisticError):
            mech.relative_activity_from_omic(_ratio_column([np.nan]), [0.0])


class TestMechanisticActivity:
    @staticmethod
    def arrays():
        return {
            "activity_score": np.array([2.0, 0.0, 1.25]),
            "inhibitor": np.array([0.0, 0.0, 0.0]),
            "ratios": _ratio_column([1.0, 1.0, 1.0]),
        }

    def test_the_control_arm_carries_no_cyp2d6_information(self, mech):
        got = mech.mechanistic_activity(self.arrays(), "none", [1.0])
        assert np.allclose(got, 1.0)

    def test_the_pgx_arm_is_the_genotype_path(self, mech):
        arrays = self.arrays()
        got = mech.mechanistic_activity(arrays, "pgx", [1.0])
        expected = mech.relative_activity_from_activity_score(
            arrays["activity_score"]
        )
        assert np.allclose(got, expected)

    def test_the_omic_arm_is_the_omic_path(self, mech):
        arrays = self.arrays()
        got = mech.mechanistic_activity(arrays, "omic", [1.0])
        assert np.allclose(got, 1.0)

    def test_weight_zero_recovers_the_genotype_arm(self, mech):
        arrays = self.arrays()
        blended = mech.mechanistic_activity(arrays, "pgx_omic", [0.0],
                                            blend_weight=0.0)
        pgx = mech.mechanistic_activity(arrays, "pgx", [0.0])
        assert np.allclose(blended, pgx, atol=1e-5)

    def test_weight_one_recovers_the_omic_arm(self, mech):
        arrays = self.arrays()
        blended = mech.mechanistic_activity(arrays, "pgx_omic", [0.0],
                                            blend_weight=1.0)
        omic = mech.mechanistic_activity(arrays, "omic", [0.0])
        assert np.allclose(blended, omic, atol=1e-5)

    def test_an_unknown_arm_is_refused(self, mech):
        with pytest.raises(mech.MechanisticError):
            mech.mechanistic_activity(self.arrays(), "whatever", [1.0])


# --- the structural backbone ---------------------------------------------
class TestSteadyStateTrough:
    def test_one_compartment_matches_an_independent_closed_form(self, mech):
        backbone = mech.BACKBONES["Chai2026_IR"]
        dose, cl, v = 47.5, 110.0, 562.0
        got = mech.steady_state_trough(dose, cl, v, backbone)[0]
        ka, tau = backbone["ka"], 12.0
        u = (12.0 - backbone["tlag"]) % tau
        k = cl / v
        expected = (dose * 1000.0 * ka / v) * (
            (np.exp(-k * u) / (1 - np.exp(-k * tau))
             - np.exp(-ka * u) / (1 - np.exp(-ka * tau))) / (ka - k)
        )
        assert got == pytest.approx(expected, rel=1e-10)

    def test_two_compartment_roots_satisfy_their_invariants(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        got = mech.steady_state_trough(47.5, 110.0, backbone["v_ref"], backbone)
        assert np.isfinite(got).all() and got[0] > 0

    @pytest.mark.parametrize("key", ["FDA2006_ER", "Chai2026_IR", "Eugene2016_IR"])
    def test_the_trough_is_proportional_to_the_dose(self, mech, key):
        backbone = mech.BACKBONES[key]
        single = mech.steady_state_trough(47.5, 110.0, backbone["v_ref"],
                                          backbone)[0]
        double = mech.steady_state_trough(95.0, 110.0, backbone["v_ref"],
                                          backbone)[0]
        assert double == pytest.approx(2.0 * single, rel=1e-10)

    @pytest.mark.parametrize("key", ["FDA2006_ER", "Chai2026_IR", "Eugene2016_IR"])
    def test_more_clearance_gives_a_lower_trough(self, mech, key):
        backbone = mech.BACKBONES[key]
        low = mech.steady_state_trough(47.5, 50.0, backbone["v_ref"], backbone)[0]
        high = mech.steady_state_trough(47.5, 200.0, backbone["v_ref"],
                                        backbone)[0]
        assert high < low

    def test_it_is_vectorised_over_patients(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        got = mech.steady_state_trough([47.5, 95.0], [110.0, 110.0],
                                       [180.0, 180.0], backbone)
        assert got.shape == (2,)
        assert got[1] == pytest.approx(2.0 * got[0], rel=1e-10)

    def test_an_absurd_clearance_is_refused_rather_than_returned(self, mech):
        backbone = mech.BACKBONES["Chai2026_IR"]
        with pytest.raises(mech.MechanisticError):
            mech.steady_state_trough(47.5, np.inf, 562.0, backbone)


class TestStructuralPrediction:
    @staticmethod
    def arrays(n=6):
        rng = np.random.default_rng(7)
        return {
            "observed": np.full(n, 50.0),
            "activity_score": np.tile([2.0, 1.25, 0.0], n // 3 + 1)[:n],
            "inhibitor": np.zeros(n),
            "weight": 70.0 + rng.normal(0, 5, n),
            "male": np.tile([0.0, 1.0], n // 2 + 1)[:n],
            "daily_dose_mg": np.full(n, 95.0),
            "ratios": np.full((n, 1), 1.4),
        }

    def test_allometry_lowers_the_trough_of_a_heavier_patient(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        light = dict(self.arrays(3), weight=np.array([50.0, 50.0, 50.0]))
        heavy = dict(self.arrays(3), weight=np.array([100.0, 100.0, 100.0]))
        a = mech.structural_prediction(light, 110.0, "pgx", backbone, [1.4])
        b = mech.structural_prediction(heavy, 110.0, "pgx", backbone, [1.4])
        assert np.all(b < a)

    def test_the_sex_effect_lowers_a_male_trough(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        base = self.arrays(2)
        female = dict(base, male=np.zeros(2))
        male = dict(base, male=np.ones(2))
        a = mech.structural_prediction(female, 110.0, "pgx", backbone, [1.4])
        b = mech.structural_prediction(male, 110.0, "pgx", backbone, [1.4])
        assert np.all(b < a)

    def test_a_poor_metaboliser_has_the_higher_trough(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        arrays = dict(self.arrays(2), activity_score=np.array([2.0, 0.0]),
                      weight=np.array([70.0, 70.0]), male=np.zeros(2))
        got = mech.structural_prediction(arrays, 110.0, "pgx", backbone, [1.4])
        assert got[1] > got[0]
        assert got[1] / got[0] == pytest.approx(5.0, rel=0.10)

    def test_a_trough_ratio_is_not_the_published_exposure_ratio(self, mech):
        """The published five-fold is an AUC statement, not a trough statement.

        Exposure is proportional to 1/CL, so the clearance multiplier gives
        exactly five. A trough is not proportional to 1/CL, so the same
        multiplier produces a slightly larger trough ratio. The difference is
        a property of the endpoint and must not be reported as a five-fold
        trough.
        """
        backbone = mech.BACKBONES["FDA2006_ER"]
        arrays = dict(TestStructuralPrediction.arrays(2),
                      activity_score=np.array([2.0, 0.0]),
                      weight=np.array([70.0, 70.0]), male=np.zeros(2))
        got = mech.structural_prediction(arrays, 110.0, "pgx", backbone, [1.4])
        clearance_ratio = (
            mech.clearance_multiplier([1.0])[0] / mech.clearance_multiplier([0.0])[0]
        )
        assert clearance_ratio == pytest.approx(5.0)
        assert got[1] / got[0] != pytest.approx(clearance_ratio, rel=1e-3)

    def test_the_control_arm_ignores_the_genotype(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        arrays = dict(self.arrays(2), activity_score=np.array([2.0, 0.0]),
                      weight=np.array([70.0, 70.0]), male=np.zeros(2))
        got = mech.structural_prediction(arrays, 110.0, "none", backbone, [1.4])
        assert got[0] == pytest.approx(got[1])


class TestBaselineLocking:
    def test_it_recovers_the_generating_baseline_without_noise(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        arrays = TestStructuralPrediction.arrays(9)
        truth = 137.0
        arrays["observed"] = mech.structural_prediction(
            arrays, truth, "pgx", backbone, [1.4]
        )
        got = mech.lock_baseline_clearance(arrays, "pgx", backbone, [1.4])
        assert got == pytest.approx(truth, rel=1e-4)

    def test_it_never_sees_the_test_patients(self, mech):
        text = _flat(MODULE)
        assert "training patients of this outer fold" in text

    def test_only_one_parameter_is_fitted(self, mech):
        import inspect

        source = inspect.getsource(mech.lock_baseline_clearance)
        assert "minimize_scalar" in source


# --- the PopPK covariate arms --------------------------------------------
class TestPopPKArms:
    def test_the_published_arm_is_the_mechanistic_arm(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        arrays = TestStructuralPrediction.arrays(6)
        design = mech.poppk_design(arrays, [1.4])
        mechanistic = mech.structural_prediction(arrays, 123.0, "pgx",
                                                 backbone, [1.4])
        poppk = mech.poppk_predict(arrays, design,
                                   {"log_cl0": float(np.log(123.0))},
                                   backbone, published=True)
        assert np.allclose(np.log(mechanistic), np.log(poppk), atol=1e-12)

    def test_the_design_carries_the_declared_columns(self, mech):
        arrays = TestStructuralPrediction.arrays(4)
        design = mech.poppk_design(arrays, [1.4])
        assert set(design) == {"log_weight", "male", "relative_activity",
                               "log_relative_activity", "omic_centred"}

    def test_the_omic_column_is_centred_on_the_training_reference(self, mech):
        arrays = dict(TestStructuralPrediction.arrays(3),
                      ratios=np.array([[1.4], [2.4], [0.4]]))
        design = mech.poppk_design(arrays, [1.4])
        assert np.allclose(design["omic_centred"], [0.0, 1.0, -1.0])

    def test_the_log_activity_offset_keeps_a_poor_metaboliser_finite(self, mech):
        arrays = dict(TestStructuralPrediction.arrays(3),
                      activity_score=np.zeros(3))
        design = mech.poppk_design(arrays, [1.4])
        assert np.all(np.isfinite(design["log_relative_activity"]))

    def test_the_fit_recovers_its_generating_coefficients(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        arrays = TestStructuralPrediction.arrays(24)
        design = mech.poppk_design(arrays, [1.4])
        truth = {"log_cl0": float(np.log(130.0)), "sex": 0.4,
                 "activity_score": 0.6}
        arrays["observed"] = mech.poppk_predict(arrays, design, truth,
                                                backbone, published=False)
        fit = mech.poppk_fit(arrays, design, backbone, with_omic=False)
        assert fit["n_free_parameters"] == 3
        for name, value in truth.items():
            assert fit["theta"][name] == pytest.approx(value, abs=5e-3)

    def test_the_omic_variant_fits_one_more_parameter(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        arrays = TestStructuralPrediction.arrays(24)
        design = mech.poppk_design(arrays, [1.4])
        fit = mech.poppk_fit(arrays, design, backbone, with_omic=True)
        assert fit["n_free_parameters"] == 4
        assert "omic" in fit["theta"]

    def test_the_fit_reports_a_convergence_cross_check(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        arrays = TestStructuralPrediction.arrays(24)
        design = mech.poppk_design(arrays, [1.4])
        fit = mech.poppk_fit(arrays, design, backbone, with_omic=False)
        assert set(fit) >= {"cross_check_max_parameter_difference",
                            "cross_check_objective_gain", "converged",
                            "solver_success_flag"}
        assert fit["cross_check_max_parameter_difference"] >= 0.0
        assert fit["converged"] is True

    def test_convergence_is_not_read_off_the_solver_flag(self, mech):
        import inspect

        source = inspect.getsource(mech.poppk_fit)
        assert '"converged": bool(discrepancy' in source
        assert "not used as evidence" in _flat(MODULE) or (
            "is a statement about its own stopping rule" in _flat(MODULE)
        )

    def test_the_refit_never_worsens_the_objective(self, mech):
        backbone = mech.BACKBONES["FDA2006_ER"]
        arrays = TestStructuralPrediction.arrays(24)
        design = mech.poppk_design(arrays, [1.4])
        fit = mech.poppk_fit(arrays, design, backbone, with_omic=True)
        assert fit["cross_check_objective_gain"] >= -1e-12

    def test_no_random_effect_and_so_no_map_bayesian_step(self, mech):
        text = _flat(MODULE)
        assert "no MAP-Bayesian step" in text
        assert "no claim of an NLME fit" in text


# --- governance text in the runner ---------------------------------------
class TestRunnerGovernance:
    def test_it_refuses_to_call_the_target_an_auc(self):
        text = _flat(RUNNER)
        assert "NOT an AUC" in text
        assert "NOT a clearance" in text

    def test_it_caps_the_status_at_exploratory(self):
        assert "exploratory" in _flat(RUNNER)
        assert "may be reported as confirmatory" in _flat(RUNNER)

    def test_it_excludes_map_be_by_construction(self):
        assert "circular" in _flat(RUNNER)

    def test_it_declares_the_dosing_assumption(self):
        assert "tau = 12 h" in _flat(RUNNER)

    def test_it_states_the_multiplicity_rule(self):
        text = _flat(RUNNER)
        assert "Bonferroni" in text
        assert "additional independent subjects" in text

    def test_the_adjusted_interval_is_computed_not_promised(self, runner):
        rng = np.random.default_rng(3)
        new = np.abs(rng.normal(10.0, 2.0, 47))
        baseline = new + 1.0
        nominal = runner.H1.paired_contrast(new, baseline, n_boot=400)
        adjusted = runner.bonferroni_interval(new, baseline, 400, 7)
        assert adjusted["family_size"] == 7
        low, high = adjusted["interval"]
        nominal_low, nominal_high = nominal["descriptive_95_interval"]
        assert low <= nominal_low and high >= nominal_high

    def test_a_family_of_one_reproduces_the_nominal_level(self, runner):
        rng = np.random.default_rng(11)
        new = np.abs(rng.normal(10.0, 2.0, 47))
        baseline = new + 0.5
        adjusted = runner.bonferroni_interval(new, baseline, 400, 1)
        nominal = runner.H1.paired_contrast(new, baseline, n_boot=400)
        assert adjusted["interval"] == pytest.approx(
            nominal["descriptive_95_interval"], abs=1e-12
        )
        assert adjusted["level"].startswith("95")

    def test_the_sensitivity_sweeps_are_declared_before_the_run(self, runner):
        assert set(runner.SENSITIVITY_SWEEPS) == {
            "fraction_metabolised_by_cyp2d6", "omic_exponent_gamma",
            "genotype_omic_weight", "activity_score_decrement_per_inhibitor",
            "salt_dose_equivalence",
        }
        for _, (keyword, values) in runner.SENSITIVITY_SWEEPS.items():
            assert keyword in {"fm", "gamma", "blend_weight", "delta_as",
                               "dose_equivalence"}
            assert len(values) >= 2

    def test_the_published_value_is_inside_every_sweep(self, runner, mech):
        published = {"fm": mech.FM_CYP2D6, "gamma": mech.OMIC_GAMMA,
                     "blend_weight": mech.BLEND_WEIGHT,
                     "delta_as": mech.DELTA_AS_PER_INHIBITOR,
                     "dose_equivalence": 1.0}
        for label, (keyword, values) in runner.SENSITIVITY_SWEEPS.items():
            assert published[keyword] in values, label

    def test_it_says_the_sweep_does_not_choose_a_value(self):
        assert "fixed before the sweep" in _flat(RUNNER)
        assert "not to choose a value" in _flat(RUNNER)


# --- the committed artefact ----------------------------------------------
class TestReportArtefact:
    def test_the_source_checksum_is_the_registered_one(self, report):
        assert report["source_sha256"] == (
            "7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"
        )

    def test_the_target_is_a_trough_not_an_exposure(self, report):
        assert report["target"] == "c(Metoprolol)[ng/ml]"
        assert "NOT an AUC" in report["target_is"]

    def test_the_prohibited_predictor_is_named_as_absent(self, report):
        assert report["prohibited_predictor_absent"] == PROHIBITED

    def test_the_status_ceiling_is_exploratory(self, report):
        assert report["status_ceiling"].startswith("exploratory")

    def test_the_folds_are_the_primary_contrast_folds(self, report):
        validation = report["validation"]
        assert validation["folds_identical_to_the_primary_contrast"] is True
        assert validation["seed"] == 20261107
        assert validation["outer_folds"] == 5
        assert validation["outer_repeats"] == 2
        assert len(validation["fold_fingerprint_sha256"]) == 64

    def test_p1_and_m4_agree_to_floating_point(self, report):
        assert report["p1_and_m4_are_the_same_model"][
            "max_absolute_log_difference"] < 1e-9

    def test_every_literature_constant_is_declared_unfitted(self, report):
        assert report["literature_constants"]
        for entry in report["literature_constants"]:
            assert entry["estimated_on_ariana"] is False
            assert entry["source"]

    def test_the_activity_score_alignment_is_reported(self, report):
        audit = report["activity_score_alignment"]
        assert audit["n_patients_interpolated"] + audit["n_patients_on_grid"] == 47
        assert audit["scores_absent_from_the_published_grid"]

    def test_the_key_secondary_carries_its_whole_loss_hierarchy(self, report):
        contrast = report["key_secondary_contrast"]
        assert "primary_loss" in contrast
        assert set(contrast["robustness_co_criteria"]) == {
            "log_mae", "trimmed_mae_10pct"}
        assert "reversal_rule" in contrast
        assert "influence" in contrast
        assert contrast["primary_loss"]["n_independent_patients"] == 47

    def test_the_bootstrap_unit_is_the_patient(self, report):
        assert "patients" in report["validation"]["bootstrap"]["unit"]
        assert report["validation"]["bootstrap"]["n"] == 2000

    def test_the_protocol_named_secondaries_are_present(self, report):
        for label in ("P3_vs_P1", "HY_vs_P1"):
            assert label in report["other_secondary_contrasts"]

    def test_the_free_parameter_counts_are_in_the_report(self, report):
        counts = report["free_parameters_fitted_on_this_cohort"]
        assert counts["P2_poppk_estimated_coefficients"] == 3
        assert counts["P3_poppk_estimated_plus_omic"] == 4

    def test_the_limits_are_listed_not_implied(self, report):
        limits = report["what_this_cannot_establish"]
        assert len(limits) >= 4
        assert any("replication" in item for item in limits)
        assert any("SSDA" in item for item in limits)

    def test_every_estimated_fit_passed_the_convergence_cross_check(self, report):
        trace = report["arms"]["poppk_covariate"]
        assert trace["fits_failing_the_convergence_cross_check"] == 0
        check = trace["convergence_cross_check"]
        assert check["max_parameter_disagreement_over_all_fits"] < check["tolerance"]

    def test_every_other_secondary_carries_an_adjusted_interval(self, report):
        family = report["other_secondary_contrasts"]
        for label, contrast in family.items():
            adjusted = contrast["multiplicity_adjusted"]
            assert adjusted["family_size"] == len(family), label
            nominal = contrast["primary_loss"]["descriptive_95_interval"]
            assert adjusted["interval"][0] <= nominal[0] + 1e-9, label
            assert adjusted["interval"][1] >= nominal[1] - 1e-9, label

    def test_the_multiplicity_block_names_what_survives(self, report):
        block = report["multiplicity"]
        assert block["other_secondary_family_size"] == len(
            report["other_secondary_contrasts"]
        )
        survive = set(block["which_other_secondaries_survive_adjustment"])
        lose = set(block["which_lose_their_interval_on_adjustment"])
        assert not (survive & lose)
        for label in survive:
            assert report["other_secondary_contrasts"][label][
                "multiplicity_adjusted"]["interval_excludes_zero"]

    def test_the_unfair_comparison_is_labelled_as_such(self, report):
        note = report["other_secondary_contrasts"]["HY_vs_M5"][
            "not_a_fair_arm_comparison"]
        assert "descriptive only" in note

    def test_the_sensitivity_sweeps_were_run(self, report):
        assert set(report["structural_sensitivity"]) == {
            "FDA2006_ER", "Chai2026_IR", "Eugene2016_IR"}
        assert set(report["omic_panel_sensitivity"]) == {
            "oh_only", "ssda_only", "oh_and_ssda"}
        assert report["transferred_constant_sensitivity"]
