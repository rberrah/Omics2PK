"""Contract tests for the H5 solanidine marker robustness analysis.

Offline: these exercise the statistics on constructed data with known truth,
never the network and never the source workbooks. They pin the prespecified
choices so a later change is a deliberate, reviewable edit, and they encode the
conclusions the analysis is NOT allowed to draw.
"""
from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

try:  # See test_public_data_evidence_registry: the `offline` CI job has no
    # dependencies, so this must skip rather than error during discovery.
    import numpy as np
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(
        f"numpy and pytest are required for this module: {exc}"
    ) from exc

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "run_h5_solanidine_marker_robustness.py"
)


def _load():
    spec = importlib.util.spec_from_file_location("run_h5", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return _load()


# --------------------------------------------------------------------------
# prespecification and provenance
# --------------------------------------------------------------------------
def test_both_source_checksums_are_pinned(mod):
    assert mod.MULLER_SHA256 == (
        "68d53071db79f961a63e433137a8a3e7dbd40daf8c8ad0f9f7f1b268bb3a8e07"
    )
    assert mod.ARIANA_SHA256 == (
        "7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"
    )


def test_checksum_guard_rejects_a_wrong_file(mod, tmp_path):
    decoy = tmp_path / "decoy.xlsx"
    decoy.write_bytes(b"not the published workbook")
    with pytest.raises(mod.AnalysisError, match="checksum mismatch"):
        mod.verify(decoy, mod.MULLER_SHA256)


def test_reference_time_is_a_protocol_choice_not_the_best_result(mod):
    """3 h is a mid-range post-meal time. The informative time proved to be 29 h."""
    assert mod.REFERENCE_TIME_H == 3.0
    assert mod.REFERENCE_TIME_H in mod.TIMEPOINTS_H
    assert mod.TIMEPOINTS_H[-1] == 29.0


def test_prespecification_claim_is_not_overstated(mod):
    """The module must not claim the statistics were frozen. They were not.

    Protocol v1.0.0 froze the H5 objectives and criterion families; Spearman,
    the ICC, Kendall's W and the reference times were chosen at analysis time,
    and the 29 h reference was chosen after seeing a result.
    """
    doc = mod.__doc__
    assert "statistical\noperationalisation added later" in doc or (
        "operationalisation added later" in doc
    )
    assert "was added after the first run" in doc
    assert "post-hoc and exploratory" in doc
    # The discredited phrasing must not come back.
    assert "Nothing here is tuned on the results" not in doc


def test_multiplicity_family_is_two_markers_by_nine_times(mod):
    assert mod.PRIMARY_FAMILY_SIZE == 18
    assert mod.CONSERVATIVE_FAMILY_SIZE == 54
    assert mod.PRIMARY_FAMILY_SIZE == 2 * len(mod.TIMEPOINTS_H)


def test_minimum_detectable_effect_is_large_at_this_sample_size(mod):
    """A null at n=19 under correction cannot exclude a moderate association."""
    detectable = mod.minimum_detectable_spearman(19, 0.05 / mod.PRIMARY_FAMILY_SIZE)
    assert 0.70 < detectable < 0.82
    # More subjects must lower the detectable threshold.
    assert mod.minimum_detectable_spearman(
        60, 0.05 / mod.PRIMARY_FAMILY_SIZE
    ) < detectable
    assert mod.minimum_detectable_spearman(3, 0.05) is None


def test_removing_the_time_effect_centres_each_timepoint(mod):
    matrix = np.array(
        [[1.0, 11.0, 21.0], [3.0, 13.0, 23.0], [5.0, 15.0, 25.0]], dtype=float
    )
    out = mod.remove_time_effect(matrix)
    np.testing.assert_allclose(np.nanmean(out, axis=0), np.zeros(3), atol=1e-12)
    # Subject ordering within a timepoint is preserved.
    assert out[0, 0] < out[1, 0] < out[2, 0]


def test_removing_the_time_effect_ignores_masked_cells(mod):
    matrix = np.array([[1.0, 2.0], [3.0, 4.0], [100.0, 5.0]], dtype=float)
    mask = np.zeros_like(matrix, dtype=bool)
    mask[2, 0] = True
    out = mod.remove_time_effect(matrix, mask)
    assert np.isnan(out[2, 0])
    # The outlier was excluded from the column mean, so the kept cells centre
    # on the mean of 1.0 and 3.0.
    assert out[0, 0] == pytest.approx(-1.0)
    assert out[1, 0] == pytest.approx(1.0)


def test_icc_attribution_separates_subject_set_from_censoring(mod):
    """The reversed finding: floor-pinned subjects, not cell-level imputation.

    Three subjects pinned far below the rest inflate between-subject variance.
    On a COMMON subject set the censoring-only effect must be small, and must
    not reproduce the large headline drop.
    """
    rng = np.random.default_rng(11)
    matrix = rng.normal(scale=0.3, size=(mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H)))
    matrix += np.linspace(0, 1.5, mod.N_PARTICIPANTS).reshape(-1, 1)
    mask = np.zeros_like(matrix, dtype=bool)
    matrix[:3, :] = -8.0
    mask[:3, :] = True

    headline_all = mod.within_subject_stability(matrix)["icc_between_share"]
    headline_censored = mod.within_subject_stability(matrix, mask)["icc_between_share"]
    attribution = mod.icc_attribution(matrix, mask)

    assert attribution["n_common_subjects"] == mod.N_PARTICIPANTS - 3
    assert attribution["n_subjects_dropped_by_censoring"] == 3
    # The headline difference is large...
    assert headline_all - headline_censored > 0.3
    # ...while the censoring-only difference on a common subject set is small.
    assert abs(attribution["censoring_only_delta_icc"]) < 0.1


def test_power_value_is_not_presented_as_an_upper_bound(mod):
    """The correlation with 80% power is not a ceiling on the true value.

    Revision 1.1 misused it as one. The script must say so explicitly, and the
    discredited phrasing must not return.
    """
    source = SCRIPT.read_text(encoding="utf-8")
    assert "NOT an upper bound" in source
    assert "not established with the available\n            precision" in source or (
        "not established with the available" in source
    )
    assert "any association there is smaller than the minimum detectable" not in source


def test_marker_icc_difference_uses_the_exact_intersection(mod):
    """Both markers must be compared on identical subjects, with uncertainty."""
    rng = np.random.default_rng(5)
    shape = (mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H))
    subject_effect = np.linspace(0.0, 3.0, mod.N_PARTICIPANTS).reshape(-1, 1)

    # Marker A: mostly within-subject noise, so a low ICC.
    a = subject_effect * 0.1 + rng.normal(scale=1.0, size=shape)
    # Marker B: strong subject separation, so a high ICC.
    b = subject_effect + rng.normal(scale=0.1, size=shape)

    mask_a = np.zeros(shape, dtype=bool)
    mask_b = np.zeros(shape, dtype=bool)
    mask_a[0, :] = True   # subject 0 unusable for A only
    mask_b[1, :] = True   # subject 1 unusable for B only

    out = mod.icc_difference_on_exact_intersection(
        {"marker_a": a, "marker_b": b},
        {"marker_a": mask_a, "marker_b": mask_b},
        n_boot=200,
        seed=1,
    )
    # The intersection drops both partially-unusable subjects.
    assert out["n_intersection"] == mod.N_PARTICIPANTS - 2
    assert out["subjects_excluded_from_either_marker"] == 2
    assert out["difference_is"] == "marker_b minus marker_a"
    # B separates subjects better, so the difference must be positive.
    assert out["point_estimates"]["raw"]["difference"] > 0
    interval = out["bootstrap_raw"]["percentile_95_interval"]
    assert interval[0] < out["bootstrap_raw"]["median_difference"] < interval[1]


def test_marker_difference_bootstrap_is_reproducible(mod):
    rng = np.random.default_rng(9)
    shape = (mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H))
    ratios = {
        "one": rng.normal(size=shape),
        "two": np.linspace(0, 2, mod.N_PARTICIPANTS).reshape(-1, 1)
        + rng.normal(scale=0.2, size=shape),
    }
    masks = {name: np.zeros(shape, dtype=bool) for name in ratios}
    first = mod.icc_difference_on_exact_intersection(ratios, masks, n_boot=150, seed=42)
    second = mod.icc_difference_on_exact_intersection(ratios, masks, n_boot=150, seed=42)
    assert (
        first["bootstrap_raw"]["percentile_95_interval"]
        == second["bootstrap_raw"]["percentile_95_interval"]
    )


def test_marker_difference_requires_exactly_two_markers(mod):
    shape = (mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H))
    ratios = {name: np.zeros(shape) for name in ("a", "b", "c")}
    masks = {name: np.zeros(shape, dtype=bool) for name in ratios}
    with pytest.raises(mod.AnalysisError, match="exactly two ratios"):
        mod.icc_difference_on_exact_intersection(ratios, masks, n_boot=10)


def test_marker_difference_abstains_on_a_tiny_intersection(mod):
    shape = (mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H))
    ratios = {"a": np.zeros(shape), "b": np.zeros(shape)}
    masks = {name: np.ones(shape, dtype=bool) for name in ratios}
    masks["a"][:2, :] = False
    masks["b"][:2, :] = False
    out = mod.icc_difference_on_exact_intersection(ratios, masks, n_boot=10)
    assert out["insufficient"] is True


def test_icc_attribution_abstains_when_too_few_subjects_remain(mod):
    matrix = np.zeros((mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H)))
    mask = np.ones_like(matrix, dtype=bool)
    mask[:2, :] = False
    assert mod.icc_attribution(matrix, mask)["insufficient"] is True


def test_ratios_are_metabolite_over_parent_as_in_ariana(mod):
    assert mod.RATIOS["oh4_solanidine_over_solanidine"] == (
        "oh4_solanidine",
        "solanidine",
    )
    assert mod.RATIOS["ssda_over_solanidine"] == ("ssda", "solanidine")


# --------------------------------------------------------------------------
# the statistics, against constructed known truth
# --------------------------------------------------------------------------
def test_icc_is_near_one_when_subjects_differ_and_time_does_not(mod):
    offsets = np.linspace(0.0, 5.0, mod.N_PARTICIPANTS).reshape(-1, 1)
    matrix = offsets + np.zeros((1, len(mod.TIMEPOINTS_H)))
    matrix = matrix + 1e-3 * np.arange(len(mod.TIMEPOINTS_H))
    out = mod.within_subject_stability(matrix)
    assert out["icc_between_share"] > 0.99
    assert out["median_within_subject_sd_ln"] < 0.01


def test_icc_is_near_zero_when_all_variance_is_within_subject(mod):
    rng = np.random.default_rng(20261107)
    matrix = rng.normal(size=(mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H)))
    out = mod.within_subject_stability(matrix)
    assert out["icc_between_share"] < 0.25


def test_censoring_mask_changes_the_icc(mod):
    """The headline finding: dropping imputed floors must be able to move the ICC.

    Floor values at one extreme inflate between-subject spread, so removing
    them has to be capable of lowering the ICC. This is the mechanism behind
    OH-solanidine falling from 0.942 to 0.400.
    """
    rng = np.random.default_rng(1)
    matrix = rng.normal(scale=0.2, size=(mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H)))
    mask = np.zeros_like(matrix, dtype=bool)
    # Five subjects pinned at a distant floor, flagged as imputed.
    matrix[:5, :] = -6.0
    mask[:5, :] = True
    with_floor = mod.within_subject_stability(matrix)
    without_floor = mod.within_subject_stability(matrix, mask)
    assert with_floor["icc_between_share"] > without_floor["icc_between_share"]
    assert with_floor["between_subject_sd_ln"] > without_floor["between_subject_sd_ln"]


def test_kendalls_w_is_one_for_identical_rankings(mod):
    column = np.arange(mod.N_PARTICIPANTS, dtype=float).reshape(-1, 1)
    matrix = np.repeat(column, len(mod.TIMEPOINTS_H), axis=1)
    out = mod.kendalls_w(matrix)
    assert out["kendalls_w"] == pytest.approx(1.0, abs=1e-9)
    assert out["n_subjects"] == mod.N_PARTICIPANTS


def test_kendalls_w_is_low_for_unrelated_rankings(mod):
    rng = np.random.default_rng(7)
    matrix = rng.normal(size=(mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H)))
    assert mod.kendalls_w(matrix)["kendalls_w"] < 0.4


def test_kendalls_w_uses_a_balanced_block(mod):
    """Subjects with any censored value are dropped, not silently partial."""
    matrix = np.tile(
        np.arange(mod.N_PARTICIPANTS, dtype=float).reshape(-1, 1),
        (1, len(mod.TIMEPOINTS_H)),
    )
    mask = np.zeros_like(matrix, dtype=bool)
    mask[:6, 0] = True
    out = mod.kendalls_w(matrix, mask)
    assert out["n_subjects"] == mod.N_PARTICIPANTS - 6


def test_agreement_between_times_excludes_the_self_comparison(mod):
    """The earlier statistic correlated a time against the mean containing it."""
    matrix = np.tile(
        np.arange(mod.N_PARTICIPANTS, dtype=float).reshape(-1, 1),
        (1, len(mod.TIMEPOINTS_H)),
    )
    out = mod.agreement_between_times(matrix, reference_h=3.0)
    assert "3h" not in out
    assert len(out) == len(mod.TIMEPOINTS_H) - 1
    for record in out.values():
        assert record["rho"] == pytest.approx(1.0)


def test_food_effect_reports_an_interval_and_flags_the_small_group(mod):
    rng = np.random.default_rng(3)
    matrix = rng.normal(scale=0.1, size=(mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H)))
    flag = np.array([1] * 15 + [0] * 5, dtype=float)
    matrix[15:, :] += 2.0  # the non-fasted group shifted upward
    out = mod.food_effect(matrix, flag)
    assert out["n_fasted"] == 15
    assert out["n_not_fasted"] == 5
    # Fasted minus not-fasted must be negative when the non-fasted are higher.
    assert out["mean_difference_ln"] < 0
    assert out["interval_excludes_zero"] is True
    assert "descriptive" in out["warning"] or "n=5" in out["warning"]


def test_food_effect_abstains_when_a_group_is_too_small(mod):
    matrix = np.zeros((mod.N_PARTICIPANTS, len(mod.TIMEPOINTS_H)))
    flag = np.array([1] * 19 + [0], dtype=float)
    assert mod.food_effect(matrix, flag)["insufficient"] is True


def test_spearman_abstains_below_five_observations(mod):
    out = mod._spearman([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0])
    assert out["insufficient"] is True
    assert out["n"] == 4


# --------------------------------------------------------------------------
# the conclusions this analysis must never draw
# --------------------------------------------------------------------------
def test_report_declares_itself_not_pharmacokinetics(mod):
    """No drug was administered, so no PK claim is admissible from this cohort."""
    text = mod.__doc__.lower()
    assert "no drug is administered" in text
    assert "nothing here is pharmacokinetics" in text


def test_prohibited_conclusions_are_enumerated_in_the_module(mod):
    """These must be stated in the script, not only in prose elsewhere."""
    source = SCRIPT.read_text(encoding="utf-8").lower()
    for phrase in (
        "between-day reproducibility",
        "poor metaboliser",
        "method bias",
        "not numerically comparable between cohorts",
    ):
        assert phrase in source, phrase
