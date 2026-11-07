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
