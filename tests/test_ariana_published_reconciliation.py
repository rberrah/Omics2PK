"""Offline tests for the ARIANA published-correlation reconciliation.

These exercise the pure transforms and the checksum guard. They never touch
the network and never require the workbook, so they run in CI. The published
reference values are pinned here so that a change to them is a deliberate,
reviewable edit.
"""
from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

try:  # See the note in test_public_data_evidence_registry: the `offline` job
    # has no dependencies, so this must skip rather than error on import.
    import numpy as np
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(
        f"numpy and pytest are required for this module: {exc}"
    ) from exc

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "reconcile_ariana_published_correlations.py"


def _load():
    spec = importlib.util.spec_from_file_location("reconcile_ariana", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return _load()


def test_published_reference_values_are_pinned(mod):
    assert mod.PUBLISHED["ssda_vs_reference_mr"] == {"r": 0.946, "r2": 0.896}
    assert mod.PUBLISHED["oh_solanidine_vs_reference_mr"] == {"r": 0.907, "r2": 0.823}
    assert mod.PUBLISHED["activity_score_vs_reference_mr"] == {"r": 0.500, "r2": 0.250}
    assert mod.EXPECTED_SHA256 == (
        "7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"
    )


def test_checksum_guard_rejects_a_wrong_file(mod, tmp_path):
    decoy = tmp_path / "not_the_workbook.xlsx"
    decoy.write_bytes(b"not the published workbook")
    with pytest.raises(mod.ReconciliationError, match="checksum mismatch"):
        mod.verify_source(decoy)


def test_rudesheim_mapping_hits_the_published_grid(mod):
    """The activity mapping must pass exactly through Rudesheim Table 3."""
    got = mod.relative_activity_from_as(list(mod.AS_GRID))
    np.testing.assert_allclose(got, mod.KCAT_REL, rtol=0, atol=1e-12)


def test_rudesheim_mapping_is_monotonic_and_nonlinear(mod):
    scores = np.linspace(0, 3, 61)
    activity = mod.relative_activity_from_as(scores)
    assert np.all(np.diff(activity) >= -1e-12), "activity must not decrease with activity score"
    # A linear activity score would make activity proportional to the score;
    # the published mapping is deliberately not proportional.
    midpoint = mod.relative_activity_from_as(1.5)
    assert midpoint < 1.5 / 2.0 * mod.relative_activity_from_as(3.0)


def test_omic_activity_is_one_at_the_reference_level(mod):
    assert mod.omic_activity(2.0, 2.0) == pytest.approx(1.0)


def test_omic_activity_cap_only_truncates_upward(mod):
    ln_mr = np.array([-5.0, 0.0, 5.0])
    uncapped = mod.omic_activity(ln_mr, 0.0, cap=None)
    capped = mod.omic_activity(ln_mr, 0.0)
    assert np.all(capped <= uncapped + 1e-12)
    assert capped.max() == pytest.approx(mod.ACTIVITY_CAP)
    # Below the cap the two agree exactly, so capping cannot alter low activity.
    below = uncapped < mod.ACTIVITY_CAP
    np.testing.assert_allclose(capped[below], uncapped[below])


def test_cap_cannot_increase_a_correlation(mod):
    """Capping compresses the upper tail; it may only lose information.

    This is the mechanism behind the one v5 number that differed from the
    publication (SSDA, 15 of 47 patients at the cap).
    """
    rng = np.random.default_rng(20261107)
    truth = rng.normal(size=400)
    ln_mr = truth * 1.5
    uncapped = np.log(mod.omic_activity(ln_mr, 0.0, cap=None))
    capped = np.log(mod.omic_activity(ln_mr, 0.0))
    r_uncapped = abs(np.corrcoef(uncapped, truth)[0, 1])
    r_capped = abs(np.corrcoef(capped, truth)[0, 1])
    assert r_capped <= r_uncapped + 1e-12


def test_phenoconversion_correction_only_lowers_activity(mod):
    scores = np.array([0.5, 1.0, 2.0, 3.0])
    with_inhibitor = np.ones_like(scores)
    corrected = np.clip(scores - mod.DELTA_AS_INHIBITOR * with_inhibitor, 0.0, None)
    assert np.all(corrected <= scores)
    assert np.all(mod.relative_activity_from_as(corrected) <= mod.relative_activity_from_as(scores))


def test_delta_as_matches_the_transferred_coefficient(mod):
    assert mod.DELTA_AS_INHIBITOR == 0.54
    assert mod.ACTIVITY_CAP == 2.13
    assert mod.NM_WINDOW == (1.25, 2.25)


def test_reference_mr_column_is_never_a_predictor_name(mod):
    """The reference ratio shares a denominator with the target concentration."""
    assert mod.COL_REFERENCE_MR == "lnMR OH-metoprolol/metoprolol"
    # It is correlated against, never transformed into an activity estimate.
    assert mod.COL_REFERENCE_MR not in (mod.COL_OH, mod.COL_SSDA, mod.COL_AS)
