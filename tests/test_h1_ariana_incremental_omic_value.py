"""Contract tests for the H1 ARIANA incremental-omic-value analysis.

Offline: constructed data with known truth, pinned provenance constants, and
the governance rules the analysis must not break. The source workbook is never
touched and no network call is made.

The rule these tests exist above all to defend: the
``lnMR OH-metoprolol/metoprolol`` column shares its denominator with the target
concentration and must never appear as a predictor in any arm.
"""
from __future__ import annotations

import importlib.util
import json
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
SCRIPT = ROOT / "scripts" / "run_h1_ariana_incremental_omic_value.py"
REPORT_DATA = ROOT / "reports" / "data"


def _load():
    spec = importlib.util.spec_from_file_location("run_h1", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return _load()


# --------------------------------------------------------------------------
# the prohibition that matters most
# --------------------------------------------------------------------------
def test_target_contaminated_ratio_is_in_no_arm(mod):
    assert mod.COL_PROHIBITED == "lnMR OH-metoprolol/metoprolol"
    for arm, features in mod.ARMS.items():
        assert mod.COL_PROHIBITED not in features, arm
    mod.assert_no_prohibited_predictor()  # must not raise


def test_the_guard_actually_fires(mod, monkeypatch):
    """A future edit that smuggles the ratio into an arm must fail loudly."""
    monkeypatch.setitem(mod.ARMS, "rogue_arm", [mod.COL_PROHIBITED])
    with pytest.raises(mod.AnalysisError, match="prohibited predictor"):
        mod.assert_no_prohibited_predictor()


def test_omic_features_do_not_contain_metoprolol(mod):
    """The omic predictors must be solanidine-derived, sharing no analyte."""
    for feature in mod.OMIC:
        assert "solanidine" in feature
        assert "metoprolol" not in feature.lower()


# --------------------------------------------------------------------------
# provenance and prespecification
# --------------------------------------------------------------------------
def test_source_checksum_is_pinned(mod):
    assert mod.EXPECTED_SHA256 == (
        "7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"
    )


def test_checksum_guard_rejects_a_wrong_file(mod, tmp_path):
    decoy = tmp_path / "decoy.xlsx"
    decoy.write_bytes(b"not the published workbook")
    with pytest.raises(mod.AnalysisError, match="checksum mismatch"):
        mod.verify_source(decoy)


def test_primary_contrast_is_the_frozen_one(mod):
    assert mod.PRIMARY_CONTRAST == ("M3_clinical_pgx_omic", "M1_clinical_pgx")
    assert mod.PRIMARY_MODEL == "ridge"
    assert mod.PRIMARY_LOSS == "mae"
    # Ridge must be the primary model: it cannot borrow literature strength.
    assert "ridge" in mod.model_factories()


def test_arms_are_nested_as_the_protocol_defines(mod):
    arms = mod.ARMS
    assert set(arms["M0_clinical"]) < set(arms["M1_clinical_pgx"])
    assert set(arms["M0_clinical"]) < set(arms["M2_clinical_omic"])
    assert set(arms["M1_clinical_pgx"]) < set(arms["M3_clinical_pgx_omic"])
    assert set(arms["M2_clinical_omic"]) < set(arms["M3_clinical_pgx_omic"])


def test_validation_scheme_is_the_frozen_one(mod):
    assert mod.OUTER_FOLDS == 5
    assert mod.OUTER_REPEATS == 2
    assert mod.INNER_FOLDS == 3
    assert mod.RELATIVE_IMPROVEMENT_THRESHOLD == 0.10
    assert mod.TRIM_FRACTION == 0.10


def test_xgboost_is_present_so_no_substitution_is_claimed(mod):
    """The protocol names XGBoost; it must really be XGBoost."""
    factories = mod.model_factories()
    assert "xgboost" in factories
    estimator = factories["xgboost"][0](n_estimators=10, max_depth=1)
    assert type(estimator).__module__.startswith("xgboost")


# --------------------------------------------------------------------------
# folds and leakage
# --------------------------------------------------------------------------
def test_folds_partition_every_patient_once_per_repeat(mod):
    n = 47
    folds = mod.make_folds(n, mod.OUTER_FOLDS, mod.OUTER_REPEATS, mod.SEED)
    assert len(folds) == mod.OUTER_FOLDS * mod.OUTER_REPEATS
    for repeat in range(mod.OUTER_REPEATS):
        block = folds[repeat * mod.OUTER_FOLDS:(repeat + 1) * mod.OUTER_FOLDS]
        covered = np.concatenate([test for _, test in block])
        assert sorted(covered) == list(range(n))
    for train, test in folds:
        assert not set(train) & set(test), "a patient appears in train and test"


def test_folds_are_identical_across_calls(mod):
    """Every arm and model must see the same folds, so the seed must bind."""
    first = mod.make_folds(47, mod.OUTER_FOLDS, mod.OUTER_REPEATS, mod.SEED)
    second = mod.make_folds(47, mod.OUTER_FOLDS, mod.OUTER_REPEATS, mod.SEED)
    for (train_a, test_a), (train_b, test_b) in zip(first, second):
        np.testing.assert_array_equal(train_a, train_b)
        np.testing.assert_array_equal(test_a, test_b)


def test_out_of_fold_recovers_a_learnable_signal(mod):
    """Sanity: a linear target in one feature must be predicted out-of-fold."""
    from sklearn.linear_model import Ridge

    rng = np.random.default_rng(3)
    n = 47
    X = rng.normal(size=(n, 3))
    y_log = 2.0 * X[:, 0] + 0.05 * rng.normal(size=n)
    folds = mod.make_folds(n, mod.OUTER_FOLDS, mod.OUTER_REPEATS, mod.SEED)
    predictions = mod.out_of_fold_predictions(
        X, y_log, folds, lambda **kw: Ridge(**kw), {"alpha": [0.01, 1.0]}
    )
    assert predictions.shape == (mod.OUTER_REPEATS, n)
    assert not np.isnan(predictions).any()
    correlation = np.corrcoef(predictions.mean(axis=0), y_log)[0, 1]
    assert correlation > 0.9


def test_out_of_fold_gives_no_skill_on_pure_noise(mod):
    """Negative control: an unlearnable target must not be predicted well.

    If preprocessing or hyperparameter choice leaked the test fold, this would
    show skill where none exists.
    """
    from sklearn.linear_model import Ridge

    rng = np.random.default_rng(17)
    n = 47
    X = rng.normal(size=(n, 5))
    y_log = rng.normal(size=n)  # independent of X by construction
    folds = mod.make_folds(n, mod.OUTER_FOLDS, mod.OUTER_REPEATS, mod.SEED)
    predictions = mod.out_of_fold_predictions(
        X, y_log, folds, lambda **kw: Ridge(**kw), {"alpha": [0.1, 1.0, 10.0]}
    )
    correlation = np.corrcoef(predictions.mean(axis=0), y_log)[0, 1]
    assert correlation < 0.5, f"suspicious out-of-fold skill on noise: {correlation}"


# --------------------------------------------------------------------------
# losses and contrasts
# --------------------------------------------------------------------------
def test_trimmed_mean_drops_the_largest_errors(mod):
    values = [1.0] * 9 + [1000.0]
    assert mod.trimmed_mean(values, fraction=0.10) == pytest.approx(1.0)
    assert mod.trimmed_mean(values, fraction=0.0) == pytest.approx(100.9)


def test_paired_contrast_sign_and_reproducibility(mod):
    better = np.full(47, 1.0)
    worse = np.full(47, 3.0)
    out = mod.paired_contrast(better, worse, n_boot=200, seed=1)
    assert out["mean_difference"] == pytest.approx(-2.0)
    assert out["direction"] == "improvement"
    assert out["relative_change_vs_baseline"] == pytest.approx(-2.0 / 3.0)
    assert out["n_independent_patients"] == 47
    again = mod.paired_contrast(better, worse, n_boot=200, seed=1)
    assert out["descriptive_95_interval"] == again["descriptive_95_interval"]


def test_paired_contrast_interval_covers_zero_when_there_is_no_effect(mod):
    rng = np.random.default_rng(2)
    base = rng.normal(loc=10.0, scale=2.0, size=47)
    out = mod.paired_contrast(base, base.copy(), n_boot=300, seed=5)
    assert out["mean_difference"] == pytest.approx(0.0)
    assert out["interval_excludes_zero"] is False


def test_influence_report_detects_concentrated_gain(mod):
    """Five patients carrying the whole gain must be surfaced, not hidden."""
    baseline = np.full(47, 10.0)
    new = baseline.copy()
    new[:5] = 0.0  # only five patients improve
    y_true = np.linspace(10.0, 500.0, 47)
    out = mod.influence_report(new, baseline, y_true)
    assert out["n_improved"] == 5
    assert out["n_worsened"] == 0
    assert out["share_of_net_gain_from_five_most_improved"] == pytest.approx(1.0)
    assert set(out["mean_gain_by_endpoint_quartile"]) <= {"Q1", "Q2", "Q3", "Q4"}


def test_reversal_rule_caps_the_status(mod):
    improved = {"mean_difference": -5.0}
    reversed_up = {"mean_difference": +0.2}
    out = mod.reversal_check(improved, reversed_up, reversed_up)
    assert out["primary_improved"] is True
    assert out["status_capped_at_partially_supported"] is True
    # One co-criterion agreeing is enough to avoid the cap.
    agreeing = {"mean_difference": -0.1}
    out = mod.reversal_check(improved, agreeing, reversed_up)
    assert out["status_capped_at_partially_supported"] is False


# --------------------------------------------------------------------------
# the committed report inputs
# --------------------------------------------------------------------------
def _report_files():
    return sorted(REPORT_DATA.glob("*.json")) if REPORT_DATA.is_dir() else []


def test_committed_report_json_carries_no_per_patient_vectors():
    """Aggregates only. A length-47 vector would be participant-level data."""
    def arrays(obj, path=""):
        found = []
        if isinstance(obj, dict):
            for key, value in obj.items():
                found += arrays(value, f"{path}.{key}")
        elif isinstance(obj, list):
            numeric = [x for x in obj if isinstance(x, (int, float))]
            if len(numeric) >= 10:
                found.append((path, len(numeric)))
            for index, value in enumerate(obj):
                if not isinstance(value, (int, float)):
                    found += arrays(value, f"{path}[{index}]")
        return found

    files = _report_files()
    assert files, "report inputs should be committed for the Rmd to render"
    for path in files:
        payload = json.loads(path.read_text(encoding="utf-8"))
        offenders = arrays(payload)
        assert not offenders, f"{path.name} carries numeric arrays: {offenders}"


def test_committed_h1_report_records_the_status_ceiling():
    path = REPORT_DATA / "h1_ariana_trough.json"
    if not path.exists():
        pytest.skip("report inputs not generated")
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["objective"] == "H1"
    assert "exploratory" in payload["status_ceiling"]
    assert payload["prohibited_predictor_absent"] == "lnMR OH-metoprolol/metoprolol"
    assert "excluded by construction" in payload["map_be"]
    # The endpoint must not be relabelled as an exposure metric.
    assert "NOT an AUC" in payload["target_is"]
