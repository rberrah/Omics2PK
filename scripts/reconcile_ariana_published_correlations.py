"""Reproduce the ARIANA source publication's correlations from its public data.

The source publication releases its data openly:
  https://github.com/katjajust/Solandine-Metoprolol-Correlations
  Solanidine_metoprolol_ratios_PGx.xlsx
  sha256 7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37

This script verifies that checksum, reproduces the published correlation
values, and quantifies three implementation questions that the Phase 1 audit
answers:

1. the effect of the 2.13 relative-activity cap on the SSDA arm;
2. how much of a genotype-only activity estimate's accuracy comes from the
   nonlinear activity mapping versus the phenoconversion correction;
3. whether dose normalisation changes the marker-endpoint association.

No patient-level value is printed: only correlation statistics and counts.
The workbook is never written into the repository (.gitignore excludes *.xlsx).

Usage
-----
    python scripts/reconcile_ariana_published_correlations.py --source PATH
    python scripts/reconcile_ariana_published_correlations.py --download DIR
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

DATA_URL = (
    "https://raw.githubusercontent.com/katjajust/"
    "Solandine-Metoprolol-Correlations/HEAD/Solanidine_metoprolol_ratios_PGx.xlsx"
)
EXPECTED_SHA256 = "7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"

#: Values printed in the source publication, used as assertion targets.
PUBLISHED = {
    "oh_solanidine_vs_reference_mr": {"r": 0.907, "r2": 0.823},
    "ssda_vs_reference_mr": {"r": 0.946, "r2": 0.896},
    "activity_score_vs_reference_mr": {"r": 0.500, "r2": 0.250},
    "oh_solanidine_vs_trough": {"r": -0.513, "r2": 0.263},
    "ssda_vs_trough": {"r": -0.533, "r2": 0.284},
    "activity_score_vs_trough": {"r": -0.288, "r2": 0.083},
}
TOLERANCE = 0.002

COL_REFERENCE_MR = "lnMR OH-metoprolol/metoprolol"
COL_OH = "lnMR OH-solanidine/solanidine"
COL_SSDA = "lnMR SSDA/solanidine"
COL_TROUGH = "c(Metoprolol)[ng/ml]"
COL_TROUGH_WEIGHTED = "c(metoprolol)weighted [ng/ml]"
COL_AS = "AS"
COL_INHIBITOR = "CYP2D6 inhibitor intake"

# Rudesheim et al. 2020, Table 3: CYP2D6 relative activity by activity score.
AS_GRID = (0.0, 0.5, 1.25, 1.5, 2.0, 3.0)
KCAT_REL = (0.0, 0.19, 0.64, 0.72, 1.00, 2.13)
ACTIVITY_CAP = 2.13
#: Mean activity-score reduction per CYP2D6 substrate/inhibitor (external n=88).
DELTA_AS_INHIBITOR = 0.54
#: Genotype-normal activity-score window used to set the omic reference level.
NM_WINDOW = (1.25, 2.25)


class ReconciliationError(RuntimeError):
    pass


def sha256_of(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_source(path):
    """Return the checksum, raising if it does not match the published file."""
    actual = sha256_of(path)
    if actual != EXPECTED_SHA256:
        raise ReconciliationError(
            f"checksum mismatch for {path}\n  expected {EXPECTED_SHA256}\n  actual   {actual}"
        )
    return actual


def download(target_dir):
    """Fetch the public workbook into *target_dir* and verify its checksum."""
    from urllib.request import urlopen

    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / "Solanidine_metoprolol_ratios_PGx.xlsx"
    with urlopen(DATA_URL, timeout=60) as response:  # noqa: S310 - fixed https URL
        target.write_bytes(response.read())
    verify_source(target)
    return target


def relative_activity_from_as(activity_score):
    """Rudesheim relative CYP2D6 activity for an activity score."""
    import numpy as np

    return np.interp(np.asarray(activity_score, dtype=float), AS_GRID, KCAT_REL)


def omic_activity(ln_mr, reference_ln_mr, *, gamma=1.0, cap=ACTIVITY_CAP):
    """Relative activity implied by a metabolic ratio, optionally capped.

    ``cap=None`` disables capping, which is what the published correlation
    corresponds to.
    """
    import numpy as np

    raw = np.exp(gamma * (np.asarray(ln_mr, dtype=float) - reference_ln_mr))
    return raw if cap is None else np.minimum(raw, cap)


def _pearson(x, y):
    import numpy as np
    from scipy import stats

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    n = int(keep.sum())
    if n < 4:
        raise ReconciliationError("fewer than four usable observations")
    r, p = stats.pearsonr(x[keep], y[keep])
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    return {
        "n": n,
        "r": float(r),
        "ci_low": float(np.tanh(z - 1.96 * se)),
        "ci_high": float(np.tanh(z + 1.96 * se)),
        "r2": float(r * r),
        "p": float(p),
    }


def reconcile(source):
    """Return a JSON-serialisable reconciliation report."""
    import numpy as np
    import pandas as pd

    checksum = verify_source(source)
    frame = pd.read_excel(source)
    frame.columns = [str(column).strip() for column in frame.columns]

    required = [
        COL_REFERENCE_MR, COL_OH, COL_SSDA, COL_TROUGH,
        COL_TROUGH_WEIGHTED, COL_AS, COL_INHIBITOR,
    ]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ReconciliationError(f"source workbook is missing columns: {missing}")

    reference = frame[COL_REFERENCE_MR]
    report = {
        "source_sha256": checksum,
        "n_rows": int(len(frame)),
        "published_reproduction": {},
        "activity_cap": {},
        "genotype_decomposition": {},
        "dose_normalisation": {},
    }

    pairs = {
        "oh_solanidine_vs_reference_mr": (frame[COL_OH], reference),
        "ssda_vs_reference_mr": (frame[COL_SSDA], reference),
        "activity_score_vs_reference_mr": (frame[COL_AS], reference),
        "oh_solanidine_vs_trough": (frame[COL_OH], frame[COL_TROUGH]),
        "ssda_vs_trough": (frame[COL_SSDA], frame[COL_TROUGH]),
        "activity_score_vs_trough": (frame[COL_AS], frame[COL_TROUGH]),
    }
    discrepancies = []
    for name, (left, right) in pairs.items():
        stats_out = _pearson(left, right)
        expected = PUBLISHED[name]
        stats_out["published_r"] = expected["r"]
        stats_out["published_r2"] = expected["r2"]
        stats_out["matches_published"] = (
            abs(stats_out["r"] - expected["r"]) <= TOLERANCE
            and abs(stats_out["r2"] - expected["r2"]) <= TOLERANCE
        )
        if not stats_out["matches_published"]:
            discrepancies.append(name)
        report["published_reproduction"][name] = stats_out
    report["all_published_values_reproduced"] = not discrepancies
    report["discrepancies"] = discrepancies

    # The reference level is the mean ratio of genotype-normal metabolisers
    # without an inhibitor, exactly as the mechanistic module defines it.
    low, high = NM_WINDOW
    normal = (
        frame[COL_AS].between(low, high) & (frame[COL_INHIBITOR] == 0)
    )
    report["reference_group_n"] = int(normal.sum())

    for label, column in (("oh_solanidine", COL_OH), ("ssda", COL_SSDA)):
        reference_level = float(frame.loc[normal, column].mean())
        uncapped = omic_activity(frame[column], reference_level, cap=None)
        capped = omic_activity(frame[column], reference_level)
        report["activity_cap"][label] = {
            "reference_ln_mr": reference_level,
            "n_at_cap": int((uncapped > ACTIVITY_CAP).sum()),
            "uncapped": _pearson(np.log(uncapped), reference),
            "capped": _pearson(np.log(capped), reference),
        }

    report["genotype_decomposition"] = {
        "raw_activity_score": _pearson(frame[COL_AS], reference),
        "nonlinear_mapping_only": _pearson(
            np.log(relative_activity_from_as(frame[COL_AS]) + 1e-6), reference
        ),
        "nonlinear_mapping_plus_phenoconversion": _pearson(
            np.log(
                relative_activity_from_as(
                    np.clip(
                        frame[COL_AS] - DELTA_AS_INHIBITOR * frame[COL_INHIBITOR],
                        0.0,
                        None,
                    )
                )
                + 1e-6
            ),
            reference,
        ),
        "n_with_inhibitor": int(frame[COL_INHIBITOR].sum()),
    }

    for label, column in (("oh_solanidine", COL_OH), ("ssda", COL_SSDA), ("activity_score", COL_AS)):
        report["dose_normalisation"][label] = {
            "unweighted_trough": _pearson(frame[column], frame[COL_TROUGH]),
            "dose_weighted_trough": _pearson(frame[column], frame[COL_TROUGH_WEIGHTED]),
        }

    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--source", help="path to the public workbook")
    group.add_argument("--download", metavar="DIR", help="download it into DIR first")
    args = parser.parse_args(argv)

    try:
        source = download(args.download) if args.download else args.source
        report = reconcile(source)
    except ReconciliationError as exc:
        parser.exit(2, f"INVALID: {exc}\n")

    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))
    return 0 if report["all_published_values_reproduced"] else 1


if __name__ == "__main__":
    sys.exit(main())
