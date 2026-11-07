"""H5 — biological, analytical and temporal robustness of the solanidine markers.

Executes the H5 items frozen in OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL
v1.0.0 §8, on the Mueller 2026 standardized-potato-meal cohort
(doi 10.3389/fphar.2026.1886427), with an optional direction-consistency check
against ARIANA (doi 10.3389/fphar.2026.1812007).

The criteria were frozen before this dataset was known to exist, which is
verifiable from the v1.0.0 commit. Nothing here is tuned on the results.

What this measures
------------------
1. Intra-individual temporal stability of the two metabolic ratios across the
   nine sampling times: within-subject geometric CV, and the intraclass
   correlation (between-subject share of total variance).
2. Rank stability: Kendall's W across timepoints, and Spearman agreement
   between each timepoint and the 3 h reference time.
3. Association with the genotype-derived CYP2D6 activity score at each time.
4. Food effect: 3-day potato fasting (15 subjects) versus not (5).
5. Analytical-quality sensitivity: the same quantities recomputed with
   below-LOQ values excluded, and with the single poor metaboliser excluded.
6. Direction consistency of the ratio-versus-activity-score association
   between the two cohorts.

What this must never claim
--------------------------
No drug is administered in this cohort, so NOTHING here is pharmacokinetics and
nothing replicates a PK association. Nine samples after ONE meal cannot show
between-day reproducibility. The cohort's only poor metaboliser is below the
limit of quantification at every timepoint for both metabolites, so no
conclusion about poor-metaboliser discrimination is drawn. Shared authors and
platform with ARIANA mean agreement between the two could be method bias.

Only aggregate statistics are printed; no individual measurement is emitted.

Usage
-----
    python scripts/run_h5_solanidine_marker_robustness.py \
        --muller PATH [--ariana PATH]
"""
from __future__ import annotations

import argparse
import json
import sys

MULLER_SHA256 = "68d53071db79f961a63e433137a8a3e7dbd40daf8c8ad0f9f7f1b268bb3a8e07"
ARIANA_SHA256 = "7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"

SHEET_PARTICIPANTS = "Participant data"
SHEET_PLASMA = "Plasma concentration data"

#: First data row (1-indexed) of each stacked analyte block.
BLOCK_FIRST_ROW = {"solanidine": 3, "oh4_solanidine": 26, "ssda": 49}
TIMEPOINTS_H = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 29.0)
N_PARTICIPANTS = 20

#: Reference time for pairwise rank agreement. 3 h is the protocol's choice of
#: a mid-range post-meal time; it is NOT selected on the results.
REFERENCE_TIME_H = 3.0

COL_PARTICIPANT = "Participant-No."
COL_ACTIVITY_SCORE = "CYP2D6 activity score"
COL_PHENOTYPE = "CYP2D6 genotype-predicted phenotype"
COL_FASTED = "Fasted on potato products for 3 days (Yes=1, No=0)"

#: The two ratios, as defined in ARIANA: metabolite over parent solanidine.
RATIOS = {
    "oh4_solanidine_over_solanidine": ("oh4_solanidine", "solanidine"),
    "ssda_over_solanidine": ("ssda", "solanidine"),
}

ARIANA_COL_AS = "AS"
ARIANA_COL_OH = "lnMR OH-solanidine/solanidine"
ARIANA_COL_SSDA = "lnMR SSDA/solanidine"


class AnalysisError(RuntimeError):
    pass


def sha256_of(path):
    import hashlib
    from pathlib import Path

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify(path, expected):
    actual = sha256_of(path)
    if actual != expected:
        raise AnalysisError(
            f"checksum mismatch for {path}\n  expected {expected}\n  actual   {actual}"
        )
    return actual


# --------------------------------------------------------------------------
# source reading
# --------------------------------------------------------------------------
def read_muller(path):
    """Return (participants frame, values dict, imputed dict).

    ``values[analyte]`` and ``imputed[analyte]`` are 20 x 9 arrays, the latter
    boolean: the source flags below-LOQ imputation with bold font only.
    """
    import numpy as np
    import pandas as pd
    from openpyxl import load_workbook

    participants = pd.read_excel(path, sheet_name=SHEET_PARTICIPANTS)
    for column in (COL_PARTICIPANT, COL_ACTIVITY_SCORE, COL_PHENOTYPE, COL_FASTED):
        if column not in participants.columns:
            raise AnalysisError(f"participant sheet is missing column {column!r}")
    if len(participants) != N_PARTICIPANTS:
        raise AnalysisError(f"expected {N_PARTICIPANTS} participants")

    worksheet = load_workbook(path, data_only=True)[SHEET_PLASMA]
    values, imputed = {}, {}
    for analyte, first_row in BLOCK_FIRST_ROW.items():
        matrix = np.full((N_PARTICIPANTS, len(TIMEPOINTS_H)), np.nan)
        flags = np.zeros_like(matrix, dtype=bool)
        ids = []
        for row_index in range(N_PARTICIPANTS):
            row = first_row + row_index
            ids.append(worksheet.cell(row=row, column=1).value)
            for time_index in range(len(TIMEPOINTS_H)):
                cell = worksheet.cell(row=row, column=2 + time_index)
                if cell.value is not None:
                    matrix[row_index, time_index] = float(cell.value)
                    flags[row_index, time_index] = bool(
                        cell.font is not None and cell.font.bold
                    )
        if ids != list(range(1, N_PARTICIPANTS + 1)):
            raise AnalysisError(f"unexpected participant order in block {analyte!r}")
        values[analyte] = matrix
        imputed[analyte] = flags
    return participants, values, imputed


def build_ratios(values, imputed):
    """Return ln-ratio matrices and a matrix flagging any imputed input."""
    import numpy as np

    ln_ratio, ratio_imputed = {}, {}
    for name, (numerator, denominator) in RATIOS.items():
        with np.errstate(divide="ignore", invalid="ignore"):
            ln_ratio[name] = np.log(values[numerator] / values[denominator])
        ratio_imputed[name] = imputed[numerator] | imputed[denominator]
    return ln_ratio, ratio_imputed


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def _finite_rows(matrix, mask=None):
    """Yield (row_index, finite values) for rows with at least 2 usable values."""
    import numpy as np

    usable = np.isfinite(matrix)
    if mask is not None:
        usable &= ~mask
    for row_index in range(matrix.shape[0]):
        row = matrix[row_index][usable[row_index]]
        if row.size >= 2:
            yield row_index, row


def within_subject_stability(ln_matrix, mask=None):
    """Within-subject dispersion of a ln-scale marker, and the ICC."""
    import numpy as np

    per_subject_sd, per_subject_n = [], []
    for _, row in _finite_rows(ln_matrix, mask):
        per_subject_sd.append(float(np.std(row, ddof=1)))
        per_subject_n.append(int(row.size))
    if len(per_subject_sd) < 3:
        return {"n_subjects": len(per_subject_sd), "insufficient": True}

    sd = np.asarray(per_subject_sd)
    # Pooled within-subject SD, weighted by degrees of freedom.
    dof = np.asarray(per_subject_n) - 1
    within_var = float(np.sum(dof * sd**2) / np.sum(dof))

    # Between-subject variance from the subject means.
    means = []
    for _, row in _finite_rows(ln_matrix, mask):
        means.append(float(np.mean(row)))
    means = np.asarray(means)
    mean_var = float(np.var(means, ddof=1))
    harmonic_n = float(len(per_subject_n) / np.sum(1.0 / np.asarray(per_subject_n)))
    between_var = max(mean_var - within_var / harmonic_n, 0.0)

    total = between_var + within_var
    return {
        "n_subjects": int(len(per_subject_sd)),
        "median_within_subject_sd_ln": float(np.median(sd)),
        "max_within_subject_sd_ln": float(np.max(sd)),
        # Geometric CV implied by the pooled within-subject ln SD.
        "within_subject_geometric_cv_percent": float(
            100.0 * np.sqrt(np.expm1(within_var))
        ),
        "between_subject_sd_ln": float(np.sqrt(between_var)),
        "icc_between_share": float(between_var / total) if total > 0 else None,
    }


def kendalls_w(ln_matrix, mask=None):
    """Concordance of subject rankings across timepoints.

    Columns are timepoints, rows subjects. Only timepoints with a complete set
    of usable subjects are included, so the statistic is computed on a balanced
    block rather than on shifting subsets.
    """
    import numpy as np
    from scipy import stats

    usable = np.isfinite(ln_matrix)
    if mask is not None:
        usable &= ~mask
    subject_ok = usable.all(axis=1)
    block = ln_matrix[subject_ok]
    n_subjects, n_times = block.shape
    if n_subjects < 3 or n_times < 2:
        return {"insufficient": True, "n_subjects": int(n_subjects)}

    ranks = np.apply_along_axis(stats.rankdata, 0, block)  # rank subjects per time
    rank_sums = ranks.sum(axis=1)
    mean_rank_sum = rank_sums.mean()
    s = float(((rank_sums - mean_rank_sum) ** 2).sum())
    denominator = n_times**2 * (n_subjects**3 - n_subjects) / 12.0
    w = s / denominator if denominator > 0 else None
    return {
        "n_subjects": int(n_subjects),
        "n_timepoints": int(n_times),
        "kendalls_w": float(w) if w is not None else None,
    }


def _spearman(x, y):
    import numpy as np
    from scipy import stats

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    n = int(keep.sum())
    if n < 5:
        return {"n": n, "insufficient": True}
    rho, p = stats.spearmanr(x[keep], y[keep])
    # Fisher z interval on the rank correlation; descriptive at this n.
    if abs(rho) < 1.0:
        z = np.arctanh(rho)
        se = 1.0 / np.sqrt(n - 3)
        low, high = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    else:
        low = high = float(rho)
    return {
        "n": n,
        "rho": float(rho),
        "ci_low": float(low),
        "ci_high": float(high),
        "p": float(p),
    }


def agreement_between_times(ln_matrix, mask=None, reference_h=REFERENCE_TIME_H):
    """Spearman agreement of subject ranking at each time with a reference time.

    This replaces an earlier statistic that correlated one timepoint against the
    mean across timepoints, which was circular because the mean contained it.
    """
    import numpy as np

    reference_index = TIMEPOINTS_H.index(reference_h)
    reference = ln_matrix[:, reference_index].copy()
    if mask is not None:
        reference[mask[:, reference_index]] = np.nan
    out = {}
    for time_index, time in enumerate(TIMEPOINTS_H):
        if time_index == reference_index:
            continue
        column = ln_matrix[:, time_index].copy()
        if mask is not None:
            column[mask[:, time_index]] = np.nan
        out[f"{time:g}h"] = _spearman(column, reference)
    return out


def association_by_time(ln_matrix, activity_score, mask=None, keep_subjects=None):
    """Spearman association with the activity score at each timepoint."""
    import numpy as np

    out = {}
    for time_index, time in enumerate(TIMEPOINTS_H):
        column = ln_matrix[:, time_index].copy()
        if mask is not None:
            column[mask[:, time_index]] = np.nan
        score = np.asarray(activity_score, dtype=float).copy()
        if keep_subjects is not None:
            column = column[keep_subjects]
            score = score[keep_subjects]
        out[f"{time:g}h"] = _spearman(column, score)
    return out


def food_effect(ln_matrix, fasted_flag, mask=None, keep_subjects=None):
    """Descriptive comparison of the subject-level mean ln ratio by fasting."""
    import numpy as np
    from scipy import stats

    usable = np.isfinite(ln_matrix)
    if mask is not None:
        usable &= ~mask
    means = np.full(ln_matrix.shape[0], np.nan)
    for row_index in range(ln_matrix.shape[0]):
        row = ln_matrix[row_index][usable[row_index]]
        if row.size:
            means[row_index] = float(np.mean(row))

    flag = np.asarray(fasted_flag, dtype=float)
    if keep_subjects is not None:
        means = means[keep_subjects]
        flag = flag[keep_subjects]

    fasted = means[(flag == 1) & np.isfinite(means)]
    not_fasted = means[(flag == 0) & np.isfinite(means)]
    if fasted.size < 2 or not_fasted.size < 2:
        return {"insufficient": True, "n_fasted": int(fasted.size),
                "n_not_fasted": int(not_fasted.size)}

    difference = float(np.mean(fasted) - np.mean(not_fasted))
    se = float(
        np.sqrt(np.var(fasted, ddof=1) / fasted.size
                + np.var(not_fasted, ddof=1) / not_fasted.size)
    )
    result = stats.ttest_ind(fasted, not_fasted, equal_var=False)
    return {
        "n_fasted": int(fasted.size),
        "n_not_fasted": int(not_fasted.size),
        "mean_difference_ln": difference,
        "fold_difference": float(np.exp(difference)),
        "descriptive_95_interval_ln": [difference - 1.96 * se, difference + 1.96 * se],
        "welch_p": float(result.pvalue),
        "interval_excludes_zero": bool(
            (difference - 1.96 * se) * (difference + 1.96 * se) > 0
        ),
        "warning": "n=5 in the non-fasted group; this is descriptive only.",
    }


# --------------------------------------------------------------------------
# ARIANA direction check
# --------------------------------------------------------------------------
def ariana_direction(path):
    """Spearman of each ln ratio against the activity score in ARIANA."""
    import pandas as pd

    verify(path, ARIANA_SHA256)
    frame = pd.read_excel(path)
    frame.columns = [str(c).strip() for c in frame.columns]
    for column in (ARIANA_COL_AS, ARIANA_COL_OH, ARIANA_COL_SSDA):
        if column not in frame.columns:
            raise AnalysisError(f"ARIANA workbook is missing column {column!r}")
    return {
        "oh4_solanidine_over_solanidine": _spearman(
            frame[ARIANA_COL_OH], frame[ARIANA_COL_AS]
        ),
        "ssda_over_solanidine": _spearman(frame[ARIANA_COL_SSDA], frame[ARIANA_COL_AS]),
        "note": (
            "ARIANA ratios are already on the ln scale in the source. The "
            "activity-score convention differs between cohorts (*41 scores 0.25 "
            "in Mueller), so only the DIRECTION and rough magnitude are compared, "
            "never the numeric value."
        ),
    }


# --------------------------------------------------------------------------
# orchestration
# --------------------------------------------------------------------------
def run(muller_path, ariana_path=None):
    import numpy as np

    checksum = verify(muller_path, MULLER_SHA256)
    participants, values, imputed = read_muller(muller_path)
    ln_ratio, ratio_imputed = build_ratios(values, imputed)

    activity_score = participants[COL_ACTIVITY_SCORE].to_numpy(dtype=float)
    fasted = participants[COL_FASTED].to_numpy(dtype=float)
    phenotype = participants[COL_PHENOTYPE].astype(str).to_numpy()
    not_poor = phenotype != "poor metabolizer"

    report = {
        "objective": "H5",
        "cohort": "muller_potato_meal_20",
        "source_sha256": checksum,
        "endpoint_class": "biomarker_only",
        "is_pharmacokinetics": False,
        "prespecification": (
            "H5 items frozen in protocol v1.0.0 before this dataset was known to "
            "exist. Reference time 3 h and Spearman as the primary association "
            "measure are protocol choices, not result-driven ones."
        ),
        "analyses": {},
    }

    for name, matrix in ln_ratio.items():
        mask = ratio_imputed[name]
        n_imputed = int(mask.sum())
        scenarios = {
            # As published: imputed floors retained, as the source supplies them.
            "all_values_including_imputed": dict(mask=None, keep=None),
            # Analytical-quality sensitivity: drop every below-LOQ input.
            "excluding_below_loq": dict(mask=mask, keep=None),
            # The single poor metaboliser is wholly at the floor.
            "excluding_poor_metaboliser": dict(mask=None, keep=not_poor),
        }
        entry = {
            "ratio_definition": "/".join(RATIOS[name]),
            "cells_with_any_imputed_input": n_imputed,
            "cells_total": int(mask.size),
            "scenarios": {},
        }
        for label, options in scenarios.items():
            keep = options["keep"]
            sub_mask = options["mask"]
            matrix_used = matrix if keep is None else matrix[keep]
            mask_used = (
                None if sub_mask is None
                else (sub_mask if keep is None else sub_mask[keep])
            )
            entry["scenarios"][label] = {
                "temporal_stability": within_subject_stability(matrix_used, mask_used),
                "rank_concordance": kendalls_w(matrix_used, mask_used),
                "agreement_with_3h_ranking": agreement_between_times(
                    matrix_used, mask_used, reference_h=REFERENCE_TIME_H
                ),
                "agreement_with_29h_ranking": agreement_between_times(
                    matrix_used, mask_used, reference_h=29.0
                ),
                "association_with_activity_score": association_by_time(
                    matrix, activity_score, mask=sub_mask, keep_subjects=keep
                ),
                "food_effect": food_effect(
                    matrix, fasted, mask=sub_mask, keep_subjects=keep
                ),
            }
        report["analyses"][name] = entry

    if ariana_path:
        report["ariana_direction_check"] = ariana_direction(ariana_path)

    report["must_not_be_concluded"] = [
        "Nothing here is pharmacokinetics; no drug was administered.",
        "Nine samples after ONE meal cannot establish between-day reproducibility.",
        "The only poor metaboliser is below the limit of quantification at every "
        "timepoint for both metabolites, so poor-metaboliser discrimination is not "
        "assessed.",
        "Agreement with ARIANA may be shared method bias: the cohorts share "
        "authors, laboratory and analytical platform.",
        "Activity-score values are not numerically comparable between cohorts.",
    ]
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--muller", required=True, help="Mueller 2026 workbook")
    parser.add_argument("--ariana", help="ARIANA workbook, for the direction check")
    args = parser.parse_args(argv)
    try:
        report = run(args.muller, args.ariana)
    except AnalysisError as exc:
        parser.exit(2, f"INVALID: {exc}\n")
    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
