#!/usr/bin/env python
"""Post-hoc diagnostic: why power peaks at moderate redundancy, and the fix.

This is **not** part of the frozen cohort sizing. It was written after that run
produced a result that looked like noise and was not: under the frozen
pipeline, the power to detect an incremental layer is **non-monotone** in the
layer's redundancy with the baseline, peaking around 0.3.

The explanation, and it is a property of the method rather than of biology:
the frozen pipeline selects candidate features by their absolute correlation
with **the outcome**. A feature that is redundant with the baseline correlates
with the outcome *through the baseline's own signal*, so redundancy makes a
feature **easier to select** while making what it carries **less incremental**.
The two effects pull in opposite directions and their product has an interior
maximum.

The fix is one line in a future protocol: select on the **residual after the
baseline**, not on the raw outcome. This script measures both criteria so the
recommendation rests on a number.

Usage:
    python -I scripts/diagnose_selection_criterion_and_redundancy.py \
        --out reports/data/selection_criterion_diagnostic.json
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.stdout.reconfigure(encoding="utf-8")

from omics2pk import cohort_size_simulation as CSS  # noqa: E402

REDUNDANCY_GRID = (0.0, 0.3, 0.6, 0.85)
POWER_SIZES = (800, 1000)
POWER_REPLICATES = 500
SELECTION_REPLICATES = 300
SELECTION_N = 1000

#: The frozen question-B reference scenario, unchanged.
REFERENCE = {
    "baseline_r2": 0.3, "target_mae_reduction": 0.1, "n_features": 200,
    "n_informative": 10, "redundancy": 0.3, "analytical_cv": 0.3,
    "missing_rate": 0.05, "selected": 10,
}


def power_by_redundancy():
    """Power at FIXED n, which is stable where a bisection on n is not."""
    out = {}
    for n in POWER_SIZES:
        row = {}
        for redundancy in REDUNDANCY_GRID:
            scenario = dict(REFERENCE)
            scenario["redundancy"] = redundancy
            row[f"redundancy_{redundancy}"] = round(
                CSS.power_at(n, CSS.simulate_incremental_layer, scenario,
                             replicates=POWER_REPLICATES), 3)
        out[f"n_{n}"] = row
    return out


def selection_quality():
    """How often each criterion picks genuinely informative features."""
    out = {}
    for redundancy in REDUNDANCY_GRID:
        rng = np.random.default_rng(5)
        on_target, on_residual = [], []
        for _ in range(SELECTION_REPLICATES):
            n = SELECTION_N
            baseline_r2 = REFERENCE["baseline_r2"]
            informative = REFERENCE["n_informative"]
            width = REFERENCE["n_features"]
            residual_share = 1.0 - baseline_r2
            layer_share = residual_share * (
                1.0 - (1.0 - REFERENCE["target_mae_reduction"]) ** 2)
            base_signal = rng.normal(size=n)
            layer_signal = (redundancy * base_signal
                            + math.sqrt(1.0 - redundancy ** 2)
                            * rng.normal(size=n))
            outcome = (math.sqrt(baseline_r2) * base_signal
                       + math.sqrt(layer_share) * layer_signal
                       + math.sqrt(max(residual_share - layer_share, 0.0))
                       * rng.normal(size=n))
            loading = 1.0 / math.sqrt(informative)
            layer = rng.standard_normal((n, width))
            layer[:, :informative] *= math.sqrt(1.0 - loading ** 2)
            layer[:, :informative] += loading * layer_signal[:, None]
            layer += REFERENCE["analytical_cv"] * rng.standard_normal(
                (n, width))
            baseline_block = (base_signal
                              + 0.10 * rng.normal(size=n)).reshape(-1, 1)

            picked = CSS._select(layer, outcome, REFERENCE["selected"])
            on_target.append(float(np.mean(picked < informative)))

            design = np.c_[np.ones(n), baseline_block]
            coefficients = np.linalg.lstsq(design, outcome, rcond=None)[0]
            residual = outcome - design @ coefficients
            picked_residual = CSS._select(layer, residual,
                                          REFERENCE["selected"])
            on_residual.append(float(np.mean(picked_residual < informative)))

        out[f"redundancy_{redundancy}"] = {
            "informative_fraction_selecting_on_the_outcome":
                round(float(np.mean(on_target)), 3),
            "informative_fraction_selecting_on_the_residual":
                round(float(np.mean(on_residual)), 3),
        }
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",
                        default=str(ROOT / "reports" / "data"
                                    / "selection_criterion_diagnostic.json"))
    args = parser.parse_args(argv)

    power = power_by_redundancy()
    selection = selection_quality()

    at_1000 = power["n_1000"]
    peak = max(at_1000, key=at_1000.get)
    on_outcome = [selection[f"redundancy_{r}"][
        "informative_fraction_selecting_on_the_outcome"]
        for r in REDUNDANCY_GRID]
    on_residual = [selection[f"redundancy_{r}"][
        "informative_fraction_selecting_on_the_residual"]
        for r in REDUNDANCY_GRID]

    payload = {
        "diagnostic": "selection_criterion_and_redundancy",
        "generated": "2026-11-10",
        "status": "POST_HOC, written after the frozen sizing run",
        "not_part_of": "docs/OMICS2PK_COHORT_SIZING_PRESPECIFICATION_20261110.md",
        "reference_scenario": REFERENCE,
        "power_at_fixed_n_by_redundancy": power,
        "power_peaks_at": peak,
        "power_is_non_monotone_in_redundancy": bool(
            peak != f"redundancy_{REDUNDANCY_GRID[0]}"
            and peak != f"redundancy_{REDUNDANCY_GRID[-1]}"),
        "selection_quality_by_criterion": selection,
        "selecting_on_the_outcome_rewards_redundancy": bool(
            on_outcome == sorted(on_outcome)),
        "selecting_on_the_residual_penalises_redundancy": bool(
            on_residual == sorted(on_residual, reverse=True)),
        "the_finding": (
            "Selecting candidate features by their absolute correlation with "
            "the OUTCOME rewards redundancy with the baseline: a redundant "
            "feature correlates with the outcome through the baseline's own "
            "signal, so it is easier to select while carrying less "
            "incremental information. Selection accuracy therefore RISES with "
            "redundancy under that criterion and FALLS under a "
            "residual-based one, which is the behaviour an incremental-value "
            "question actually wants."),
        "the_recommendation": (
            "A protocol that tests whether a layer ADDS to a baseline should "
            "select the layer's features on the RESIDUAL after the baseline, "
            "inside the training folds, not on the raw outcome. Otherwise the "
            "selection step is confounded with the question being asked, and "
            "a layer can look most useful precisely where it is most "
            "redundant."),
        "what_this_is_not": (
            "not a biological result. It is a property of a feature-selection "
            "criterion, and it would apply to any wide block in any domain."),
        "no_individual_level_data_in_this_artefact": True,
    }

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")

    print("power at fixed n, by redundancy with the baseline")
    for size, row in power.items():
        print(f"  {size}: " + "  ".join(f"{k.split('_')[1]}: {v}"
                                        for k, v in row.items()))
    print(f"  power peaks at {peak}, non-monotone: "
          f"{payload['power_is_non_monotone_in_redundancy']}")
    print("\nfraction of selected features that are genuinely informative")
    print(f"  {'redundancy':>12} {'on the outcome':>16} {'on the residual':>18}")
    for redundancy in REDUNDANCY_GRID:
        cell = selection[f"redundancy_{redundancy}"]
        print(f"  {redundancy:>12} "
              f"{cell['informative_fraction_selecting_on_the_outcome']:>16} "
              f"{cell['informative_fraction_selecting_on_the_residual']:>18}")
    print(f"\nartefact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
