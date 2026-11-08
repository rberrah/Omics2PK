"""Hybrid selection: knowledge proposes a shortlist, training data disposes.

Why this exists, and where the idea came from
---------------------------------------------
The H2a run on one endpoint showed that a purely mechanistic ordering keeps a
biologically plausible variable even after its information has been destroyed.
With the target transcript's donor order permuted, the outcome-blind arm went on
choosing it and lost its whole advantage, landing on the clinical baseline,
while a selector that reads the endpoint chose something else and recovered half
the gap. Mechanistic relevance does not imply predictive relevance.

This module is the answer to that, and its provenance must travel with it: **the
rule was designed after seeing that result, on one endpoint.** It is therefore
hypothesis-generating with respect to that endpoint and can only be confirmatory
on endpoints that did not generate it. The prespecification document states
which, and the rule must not be presented as having been frozen before the
observation that motivated it.

The rule
--------
1. **Knowledge proposes.** Take the top `oversample * k` genes of the
   outcome-blind mechanistic ordering. Nothing here reads the endpoint.
2. **Data disposes, in training folds only.** For each shortlisted gene,
   measure the information it adds over the clinical baseline as the absolute
   partial rank correlation with the endpoint, after both have been residualised
   on the baseline covariates.
3. **Calibrate.** Build a null by permuting the residualised endpoint among
   training donors, and retain genes with the permutation false-discovery
   estimator already used elsewhere in this project -- the one with no
   resolution floor, since an empirical p-value from `B` permutations cannot
   fall below `1/(B+1)`.
4. **Keep, demote or abstain.** Retain at most `k`. If nothing survives, retain
   nothing and let the model fall back to the clinical covariates. Abstention is
   an allowed answer, not a failure to produce one.

What this changes about the comparison
--------------------------------------
The hybrid reads the endpoint, so its selection varies between folds. Its
Jaccard stability is no longer 1.0 by construction, which is the first time the
stability comparison against a statistical selector means anything. That is a
gain in interpretability and a loss of the structural advantage the blind
ordering had on that one metric; both must be reported.

The most likely negative outcome is stated in advance: if the shortlist is wide
enough, the data step does all the work and the hybrid collapses onto the
univariate selector, which would mean the knowledge contributed nothing beyond
shortlisting. The overlap with the univariate choice is therefore reported
alongside every result, and the oversampling factor is kept small.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from omics2pk.mechanistic_graph_selection import (
    MechanisticSelectionError,
    permutation_fdr_cutoff,
)

FROZEN = {
    "oversample": 3,
    "null_permutations": 200,
    "target_fdr": 0.05,
    "permutation_seed": 20261108,
}


class HybridSelectionError(ValueError):
    """Raised when a shortlist, a design or a parameter is unusable."""


@dataclass(frozen=True)
class HybridSelection:
    """What survived, what was demoted, and the evidence for both."""

    shortlist: tuple[str, ...]
    retained: tuple[str, ...]
    demoted: tuple[str, ...]
    statistics: Mapping[str, float]
    abstained: bool
    budget: int
    parameters: Mapping[str, object]
    how_to_read: str = (
        "The shortlist is proposed by knowledge alone. Retention is decided on "
        "training donors only. An empty retained set is an abstention, which is "
        "an allowed answer; it is not a failure to select."
    )


def _rank(values: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    return rankdata(values, axis=0).astype(np.float64)


def _residualise(target: np.ndarray, design: np.ndarray) -> np.ndarray:
    """Remove the baseline's linear part, by least squares with an intercept."""
    augmented = np.column_stack([np.ones(design.shape[0]), design])
    coefficients, *_ = np.linalg.lstsq(augmented, target, rcond=None)
    return target - augmented @ coefficients


def partial_rank_strength(matrix: np.ndarray, endpoint: np.ndarray,
                          baseline: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Absolute partial rank correlation of each column with the endpoint.

    Both sides are rank-transformed and then residualised on the baseline, so
    the statistic measures information the clinical covariates do not already
    carry. Returns the strengths and the residualised endpoint, the latter so a
    permutation null can reuse it instead of destroying the adjustment.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    endpoint = np.asarray(endpoint, dtype=np.float64)
    baseline = np.asarray(baseline, dtype=np.float64)
    if matrix.ndim != 2:
        raise HybridSelectionError("The candidate block must be two-dimensional")
    if endpoint.shape[0] != matrix.shape[0] or baseline.shape[0] != matrix.shape[0]:
        raise HybridSelectionError("Endpoint, baseline and candidates must align")
    if matrix.shape[1] == 0:
        return np.zeros(0, dtype=np.float64), np.zeros(endpoint.shape[0])

    # The adjustment happens in the same space as what it adjusts. Residualising
    # a rank-transformed variable on raw covariates leaves curvature behind, and
    # a variable that merely repeats a baseline covariate then keeps a visible
    # partial correlation it has not earned.
    design = _rank(baseline)
    ranked = _rank(matrix)
    residual_columns = np.column_stack(
        [_residualise(ranked[:, index], design) for index in range(ranked.shape[1])])
    residual_endpoint = _residualise(_rank(endpoint), design)
    return _correlate(residual_columns, residual_endpoint), residual_endpoint


def _correlate(columns: np.ndarray, target: np.ndarray) -> np.ndarray:
    centred = columns - columns.mean(axis=0)
    aim = target - target.mean()
    denominator = np.sqrt((centred ** 2).sum(axis=0) * (aim ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = (centred * aim[:, None]).sum(axis=0) / denominator
    return np.nan_to_num(np.abs(correlation))


def select_hybrid(ordering: Sequence[str], available: Mapping[str, int],
                  matrix: np.ndarray, endpoint: np.ndarray, baseline: np.ndarray,
                  *, budget: int,
                  oversample: int | None = None,
                  null_permutations: int | None = None,
                  q: float | None = None,
                  random_seed: int | None = None) -> HybridSelection:
    """Apply the rule. `matrix`, `endpoint` and `baseline` hold training donors.

    `ordering` is the outcome-blind mechanistic ranking; `available` maps a gene
    name to its column in `matrix`. Genes absent from `available` are skipped
    without comment, exactly as in the blind arm.
    """
    factor = FROZEN["oversample"] if oversample is None else int(oversample)
    permutations = (FROZEN["null_permutations"] if null_permutations is None
                    else int(null_permutations))
    rate = FROZEN["target_fdr"] if q is None else float(q)
    seed = FROZEN["permutation_seed"] if random_seed is None else int(random_seed)
    budget = int(budget)
    if budget < 0:
        raise HybridSelectionError("A budget cannot be negative")
    if factor < 1:
        raise HybridSelectionError("The oversampling factor must be at least one")
    if permutations < 1:
        raise HybridSelectionError("At least one permutation is required")

    shortlist = [name for name in ordering if name in available][: factor * budget]
    parameters = {
        "oversample": factor,
        "null_permutations": permutations,
        "target_fdr": rate,
        "permutation_seed": seed,
        "budget": budget,
        "shortlist_size": len(shortlist),
    }
    if not shortlist or budget == 0:
        return HybridSelection(shortlist=tuple(shortlist), retained=(),
                               demoted=tuple(shortlist), statistics={},
                               abstained=True, budget=budget,
                               parameters=parameters)

    columns = [available[name] for name in shortlist]
    block = matrix[:, columns]
    strengths, residual_endpoint = partial_rank_strength(block, endpoint, baseline)

    rng = np.random.default_rng(seed)
    design = _rank(baseline)
    ranked = _rank(block)
    residual_columns = np.column_stack(
        [_residualise(ranked[:, index], design)
         for index in range(ranked.shape[1])])
    null = np.empty((len(shortlist), permutations), dtype=np.float64)
    for draw in range(permutations):
        # Permuting the residualised endpoint reintroduces a baseline component
        # that the residualised columns are orthogonal to, which would shrink
        # the null and make the test anti-conservative. Re-residualising gives
        # the null statistic the same geometry as the observed one.
        shuffled = _residualise(
            residual_endpoint[rng.permutation(residual_endpoint.size)], design)
        null[:, draw] = _correlate(residual_columns, shuffled)

    retained_count, _ = permutation_fdr_cutoff(strengths, null, rate)
    retained_count = min(retained_count, budget)
    order = np.argsort(-strengths, kind="stable")
    retained = tuple(shortlist[int(index)] for index in order[:retained_count])
    demoted = tuple(name for name in shortlist if name not in set(retained))
    return HybridSelection(
        shortlist=tuple(shortlist),
        retained=retained,
        demoted=demoted,
        statistics={name: float(strengths[position])
                    for position, name in enumerate(shortlist)},
        abstained=bool(retained_count == 0),
        budget=budget,
        parameters=parameters,
    )


def collapse_onto(selection: HybridSelection, statistical_choice: Sequence[str]
                  ) -> float:
    """Fraction of the hybrid's retained set that a statistical arm also picked.

    Reported with every result because it measures the stated failure mode: a
    hybrid whose choices coincide with the statistical selector's has shown that
    the knowledge contributed nothing beyond shortlisting.
    """
    retained = set(selection.retained)
    if not retained:
        return float("nan")
    return len(retained & set(statistical_choice)) / len(retained)


__all__ = [
    "FROZEN",
    "HybridSelection",
    "HybridSelectionError",
    "MechanisticSelectionError",
    "collapse_onto",
    "partial_rank_strength",
    "select_hybrid",
]
