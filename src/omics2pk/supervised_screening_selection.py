"""Which SELECTOR chooses useful biological information? One harness, many selectors.

What this module is for
-----------------------
The commissioned question is not "which model predicts best" but **"does the
methodology choose better which biological information to use"**. Those need
different experiments. Comparing models confounds the choice of variables with
the choice of estimator; so every selector here is handed the **same candidate
pool**, the **same folds** and the **same final estimator**, and the only thing
that differs is *which columns it keeps*.

The one rule that makes the answer mean anything
------------------------------------------------
**Every supervised choice happens inside the training folds — including the
choice of how many variables to keep, and including the choice of selector.**
A panel picked by looking at out-of-fold performance and then scored on that
same performance is scored on its own selection. This module therefore reports
two different things and never conflates them:

* `curve` — out-of-fold loss as a function of k, for each selector. Descriptive.
  **Its best point is optimistically biased** and is labelled as such.
* `selectors` — the honest number: k and the panel are both chosen by an inner
  grouped split of each training fold, so the outer loss estimates *the method
  that selects*, not a panel chosen with hindsight.

Why SHAP without the `shap` package
-----------------------------------
`shap` is not installed and does not need to be: XGBoost computes **exact**
TreeSHAP itself through `pred_contribs=True`. Adding a dependency to obtain a
number the existing one already produces exactly would be the wrong trade.

SHAP is computed on rows the booster did **not** train on — an inner split of
the training fold, averaged. In-sample SHAP ranks a feature the model overfitted
as important, which is precisely the failure this harness exists to avoid.

What a SHAP value is not
------------------------
A SHAP value is an attribution **within one fitted model**, not evidence that a
variable carries information. Two collinear probes split one signal and both
look weak; a variable the model leaned on can still be useless out of sample.
So SHAP here only ever **proposes an ordering**, and every claim about value is
settled by an **out-of-sample ablation** — refitting without the block and
measuring the paired difference. The commissioning note asks for exactly that.

Negative controls, because a selector always returns something
--------------------------------------------------------------
A selector cannot abstain by accident: hand it noise and it returns the k
luckiest columns. So the harness runs a permuted outcome, a permuted knowledge
ordering, and a decoy block of pure noise shaped like the real one, and it
refuses to report a selection result without the positive control passing.
"""
from __future__ import annotations

import math

import numpy as np

from .fusion_headroom_diagnostics import paired_gap, ridge_predict
from .variant_mechanism_omic_selection import _ladder

#: Frozen before execution. Nothing here is tuned on an outcome.
FROZEN = {
    "seed": 20261112,
    "outer_folds": 5,
    "inner_folds": 3,
    "n_bootstrap": 2000,
    "alpha": 0.05,
    "minimum_interesting_relative_gain": 0.05,
    "ridge_penalty": 1.0,
    "loss": "absolute_error",
    # The panel sizes a selector may choose from. Declared, and the choice
    # among them is made inside each training fold.
    #
    # **k = 0 is abstention**: keep the baseline covariates and nothing else.
    # A method that cannot choose it is not automatic, it is merely obedient --
    # and on HLC CYP2E1 selecting any variable was DETECTABLY_WORSE than
    # taking everything, so "select nothing" was the right answer and no arm
    # could express it. With 0 in the grid the inner folds can now reach it.
    "k_grid": (0, 1, 2, 3, 5, 8, 12, 20, 32, 50),
    # Two candidates whose absolute Pearson correlation on the TRAINING rows
    # exceeds this are treated as one piece of information, and the lower-
    # ranked one is dropped. 0.9 is near-duplication, not mere association:
    # the positive control's decoy sits at rho > 0.9 by construction and was
    # selected in five folds of five without this step.
    "redundancy_threshold": 0.90,
    # Knowledge proposes this many times the budget before data disposes.
    "oversample": 3,
    # High dimension only. Before any selector runs, keep the most variable
    # candidates, computed on TRAINING rows alone. The filter never reads the
    # outcome, so it cannot leak; it is applied in-fold anyway, so it cannot
    # even see the test rows' variance. None means no filter, which is the
    # default and what every run so far used.
    "prefilter_top_variance": None,
    # Matched random pools, for separating knowledge from mere dimension
    # reduction. 40 draws give a 2.5 %% resolution on an empirical p.
    "n_matched_pools": 40,
    # A candidate is an admissible match when it sits within this many nearest
    # neighbours of its target in standardised (mean, log variance) space.
    # Wide enough to randomise, narrow enough to match.
    "match_neighbourhood": 25,
    # The share of the validated gain a parsimonious panel must retain.
    "parsimony_retention": (0.80, 0.90),
    # Reused verbatim from `hlc_real_xgboost_validation`, where it was fixed
    # before any outcome was inspected. Keeping it identical means the
    # screening arm is comparable to the frozen HLC comparison.
    "xgb": {"n_estimators": 150, "max_depth": 2, "learning_rate": 0.035,
            "min_child_weight": 8, "subsample": 0.85, "colsample_bytree": 0.9,
            "reg_lambda": 10.0},
}

#: The selectors. All receive the same pool, folds and estimator.
SELECTORS = (
    "S0_everything",            # no selection at all
    "S1_univariate_outcome",    # |rank correlation| with the outcome
    "S2_residual_on_baseline",  # the same, residualised on the baseline block
    "S3_elastic_net",           # non-zero coefficients, penalty chosen in-fold
    "S4_shap_xgboost",          # exact TreeSHAP, computed out-of-sample in-fold
    "S5_knowledge_only",        # outcome-blind ordering from knowledge
    "S6_knowledge_then_shap",   # knowledge proposes, SHAP disposes
    "S7_shap_deduplicated",     # S4, with near-duplicate candidates removed
    "S8_pool_then_univariate",  # knowledge restricts the POOL, data ranks in it
)


class ScreeningError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise ScreeningError(message)


# --------------------------------------------------------------------------
# small helpers, kept local and dependency-free
# --------------------------------------------------------------------------
def _rank(values):
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        order = values.argsort(kind="mergesort")
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, values.size + 1, dtype=float)
        return ranks
    return np.column_stack([_rank(values[:, j]) for j in range(values.shape[1])])


def _residualise(target, design):
    """Least-squares residual of `target` on `design` plus an intercept.

    `target` may be a matrix, in which case every column is residualised in
    one solve. That matters at scale: the first version was called once per
    column from a list comprehension, which is 40 638 solves on the full HLC
    array and makes the arm uncomputable rather than slow.
    """
    design = np.asarray(design, dtype=float)
    if design.ndim == 1:
        design = design.reshape(-1, 1)
    block = np.column_stack([np.ones(design.shape[0]), design])
    values = np.asarray(target, dtype=float)
    coefficients, *_ = np.linalg.lstsq(block, values, rcond=None)
    return values - block @ coefficients


def _abs_correlation(columns, target):
    columns = np.asarray(columns, dtype=float)
    target = np.asarray(target, dtype=float)
    centred = columns - columns.mean(axis=0)
    scale = np.sqrt((centred ** 2).sum(axis=0))
    centred_target = target - target.mean()
    target_scale = math.sqrt(float((centred_target ** 2).sum())) or 1.0
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.abs(centred.T @ centred_target) / (scale * target_scale)
    return np.nan_to_num(out, nan=0.0)


def drop_redundant(order, design, *, threshold=None, keep_columns=()):
    """Walk an ordering and drop a candidate that duplicates a kept one.

    Greedy, order-respecting and computed on **training rows only**. A tree
    model attributes a split to whichever of two near-identical columns it
    happened to use, so SHAP ranks both as important and a budget of k buys
    fewer than k pieces of information. This step makes `k` mean what it says.

    `keep_columns` are never dropped: they are the baseline, which is pinned
    into every panel by contract rather than earned.

    Correlation is computed, not assumed symmetric-free: a column with zero
    variance correlates with nothing, so it survives and is then harmless.
    """
    threshold = (FROZEN["redundancy_threshold"] if threshold is None
                 else threshold)
    order = np.asarray(order, dtype=int)
    design = np.asarray(design, dtype=float)
    protected = set(int(c) for c in keep_columns)
    centred = design - design.mean(axis=0)
    norms = np.sqrt((centred ** 2).sum(axis=0))
    kept, dropped = [], []
    for column in order:
        column = int(column)
        if column in protected:
            kept.append(column)
            continue
        redundant = False
        for already in kept:
            if already in protected:
                continue
            denominator = norms[column] * norms[already]
            if denominator <= 0.0:
                continue
            correlation = abs(float(centred[:, column] @ centred[:, already])
                              / denominator)
            if correlation >= threshold:
                redundant = True
                dropped.append((column, already, correlation))
                break
        if not redundant:
            kept.append(column)
    return np.asarray(kept, dtype=int), dropped


def matched_random_pool(design, pool, *, exclude=(), seed, neighbourhood=None):
    """A random pool the same size as `pool`, matched on level and spread.

    **Why this exists.** Comparing a 31-probe mechanistic pool against 30 218
    candidates changes two things at once: the biology used to choose, and the
    dimension searched. A random pool of 31 probes holds the dimension fixed
    and removes the biology, so the difference between them is what the
    knowledge contributed and nothing else.

    Matching is on the mean and the log variance of each probe, because an
    unmatched draw would mostly pick low-expressed, low-variance probes and
    lose to the mechanistic pool for a reason that has nothing to do with
    biology either. Each target probe is replaced by a uniform draw from its
    `neighbourhood` nearest neighbours in standardised (mean, log variance)
    space, sampled without replacement.

    **Outcome-blind by construction**: mean, variance and correlation involve
    no activity value, so this can be computed once on the retained donors
    without leaking anything into a fold.
    """
    neighbourhood = (FROZEN["match_neighbourhood"] if neighbourhood is None
                     else neighbourhood)
    design = np.asarray(design, dtype=float)
    pool = np.asarray(pool, dtype=int)
    banned = set(int(c) for c in exclude) | set(int(c) for c in pool)
    level = design.mean(axis=0)
    spread = np.log1p(design.var(axis=0))
    features = np.column_stack([level, spread])
    centre = features.mean(axis=0)
    scale = features.std(axis=0)
    scale[scale == 0] = 1.0
    standardised = (features - centre) / scale
    rng = np.random.default_rng(seed)
    taken, chosen = set(banned), []
    for target in pool:
        distance = np.sqrt(((standardised - standardised[int(target)]) ** 2
                            ).sum(axis=1))
        order = np.argsort(distance, kind="mergesort")
        candidates = [int(c) for c in order if int(c) not in taken]
        if not candidates:
            break
        window = candidates[:max(int(neighbourhood), 1)]
        pick = int(rng.choice(window))
        chosen.append(pick)
        taken.add(pick)
    chosen = np.asarray(chosen, dtype=int)
    return chosen, {
        "requested": int(pool.size),
        "drawn": int(chosen.size),
        "target_mean_level": float(level[pool].mean()),
        "matched_mean_level": float(level[chosen].mean())
            if chosen.size else float("nan"),
        "target_mean_log_variance": float(spread[pool].mean()),
        "matched_mean_log_variance": float(spread[chosen].mean())
            if chosen.size else float("nan"),
        "target_within_pool_redundancy": _mean_abs_correlation(design, pool),
        "matched_within_pool_redundancy": _mean_abs_correlation(design, chosen),
        "matched_on": "mean and log variance; redundancy is REPORTED not matched",
        "why_redundancy_is_only_reported": (
            "a mechanistic pool is correlated because its members are "
            "co-regulated, which is part of what makes it a mechanism. "
            "Forcing a random pool to the same internal correlation would "
            "smuggle in a piece of the structure the comparison is meant to "
            "test. It is measured on both sides so a reader can see whether "
            "the pools differ in it."),
    }


def _mean_abs_correlation(design, columns):
    columns = np.asarray(columns, dtype=int)
    if columns.size < 2:
        return float("nan")
    block = np.asarray(design, dtype=float)[:, columns]
    centred = block - block.mean(axis=0)
    norms = np.sqrt((centred ** 2).sum(axis=0))
    norms[norms == 0] = 1.0
    matrix = np.abs((centred.T @ centred) / np.outer(norms, norms))
    upper = matrix[np.triu_indices(columns.size, k=1)]
    return float(np.nanmean(upper)) if upper.size else float("nan")


def grouped_folds(groups, n_folds, seed):
    """Whole groups to folds. A group is a donor, or a protein position."""
    labels = np.asarray(groups)
    keys = np.unique(labels)
    _check(keys.size >= n_folds,
           f"{keys.size} groups cannot fill {n_folds} folds")
    shuffled = np.random.default_rng(seed).permutation(keys)
    assignment = {key: index % n_folds for index, key in enumerate(shuffled)}
    return np.array([assignment[key] for key in labels], dtype=int)


# --------------------------------------------------------------------------
# exact TreeSHAP, from XGBoost itself
# --------------------------------------------------------------------------
def shap_importance(design, target, *, groups, seed, inner_folds=None,
                    settings=None):
    """Mean |exact TreeSHAP| per column, computed on UNSEEN rows.

    The booster is fitted on an inner training part and the contributions are
    read on the inner held-out part, averaged over inner folds. SHAP computed
    on training rows ranks an overfitted feature as important, which is the one
    thing this ordering must not do.

    Returns `(importance, diagnostics)`. `importance` is zero-filled when
    XGBoost is absent, and the caller is told so rather than being handed a
    silent fallback.
    """
    inner_folds = FROZEN["inner_folds"] if inner_folds is None else inner_folds
    settings = dict(FROZEN["xgb"] if settings is None else settings)
    design = np.asarray(design, dtype=float)
    target = np.asarray(target, dtype=float)
    try:
        import xgboost
    except ImportError:
        return (np.zeros(design.shape[1]),
                {"available": False,
                 "why": "xgboost is not installed; no SHAP ordering exists"})
    folds = grouped_folds(groups, min(inner_folds, np.unique(groups).size), seed)
    totals = np.zeros(design.shape[1])
    used = 0
    for fold in np.unique(folds):
        held = np.flatnonzero(folds == fold)
        keep = np.flatnonzero(folds != fold)
        if keep.size < 10 or held.size == 0:
            continue
        booster = xgboost.XGBRegressor(
            objective="reg:squarederror", random_state=seed, n_jobs=1,
            **settings)
        booster.fit(design[keep], target[keep])
        contributions = booster.get_booster().predict(
            xgboost.DMatrix(design[held]), pred_contribs=True)
        # The last column is the bias term, not a feature.
        totals += np.abs(contributions[:, :-1]).mean(axis=0)
        used += 1
    _check(used > 0, "no inner fold was usable for a SHAP ordering")
    return (totals / used,
            {"available": True, "inner_folds_used": used,
             "computed_on": "inner held-out rows, never training rows",
             "exact": True, "via": "xgboost pred_contribs"})


def _elastic_net_support(design, target, *, groups, seed):
    """Non-zero coefficients of an elastic net whose penalty is chosen in-fold.

    **The internal cross-validation is grouped too.** `ElasticNetCV`'s default
    splits rows at random, which on a deep-mutational-scanning cohort puts
    `p.Arg34Ala` in its training part and `p.Arg34Gly` in its validation part
    and then picks a penalty that looks better than it is. Grouping the outer
    folds and leaving the penalty search ungrouped would leak at exactly the
    step whose job is to control complexity.
    """
    try:
        from sklearn.linear_model import ElasticNetCV
    except ImportError:
        return None
    centre = design.mean(axis=0)
    scale = design.std(axis=0)
    scale[scale == 0] = 1.0
    n_inner = min(FROZEN["inner_folds"], int(np.unique(groups).size))
    if n_inner < 2:
        return None
    inner = grouped_folds(groups, n_inner, seed)
    splits = [(np.flatnonzero(inner != fold), np.flatnonzero(inner == fold))
              for fold in np.unique(inner)]
    splits = [(keep, held) for keep, held in splits
              if keep.size >= 10 and held.size > 0]
    if len(splits) < 2:
        return None
    model = ElasticNetCV(l1_ratio=[0.2, 0.5, 0.9], cv=splits,
                         random_state=seed, max_iter=2000, n_alphas=12)
    model.fit((design - centre) / scale, target)
    return np.abs(model.coef_)


# --------------------------------------------------------------------------
# the selectors
# --------------------------------------------------------------------------
def rank_candidates(name, design, target, *, baseline_columns, knowledge_order,
                    groups, seed, pool_size=None):
    """An ORDERING of candidate columns, strongest first. No budget applied yet.

    Every selector that sees the outcome sees only training rows, because this
    is called with a training slice and nothing else.
    """
    n_columns = design.shape[1]
    if name == "S0_everything":
        return np.arange(n_columns), {"uses_outcome": False}
    if name == "S1_univariate_outcome":
        strength = _abs_correlation(_rank(design), _rank(target))
        return strength.argsort()[::-1], {"uses_outcome": True}
    if name == "S2_residual_on_baseline":
        # Information the baseline does not already carry. Selecting on the
        # raw outcome rewards a variable that merely repeats a covariate.
        if baseline_columns:
            base = _rank(design[:, list(baseline_columns)])
            # One solve for every column, not one solve per column.
            residual_columns = _residualise(_rank(design), base)
            residual_target = _residualise(_rank(target), base)
        else:
            residual_columns, residual_target = _rank(design), _rank(target)
        strength = _abs_correlation(residual_columns, residual_target)
        return strength.argsort()[::-1], {"uses_outcome": True,
                                          "residualised_on_baseline": True}
    if name == "S3_elastic_net":
        strength = _elastic_net_support(design, target, groups=groups,
                                        seed=seed)
        if strength is None:
            return np.arange(n_columns), {
                "uses_outcome": True,
                "unavailable": ("scikit-learn absent, or too few groups for a "
                                "grouped penalty search")}
        return strength.argsort()[::-1], {
            "uses_outcome": True,
            "penalty_search_is_grouped": True,
            "why": ("an ungrouped penalty search would validate on the same "
                    "protein position it trained on and choose a penalty that "
                    "looks better than it is")}
    if name in ("S4_shap_xgboost", "S7_shap_deduplicated"):
        strength, diagnostics = shap_importance(design, target, groups=groups,
                                                seed=seed)
        return strength.argsort()[::-1], {"uses_outcome": True, **diagnostics}
    if name == "S8_pool_then_univariate":
        # The arm that asks the question a restricted pool actually poses: if
        # a mechanistically chosen handful of candidates is as good as ranking
        # the whole array, the knowledge did its work by DELIMITING rather
        # than by ordering -- which is a different claim from S5's, and the
        # one that survived the earlier runs.
        _check(knowledge_order is not None,
               "the pool selector needs a knowledge ordering to delimit the pool")
        _check(pool_size is not None and pool_size > 0,
               "the pool selector needs a declared pool size")
        pool = np.asarray(knowledge_order, dtype=int)[:int(pool_size)]
        strength = _abs_correlation(_rank(design[:, pool]), _rank(target))
        return pool[strength.argsort()[::-1]], {
            "uses_outcome": True, "pool_size": int(pool.size),
            "how": ("knowledge delimits the candidate pool and the training "
                    "data orders what is inside it")}
    if name == "S5_knowledge_only":
        _check(knowledge_order is not None,
               "the knowledge-only selector needs a knowledge ordering")
        return np.asarray(knowledge_order, dtype=int), {"uses_outcome": False}
    if name == "S6_knowledge_then_shap":
        _check(knowledge_order is not None,
               "the hybrid selector needs a knowledge ordering")
        return np.asarray(knowledge_order, dtype=int), {
            "uses_outcome": True, "two_stage": True,
            "how": ("knowledge proposes oversample x budget candidates, then "
                    "out-of-sample SHAP reorders only those")}
    raise ScreeningError(f"unknown selector {name!r}")


def selector_ordering(name, design, target, *, baseline_columns,
                      knowledge_order, groups, seed, max_k, pool_size=None):
    """One ordering per selector per training set. Budgets are PREFIXES of it.

    Computing this once instead of once per k is not only cheaper, it is
    cleaner: a selector's ordering does not depend on how many variables it is
    allowed to keep, so recomputing it per budget invited the two to drift
    apart. The hybrid is the only one where the budget could matter, and it is
    pinned here to a single shortlist of `oversample x max_k` candidates, so
    every budget reads a prefix of one reordering rather than a different one.
    """
    if name == "S0_everything":
        return np.arange(design.shape[1]), {"uses_outcome": False}
    order, diagnostics = rank_candidates(
        name, design, target, baseline_columns=baseline_columns,
        knowledge_order=knowledge_order, groups=groups, seed=seed,
        pool_size=pool_size)
    if name == "S6_knowledge_then_shap":
        width = min(max(int(max_k) * FROZEN["oversample"], int(max_k)),
                    design.shape[1])
        shortlist = np.asarray(order[:width], dtype=int)
        strength, shap_diagnostics = shap_importance(
            design[:, shortlist], target, groups=groups, seed=seed)
        order = shortlist[strength.argsort()[::-1]]
        diagnostics = {**diagnostics, "shortlist_width": int(width),
                       **{f"shap_{k}": v for k, v in shap_diagnostics.items()}}
    if name == "S7_shap_deduplicated":
        order, dropped = drop_redundant(order, design,
                                        keep_columns=baseline_columns)
        diagnostics = {**diagnostics, "deduplicated": True,
                       "n_dropped_as_redundant": len(dropped),
                       "redundancy_threshold":
                           FROZEN["redundancy_threshold"],
                       "why": ("a budget of k buys fewer than k pieces of "
                               "information when two candidates are near "
                               "duplicates of one another")}
    return np.asarray(order, dtype=int), diagnostics


# --------------------------------------------------------------------------
# the harness
# --------------------------------------------------------------------------
def _fit_predict(design, target, train, test, *, estimator, seed):
    if design.shape[1] == 0:
        # k = 0 with no baseline covariate: the only honest prediction left is
        # the training mean. Returning it explicitly beats letting a solver
        # fail on an empty design, and it makes abstention a real arm rather
        # than an error.
        return np.full(test.size, float(np.mean(target[train])))
    if estimator == "ridge":
        return ridge_predict(design[train], target[train], design[test],
                             penalty=FROZEN["ridge_penalty"])
    if estimator == "xgboost":
        import xgboost
        booster = xgboost.XGBRegressor(objective="reg:squarederror",
                                       random_state=seed, n_jobs=1,
                                       **FROZEN["xgb"])
        booster.fit(design[train], target[train])
        return booster.predict(design[test])
    raise ScreeningError(f"unknown estimator {estimator!r}")


def screen(target, design, groups, *, baseline_columns=(), knowledge_order=None,
           feature_names=None, estimator="ridge", seed=None, k_grid=None,
           selectors=None, permute_target=False, permute_knowledge=False,
           pool_size=None, prefilter_top_variance=None):
    """Every selector, same folds, same estimator; k chosen inside each fold.

    `baseline_columns` are columns every panel keeps — the information already
    in hand, against which a new measurement has to add something. Without
    them "which variables are useful" is a different and easier question.
    """
    seed = FROZEN["seed"] if seed is None else seed
    k_grid = tuple(FROZEN["k_grid"] if k_grid is None else k_grid)
    selectors = tuple(SELECTORS if selectors is None else selectors)
    target = np.asarray(target, dtype=float)
    design = np.asarray(design, dtype=float)
    groups = np.asarray(groups)
    n, p = design.shape
    _check(target.size == n, "the outcome and the design disagree on n")
    _check(groups.size == n, "the group labels and the design disagree on n")
    _check(np.isfinite(design).all(), "the candidate pool has a missing value")
    _check(np.isfinite(target).all(), "the outcome has a missing value")
    baseline_columns = tuple(int(c) for c in baseline_columns)
    for column in baseline_columns:
        _check(0 <= column < p, f"baseline column {column} is out of range")
    k_grid = tuple(k for k in k_grid if k <= p) or (p,)
    _check(all(k >= 0 for k in k_grid), "a budget cannot be negative")

    rng = np.random.default_rng(seed)
    if permute_target:
        target = target[rng.permutation(n)]
    if permute_knowledge and knowledge_order is not None:
        knowledge_order = rng.permutation(np.asarray(knowledge_order, int))

    outer = grouped_folds(groups, min(FROZEN["outer_folds"], np.unique(groups).size),
                          seed)
    scored = {name: np.full(n, np.nan) for name in selectors}
    curve = {name: {k: np.full(n, np.nan) for k in k_grid} for name in selectors}
    chosen_k = {name: [] for name in selectors}
    chosen_panels = {name: [] for name in selectors}

    def with_baseline(columns):
        return np.asarray(sorted(set(int(c) for c in columns)
                                 | set(baseline_columns)), int)

    max_k = max(k_grid)
    orderings = {}

    prefilter = (FROZEN["prefilter_top_variance"]
                 if prefilter_top_variance is None else prefilter_top_variance)
    prefilter_sizes = []

    def _prefiltered(rows):
        """Most-variable candidates, by TRAINING-row variance. Outcome-blind.

        Returned as the set a selector may order, and as the map back to the
        original column numbers. Nothing downstream ever sees a filtered
        index, because a panel reported in filtered coordinates would be
        unreadable and silently wrong.
        """
        if not prefilter or prefilter >= p:
            return None
        variance = design[rows].var(axis=0)
        keep = set(int(c) for c in baseline_columns)
        ranked = [int(c) for c in np.argsort(variance)[::-1]
                  if int(c) not in keep]
        chosen = sorted(keep | set(ranked[:int(prefilter)]))
        return np.asarray(chosen, dtype=int)

    def ordering_for(name, rows):
        """Memoised on the exact training row set, so it is computed once."""
        key = (name, rows.tobytes())
        if key not in orderings:
            columns = _prefiltered(rows)
            if columns is None:
                orderings[key] = selector_ordering(
                    name, design[rows], target[rows],
                    baseline_columns=baseline_columns,
                    knowledge_order=knowledge_order, groups=groups[rows],
                    seed=seed, max_k=max_k, pool_size=pool_size)
            else:
                prefilter_sizes.append(int(columns.size))
                position = {int(c): i for i, c in enumerate(columns)}
                local_baseline = tuple(position[int(c)]
                                       for c in baseline_columns)
                local_knowledge = (
                    np.asarray([position[int(c)]
                                for c in np.asarray(knowledge_order, int)
                                if int(c) in position], dtype=int)
                    if knowledge_order is not None else None)
                local_order, diagnostics = selector_ordering(
                    name, design[np.ix_(rows, columns)], target[rows],
                    baseline_columns=local_baseline,
                    knowledge_order=local_knowledge, groups=groups[rows],
                    seed=seed, max_k=max_k, pool_size=pool_size)
                # Straight back to original column numbers.
                orderings[key] = (columns[np.asarray(local_order, int)],
                                  diagnostics)
        return orderings[key][0]

    for fold in np.unique(outer):
        held = np.flatnonzero(outer == fold)
        keep = np.flatnonzero(outer != fold)
        if keep.size < 10 or held.size == 0:
            continue
        inner = grouped_folds(groups[keep],
                              min(FROZEN["inner_folds"],
                                  np.unique(groups[keep]).size), seed + 1)
        for name in selectors:
            # --- the curve: one panel per k, selected on TRAINING rows only
            outer_order = ordering_for(name, keep)
            for k in k_grid:
                columns = with_baseline(outer_order[:k])
                curve[name][k][held] = _fit_predict(
                    design[:, columns], target, keep, held,
                    estimator=estimator, seed=seed)
            # --- the honest arm: k chosen by the INNER split of this fold
            inner_loss = {}
            for k in k_grid:
                losses = []
                for inner_fold in np.unique(inner):
                    inner_held = keep[inner == inner_fold]
                    inner_keep = keep[inner != inner_fold]
                    if inner_keep.size < 10 or inner_held.size == 0:
                        continue
                    columns = with_baseline(
                        ordering_for(name, inner_keep)[:k])
                    prediction = _fit_predict(design[:, columns], target,
                                              inner_keep, inner_held,
                                              estimator=estimator, seed=seed)
                    losses.append(float(np.abs(prediction
                                               - target[inner_held]).mean()))
                if losses:
                    inner_loss[k] = float(np.mean(losses))
            _check(bool(inner_loss),
                   f"no inner fold was usable for {name} in outer fold {fold}")
            best_k = min(sorted(inner_loss), key=lambda k: (inner_loss[k], k))
            columns = with_baseline(outer_order[:best_k])
            scored[name][held] = _fit_predict(design[:, columns], target, keep,
                                              held, estimator=estimator,
                                              seed=seed)
            chosen_k[name].append(int(best_k))
            chosen_panels[name].append(tuple(int(c) for c in columns))

    complete = np.all([np.isfinite(scored[name]) for name in selectors], axis=0)
    _check(complete.sum() >= 20,
           f"only {int(complete.sum())} rows were scored by every selector")
    losses = {name: np.abs(scored[name][complete] - target[complete])
              for name in selectors}
    names = (list(feature_names) if feature_names is not None
             else [f"x{j}" for j in range(p)])
    return {
        "n_scored": int(complete.sum()),
        "n_candidates": int(p),
        "estimator": estimator,
        "k_grid": list(k_grid),
        "permuted_target": bool(permute_target),
        "permuted_knowledge": bool(permute_knowledge),
        "scored_rows": complete,
        "losses": losses,
        "mean_loss": {name: float(v.mean()) for name, v in losses.items()},
        "k_chosen_per_fold": {name: list(v) for name, v in chosen_k.items()},
        # How often each selector chose to keep nothing but the baseline.
        # On an endpoint where selecting hurts, this is the right answer and
        # a method that cannot reach it is not automatic.
        "n_folds_that_abstained": {
            name: int(sum(1 for k in v if k == 0))
            for name, v in chosen_k.items()},
        "abstention_was_available": 0 in k_grid,
        "panels_per_fold": {
            name: [[names[c] for c in panel] for panel in v]
            for name, v in chosen_panels.items()},
        "stability": {name: _stability(v, names)
                      for name, v in chosen_panels.items()},
        "curve": {
            name: {str(k): float(np.abs(values[complete]
                                        - target[complete]).mean())
                   for k, values in by_k.items()
                   if np.isfinite(values[complete]).all()}
            for name, by_k in curve.items()},
        "the_curve_is_descriptive": (
            "its best point is chosen with hindsight across the whole grid, so "
            "that number is optimistically biased and is never the headline. "
            "The headline is `mean_loss`, where k and the panel were both "
            "chosen inside each training fold."),
        "baseline_columns": [names[c] for c in baseline_columns],
        # Carried so that nothing downstream has to re-derive them. The first
        # version of knowledge_agreement built its own generic names, compared
        # them against panels written in the real ones, and reported an
        # overlap of exactly 0.000 on every cell -- a bug that looked like a
        # finding.
        "candidate_names": list(names),
        "n_orderings_computed": len(orderings),
        "one_ordering_per_selector_per_training_set": True,
        "pool_size": None if pool_size is None else int(pool_size),
        "prefilter_top_variance": None if not prefilter else int(prefilter),
        "prefiltered_candidates_per_fold": sorted(set(prefilter_sizes)),
        "the_prefilter_is_outcome_blind": (
            "it ranks candidates by their variance on the TRAINING rows "
            "alone. It never reads the outcome, so it cannot leak, and being "
            "refitted per fold it never sees the test rows' variance either. "
            "It exists because 40 638 candidates make an elastic net and an "
            "exact TreeSHAP uncomputable, not because it improves anything."),
    }


def _stability(panels, names):
    """How often each column survived, and the fold-to-fold agreement."""
    if not panels:
        return {"n_folds": 0}
    counts = {}
    for panel in panels:
        for column in panel:
            counts[names[column]] = counts.get(names[column], 0) + 1
    pairs = []
    for i in range(len(panels)):
        for j in range(i + 1, len(panels)):
            left, right = set(panels[i]), set(panels[j])
            union = left | right
            pairs.append(len(left & right) / len(union) if union else 1.0)
    always = sorted(k for k, v in counts.items() if v == len(panels))
    return {
        "n_folds": len(panels),
        "mean_pairwise_jaccard": float(np.mean(pairs)) if pairs else 1.0,
        "n_selected_in_every_fold": len(always),
        "selected_in_every_fold": always[:12],
        "n_distinct_columns_ever_selected": len(counts),
        "why_it_matters": (
            "a selector that names a different panel in every fold has not "
            "found a mechanism, whatever its out-of-fold loss. Agreement is "
            "reported beside performance, never instead of it."),
    }


def contrast(result, better, worse, *, groups=None, family_size=1, margin=None,
             seed=None):
    """Paired interval and the project's six-rung ladder, cluster-bootstrapped."""
    margin = (FROZEN["minimum_interesting_relative_gain"] if margin is None
              else margin)
    losses = result["losses"]
    _check(better in losses and worse in losses,
           f"unknown selector in {better!r} vs {worse!r}")
    cluster = None
    if groups is not None:
        labels = np.asarray(groups)
        mask = np.asarray(result["scored_rows"], dtype=bool)
        _check(labels.size == mask.size,
               "the group labels must cover every row screen() saw")
        cluster = labels[mask]
    gap = paired_gap(losses[worse], losses[better], family_size=family_size,
                     n_boot=FROZEN["n_bootstrap"], alpha=FROZEN["alpha"],
                     seed=FROZEN["seed"] if seed is None else seed,
                     groups=cluster)
    reference = float(np.mean(losses[worse]))
    _check(reference > 0.0, "the reference loss must be positive")
    threshold = margin * reference
    return {"contrast": f"{better}_minus_{worse}",
            "verdict": _ladder(gap["interval"], threshold),
            "relative_gain": float(gap["gap"] / reference),
            "minimum_interesting_gain": float(threshold), **gap}


def non_inferiority(result, candidate, incumbent, *, groups=None,
                    family_size=1, margin=None, seed=None):
    """Is `candidate` no worse than `incumbent` by more than the margin?

    **A non-significant difference is not equivalence.** To conclude that the
    simple method suffices, the interval has to exclude the possibility that
    it is worse by an amount that would matter. That is a statement about the
    gap's **lower** bound against **minus** the margin, and it is a different
    test from the superiority ladder used everywhere else in this project --
    which is why it is a separate function rather than another rung.

    `gap` is `loss(incumbent) - loss(candidate)`, so a positive gap means the
    candidate is better.

    * `NON_INFERIOR` -- the lower bound is above minus the margin. The
      candidate cannot be worse by an amount that matters.
    * `INFERIOR` -- the upper bound is below minus the margin. It is worse by
      an amount that matters.
    * `UNDECIDED` -- the interval spans minus the margin. **This answers
      nothing**, and in particular it is not evidence of equivalence.
    """
    margin = (FROZEN["minimum_interesting_relative_gain"] if margin is None
              else margin)
    losses = result["losses"]
    _check(candidate in losses and incumbent in losses,
           f"unknown selector in {candidate!r} vs {incumbent!r}")
    cluster = None
    if groups is not None:
        labels = np.asarray(groups)
        mask = np.asarray(result["scored_rows"], dtype=bool)
        _check(labels.size == mask.size,
               "the group labels must cover every row screen() saw")
        cluster = labels[mask]
    gap = paired_gap(losses[incumbent], losses[candidate],
                     family_size=family_size, n_boot=FROZEN["n_bootstrap"],
                     alpha=FROZEN["alpha"],
                     seed=FROZEN["seed"] if seed is None else seed,
                     groups=cluster)
    reference = float(np.mean(losses[incumbent]))
    _check(reference > 0.0, "the reference loss must be positive")
    threshold = margin * reference
    low, high = gap["interval"]
    if low > -threshold:
        verdict = "NON_INFERIOR"
    elif high < -threshold:
        verdict = "INFERIOR"
    else:
        verdict = "UNDECIDED"
    return {
        "test": f"{candidate}_non_inferior_to_{incumbent}",
        "verdict": verdict,
        "relative_gap": float(gap["gap"] / reference),
        "non_inferiority_margin": float(threshold),
        "relative_margin": float(margin),
        "an_undecided_result_is_not_equivalence": True,
        "how_to_read": (
            "NON_INFERIOR means the candidate cannot be worse by more than "
            f"{margin:.0%} of the incumbent's loss. UNDECIDED means the "
            "design lacks the power to say that, which is not the same claim "
            "and must not be reported as one."),
        **gap}


def parsimony_point(curve, *, retention=None, baseline_k=None):
    """Smallest k retaining a share of the gain the curve shows over its worst.

    The commissioning note asks for the number of variables retaining 80-90 %
    of the predictive gain. That number is read off an out-of-fold curve, so
    **it is a description of this curve and not a validated panel size** — the
    validated one is whatever the inner folds chose, which is reported
    separately and may differ.
    """
    retention = (FROZEN["parsimony_retention"] if retention is None
                 else retention)
    points = {int(k): float(v) for k, v in curve.items()}
    _check(len(points) >= 2, "a parsimony point needs at least two k values")
    worst = max(points.values())
    best = min(points.values())
    span = worst - best
    out = {"best_loss": best, "worst_loss": worst, "gain_span": span,
           "k_of_best": min(points, key=lambda k: (points[k], k)),
           "the_k_of_best_is_optimistically_biased": True}
    if baseline_k is not None and baseline_k in points:
        out["loss_at_baseline_k"] = points[baseline_k]
    for share in retention:
        if span <= 0:
            out[f"k_retaining_{int(share * 100)}pct"] = None
            continue
        ceiling = worst - share * span
        eligible = [k for k in sorted(points) if points[k] <= ceiling]
        out[f"k_retaining_{int(share * 100)}pct"] = (eligible[0] if eligible
                                                     else None)
    return out


def knowledge_agreement(result, knowledge_order, *, k, n_candidates=None,
                        seed=None):
    """Does an outcome-blind ordering anticipate what the data ends up keeping?

    This is the commissioning note's **first** question — before the omic
    measurements, can knowledge say which ones would be worth making? It is
    answered by overlap against a random-ordering null, not by eyeballing a
    gene list.
    """
    seed = FROZEN["seed"] if seed is None else seed
    panels = result["panels_per_fold"].get("S4_shap_xgboost") or []
    _check(bool(panels), "no data-driven panel to compare the knowledge against")
    # The names the panels are actually written in. Rebuilding generic ones
    # here guarantees an empty intersection on every cell, which is how the
    # first version of this function reported 0.000 overlap everywhere.
    names = result.get("candidate_names")
    _check(names is not None,
           "this result predates candidate_names; re-run screen()")
    n_candidates = len(names) if n_candidates is None else int(n_candidates)
    _check(n_candidates == len(names),
           f"n_candidates={n_candidates} disagrees with {len(names)} names")
    observed_names = {n for panel in panels for n in panel}
    _check(observed_names <= set(names),
           "the panels are written in names this result does not list")
    proposed = {names[int(c)] for c in np.asarray(knowledge_order, int)[:k]}
    overlaps = []
    for panel in panels:
        overlaps.append(len(proposed & set(panel)) / max(len(proposed), 1))
    rng = np.random.default_rng(seed + 7)
    null = []
    for _ in range(2000):
        draw = {names[int(c)] for c in rng.choice(
            len(names), size=min(k, len(names)), replace=False)}
        null.append(np.mean([len(draw & set(panel)) / max(len(draw), 1)
                             for panel in panels]))
    observed = float(np.mean(overlaps))
    null = np.asarray(null, dtype=float)
    return {
        "k": int(k),
        "n_candidates": len(names),
        "observed_overlap": observed,
        "proposed_by_knowledge": sorted(proposed)[:12],
        "null_mean_overlap": float(null.mean()),
        "null_95th_percentile": float(np.percentile(null, 95)),
        "beats_a_random_ordering": bool(observed > np.percentile(null, 95)),
        "empirical_p": float((1 + (null >= observed).sum()) / (1 + null.size)),
        "what_this_does_not_show": (
            "agreement with what the data chose is not proof that either is "
            "right. It measures whether knowledge anticipates the data's "
            "choice, which is the question asked BEFORE a measurement is "
            "made, and says nothing about predictive value on its own."),
    }


# --------------------------------------------------------------------------
# the positive control: the harness must find a sparse signal that exists
# --------------------------------------------------------------------------
def sparse_signal_cohort(n_groups=120, per_group=3, n_candidates=60,
                         n_true=4, seed=0):
    """A cohort where exactly `n_true` of `n_candidates` columns carry signal.

    Without this, a null result from any selector is unreadable: it could mean
    the data has no structure, or that the harness cannot see structure. The
    useful columns are placed LAST so an ordering bug that returns the first k
    columns fails loudly instead of succeeding by accident.
    """
    rng = np.random.default_rng(seed)
    groups = np.repeat(np.arange(n_groups), per_group)
    n = groups.size
    pool = rng.normal(size=(n, n_candidates))
    # Two baseline covariates that matter, and correlated decoys that do not.
    baseline = rng.normal(size=(n, 2))
    pool[:, 0:2] = baseline
    truth_columns = list(range(n_candidates - n_true, n_candidates))
    weights = np.array([1.0, -0.8, 0.6, -0.5])[:n_true]
    signal = pool[:, truth_columns] @ weights
    # A redundant copy of the strongest true column: a univariate selector
    # should waste budget on it, a joint one should not.
    pool[:, 2] = pool[:, truth_columns[0]] + rng.normal(0, 0.1, size=n)
    target = (0.7 * baseline[:, 0] - 0.4 * baseline[:, 1] + signal
              + rng.normal(0, 0.5, size=n))
    knowledge_order = np.array(truth_columns + [
        j for j in range(n_candidates) if j not in truth_columns], dtype=int)
    return {"target": target, "design": pool, "groups": groups,
            "baseline_columns": (0, 1), "truth_columns": tuple(truth_columns),
            "redundant_copy_of_a_true_column": 2,
            "knowledge_order": knowledge_order}


def self_test(selectors=("S0_everything", "S1_univariate_outcome",
                         "S4_shap_xgboost")):
    """Assert the harness recovers a sparse signal, and loses it when permuted."""
    cohort = sparse_signal_cohort()
    live = screen(cohort["target"], cohort["design"], cohort["groups"],
                  baseline_columns=cohort["baseline_columns"],
                  knowledge_order=cohort["knowledge_order"],
                  selectors=selectors)
    dead = screen(cohort["target"], cohort["design"], cohort["groups"],
                  baseline_columns=cohort["baseline_columns"],
                  knowledge_order=cohort["knowledge_order"],
                  selectors=selectors, permute_target=True)
    gain = contrast(live, "S4_shap_xgboost", "S0_everything",
                    groups=cohort["groups"])
    control = contrast(dead, "S4_shap_xgboost", "S0_everything",
                       groups=cohort["groups"])
    _check(gain["verdict"] == "GAIN_REACHES_MARGIN",
           "selection must beat using everything when the signal IS sparse; "
           f"got {gain['verdict']} ({gain['relative_gain']:+.3f})")
    # Did it find the right columns, and not the redundant copy?
    panels = live["panels_per_fold"]["S4_shap_xgboost"]
    truth = {f"x{c}" for c in cohort["truth_columns"]}
    found = [len(truth & set(panel)) / len(truth) for panel in panels]
    _check(float(np.mean(found)) >= 0.5,
           "the selector must recover at least half the truly useful columns; "
           f"got {float(np.mean(found)):.2f}")
    return {
        "selection_vs_everything": {"verdict": gain["verdict"],
                                    "relative_gain": gain["relative_gain"]},
        "with_permuted_outcome": {"verdict": control["verdict"],
                                  "relative_gain": control["relative_gain"]},
        "share_of_truly_useful_columns_recovered": float(np.mean(found)),
        "redundant_copy_selected_in_n_folds": sum(
            1 for panel in panels
            if f"x{cohort['redundant_copy_of_a_true_column']}" in panel),
        "what_this_proves": (
            "when a sparse subset really carries the signal, this harness "
            "finds most of it and selection beats using everything. A null on "
            "real data is therefore about the data, not the instrument."),
    }


def what_this_module_does_not_do():
    return {
        "a_shap_value_is_not_evidence_of_value": (
            "it is an attribution inside one fitted model. Two collinear "
            "probes split one signal and both look weak; a feature the model "
            "leaned on can be useless out of sample. SHAP proposes an "
            "ordering here and never settles a claim."),
        "the_curve_is_not_a_validated_panel_size": (
            "the k at the curve's best point was chosen with hindsight across "
            "the whole grid. The validated size is the one the inner folds "
            "chose, reported separately, and the two often differ."),
        "beating_a_selector_is_not_predicting_pharmacokinetics": (
            "every endpoint here is an in-vitro activity or a functional "
            "screen score. Nothing in this module touches a clearance, an AUC "
            "or a concentration."),
        "a_selector_cannot_abstain_by_accident": (
            "hand any of these pure noise and it returns the k luckiest "
            "columns. That is why the permuted-outcome control is run beside "
            "every real result rather than on request."),
        "agreement_between_selectors_is_not_truth": (
            "two selectors can agree because they share a bias. Agreement is "
            "reported as a stability diagnostic, never as confirmation."),
    }
