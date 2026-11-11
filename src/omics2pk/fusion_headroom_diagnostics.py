"""How much headroom a fusion method could have, and whether any of it is reachable.

This module computes ceilings and diagnostics from predictions that **already
exist**. It implements no fusion method, fits no gate and introduces no
predictive architecture. Its only purpose is to decide whether building one
could pay, and to be able to say *no* on evidence.

The correction this module exists to carry
------------------------------------------
An earlier version of the refutation protocol used a single ceiling,

    L_const = min over a CONSTANT weight vector w in the simplex
              of the mean loss of sum_m w_m * yhat_im,

and proposed to stop the programme when `best single arm - L_const` was not
detectably positive. **That rule is wrong, and it is wrong in the direction
that would have eliminated the hypothesis it was meant to test.** A constant
weight vector cannot express a gate that chooses differently for different
individuals, so `L_const` does not bound contextual fusion.

The counterexample is two arms, each with mean absolute error 0.5, whose
predictions bracket the truth on opposite halves of the cohort: every constant
weight keeps the error at exactly 0.5, and a perfect contextual gate reaches 0.
``counterexample_constant_weights_cannot_see_a_gate`` builds it, and
``self_test`` asserts that this module's ceilings separate the two cases.

So four quantities are computed, not one, and they answer different questions:

=========================  ====================================================
`best_single_oracle`       the best arm, chosen knowing the answers
`constant_weight_oracle`   best fixed mixture - bounds **stacking only**
`individual_selection_oracle`  best arm **per individual** - bounds per-individual
                           *selection*
`individual_convex_oracle` best convex combination **per individual** - bounds
                           **every** method in the family, gates included
=========================  ====================================================

For any loss that is monotone in |prediction - truth|, the last one has a
closed form: the reachable predictions for individual *i* are exactly the
interval between that individual's smallest and largest arm prediction, so the
per-individual convex oracle error is the distance from the truth to that
interval, and is **zero whenever the arms bracket the truth**.

That closed form is also the reason a ceiling alone cannot guide anything. With
several arms the bracketing rate is often high and `individual_convex_oracle`
is often near zero, which makes the ceiling permissive almost always. A
permissive ceiling refutes nothing and licenses nothing. So the module also
computes what is **recoverable**:

* ``nested_constant_stacking`` fits the weights on one partition's
  out-of-fold predictions and scores them on another's, because fitting weights
  on the same out-of-fold predictions that then measure them is leakage even
  when the base predictions are themselves honestly out-of-fold;
* ``context_informativeness`` asks whether any available context variable
  predicts, out of sample, *which arm will win* for an individual. If it cannot,
  a contextual gate has nothing to condition on, and that refutes contextual
  fusion on this cohort **without** depending on the ceiling being tight.

And the verdict rule no longer treats an interval containing zero as a
refutation. Refuting requires an equivalence-style statement: the ceiling's
**upper** confidence bound must fall below a declared minimum interesting gain.
Anything else is a non-detection, which at these sample sizes usually means
low power.
"""
from __future__ import annotations

import math

import numpy as np

#: Frozen before execution.
FROZEN = {
    "seed": 20261111,
    "n_bootstrap": 2000,
    "alpha": 0.05,
    # The smallest relative reduction in the primary loss that would make a
    # fusion method worth building. 5 % matches the non-inferiority margin H2
    # already uses; 10 % matches H1's frozen incremental criterion and is
    # reported as a sensitivity.
    "minimum_interesting_relative_gain": 0.05,
    "sensitivity_relative_gain": 0.10,
    # The smallest out-of-sample rank association between a context variable
    # and the per-individual arm advantage that would make a gate worth
    # fitting. Below this a gate has nothing to condition on.
    "minimum_interesting_rank_association": 0.30,
    "frank_wolfe_iterations": 500,
    "context_probe_folds": 5,
    "loss": "absolute_error",
}


class FusionDiagnosticsError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise FusionDiagnosticsError(message)


def _as_matrix(predictions, truth):
    matrix = np.asarray(predictions, dtype=float)
    target = np.asarray(truth, dtype=float)
    _check(matrix.ndim == 2, "predictions must be (n_individuals, n_arms)")
    _check(matrix.shape[0] == target.size,
           f"{matrix.shape[0]} prediction rows against {target.size} outcomes")
    _check(matrix.shape[1] >= 2, "a ceiling needs at least two arms")
    _check(np.isfinite(matrix).all(), "a prediction is missing or non-finite")
    _check(np.isfinite(target).all(), "an outcome is missing or non-finite")
    return matrix, target


def losses(predictions, truth):
    """Per-individual absolute error, the module's declared loss."""
    matrix, target = _as_matrix(predictions, truth)
    return np.abs(matrix - target[:, None])


# --------------------------------------------------------------------------
# the four ceilings
# --------------------------------------------------------------------------
def best_single_oracle(predictions, truth):
    """The best arm's mean loss, the arm chosen knowing the answers.

    Optimistic by one choice among M. It is the reference the ceilings are
    measured against, and it is NOT what a method would achieve, because a
    method must choose its arm on training data. ``best_single_honest``
    measures that.
    """
    per_arm = losses(predictions, truth).mean(axis=0)
    index = int(np.argmin(per_arm))
    return {"mean_loss": float(per_arm[index]), "arm_index": index,
             "per_arm_mean_loss": [float(v) for v in per_arm]}


def constant_weight_oracle(predictions, truth, *, iterations=None):
    """Best FIXED mixture, by Frank-Wolfe on a convex objective.

    **This bounds stacking and uniform averaging. It does NOT bound a
    contextual gate**, which is the correction this module carries.

    The objective mean_i |sum_m w_m yhat_im - y_i| is convex in w, so a
    conditional-gradient method converges to the optimum over the simplex and
    needs no optimiser to debug: each step moves towards the vertex with the
    smallest directional derivative. The vertices are also evaluated directly,
    so the result can never be worse than the best single arm.
    """
    iterations = (FROZEN["frank_wolfe_iterations"] if iterations is None
                  else iterations)
    matrix, target = _as_matrix(predictions, truth)
    n, arms = matrix.shape

    def objective(weights):
        return float(np.abs(matrix @ weights - target).mean())

    # Start at the best vertex, so the answer is never worse than best_single.
    vertex_losses = np.abs(matrix - target[:, None]).mean(axis=0)
    weights = np.zeros(arms)
    weights[int(np.argmin(vertex_losses))] = 1.0
    best = objective(weights)

    for step in range(1, iterations + 1):
        residual = matrix @ weights - target
        # Subgradient of the mean absolute error with respect to w.
        gradient = (np.sign(residual) @ matrix) / n
        direction = int(np.argmin(gradient))
        candidate = np.zeros(arms)
        candidate[direction] = 1.0
        gamma = 2.0 / (step + 2.0)
        moved = (1.0 - gamma) * weights + gamma * candidate
        value = objective(moved)
        if value <= best:
            weights, best = moved, value

    return {"mean_loss": best,
            "weights": {f"arm_{index}": float(value)
                        for index, value in enumerate(weights)},
            "bounds": ("constant-weight mixtures only: stacking and uniform "
                       "averaging. It does NOT bound a contextual gate."),
            "iterations": int(iterations)}


def individual_selection_oracle(predictions, truth):
    """Best arm PER INDIVIDUAL. Bounds per-individual selection, not averaging.

    Averaging can beat it: arms erring +1 and -1 give a per-individual best of
    1 and an average of 0.
    """
    per_individual = losses(predictions, truth)
    chosen = per_individual.argmin(axis=1)
    shares = np.bincount(chosen, minlength=per_individual.shape[1]) / chosen.size
    return {"mean_loss": float(per_individual.min(axis=1).mean()),
            "best_arm_shares": {f"arm_{index}": float(value)
                                for index, value in enumerate(shares)},
            "bounds": "per-individual arm selection"}


def individual_convex_oracle(predictions, truth):
    """Best convex combination PER INDIVIDUAL — the ceiling that bounds gates.

    For a loss monotone in |prediction - truth|, individual *i*'s reachable
    predictions are exactly [min_m yhat_im, max_m yhat_im], so the achievable
    error is the distance from y_i to that interval. It is **zero whenever the
    arms bracket the truth**, which is why a high bracketing rate makes this
    ceiling permissive and why it can refute but never guide.
    """
    matrix, target = _as_matrix(predictions, truth)
    low = matrix.min(axis=1)
    high = matrix.max(axis=1)
    below = np.maximum(low - target, 0.0)
    above = np.maximum(target - high, 0.0)
    per_individual = below + above
    bracketed = (per_individual <= 0.0)
    return {"mean_loss": float(per_individual.mean()),
            "bracketing_rate": float(bracketed.mean()),
            "bounds": ("EVERY method in the family, contextual gates "
                       "included. A gap of zero here refutes the family."),
            "why_it_cannot_guide": (
                "it is zero wherever the arms bracket the truth, so with "
                "several arms it is usually small and usually permissive. "
                "Use it to refute, never to promise."),
            "per_individual_losses": per_individual}


# --------------------------------------------------------------------------
# what is recoverable, which is a different question from what is possible
# --------------------------------------------------------------------------
def _partitions(n, folds, rng):
    order = rng.permutation(n)
    return [order[index::folds] for index in range(folds)]


def best_single_honest(predictions, truth, *, folds=None, seed=None):
    """The best arm chosen on a training partition, scored on the held-out one.

    This is what a method actually achieves by selecting, and it is the only
    honest baseline for a fusion gain. It can be worse than
    ``best_single_oracle``, and the difference is the cost of not knowing which
    arm is best.
    """
    folds = FROZEN["context_probe_folds"] if folds is None else folds
    seed = FROZEN["seed"] if seed is None else seed
    matrix, target = _as_matrix(predictions, truth)
    per_individual = np.abs(matrix - target[:, None])
    rng = np.random.default_rng(seed)
    scored = np.zeros(matrix.shape[0])
    picked = []
    for held in _partitions(matrix.shape[0], folds, rng):
        keep = np.setdiff1d(np.arange(matrix.shape[0]), held)
        choice = int(per_individual[keep].mean(axis=0).argmin())
        picked.append(choice)
        scored[held] = per_individual[held, choice]
    return {"mean_loss": float(scored.mean()),
            "arm_chosen_per_partition": picked,
            "per_individual_losses": scored}


def nested_constant_stacking(predictions, truth, *, folds=None, seed=None):
    """Constant weights fitted on one partition, scored on another.

    Fitting weights on the same out-of-fold predictions that then measure them
    is leakage even when the base predictions are honestly out-of-fold: the
    weights see every individual's realised error. This is the check the
    commissioning note asked for.
    """
    folds = FROZEN["context_probe_folds"] if folds is None else folds
    seed = FROZEN["seed"] if seed is None else seed
    matrix, target = _as_matrix(predictions, truth)
    rng = np.random.default_rng(seed + 1)
    scored = np.zeros(matrix.shape[0])
    for held in _partitions(matrix.shape[0], folds, rng):
        keep = np.setdiff1d(np.arange(matrix.shape[0]), held)
        fitted = constant_weight_oracle(matrix[keep], target[keep])
        weights = np.array([fitted["weights"][f"arm_{index}"]
                            for index in range(matrix.shape[1])])
        scored[held] = np.abs(matrix[held] @ weights - target[held])
    return {"mean_loss": float(scored.mean()),
            "weights_refitted_per_partition": True,
            "per_individual_losses": scored,
            "what_this_is": ("the gain constant-weight stacking actually "
                             "recovers out of sample, as opposed to the "
                             "ceiling it could reach in hindsight")}


def ridge_predict(train_design, train_target, test_design, *, penalty=1.0):
    """Standardise on the training rows, solve, apply. No sklearn dependency.

    Standardisation uses the **training** mean and spread only, so no test row
    influences the fit through the scaling either.
    """
    train_design = np.asarray(train_design, dtype=float)
    test_design = np.asarray(test_design, dtype=float)
    train_target = np.asarray(train_target, dtype=float)
    mean = train_design.mean(axis=0)
    sd = train_design.std(axis=0, ddof=0)
    safe = np.where(sd > 0.0, sd, 1.0)
    train = (train_design - mean) / safe
    test = (test_design - mean) / safe
    centre = float(train_target.mean())
    gram = train.T @ train + penalty * np.eye(train.shape[1])
    coefficients = np.linalg.solve(gram, train.T @ (train_target - centre))
    return test @ coefficients + centre


def expand_basis(block, basis):
    """`linear` returns the columns; `quadratic` adds squares and products.

    The commissioning review's point: *the failure of a linear ridge probe on
    two arms does not refute contextual fusion, which may be non-linear.* A
    second-order expansion is the cheapest honest answer — it stays a linear
    solve, needs no new dependency, and is declared as its own member of the
    probe family rather than tried silently until something clears.
    """
    block = np.asarray(block, dtype=float)
    if block.ndim == 1:
        block = block.reshape(-1, 1)
    if basis == "linear":
        return block
    _check(basis == "quadratic", f"unknown basis {basis!r}")
    columns = [block]
    for i in range(block.shape[1]):
        for j in range(i, block.shape[1]):
            columns.append((block[:, i] * block[:, j]).reshape(-1, 1))
    return np.column_stack(columns)


def context_informativeness(predictions, truth, context, *, folds=None,
                            seed=None, family_size=1, basis="linear",
                            arm_pair=None):
    """Can any context variable predict WHICH ARM WINS, out of sample?

    This is the test that bears on contextual fusion directly, and it does not
    depend on the ceiling being tight. The target is the signed advantage of
    two arms for each individual; a ridge on the context is fitted on a
    training partition and its out-of-sample rank agreement with the realised
    advantage is reported.

    **If this has no out-of-sample skill, a contextual gate has nothing to
    condition on**, and contextual fusion is refuted on this cohort whatever
    the ceiling says.

    The leak this version closes
    ----------------------------
    Version 2.0.0 picked the two compared arms from the mean losses over
    **every** individual, before the train/test split, and only then
    cross-fitted the ridge. The ridge was honest; the choice of what it
    predicted was not. The bias has **no reliable sign** — on synthetic null
    cohorts it runs pessimistic, because the globally-best arm is selected to
    have the lowest total loss and the advantage it defines is then
    systematically negative — so it has to be measured, not assumed. Measured
    on the HLC predictions it ran optimistic. The pair
    is now chosen **inside each fold from the training rows only**, so the
    realised advantage is itself fold-dependent — which is exactly what a gate
    faces in deployment. ``arm_pair_stability`` reports how often the folds
    agreed, because an unstable pair is itself a reason a gate would struggle.
    """
    folds = FROZEN["context_probe_folds"] if folds is None else folds
    seed = FROZEN["seed"] if seed is None else seed
    matrix, target = _as_matrix(predictions, truth)
    block = np.asarray(context, dtype=float)
    if block.ndim == 1:
        block = block.reshape(-1, 1)
    _check(block.shape[0] == matrix.shape[0],
           "the context block and the predictions disagree on n")
    _check(np.isfinite(block).all(), "a context value is missing")

    per_individual = np.abs(matrix - target[:, None])
    design = expand_basis(block, basis)

    rng = np.random.default_rng(seed + 2)
    predicted = np.zeros(matrix.shape[0])
    # Fold-dependent, because the compared pair is now a training-set decision.
    advantage = np.zeros(matrix.shape[0])
    pairs = []
    for held in _partitions(matrix.shape[0], folds, rng):
        keep = np.setdiff1d(np.arange(matrix.shape[0]), held)
        if keep.size < 3:
            continue
        if arm_pair is None:
            # TRAINING ROWS ONLY. This line is the correction.
            order = per_individual[keep].mean(axis=0).argsort()
            first, second = int(order[0]), int(order[1])
        else:
            first, second = int(arm_pair[0]), int(arm_pair[1])
        pairs.append((first, second))
        # Positive advantage means the runner-up wins for this individual,
        # which is exactly what a gate would have to detect.
        train_advantage = per_individual[keep, first] - per_individual[keep, second]
        advantage[held] = per_individual[held, first] - per_individual[held, second]
        predicted[held] = ridge_predict(design[keep], train_advantage,
                                        design[held])

    observed_rho = _spearman(predicted, advantage)
    # The probe gets the same three-valued treatment as the ceilings: a null
    # probe on 15 or 47 individuals is as underpowered as a null gap, and
    # calling it a refutation would repeat the error this module corrects.
    rng_boot = np.random.default_rng(seed + 4)
    draws = []
    for _ in range(FROZEN["n_bootstrap"]):
        picked = rng_boot.integers(0, advantage.size, size=advantage.size)
        # Two distinct values are enough for a rank correlation; requiring
        # three would discard every draw whenever the advantage is binary,
        # which is exactly the counterexample's shape.
        if np.unique(advantage[picked]).size < 2:
            continue
        draws.append(_spearman(predicted[picked], advantage[picked]))
    half = 100.0 * FROZEN["alpha"] / 2.0
    adjusted_half = half / float(family_size)
    interval = ([float(np.percentile(draws, half)),
                 float(np.percentile(draws, 100.0 - half))]
                if len(draws) >= 100 else [float("nan"), float("nan")])
    adjusted = ([float(np.percentile(draws, adjusted_half)),
                 float(np.percentile(draws, 100.0 - adjusted_half))]
                if len(draws) >= 100 else [float("nan"), float("nan")])
    floor = FROZEN["minimum_interesting_rank_association"]
    name = probe_verdict(interval, adjusted, floor=floor)

    modal = max(set(pairs), key=pairs.count) if pairs else (-1, -1)
    return {
        "compared_arms": list(modal),
        "compared_arms_per_fold": [list(pair) for pair in pairs],
        "arm_pair_stability": (float(pairs.count(modal)) / len(pairs)
                               if pairs else float("nan")),
        "arm_pair_chosen_inside_each_fold": arm_pair is None,
        "basis": basis,
        "spearman_out_of_sample": observed_rho,
        "spearman_interval": interval,
        "spearman_interval_bonferroni": adjusted,
        "family_size": int(family_size),
        "verdict": name,
        "minimum_interesting_rank_association": float(floor),
        "reaches_minimum_interesting_association": bool(
            np.isfinite(adjusted[0]) and adjusted[0] >= floor),
        "sign_agreement": float(np.mean(
            np.sign(predicted) == np.sign(advantage))),
        "share_where_the_runner_up_wins": float(np.mean(advantage > 0.0)),
        "how_to_read": (
            "A rank correlation near zero and a sign agreement near the share "
            "of the majority class mean the context carries no information "
            "about which arm will win, so a gate would have nothing to learn. "
            "But that is only a REFUTATION when the interval's upper bound "
            "falls below the minimum interesting association; otherwise it is "
            "a non-detection, and at these cohort sizes it usually is. And an "
            "interval that excludes zero while sitting below the floor is a "
            "real but TOO SMALL association, which is not the same claim as "
            "reaching the declared margin."),
    }


def probe_verdict(interval, adjusted, *, floor=None):
    """Six outcomes, because *detectably positive* and *big enough* differ.

    Version 2.0.0 declared a minimum interesting rank association of 0.30 and
    then awarded `CONTEXT_CARRIES_SIGNAL` to any probe whose adjusted interval
    excluded zero. That let an association of 0.17 be reported under a name
    implying the declared margin had been met. The floor is now applied to
    **both** sides:

    * reaching it requires the adjusted **lower** bound at or above the floor;
    * an interval that excludes zero but whose **upper** bound is below the
      floor is a real association that is simultaneously **too small to be
      worth building on** — a distinct, informative outcome;
    * refuting still requires the upper bound below the floor with no
      detectable positive association.
    """
    floor = (FROZEN["minimum_interesting_rank_association"] if floor is None
             else floor)
    low_nominal, high_nominal = interval
    low_adjusted, high_adjusted = adjusted
    finite = np.isfinite(high_adjusted)
    if np.isfinite(low_adjusted) and low_adjusted >= floor:
        return "SIGNAL_REACHES_MARGIN"
    if np.isfinite(low_adjusted) and low_adjusted > 0.0:
        return ("SIGNAL_POSITIVE_BUT_BELOW_MARGIN" if finite
                and high_adjusted < floor
                else "SIGNAL_POSITIVE_MARGIN_UNDECIDED")
    if np.isfinite(low_nominal) and low_nominal > 0.0:
        return ("SIGNAL_POSITIVE_UNADJUSTED_ONLY_BELOW_MARGIN"
                if np.isfinite(high_nominal) and high_nominal < floor
                else "SIGNAL_POSITIVE_UNADJUSTED_ONLY_MARGIN_UNDECIDED")
    if np.isfinite(high_nominal) and high_nominal < floor:
        return "REFUTED_CONTEXT_UNINFORMATIVE"
    return "NO_DETECTABLE_SIGNAL_UNDERPOWERED"


def _spearman(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.size < 3:
        return float("nan")

    def rank(values):
        order = values.argsort()
        ranks = np.empty(values.size, dtype=float)
        ranks[order] = np.arange(1, values.size + 1, dtype=float)
        return ranks

    a, b = rank(left), rank(right)
    a = a - a.mean()
    b = b - b.mean()
    denominator = math.sqrt(float((a ** 2).sum() * (b ** 2).sum()))
    return float((a * b).sum() / denominator) if denominator > 0 else float("nan")


# --------------------------------------------------------------------------
# intervals, and a verdict rule that does not mistake non-detection for proof
# --------------------------------------------------------------------------
def paired_gap(reference_losses, ceiling_losses, *, n_boot=None, seed=None,
               alpha=None, family_size=1, groups=None):
    """Interval for the per-row gap `reference - ceiling`, resampled.

    The resampling unit must match the unit of independence, and those are not
    always the row. In the MaveDB benchmark a row is one *substitution* and up
    to nineteen of them share a protein *position*: they see the same residue,
    the same structural environment and largely the same assay artefacts. An
    i.i.d. row bootstrap then counts nineteen correlated observations as
    nineteen independent ones and returns an interval that is too narrow — in
    the wrong direction, because every verdict in this project is decided by
    whether the interval's **lower** bound clears a margin.

    Pass `groups` (one label per row) to resample **clusters** with
    replacement instead. The estimator is unchanged: the mean over the
    resampled rows, computed exactly as a ratio of summed group totals to
    summed group sizes, so no row is materialised twice.

    `groups=None` keeps the row bootstrap bit-for-bit, because the frozen HLC
    analyses were published with it and a row there is a donor, which *is* the
    unit of independence.
    """
    n_boot = FROZEN["n_bootstrap"] if n_boot is None else n_boot
    seed = FROZEN["seed"] if seed is None else seed
    alpha = FROZEN["alpha"] if alpha is None else alpha
    reference = np.asarray(reference_losses, dtype=float)
    ceiling = np.asarray(ceiling_losses, dtype=float)
    _check(reference.shape == ceiling.shape,
           "the two loss vectors must cover the same individuals")
    difference = reference - ceiling
    rng = np.random.default_rng(seed + 3)
    if groups is None:
        draws = difference[rng.integers(
            0, difference.size, size=(n_boot, difference.size))].mean(axis=1)
        unit, n_clusters = "individual", int(difference.size)
    else:
        labels = np.asarray(groups)
        _check(labels.shape == difference.shape,
               "one group label per row is required")
        keys, codes = np.unique(labels, return_inverse=True)
        _check(keys.size >= 3,
               f"a cluster bootstrap needs 3 clusters, got {keys.size}")
        # mean over a concatenation of clusters == sum of sums / sum of sizes,
        # so the cluster bootstrap is exact without building the index lists.
        sums = np.bincount(codes, weights=difference, minlength=keys.size)
        sizes = np.bincount(codes, minlength=keys.size).astype(float)
        picked = rng.integers(0, keys.size, size=(n_boot, keys.size))
        draws = sums[picked].sum(axis=1) / sizes[picked].sum(axis=1)
        unit, n_clusters = "cluster (group resampled with replacement)",             int(keys.size)
    half = 100.0 * alpha / (2.0 * family_size)
    low, high = (float(v) for v in np.percentile(draws, [half, 100.0 - half]))
    return {"gap": float(difference.mean()),
            "interval": [low, high],
            "excludes_zero": bool(low > 0.0),
            "level": f"{100.0 * (1.0 - alpha / family_size):.2f} %",
            "family_size": int(family_size),
            "n_bootstrap": int(n_boot),
            "n_resampling_units": n_clusters,
            "resampling_unit": unit}


def verdict(gap, reference_mean_loss, *, margin=None):
    """Three outcomes, and only one of them is a refutation.

    An interval containing zero is **not** a refutation; at the cohort sizes
    this project has it usually means low power. Refuting requires an
    equivalence statement: the gap's **upper** bound must fall below the
    smallest gain that would make the method worth building.
    """
    margin = (FROZEN["minimum_interesting_relative_gain"] if margin is None
              else margin)
    _check(reference_mean_loss > 0.0, "the reference loss must be positive")
    threshold = margin * reference_mean_loss
    low, high = gap["interval"]
    if low > 0.0:
        name = "HEADROOM_PRESENT"
        reading = ("the ceiling is detectably above the reference. A method "
                   "could in principle recover some of it; whether one can is "
                   "a different question.")
    elif high < threshold:
        name = "REFUTED_CEILING_BELOW_MARGIN"
        reading = (f"the ceiling's upper bound ({high:+.4f}) is below the "
                   f"smallest interesting gain ({threshold:+.4f}). Even in "
                   "hindsight there is less to win than would be worth "
                   "winning. This IS a refutation.")
    else:
        name = "NO_DETECTABLE_HEADROOM_UNDERPOWERED"
        reading = (f"the interval contains zero but reaches {high:+.4f}, above "
                   f"the {threshold:+.4f} that would matter. Nothing is "
                   "detected and nothing is refuted; the design lacks the "
                   "power to tell those apart.")
    return {"verdict": name, "how_to_read": reading,
            "minimum_interesting_gain": float(threshold),
            "relative_margin": float(margin),
            "minimum_detectable_gap": float(
                (gap["interval"][1] - gap["gap"]) or float("nan")),
            "an_interval_containing_zero_is_not_a_refutation": True}


# --------------------------------------------------------------------------
# the counterexample, kept in code because it is the reason this module exists
# --------------------------------------------------------------------------
def counterexample_constant_weights_cannot_see_a_gate(n=200):
    """Two arms, each MAE 0.5; every constant weight keeps 0.5; a gate gets 0.

    Half the cohort has the truth at one arm's prediction and half at the
    other's, and the two predictions bracket the truth everywhere. A constant
    weight must trade one half against the other; a gate does not.
    """
    _check(n % 2 == 0, "use an even cohort so the halves are equal")
    truth = np.zeros(n)
    arm_a = np.concatenate([np.zeros(n // 2), np.ones(n // 2)])
    arm_b = np.concatenate([np.ones(n // 2), np.zeros(n // 2)])
    context = np.concatenate([np.zeros(n // 2), np.ones(n // 2)])
    return np.column_stack([arm_a, arm_b]), truth, context


def self_test():
    """Assert that the ceilings separate the case the old rule got wrong."""
    predictions, truth, context = \
        counterexample_constant_weights_cannot_see_a_gate()
    single = best_single_oracle(predictions, truth)["mean_loss"]
    constant = constant_weight_oracle(predictions, truth)["mean_loss"]
    selection = individual_selection_oracle(predictions, truth)["mean_loss"]
    convex = individual_convex_oracle(predictions, truth)["mean_loss"]
    probe = context_informativeness(predictions, truth, context)

    _check(abs(single - 0.5) < 1e-9,
           f"each arm should have MAE 0.5, got {single}")
    _check(abs(constant - 0.5) < 1e-6,
           "the constant-weight ceiling must show NO gain here; if it shows "
           f"one the counterexample is broken (got {constant})")
    _check(selection < 1e-9, f"per-individual selection should reach 0, got {selection}")
    _check(convex < 1e-9, f"the per-individual convex ceiling should reach 0, got {convex}")
    _check(abs(probe["sign_agreement"] - 1.0) < 1e-9,
           "the context should perfectly predict which arm wins here, got "
           f"{probe['sign_agreement']}")
    # Positive control for the probe itself, now that the pair is chosen inside
    # each fold: if the ladder cannot award the margin on a cohort where the
    # context is a perfect gate, a null elsewhere would be uninterpretable.
    _check(probe["verdict"] == "SIGNAL_REACHES_MARGIN",
           "the probe must reach the declared margin on a perfect gate, got "
           f"{probe['verdict']} (rho {probe['spearman_out_of_sample']})")
    return {
        "best_single_oracle": single,
        "constant_weight_oracle": constant,
        "individual_selection_oracle": selection,
        "individual_convex_oracle": convex,
        "context_sign_agreement": probe["sign_agreement"],
        "context_probe_verdict": probe["verdict"],
        "arm_pair_chosen_inside_each_fold": probe[
            "arm_pair_chosen_inside_each_fold"],
        "what_this_proves": (
            "A constant-weight ceiling can show exactly zero headroom on a "
            "cohort where a contextual gate is perfect. Stopping the "
            "programme on the constant-weight ceiling alone would have "
            "discarded the hypothesis it was meant to test."),
    }


def what_this_module_does_not_do():
    return {
        "implements_no_fusion_method": True,
        "fits_no_gate": True,
        "ceilings_are_not_performances": (
            "every ceiling here is computed knowing the answers and is an "
            "upper bound on what any honest method could reach. None may be "
            "quoted as a result."),
        "a_permissive_ceiling_licenses_nothing": (
            "the per-individual convex ceiling is zero wherever the arms "
            "bracket the truth. It can refute; it cannot promise."),
        "non_detection_is_not_refutation": (
            "an interval containing zero yields "
            "NO_DETECTABLE_HEADROOM_UNDERPOWERED, never a refutation."),
    }
