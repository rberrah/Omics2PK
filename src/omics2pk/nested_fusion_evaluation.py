"""Comparators `S`, `U`, `K` and `X` under a genuinely nested evaluation.

What this module exists to fix
------------------------------
The headroom diagnostics reused **one** layer of out-of-fold predictions
computed over the whole cohort, and then split *those predictions* to fit and
score the stacking weights. The base predictions were honestly out-of-fold for
their own individual, so this looked nested. It is not. Writing `keep` for the
weight-fitting rows and `held` for the scoring rows:

* individual *i* in `held` has a prediction from a base model trained without
  *i*, so `yhat_i` is honest with respect to `y_i`;
* but every individual *j* in `keep` has a prediction from a base model whose
  training set **contained** *i*, so `yhat_j` depends on `y_i`;
* the weights are fitted on `{yhat_j, y_j}` for *j* in `keep`, and therefore
  depend on `y_i` — the outcome of an individual they are then scored on.

The leak is through the base models' training sets, which no amount of
splitting the prediction layer can undo. The only repair is to refit the base
models inside the outer loop:

    for each outer fold:
        inner out-of-fold predictions on the outer-TRAIN rows only,
            each from arms refit on a subset of the outer-train rows
        -> fit every comparator on those, and nothing else
        arms refit on ALL outer-train rows -> predict the outer-TEST rows
        -> score every comparator there

No outer-test outcome enters any fit, at any level, through any path. Every
comparator sees the same outer folds, so the comparisons are paired at the
individual level, which is what the commissioning review asked for in asking
for harmonised folds.

The comparators, as the protocol declares them
----------------------------------------------
=====  =====================================================================
`S`    the single arm with the lowest inner out-of-fold loss. The arm to beat.
`U`    the uniform average. The control known to lose.
`K`    constant weights fitted by Frank-Wolfe on the inner out-of-fold
       predictions. **This is stacking.**
`X`    a contextual gate with **no mechanistic prior**: per arm, a ridge
       predicts that arm's loss from the context, and the weights are a
       softmax of minus the predicted losses. Reduces to a constant mixture
       when the ridge coefficients vanish, and to hard per-individual
       selection as the temperature goes to zero.
=====  =====================================================================

`X` is deliberately the smallest contextual gate that can express a gate at
all. Making it elaborate before knowing whether context helps would confound
*does context help* with *is this particular architecture good*.

What this module does NOT do
----------------------------
It fits no mechanistic prior, so it cannot answer Q2 (`F` vs `X`). It is the
machinery for Q1 (`X` vs `K`) and for the controls that make a Q1 answer
readable. `X` is permitted to lose, and `contrast_verdict` has a name for
losing detectably, because that is the outcome the project's own history makes
most likely.
"""
from __future__ import annotations

import numpy as np

from .fusion_headroom_diagnostics import (
    FROZEN as HEADROOM_FROZEN,
    constant_weight_oracle,
    expand_basis,
    individual_convex_oracle,
    paired_gap,
    ridge_predict,
)

#: Frozen before execution.
FROZEN = {
    "seed": 20261112,
    "outer_folds": 5,
    "inner_folds": 5,
    "n_bootstrap": 2000,
    "alpha": 0.05,
    # Same margin as the headroom diagnostics, so a gain and a ceiling are
    # judged against the same yardstick.
    "minimum_interesting_relative_gain": 0.05,
    "sensitivity_relative_gain": 0.10,
    # Multiplies the median inner loss, so the grid is scale-free and the same
    # numbers mean the same thing on activity units and on ng/mL.
    "temperature_grid": (0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0),
    "ridge_penalty": 1.0,
    "gate_basis": "linear",
    "loss": "absolute_error",
    "comparators": ("S", "U", "K", "X"),
}

COMPARATORS = ("S", "U", "K", "X")


class NestedEvaluationError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise NestedEvaluationError(message)


def _folds(n, count, rng):
    order = rng.permutation(n)
    return [np.sort(order[index::count]) for index in range(count)]


def _softmax(scores):
    shifted = scores - scores.max(axis=1, keepdims=True)
    exponentiated = np.exp(shifted)
    return exponentiated / exponentiated.sum(axis=1, keepdims=True)


def _augment(rows, predictions, mode):
    """Build the gate's design: clinical columns, the arms' own spread, or both.

    The arms' predictions and their disagreement are legitimate context: they
    are computable for an individual **without** that individual's outcome,
    which is the only requirement. They are also where the context probe found
    its strongest association, so excluding them would let a null result rest
    on having looked in one place. They are taken from the fold's **own**
    predictions — inner out-of-fold for the training rows, outer for the test
    rows — so nothing crosses the split.
    """
    spread = predictions.max(axis=1) - predictions.min(axis=1)
    endogenous = np.column_stack([predictions, spread,
                                  predictions.mean(axis=1),
                                  predictions.std(axis=1)])
    if mode == "clinical":
        return rows
    if mode == "predictions":
        return endogenous
    _check(mode == "both", f"unknown context mode {mode!r}")
    return np.column_stack([rows, endogenous])


def _gate_weights(train_design, train_losses, test_design, temperature, *,
                  penalty):
    """Softmax of minus each arm's predicted loss. Convex by construction.

    Being convex matters: it keeps `X` inside the family bounded by the
    per-individual convex ceiling, so the ceiling remains a ceiling for it.
    """
    predicted = np.column_stack([
        ridge_predict(train_design, train_losses[:, arm], test_design,
                      penalty=penalty)
        for arm in range(train_losses.shape[1])])
    return _softmax(-predicted / temperature), predicted


def nested_fusion_evaluation(fit_predict, truth, context, *, arm_names,
                             seed=None, outer_folds=None, inner_folds=None,
                             temperature_grid=None, gate_basis=None,
                             ridge_penalty=None, permute_context=False,
                             context_mode="clinical"):
    """Run `S`, `U`, `K` and `X` on shared outer folds with nothing leaked.

    `fit_predict(train_index, test_index)` must refit every arm on
    `train_index` and return `{arm_name: predictions_for_test_index}`. It is
    called `outer_folds * (inner_folds + 1)` times. Refitting is the whole
    point: a cached prediction layer is what made the previous version leak.

    `permute_context=True` is the protocol's **permuted context** negative
    control. It shuffles the gate's design rows **within each side of each
    outer split**, which leaves every arm, every fold and `S`/`U`/`K`
    untouched and destroys only the correspondence `X` is supposed to use. If
    `X` still beats `K` under it, the gain is an artefact of the extra
    flexibility and not of the context. Permuting inside the fold rather than
    once globally is what makes the control valid for `context_mode` values
    that include the arms' own predictions: the gate then receives another
    individual's disagreement pattern, while the predictions actually being
    fused stay attached to their own individual.
    """
    seed = FROZEN["seed"] if seed is None else seed
    outer_folds = FROZEN["outer_folds"] if outer_folds is None else outer_folds
    inner_folds = FROZEN["inner_folds"] if inner_folds is None else inner_folds
    grid = tuple(FROZEN["temperature_grid"] if temperature_grid is None
                 else temperature_grid)
    gate_basis = FROZEN["gate_basis"] if gate_basis is None else gate_basis
    penalty = FROZEN["ridge_penalty"] if ridge_penalty is None else ridge_penalty

    target = np.asarray(truth, dtype=float)
    n = target.size
    block = np.asarray(context, dtype=float)
    if block.ndim == 1:
        block = block.reshape(-1, 1)
    _check(block.shape[0] == n, "the context block and the outcomes disagree on n")
    _check(np.isfinite(target).all(), "an outcome is missing")
    _check(np.isfinite(block).all(), "a context value is missing")
    _check(len(arm_names) >= 2, "a fusion comparison needs at least two arms")
    _check(n >= outer_folds * 4,
           f"n={n} is too small for {outer_folds} outer folds with an inner "
           f"{inner_folds}-fold inside each")

    rng = np.random.default_rng(seed)
    permute_rng = np.random.default_rng(seed + 7)
    arms = len(arm_names)
    scored = {name: np.full(n, np.nan) for name in COMPARATORS}
    chosen_arm, chosen_temperature = [], []
    constant_weights, gate_weight_sums, gate_weight_rows = [], np.zeros(arms), 0
    disagreement_with_selection = []
    outer_predictions = np.full((n, arms), np.nan)

    for held in _folds(n, outer_folds, rng):
        keep = np.setdiff1d(np.arange(n), held)

        # --- inner out-of-fold predictions, arms refit without the scored rows
        inner = np.full((keep.size, arms), np.nan)
        inner_rng = np.random.default_rng(seed + 1)
        for inner_held in _folds(keep.size, inner_folds, inner_rng):
            inner_keep = np.setdiff1d(np.arange(keep.size), inner_held)
            got = fit_predict(keep[inner_keep], keep[inner_held])
            inner[inner_held] = np.column_stack(
                [np.asarray(got[name], dtype=float) for name in arm_names])
        _check(np.isfinite(inner).all(),
               "fit_predict left an inner out-of-fold prediction missing")

        # --- arms refit on ALL the outer-train rows, predicting the held rows
        got = fit_predict(keep, held)
        outer = np.column_stack(
            [np.asarray(got[name], dtype=float) for name in arm_names])
        _check(np.isfinite(outer).all(),
               "fit_predict left an outer prediction missing")
        outer_predictions[held] = outer

        inner_target = target[keep]
        inner_losses = np.abs(inner - inner_target[:, None])
        held_losses = np.abs(outer - target[held][:, None])

        # --- S: the arm with the lowest inner out-of-fold loss
        index = int(inner_losses.mean(axis=0).argmin())
        chosen_arm.append(arm_names[index])
        scored["S"][held] = held_losses[:, index]

        # --- U: uniform
        scored["U"][held] = np.abs(outer.mean(axis=1) - target[held])

        # --- K: constant weights fitted on the inner out-of-fold predictions
        fitted = constant_weight_oracle(inner, inner_target)
        weights = np.array([fitted["weights"][f"arm_{arm}"]
                            for arm in range(arms)])
        constant_weights.append(weights)
        scored["K"][held] = np.abs(outer @ weights - target[held])

        # --- X: contextual gate. The temperature is chosen on a split of the
        # outer-train rows, so no outer-test row sets it either.
        design_keep = expand_basis(
            _augment(block[keep], inner, context_mode), gate_basis)
        design_held = expand_basis(
            _augment(block[held], outer, context_mode), gate_basis)
        if permute_context:
            design_keep = design_keep[permute_rng.permutation(keep.size)]
            design_held = design_held[permute_rng.permutation(held.size)]
        scale = float(np.median(inner_losses))
        scale = scale if scale > 0.0 else 1.0
        split = np.random.default_rng(seed + 2).permutation(keep.size)
        first, second = split[::2], split[1::2]
        best_temperature, best_value = grid[0], np.inf
        for candidate in grid:
            trial, _ = _gate_weights(design_keep[first], inner_losses[first],
                                     design_keep[second],
                                     candidate * scale, penalty=penalty)
            value = float(np.abs((trial * inner[second]).sum(axis=1)
                                 - inner_target[second]).mean())
            if value < best_value:
                best_temperature, best_value = candidate, value
        chosen_temperature.append(float(best_temperature))
        gate, _ = _gate_weights(design_keep, inner_losses, design_held,
                                best_temperature * scale, penalty=penalty)
        scored["X"][held] = np.abs((gate * outer).sum(axis=1) - target[held])
        gate_weight_sums += gate.sum(axis=0)
        gate_weight_rows += gate.shape[0]
        disagreement_with_selection.append(
            float(np.mean(gate.argmax(axis=1) != index)))

    for name, values in scored.items():
        _check(np.isfinite(values).all(), f"{name} left an individual unscored")

    convex = individual_convex_oracle(outer_predictions, target)
    mean_constant = np.mean(np.array(constant_weights), axis=0)
    return {
        "n_individuals": int(n),
        "arms": list(arm_names),
        "outer_folds": int(outer_folds),
        "inner_folds": int(inner_folds),
        "permuted_context": bool(permute_context),
        "context_mode": context_mode,
        "gate_basis": gate_basis,
        "losses": {name: scored[name] for name in COMPARATORS},
        # Individual-level. Kept in memory for the context probe and the
        # ceiling; it must never reach a committed artefact.
        "outer_predictions": outer_predictions,
        "mean_loss": {name: float(scored[name].mean())
                      for name in COMPARATORS},
        "arm_chosen_per_outer_fold": chosen_arm,
        "gate_temperature_per_outer_fold": chosen_temperature,
        "mean_constant_weights": {arm_names[i]: float(mean_constant[i])
                                  for i in range(arms)},
        "mean_gate_weights": {
            arm_names[i]: float(gate_weight_sums[i] / gate_weight_rows)
            for i in range(arms)},
        "share_where_the_gate_disagrees_with_selection": float(
            np.mean(disagreement_with_selection)),
        "ceiling_on_these_outer_predictions": {
            "individual_convex_oracle": convex["mean_loss"],
            "bracketing_rate": convex["bracketing_rate"],
        },
        "nothing_leaked": (
            "the arms are refit inside each outer fold; every comparator is "
            "fitted only on inner out-of-fold predictions from arms that "
            "never saw an outer-test row; the gate temperature is chosen on a "
            "split of the outer-train rows."),
    }


def contrast(result, better, worse, *, family_size=1, margin=None,
             n_boot=None, alpha=None, seed=None):
    """Paired interval and verdict for `worse - better`, positive = `better` wins.

    Four outcomes, and only two of them are conclusions. A verdict of
    `DETECTABLY_WORSE` is a real scientific result and is reported with the
    same prominence as a gain.
    """
    margin = (FROZEN["minimum_interesting_relative_gain"] if margin is None
              else margin)
    losses = result["losses"]
    _check(better in losses and worse in losses,
           f"unknown comparator in {better!r} vs {worse!r}")
    gap = paired_gap(losses[worse], losses[better],
                     family_size=family_size,
                     n_boot=FROZEN["n_bootstrap"] if n_boot is None else n_boot,
                     alpha=FROZEN["alpha"] if alpha is None else alpha,
                     seed=FROZEN["seed"] if seed is None else seed)
    reference = float(np.mean(losses[worse]))
    _check(reference > 0.0, "the reference loss must be positive")
    threshold = margin * reference
    low, high = gap["interval"]
    if low > 0.0:
        # The same correction the commissioning review forced on the context
        # probe, applied here before it could mislead: detectably positive and
        # large enough to matter are different claims, and a gain of 0.6 % does
        # not become interesting by having an interval that excludes zero. The
        # three rungs mirror the probe ladder exactly, on purpose.
        if low >= threshold:
            name = "GAIN_REACHES_MARGIN"
            reading = (f"{better} beats {worse} by at least the smallest "
                       f"interesting gain ({threshold:+.4f}). Compare it to "
                       "the ceiling before calling it large.")
        elif high < threshold:
            name = "GAIN_PRESENT_BUT_BELOW_MARGIN"
            reading = (f"{better} beats {worse} out of sample, and the whole "
                       f"interval sits below the {threshold:+.4f} that would "
                       "make it worth building. A real effect, too small to "
                       "act on.")
        else:
            name = "GAIN_PRESENT_MARGIN_UNDECIDED"
            reading = (f"{better} beats {worse} out of sample, but the "
                       "interval straddles the smallest interesting gain "
                       f"({threshold:+.4f}), so whether it is worth building "
                       "is undecided.")
    elif high < 0.0:
        name = "DETECTABLY_WORSE"
        reading = (f"{better} is detectably WORSE than {worse}. This is a "
                   "result, not a failed run.")
    elif high < threshold:
        name = "REFUTED_GAIN_BELOW_MARGIN"
        reading = (f"the gain's upper bound ({high:+.4f}) is below the "
                   f"smallest interesting gain ({threshold:+.4f}). Even the "
                   "optimistic end of the interval is not worth building. "
                   "This IS a refutation.")
    else:
        name = "NO_DETECTABLE_GAIN_UNDERPOWERED"
        reading = (f"the interval contains zero but reaches {high:+.4f}, above "
                   f"the {threshold:+.4f} that would matter. Nothing is shown "
                   "and nothing is refuted.")
    return {
        "contrast": f"{better}_minus_{worse}",
        "relative_gain": float(gap["gap"] / reference),
        "verdict": name,
        "how_to_read": reading,
        "minimum_interesting_gain": float(threshold),
        "relative_margin": float(margin),
        "an_interval_containing_zero_is_not_a_refutation": True,
        **gap,
    }


# --------------------------------------------------------------------------
# positive control: the harness must be able to SEE a gate that exists
# --------------------------------------------------------------------------
def gateable_cohort(n=400, seed=0):
    """Two arms that bracket the truth on opposite halves; context says which.

    The same shape as the counterexample that refuted the constant-weight
    ceiling, but with noise and a real split so the nested harness has
    something to learn rather than something to memorise. `K` cannot help
    here; a working gate nearly eliminates the loss.
    """
    rng = np.random.default_rng(seed)
    side = rng.integers(0, 2, size=n)
    truth = rng.normal(0.0, 1.0, size=n)
    noise = rng.normal(0.0, 0.05, size=n)
    arm_a = np.where(side == 0, truth + noise, truth + 1.0)
    arm_b = np.where(side == 0, truth + 1.0, truth + noise)
    context = np.column_stack([side.astype(float),
                               rng.normal(0.0, 1.0, size=n)])
    return arm_a, arm_b, truth, context


def self_test():
    """Assert the harness detects a real gate, and reports none when absent.

    Without this, a null `X` vs `K` on real data would be uninterpretable: it
    could mean the context carries nothing, or that the gate is broken. The
    second case is excluded here.
    """
    arm_a, arm_b, truth, context = gateable_cohort()
    arms = {"arm_a": arm_a, "arm_b": arm_b}

    def fit_predict(train_index, test_index):
        # Fixed arms, so this isolates the fusion layer from base-model noise.
        del train_index
        return {name: values[test_index] for name, values in arms.items()}

    live = nested_fusion_evaluation(fit_predict, truth, context,
                                    arm_names=list(arms))
    dead = nested_fusion_evaluation(fit_predict, truth, context,
                                    arm_names=list(arms), permute_context=True)
    gain = contrast(live, "X", "K")
    control = contrast(dead, "X", "K")

    _check(gain["verdict"] == "GAIN_REACHES_MARGIN",
           "the harness must detect a gate that exists; got "
           f"{gain['verdict']} ({gain['gap']:+.4f})")
    _check(live["mean_loss"]["X"] < 0.5 * live["mean_loss"]["K"],
           "a near-perfect gate should roughly halve the stacking loss, got "
           f"X={live['mean_loss']['X']:.4f} K={live['mean_loss']['K']:.4f}")
    _check(not control["verdict"].startswith("GAIN_"),
           "with the context permuted the gate must lose its advantage, got "
           f"{control['verdict']} ({control['gap']:+.4f})")
    return {
        "with_real_context": {
            "mean_loss": live["mean_loss"],
            "X_vs_K": gain["verdict"],
            "relative_gain": gain["relative_gain"],
        },
        "with_permuted_context": {
            "mean_loss": dead["mean_loss"],
            "X_vs_K": control["verdict"],
            "relative_gain": control["relative_gain"],
        },
        "what_this_proves": (
            "the nested harness recovers a contextual gain that is really "
            "there, and loses it when the context is permuted. A null X vs K "
            "on a real cohort is therefore a statement about that cohort and "
            "not about this implementation."),
        "headroom_margin_shared_with_the_ceilings": (
            HEADROOM_FROZEN["minimum_interesting_relative_gain"]),
    }
