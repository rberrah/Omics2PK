"""Simulation-based cohort sizing for the three open questions.

Why a simulation rather than a formula. A closed-form sample size needs an
effect expressed as a standardised mean difference on a quantity that is
measured directly. The quantity this project actually reports is the
**paired difference in out-of-sample absolute error between two arms**, after
a penalty chosen by inner cross-validation and, for a wide block, after feature
selection inside the training folds. Shrinkage, selection and the pairing all
change the distribution of that difference in ways no closed form captures --
the Axis A run showed the effect concretely, where a wider arm scored better
under a permutation null purely because it shrank harder. Simulating the whole
pipeline is the only way to size the quantity that will actually be reported.

What this module does **not** do. It builds no new predictive architecture: the
arms are the project's existing ridge-with-inner-CV arms. It fuses no
individuals across cohorts -- real data enters only as **declared scalar
parameters** (a dispersion, a correlation) read off published aggregates, never
as rows joined together. And it never converts a published rank correlation
into a fraction of explained variance by squaring it; the one place a rank
correlation is used, it passes through an explicitly declared Gaussian-copula
mapping and the assumption is recorded with the number.

The two-level rule this follows
-------------------------------
* **Level 1, transferable.** Where the public record contains an effect of the
  same kind on the same sort of endpoint, it is used and its source is named.
* **Level 2, hypothetical.** Where it does not, a declared grid of effects is
  used and every figure derived from it is labelled a **scenario**, not an
  expected effect.

Each sized question states which level it is on. Question A is level 1.
Questions B and C are level 2, and no in-vitro effect is carried into either.
"""
from __future__ import annotations

import math

import numpy as np

#: Frozen before execution.
FROZEN = {
    "seed": 20261110,
    "replicates": 200,
    "outer_folds": 5,
    "inner_folds": 3,
    "ridge_alphas": (0.001, 0.01, 0.1, 1.0, 10.0, 100.0),
    "target_power": 0.80,
    "alpha": 0.05,
    "search_low": 8,
    "search_high": 1200,
    "evidence_levels": ("transferable_literature", "declared_scenario"),
}

#: Validation is 5-fold here and leave-one-donor-out in the executed pilots.
#: The reason is declared rather than quietly chosen: leave-one-out at n = 600
#: costs n^2 fits for no statistical gain, and the sizing must use ONE scheme
#: across the whole n range or the required-n curve would mix two estimators.
VALIDATION_NOTE = (
    "5-fold outer, 3-fold inner. The executed pilots use leave-one-donor-out "
    "because they have fifteen donors. Mixing the two across an n range from "
    "10 to 2000 would make the required-n curve an artefact of the scheme, so "
    "one scheme is used throughout and the difference is recorded here."
)


class CohortSizeSimulationError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise CohortSizeSimulationError(message)


# --------------------------------------------------------------------------
# the one place a published rank correlation is used
# --------------------------------------------------------------------------
def spearman_to_pearson(rho_spearman):
    """Gaussian-copula mapping, with its assumption attached.

    rho_P = 2 sin(pi * rho_S / 6). This is exact only if the pair is a
    Gaussian copula, i.e. if both variables become jointly normal after a
    monotone transform. On log-transformed enzyme abundances and log activity
    rates that is a reasonable working assumption and it is an ASSUMPTION.

    This is the project's permitted route from a published Spearman to a
    simulation parameter. Squaring the Spearman and calling the result a
    fraction of explained variance is not permitted and is not done anywhere
    in this module.
    """
    rho = float(rho_spearman)
    _check(-1.0 <= rho <= 1.0, "a rank correlation must lie in [-1, 1]")
    return 2.0 * math.sin(math.pi * rho / 6.0)


# --------------------------------------------------------------------------
# the model, imported in spirit from the frozen pilot and kept small
# --------------------------------------------------------------------------
def _standardise(train, test):
    mean = train.mean(axis=0)
    sd = train.std(axis=0, ddof=0)
    safe = np.where(sd > 0.0, sd, 1.0)
    return (train - mean) / safe, (test - mean) / safe


def _ridge(design_train, target_train, design_test, alpha):
    centre = float(target_train.mean())
    if design_train.shape[1] == 0:
        return np.full(design_test.shape[0], centre)
    z_train, z_test = _standardise(design_train, design_test)
    gram = z_train.T @ z_train + alpha * np.eye(z_train.shape[1])
    coefficients = np.linalg.solve(gram, z_train.T @ (target_train - centre))
    return z_test @ coefficients + centre


def _folds(n, k, rng):
    order = rng.permutation(n)
    return [order[i::k] for i in range(k)]


def _select(block, target, q):
    """Top-q features by absolute correlation, computed on TRAINING rows only."""
    if block.shape[1] <= q:
        return np.arange(block.shape[1])
    centred = block - block.mean(axis=0)
    sd = centred.std(axis=0, ddof=0)
    safe = np.where(sd > 0.0, sd, np.inf)
    scores = np.abs((centred / safe).T @ (target - target.mean()))
    return np.argsort(scores)[::-1][:q]


def _inner_alpha(block, target, alphas, inner_k, rng):
    if block.shape[1] == 0:
        return float(alphas[0])
    folds = _folds(block.shape[0], min(inner_k, block.shape[0]), rng)
    errors = np.zeros(len(alphas))
    for held in folds:
        keep = np.setdiff1d(np.arange(block.shape[0]), held,
                            assume_unique=False)
        if keep.size < 2 or held.size == 0:
            continue
        for index, alpha in enumerate(alphas):
            predicted = _ridge(block[keep], target[keep], block[held], alpha)
            errors[index] += np.abs(predicted - target[held]).sum()
    return float(alphas[int(np.argmin(errors))])


def cross_validated_absolute_errors(block, target, *, q, rng,
                                    outer_k=None, inner_k=None, alphas=None):
    """Out-of-fold |error| per individual, with selection inside the folds.

    Both the feature selection and the penalty are re-estimated on the
    training rows of every outer fold, so neither sees the held-out rows.
    """
    outer_k = FROZEN["outer_folds"] if outer_k is None else outer_k
    inner_k = FROZEN["inner_folds"] if inner_k is None else inner_k
    alphas = FROZEN["ridge_alphas"] if alphas is None else alphas
    n = target.size
    errors = np.zeros(n)
    for held in _folds(n, min(outer_k, n), rng):
        keep = np.setdiff1d(np.arange(n), held)
        if keep.size < 3 or held.size == 0:
            errors[held] = np.abs(target[held] - target.mean())
            continue
        if block.shape[1] == 0:
            errors[held] = np.abs(float(target[keep].mean()) - target[held])
            continue
        chosen = _select(block[keep], target[keep], q)
        train, test = block[keep][:, chosen], block[held][:, chosen]
        alpha = _inner_alpha(train, target[keep], alphas, inner_k, rng)
        errors[held] = np.abs(_ridge(train, target[keep], test, alpha)
                              - target[held])
    return errors


# --------------------------------------------------------------------------
# paired test on the quantity that gets reported
# --------------------------------------------------------------------------
def paired_detection(errors_arm, errors_baseline, *, alpha, family_size):
    """Did this replicate detect an improvement, at the adjusted level?

    A paired two-sided t statistic on the per-individual differences in
    absolute error. Detection requires BOTH significance at alpha/family and
    the right sign, so a significant degradation never counts as a detection.
    """
    difference = np.asarray(errors_arm) - np.asarray(errors_baseline)
    n = difference.size
    if n < 3:
        return False
    sd = difference.std(ddof=1)
    if not np.isfinite(sd) or sd <= 0.0:
        return bool(difference.mean() < 0.0)
    statistic = difference.mean() / (sd / math.sqrt(n))
    # Normal approximation to the two-sided t quantile; at the n this module
    # searches over (>= 8) the difference from the exact t quantile changes no
    # reported required-n by more than one step of the search.
    critical = _two_sided_normal_quantile(alpha / family_size)
    return bool(statistic < -critical)


def _two_sided_normal_quantile(alpha):
    """Inverse normal at 1 - alpha/2, by bisection; no scipy dependency."""
    target = 1.0 - alpha / 2.0
    low, high = 0.0, 12.0
    for _ in range(200):
        middle = (low + high) / 2.0
        if 0.5 * (1.0 + math.erf(middle / math.sqrt(2.0))) < target:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


# --------------------------------------------------------------------------
# question A -- does the enzyme's abundance predict its activity? (level 1)
# --------------------------------------------------------------------------
def simulate_enzyme_activity(n, *, rho_pearson_log, analytical_cv_log,
                             missing_rate, rng, alpha, family_size):
    """One replicate of the A design: one predictor against the mean.

    Both the abundance and the activity are generated on the log10 scale,
    because that is the scale on which the project's endpoints are modelled
    and on which these quantities are roughly symmetric.
    """
    latent = rng.normal(size=n)
    activity = (rho_pearson_log * latent
                + math.sqrt(max(1.0 - rho_pearson_log ** 2, 0.0))
                * rng.normal(size=n))
    measured = latent + analytical_cv_log * rng.normal(size=n)
    if missing_rate > 0.0:
        lost = rng.random(n) < missing_rate
        # A donor with no abundance cannot be predicted by the abundance arm
        # and is carried at the baseline's error, which is what dropping the
        # row would NOT do: dropping would quietly shrink n for one arm only.
        measured = np.where(lost, np.nan, measured)
    keep = np.isfinite(measured)
    block = np.where(keep, measured, 0.0).reshape(-1, 1)
    baseline = np.zeros((n, 0))
    errors_arm = cross_validated_absolute_errors(
        block, activity, q=1, rng=rng)
    errors_base = cross_validated_absolute_errors(
        baseline, activity, q=0, rng=rng)
    errors_arm = np.where(keep, errors_arm, errors_base)
    return paired_detection(errors_arm, errors_base, alpha=alpha,
                            family_size=family_size)


# --------------------------------------------------------------------------
# questions B and C -- does a layer add beyond what is already known? (level 2)
# --------------------------------------------------------------------------
def simulate_incremental_layer(n, *, baseline_r2, target_mae_reduction,
                               n_features, n_informative, redundancy,
                               analytical_cv, missing_rate, selected,
                               rng, alpha, family_size):
    """One replicate of the B/C design: baseline+layer against baseline.

    ``target_mae_reduction`` is the fractional reduction in mean absolute
    error the layer would deliver **if it were measured without error and
    carried no redundancy with the baseline**. It is a declared scenario, not
    an estimate, and the realised reduction is lower once analytical noise,
    redundancy, selection among ``n_features`` and missingness are applied --
    which is the whole reason for simulating rather than solving.
    """
    _check(0.0 <= baseline_r2 < 1.0, "baseline_r2 must lie in [0, 1)")
    _check(0.0 <= target_mae_reduction < 1.0,
           "a reduction must be a fraction below 1")
    _check(1 <= n_informative <= n_features,
           "the informative features must be a non-empty subset")

    # A reduction f in MAE corresponds, for a Gaussian residual, to a residual
    # standard deviation multiplied by (1 - f): MAE is proportional to sigma.
    # So the layer must explain the fraction 1 - (1 - f)^2 of what the
    # baseline leaves unexplained.
    residual_share = 1.0 - baseline_r2
    layer_share = residual_share * (1.0 - (1.0 - target_mae_reduction) ** 2)

    baseline_signal = rng.normal(size=n)
    layer_signal = rng.normal(size=n)
    # Redundancy: the part of the layer's signal that the baseline already
    # carries. At redundancy 1 the layer can add nothing however well measured.
    layer_signal = (redundancy * baseline_signal
                    + math.sqrt(max(1.0 - redundancy ** 2, 0.0))
                    * layer_signal)

    outcome = (math.sqrt(baseline_r2) * baseline_signal
               + math.sqrt(layer_share) * layer_signal
               + math.sqrt(max(residual_share - layer_share, 0.0))
               * rng.normal(size=n))

    baseline_block = (baseline_signal
                      + 0.10 * rng.normal(size=n)).reshape(-1, 1)
    # Each informative feature carries the shared layer signal with loading
    # 1/sqrt(n_informative), so that the informative set jointly reconstructs
    # the signal and no single feature does. Uninformative features are pure
    # noise and exist to compete in the selection step.
    loading = 1.0 / math.sqrt(n_informative)
    layer = rng.standard_normal((n, n_features), dtype=np.float32)
    layer[:, :n_informative] *= np.float32(
        math.sqrt(max(1.0 - loading ** 2, 0.0)))
    layer[:, :n_informative] += np.float32(loading) * layer_signal.astype(
        np.float32)[:, None]
    layer += np.float32(analytical_cv) * rng.standard_normal(
        layer.shape, dtype=np.float32)
    if missing_rate > 0.0:
        mask = rng.random(layer.shape) < missing_rate
        # Mean imputation computed per column on the whole column is a
        # deliberate simplification and is declared: a fold-safe imputation
        # would only reduce the layer's apparent value, so this choice is
        # OPTIMISTIC for the layer and the required n it yields is a lower
        # bound, not a best guess.
        column_means = np.nanmean(np.where(mask, np.nan, layer), axis=0)
        layer = np.where(mask, column_means, layer)

    combined = np.hstack([baseline_block, layer])
    errors_base = cross_validated_absolute_errors(
        baseline_block, outcome, q=1, rng=rng)
    errors_arm = cross_validated_absolute_errors(
        combined, outcome, q=1 + selected, rng=rng)
    return paired_detection(errors_arm, errors_base, alpha=alpha,
                            family_size=family_size)


# --------------------------------------------------------------------------
# power and the search over n
# --------------------------------------------------------------------------
def power_at(n, simulator, scenario, *, replicates=None, seed=None,
             alpha=None, family_size=1):
    replicates = FROZEN["replicates"] if replicates is None else replicates
    seed = FROZEN["seed"] if seed is None else seed
    alpha = FROZEN["alpha"] if alpha is None else alpha
    rng = np.random.default_rng(seed + int(n))
    detected = sum(
        simulator(int(n), rng=rng, alpha=alpha, family_size=family_size,
                  **scenario)
        for _ in range(replicates))
    return detected / replicates


def required_sample_size(simulator, scenario, *, target_power=None,
                         replicates=None, seed=None, alpha=None,
                         family_size=1, low=None, high=None):
    """Smallest n reaching the target power, by bisection on a monotone curve.

    Power is monotone in n only in expectation; with 200 replicates the curve
    is noisy, so the bisection can land one step off. ``searched`` records
    every n evaluated and its power, so a reader can see the curve rather
    than trusting the single returned number.
    """
    target_power = (FROZEN["target_power"] if target_power is None
                    else target_power)
    low = FROZEN["search_low"] if low is None else low
    high = FROZEN["search_high"] if high is None else high
    searched = {}

    def power(n):
        if n not in searched:
            searched[n] = power_at(n, simulator, scenario,
                                   replicates=replicates, seed=seed,
                                   alpha=alpha, family_size=family_size)
        return searched[n]

    if power(high) < target_power:
        return {"required_n": None,
                "reached_at_ceiling": power(high),
                "ceiling": high,
                "searched": dict(sorted(searched.items())),
                "verdict": (f"the target power is not reached even at "
                            f"n = {high}; under this scenario the question is "
                            "not answerable by enlarging the cohort alone")}
    if power(low) >= target_power:
        return {"required_n": low, "searched": dict(sorted(searched.items())),
                "verdict": (f"the target power is already met at the search "
                            f"floor n = {low}; the required size is at or "
                            "below it and this design is not the binding "
                            "constraint")}
    bottom, top = low, high
    while top - bottom > max(2, int(0.02 * top)):
        middle = (bottom + top) // 2
        if power(middle) >= target_power:
            top = middle
        else:
            bottom = middle
    return {"required_n": top,
            "power_at_required_n": power(top),
            "searched": dict(sorted(searched.items())),
            "verdict": f"about {top} individuals for {target_power:.0%} power"}


def one_at_a_time_sensitivity(simulator, reference, variations, **kwargs):
    """Required n for the reference, then for one changed parameter at a time.

    This is what answers "which assumptions decide the answer". A parameter
    that moves the required n by a factor of two matters; one that moves it by
    ten per cent is noise at 200 replicates and must not be presented as a
    finding.
    """
    base = required_sample_size(simulator, reference, **kwargs)
    rows = {}
    for name, values in variations.items():
        for value in values:
            scenario = dict(reference)
            scenario[name] = value
            result = required_sample_size(simulator, scenario, **kwargs)
            rows[f"{name}={value}"] = {
                "required_n": result["required_n"],
                "verdict": result["verdict"],
            }
    ranked = sorted(
        ((key, row["required_n"]) for key, row in rows.items()),
        key=lambda pair: (pair[1] is None, -(pair[1] or 0)))
    return {
        "reference": {"scenario": reference,
                      "required_n": base["required_n"],
                      "verdict": base["verdict"]},
        "variations": rows,
        "ordered_by_required_n_descending": [k for k, _ in ranked],
        "how_to_read": (
            "The parameter whose variations span the widest range of required "
            "n is the one the sizing actually depends on. A variation that "
            "returns None is a scenario in which no cohort size inside the "
            "search ceiling answers the question."),
    }


# --------------------------------------------------------------------------
# independent validation is a different question from discovery
# --------------------------------------------------------------------------
def validation_cohort_size(effect_log10, dispersion_log10, *, power=0.80,
                           alpha=0.05, family_size=1):
    """How many individuals a CONFIRMATORY cohort needs for one stated effect.

    Discovery sizing asks how many are needed to find an effect among several
    candidates. Confirmation asks how many are needed to reject the null for
    ONE prespecified effect, with no selection and a family of one unless
    stated. The second is always smaller, and conflating them is how a funding
    case ends up asking for the wrong cohort.
    """
    _check(dispersion_log10 > 0.0, "a dispersion must be positive")
    _check(0.0 < power < 1.0, "power must lie in (0, 1)")
    standardised = abs(float(effect_log10)) / float(dispersion_log10)
    _check(standardised > 0.0, "a zero effect cannot be sized")
    z_alpha = _two_sided_normal_quantile(alpha / family_size)
    z_power = _two_sided_normal_quantile(2.0 * (1.0 - power))
    n = ((z_alpha + z_power) / standardised) ** 2
    return {
        "effect_log10": float(effect_log10),
        "dispersion_log10": float(dispersion_log10),
        "standardised_effect": standardised,
        "power": power,
        "alpha": alpha,
        "family_size": int(family_size),
        "required_n": int(math.ceil(n)),
        "design": "paired, two-sided, one prespecified contrast",
        "what_this_assumes": (
            "that the effect and the dispersion are the ones that will hold in "
            "the confirmatory cohort. Both come from a cohort of fifteen, so "
            "both carry wide uncertainty, and the dispersion from a small "
            "sample is biased low more often than high -- which makes this n "
            "an optimistic figure."),
        "what_this_is_not": (
            "not a discovery size. It assumes the contrast is fixed in "
            "advance, with no selection among candidates and no search over "
            "endpoints."),
    }


def how_to_read_any_required_n():
    return (
        "A required n from a simulation is a statement about the assumed "
        "generating process, not about biology. It is only as good as the "
        "effect size, the redundancy and the noise that were assumed, and the "
        "sensitivity table says which of those three decides the answer. "
        "Where the effect comes from a declared scenario rather than from a "
        "measured estimate, the resulting n is a scenario too and must be "
        "quoted with its assumption attached."
    )
