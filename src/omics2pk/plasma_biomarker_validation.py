"""Independent replication, and one genuinely new contrast, on the Lee / Jackson
CYP3A donor panel.

This module adds **no new predictive architecture**. Every fitting primitive --
standardisation, ridge, the inner leave-one-out penalty, the loss triple, the
permutation rank and the sample-size projection -- is imported from the frozen
``protein_clint_pilot``. What is new here is the arm set, and one interval the
Wegler pilot did not provide.

Three questions, in the order of what they are worth:

1. **Replication.** Does the target enzyme's own abundance predict its probe
   reaction out of sample, in donors and a laboratory independent of Wegler?
   That is ``V1`` against ``V0``, and it is the project's one positive
   functional claim put at risk.
2. **Incremental value, new.** Does a **plasma** metabolite measured in the
   same donors add to the **hepatic** protein? That is ``V3`` against ``V1``.
   No public dataset has let this project ask that question before, and the
   two predictors correlate only weakly with each other, so the answer is not
   arithmetically forced either way.
3. **Substitution, secondary.** Does the plasma metabolite alone beat the mean?
   That is ``V2`` against ``V0``, and it speaks to whether a non-invasive
   measurement could stand in for a liver biopsy.

The endpoint is an activity RATE at a single substrate concentration. It is not
an intrinsic clearance, not an AUC, not a clearance and not a concentration in
a patient. Nothing here is a pharmacokinetic result.

The one methodological addition
-------------------------------
The Wegler pilot's intervals resampled stored out-of-fold errors without
refitting, which prices the variability of the donors used for scoring and not
the variability of learning on fourteen donors. Here the model is small enough
to refit inside every draw, so ``out_of_bag_bootstrap`` does exactly that:
it draws a bootstrap sample of donors, refits each arm on it -- penalty
selection included -- and scores the donors left out of the draw. Both
intervals are reported, and the refitting one is the one that answers the
objection.
"""
from __future__ import annotations

import numpy as np

from . import protein_clint_pilot as PCP

#: Everything frozen before execution. Section 9 of the prespecification.
FROZEN = {
    "seed": 20261110,
    "validation": "leave_one_donor_out",
    "ridge_alphas": PCP.FROZEN["ridge_alphas"],
    "permutation_replicates": 200,
    "n_bootstrap": 2000,
    "n_out_of_bag_bootstrap": 2000,
    "primary_endpoint": "midazolam_1_hydroxylation",
    "secondary_endpoint": "ibrutinib_m37",
    "primary_arm": "V3_protein_plus_plasma",
    "primary_comparator": "V1_hepatic_protein",
    "replication_arm": "V1_hepatic_protein",
    "replication_comparator": "V0_mean",
    "substitution_arm": "V2_plasma_metabolite",
    "target_scale": "log10",
    "declared_family_size": 6,
    "status_ceiling": "exploratory_independent_replication",
}

#: The four arms, as predictor column lists. ``V0`` is deliberately empty: the
#: frozen ridge returns the training mean when handed no columns, so the
#: no-information baseline needs no special case.
ARMS = {
    "V0_mean": (),
    "V1_hepatic_protein": ("cyp3a4_abundance",),
    "V2_plasma_metabolite": ("plasma_4b_hc_cholesterol_ratio",),
    "V3_protein_plus_plasma": ("cyp3a4_abundance",
                               "plasma_4b_hc_cholesterol_ratio"),
}

#: A declared sensitivity, not a second primary: CYP3A5 contributes to midazolam
#: 1'-hydroxylation, so the sum is the mechanistically defensible alternative
#: to CYP3A4 alone.
SENSITIVITY_ARMS = {
    "V1s_cyp3a_total": ("cyp3a_total_abundance",),
    "V3s_cyp3a_total_plus_plasma": ("cyp3a_total_abundance",
                                    "plasma_4b_hc_cholesterol_ratio"),
}

#: The three contrasts, each as (arm, comparator, the block permuted under the
#: null). ``None`` means the null permutes the arm's own only predictor.
CONTRASTS = {
    "incremental_value_V3_vs_V1": {
        "arm": "V3_protein_plus_plasma",
        "comparator": "V1_hepatic_protein",
        "permuted_block": ("plasma_4b_hc_cholesterol_ratio",),
        "question": ("does a plasma metabolite add to the hepatic protein in "
                     "the same donors?"),
    },
    "replication_V1_vs_V0": {
        "arm": "V1_hepatic_protein",
        "comparator": "V0_mean",
        "permuted_block": ("cyp3a4_abundance",),
        "question": ("does the target enzyme's abundance predict its probe "
                     "reaction out of sample, independently of Wegler?"),
    },
    "substitution_V2_vs_V0": {
        "arm": "V2_plasma_metabolite",
        "comparator": "V0_mean",
        "permuted_block": ("plasma_4b_hc_cholesterol_ratio",),
        "question": ("does the plasma metabolite alone beat the "
                     "no-information baseline?"),
    },
}

#: Columns that must never be an endpoint, because the endpoint would then
#: contain its own predictor. This source has no CYP-normalised quantity, so
#: the list guards against a future one being added by hand.
FORBIDDEN_ENDPOINTS = (
    "cyp3a4_abundance",
    "cyp3a5_abundance",
    "cyp3a_total_abundance",
    "plasma_4b_hc_cholesterol_ratio",
)


class PlasmaBiomarkerValidationError(ValueError):
    """Raised when the design would not test what it claims to test."""


def _check(condition, message):
    if not condition:
        raise PlasmaBiomarkerValidationError(message)


def assert_endpoint_is_admissible(name):
    """An endpoint may not be one of its own predictors.

    Wegler's equations 8 and 9 divided the clearance by the probed CYP and
    scaled it back by the homogenate CYP concentration, which made a
    CYP-abundance predictor circular. This source has no such quantity, and
    this guard exists so that the absence stays checked rather than assumed.
    """
    _check(name not in FORBIDDEN_ENDPOINTS,
           f"{name!r} is a predictor in this design and must never also be "
           "the endpoint: the model would then be predicting its own input")
    return name


# --------------------------------------------------------------------------
# the design matrix
# --------------------------------------------------------------------------
def design(rows, donors, columns):
    """Rows x columns, in the donor order given. Missing cells are refused."""
    if not columns:
        return np.zeros((len(donors), 0))
    block = np.array([[rows[d][c] for c in columns] for d in donors],
                     dtype=float)
    _check(np.isfinite(block).all(),
           f"the design for {columns} has a missing or non-finite cell")
    return block


def log_target(rows, donors, endpoint):
    """log10 of the endpoint, which must be strictly positive."""
    assert_endpoint_is_admissible(endpoint)
    values = np.array([rows[d][endpoint] for d in donors], dtype=float)
    _check(np.isfinite(values).all(), f"{endpoint} has a missing value")
    _check((values > 0.0).all(),
           f"{endpoint} has a non-positive value, which log10 cannot take; "
           "a zero or a non-detect must be handled by a declared rule, never "
           "by silently dropping the donor")
    return np.log10(values)


# --------------------------------------------------------------------------
# leave-one-donor-out, with the penalty chosen inside the training donors
# --------------------------------------------------------------------------
def leave_one_donor_out(predictors, target, *, alphas=None):
    """Out-of-sample prediction for every donor, one at a time.

    Returns the predictions, the per-donor absolute log10 errors and the
    penalty chosen in each fold. The penalty is selected by an inner
    leave-one-out among the training donors only, so the held-out donor
    touches neither the standardisation constants nor the penalty.
    """
    alphas = FROZEN["ridge_alphas"] if alphas is None else alphas
    predictors = np.asarray(predictors, dtype=float)
    target = np.asarray(target, dtype=float)
    n = target.size
    _check(predictors.shape[0] == n, "predictors and target disagree on n")
    _check(n >= 4, f"leave-one-donor-out with an inner leave-one-out needs at "
                   f"least four donors, got {n}")
    predicted = np.zeros(n)
    penalties = np.zeros(n)
    for hold in range(n):
        keep = np.arange(n) != hold
        alpha = PCP.select_penalty(predictors[keep], target[keep], alphas)
        penalties[hold] = alpha
        predicted[hold] = PCP._ridge_fit_predict(
            predictors[keep], target[keep],
            predictors[hold:hold + 1], alpha)[0]
    return predicted, np.abs(predicted - target), penalties


def score_arm(rows, donors, endpoint, columns, *, alphas=None):
    """One arm, scored out of sample, in log and in native units."""
    target = log_target(rows, donors, endpoint)
    predicted, absolute_log_errors, penalties = leave_one_donor_out(
        design(rows, donors, columns), target, alphas=alphas)
    native = np.abs(10.0 ** predicted - 10.0 ** target)
    out = PCP.losses(absolute_log_errors, native)
    out["penalty_by_fold"] = {d: float(p) for d, p in zip(donors, penalties)}
    out["n_predictors"] = len(columns)
    out["predictors"] = list(columns)
    out["native_units"] = "pmol per million cells per minute"
    return out, absolute_log_errors


# --------------------------------------------------------------------------
# the interval that refits, which is the point of this module
# --------------------------------------------------------------------------
def out_of_bag_bootstrap(rows, donors, endpoint, arm_columns,
                         comparator_columns, *, n_boot=None, seed=None,
                         alpha=0.05, family_size=1):
    """A bootstrap that re-learns both arms inside every draw.

    Each draw samples donors with replacement to form a training set, refits
    both arms on it -- including the inner leave-one-out that picks the
    penalty -- and scores the donors the draw left out. The paired difference
    is taken on those out-of-bag donors only, so no donor is ever both
    learned from and scored.

    This prices the variability of LEARNING as well as of scoring, which the
    stored-error bootstrap cannot. The cost is that each draw scores about a
    third of the cohort, so the interval is wide; that width is the honest
    answer at fifteen donors, not a defect of the estimator.
    """
    n_boot = FROZEN["n_out_of_bag_bootstrap"] if n_boot is None else n_boot
    seed = FROZEN["seed"] if seed is None else seed
    target = log_target(rows, donors, endpoint)
    arm = design(rows, donors, arm_columns)
    comparator = design(rows, donors, comparator_columns)
    n = len(donors)
    rng = np.random.default_rng(seed)
    differences = []
    skipped_too_small = 0
    for _ in range(n_boot):
        drawn = rng.integers(0, n, size=n)
        held = np.array(sorted(set(range(n)) - set(drawn.tolist())), dtype=int)
        # Three training donors is the floor the frozen penalty selector
        # enforces, and a single out-of-bag donor gives no paired difference
        # worth the name.
        if held.size < 2 or (n - held.size) < 3:
            skipped_too_small += 1
            continue
        kept = drawn
        pair = []
        for columns, block in ((arm_columns, arm),
                               (comparator_columns, comparator)):
            if columns:
                chosen = PCP.select_penalty(block[kept], target[kept],
                                            FROZEN["ridge_alphas"])
                predicted = PCP._ridge_fit_predict(
                    block[kept], target[kept], block[held], chosen)
            else:
                predicted = np.full(held.size, float(target[kept].mean()))
            pair.append(np.abs(predicted - target[held]).mean())
        differences.append(pair[0] - pair[1])
    _check(len(differences) >= 100,
           f"only {len(differences)} of {n_boot} draws left enough donors out "
           "to score; the interval would not mean anything")
    draws = np.asarray(differences, dtype=float)
    half = 100.0 * alpha / (2.0 * family_size)
    low, high = (float(v) for v in np.percentile(draws, [half, 100.0 - half]))
    return {
        "mean_difference_log10": float(draws.mean()),
        "median_difference_log10": float(np.median(draws)),
        "interval": [low, high],
        "interval_excludes_zero": bool(low * high > 0.0),
        "level": f"{100.0 * (1.0 - alpha / family_size):.2f} %",
        "family_size": int(family_size),
        "n_draws_used": int(draws.size),
        "n_draws_skipped_too_small": int(skipped_too_small),
        "refits_the_model_inside_each_draw": True,
        "scoring_donors_per_draw": "the out-of-bag donors, about a third of 15",
        "what_this_interval_prices": (
            "both the choice of donors learned from and the choice of donors "
            "scored on, because each draw re-estimates the standardisation "
            "constants, the ridge penalty and the coefficients. It is the "
            "answer to the objection that the stored-error interval is "
            "narrower than a nested one."
        ),
        "what_it_still_does_not_price": (
            "the choice of this cohort, this endpoint and these two "
            "predictors. No resampling of fifteen donors can price the fact "
            "that they are fifteen purchased hepatocyte lots of ages 1 to 63."
        ),
        "exploratory": True,
    }


def contrast(rows, donors, endpoint, name, *, family_size=1,
             n_boot=None, n_oob=None, seed=None):
    """One contrast, with both intervals and the arm scores."""
    _check(name in CONTRASTS, f"{name!r} is not a frozen contrast")
    spec = CONTRASTS[name]
    arm_columns = ARMS[spec["arm"]]
    comparator_columns = ARMS[spec["comparator"]]
    arm_score, arm_errors = score_arm(rows, donors, endpoint, arm_columns)
    base_score, base_errors = score_arm(rows, donors, endpoint,
                                        comparator_columns)
    stored = PCP.paired_bootstrap(arm_errors, base_errors, n_boot=n_boot,
                                  seed=seed, family_size=family_size)
    refitting = out_of_bag_bootstrap(rows, donors, endpoint, arm_columns,
                                     comparator_columns, n_boot=n_oob,
                                     seed=seed, family_size=family_size)
    return {
        "question": spec["question"],
        "arm": spec["arm"],
        "comparator": spec["comparator"],
        "arm_score": arm_score,
        "comparator_score": base_score,
        "stored_error_bootstrap": stored,
        "refitting_out_of_bag_bootstrap": refitting,
        "the_two_intervals_disagree": bool(
            stored["interval_excludes_zero"]
            != refitting["interval_excludes_zero"]),
        "which_interval_governs": (
            "the refitting one. Where they disagree, the stored-error "
            "interval is the optimistic reading and must not be quoted alone."
        ),
    }, arm_errors, base_errors


# --------------------------------------------------------------------------
# the null
# --------------------------------------------------------------------------
def permutation_null(rows, donors, endpoint, name, *, replicates=None,
                     seed=None):
    """Permute the ADDED block's donor labels and refit everything.

    The comparator's predictors and the endpoint stay paired, so the null asks
    exactly one question: how much does this contrast improve when the added
    information is destroyed but its dimensionality is kept? That difference
    is not zero in expectation -- a wider arm shrinks harder under an
    inner-CV penalty and so scores better with no information at all -- which
    is why every contrast reported here carries its own null and not only the
    ones that look significant.
    """
    _check(name in CONTRASTS, f"{name!r} is not a frozen contrast")
    replicates = (FROZEN["permutation_replicates"] if replicates is None
                  else replicates)
    seed = FROZEN["seed"] if seed is None else seed
    spec = CONTRASTS[name]
    permuted = set(spec["permuted_block"])
    arm_columns = ARMS[spec["arm"]]
    comparator_columns = ARMS[spec["comparator"]]
    rng = np.random.default_rng(seed + 1)
    values = []
    for _ in range(replicates):
        order = rng.permutation(len(donors))
        shuffled = {}
        for index, donor in enumerate(donors):
            record = dict(rows[donor])
            for column in permuted:
                record[column] = rows[donors[order[index]]][column]
            shuffled[donor] = record
        _, arm_errors = score_arm(shuffled, donors, endpoint, arm_columns)
        _, base_errors = score_arm(shuffled, donors, endpoint,
                                   comparator_columns)
        values.append(float((arm_errors - base_errors).mean()))
    return values


def deciles(values):
    """Ten numbers instead of two hundred, so no artefact carries a long list."""
    array = np.asarray(sorted(float(v) for v in values), dtype=float)
    return {f"p{int(q)}": float(np.percentile(array, q))
            for q in range(0, 101, 10)}


# --------------------------------------------------------------------------
# the statements this module refuses to let a reader infer
# --------------------------------------------------------------------------
def what_this_cannot_claim():
    return {
        "not_an_intrinsic_clearance": (
            "the endpoint is a formation rate at one substrate concentration. "
            "It is not Vmax/Km, not a fitted CLint, and not comparable to "
            "Wegler's equation 7 on the same scale."),
        "not_a_pharmacokinetic_result": (
            "no dose, no plasma drug concentration, no exposure. A gain here "
            "says nothing about predicting a concentration, an AUC or a "
            "clearance in a patient."),
        "not_a_clinical_biomarker_validation": (
            "the plasma came from deceased organ donors at procurement. A "
            "4-beta-hydroxycholesterol ratio measured there is not "
            "interchangeable with one measured in a living patient."),
        "two_substrates_are_not_two_cohorts": (
            "midazolam and ibrutinib were measured in the same fifteen "
            "donors. The second endpoint is a second reading of one cohort."),
        "absence_of_a_gain_is_not_equivalence": (
            "at fifteen donors an interval that contains zero means no "
            "difference was detected. The estimate and its width are the "
            "result; a verdict of 'no difference' is not available."),
    }
