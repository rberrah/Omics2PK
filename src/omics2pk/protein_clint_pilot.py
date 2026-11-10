"""Does the microsomal proteome add to the target enzyme's own abundance?

Axis B, Wegler 2021, 15 donors, leave-one-donor-out. The design is frozen in
``docs/OMICS2PK_AXIS_B_WEGLER_PROTEIN_CLINT_PILOT_PRESPECIFICATION_20261110.md``
and this module's :data:`FROZEN` is read back against that document by the
runner, which stops before loading any data if the two disagree.

Three properties the primitives here must have, because each one is a way the
experiment could be silently wrong rather than merely weak:

* **The endpoint must not contain the predictor.** Two of the article's three
  published clearance equations scale by the specific CYP amount, so they are
  circular against a CYP-abundance predictor. Only ``eq7`` is admissible and
  :func:`assert_endpoint_is_admissible` refuses the others by name.
* **The two halves of the source number their donors differently.** Table S3
  uses ``UU1``-``UU15`` and Data S1 uses ``1``-``15``; they agree for donors
  1-9 and are permuted for 10-15. :func:`resolve_donor_mapping` recovers the
  bijection from the CYP concentrations both files publish and refuses
  anything that is not unanimous, rather than trusting a stored constant.
* **Nothing crosses from the held-out donor.** Column means, SDs,
  principal-component loadings, the target mean and the ridge penalty are all
  estimated inside the training donors, including inside the inner
  leave-one-out that selects the penalty.

Nothing computed here is a clinical-pharmacokinetic result. The endpoint is an
in-vitro intrinsic clearance.
"""
from __future__ import annotations

import math
import re
from statistics import NormalDist

import numpy as np

#: The freeze table of the prespecification, section 9.
FROZEN = {
    "seed": 20261110,
    "validation": "leave_one_donor_out",
    "ridge_alphas": (0.001, 0.01, 0.1, 1.0, 10.0, 100.0),
    "components": 2,
    "component_sensitivity": 1,
    "permutation_replicates": 200,
    "n_bootstrap": 2000,
    "primary_endpoint": "clint_mic_eq7",
    "paired_sensitivity_endpoint": "clint_hep_eq7",
    "forbidden_endpoints": ("clint_mic_eq8", "clint_mic_eq9", "clint_hep_eq8"),
    "primary_arm": "W2_target_cyp_plus_proteome",
    "primary_comparator": "W1_target_cyp",
    "secondary_comparator": "W0_mean",
    "n_substrates": 5,
    "target_scale": "log10",
    "status_ceiling": "exploratory",
}

ARMS = ("W0_mean", "W1_target_cyp", "W2_target_cyp_plus_proteome")
PRIMARY_ARM = FROZEN["primary_arm"]
PRIMARY_COMPARATOR = FROZEN["primary_comparator"]
SECONDARY_COMPARATOR = FROZEN["secondary_comparator"]

#: Substrate to the CYP whose abundance is its baseline predictor.
SUBSTRATE_CYP = {
    "Midazolam": "CYP3A4",
    "Omeprazole": "CYP2C19",
    "Diclofenac": "CYP2C9",
    "Bupropion": "CYP2B6",
    "Bufuralol": "CYP2D6",
}

#: Table S3 columns carrying the CYP concentration, by system.
BASELINE_COLUMN = {
    "HLM": "cyp_pmol_per_mg_HLM",
    "HH": "cyp_pmol_per_mg_HH",
    "HL": "cyp_pmol_per_mg_HL",
}

#: Endpoint column by system, and the Data S1 sheet each pairs with.
ENDPOINT_COLUMN = {"HLM": "clint_mic_eq7", "HH": "clint_hep_eq7"}
PROTEOME_SHEET = {"HLM": "Microsomes", "HH": "Hepatocytes"}

#: Cytochromes and their redox partners, for the declared sensitivity that
#: asks whether W2's addition is really other cytochromes.
CYTOCHROME_PATTERN = re.compile(
    r"^(CYP[0-9]+[A-Z]*[0-9]*(P[0-9]+)?|POR|CYB5A|CYB5R[0-9]*|MTRR)$", re.IGNORECASE
)

#: Columns that are the baseline predictor in another unit and must never be
#: an endpoint or an extra predictor.
NEVER_A_PREDICTOR_OR_ENDPOINT = (
    "cyp_pmol_in_incubation_HLM",
    "cyp_pmol_in_incubation_HH",
)


class ProteinClintPilotError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise ProteinClintPilotError(message)


# --- the circularity guard ------------------------------------------------
def assert_endpoint_is_admissible(name):
    """Refuse an endpoint that is scaled by the CYP abundance.

    The article's Figure S11 caption states that ``eq. 8`` scales clearance by
    the specific CYP amount and the homogenate CYP concentration, and that
    ``eq. 9`` adds Kpuu on top of that scaling. Both therefore contain the
    predictor, and predicting them from the CYP abundance is circular. Only
    ``eq. 7`` -- scaled by hepatocellularity or literature MPPGL -- is free of
    it.
    """
    _check(
        isinstance(name, str) and name,
        "an endpoint column name is required",
    )
    _check(
        name not in FROZEN["forbidden_endpoints"],
        f"endpoint {name!r} is scaled by the CYP abundance that serves as the "
        "predictor, so predicting it from that abundance would be circular. "
        f"Forbidden: {list(FROZEN['forbidden_endpoints'])}",
    )
    _check(
        name not in NEVER_A_PREDICTOR_OR_ENDPOINT,
        f"{name!r} is the CYP concentration multiplied by the incubation "
        "protein, so it is the baseline predictor in another unit",
    )
    _check(
        name.endswith("eq7"),
        f"endpoint {name!r} is not one of the admissible eq. 7 columns "
        f"({sorted(ENDPOINT_COLUMN.values())})",
    )
    return True


# --- the donor mapping ----------------------------------------------------
def resolve_donor_mapping(table_s3, data_s1, *, tolerance=0.051):
    """Recover the Data S1 index -> Table S3 batch bijection from the data.

    ``table_s3`` maps ``(gene, batch_label, system)`` to a concentration or
    None; ``data_s1`` maps ``(gene, index, system)`` the same way. Every
    (gene, system) pair that both files carry votes, and the result is accepted
    only if every index resolves to one batch unanimously and the whole is a
    bijection. A stored constant would be trusted; this is measured.
    """
    genes = sorted({gene for gene, _, _ in table_s3} & {gene for gene, _, _ in data_s1})
    systems = sorted({system for _, _, system in table_s3}
                     & {system for _, _, system in data_s1})
    _check(genes, "no gene is carried by both files, so no mapping can be proven")
    indices = sorted({index for _, index, _ in data_s1})
    batches = sorted({batch for _, batch, _ in table_s3})
    _check(
        len(indices) == len(batches),
        f"{len(indices)} Data S1 donors against {len(batches)} Table S3 "
        "batches: the files do not describe the same donor set",
    )

    vectors = 0
    votes = {index: {} for index in indices}
    for gene in genes:
        for system in systems:
            column = [(index, data_s1.get((gene, index, system))) for index in indices]
            reference = {
                batch: table_s3.get((gene, batch, system)) for batch in batches
            }
            if all(value is None for _, value in column) or all(
                value is None for value in reference.values()
            ):
                continue
            vectors += 1
            for index, value in column:
                if value is None:
                    continue
                for batch, target in reference.items():
                    if target is None:
                        continue
                    if abs(float(target) - float(value)) <= tolerance:
                        votes[index][batch] = votes[index].get(batch, 0) + 1

    _check(vectors >= 2, f"only {vectors} comparable vectors: too few to prove a mapping")

    mapping = {}
    detail = {}
    for index in indices:
        ranked = sorted(votes[index].items(), key=lambda kv: -kv[1])
        detail[index] = ranked[:3]
        _check(
            ranked, f"Data S1 donor {index} matches no Table S3 batch at all"
        )
        best, best_votes = ranked[0]
        runner_up = ranked[1][1] if len(ranked) > 1 else 0
        _check(
            best_votes == vectors,
            f"Data S1 donor {index} matches {best} on only {best_votes} of "
            f"{vectors} vectors: the mapping is not unanimous and will not be used",
        )
        _check(
            runner_up < best_votes,
            f"Data S1 donor {index} matches more than one batch equally well",
        )
        mapping[index] = best

    _check(
        len(set(mapping.values())) == len(indices),
        "the resolved mapping is not a bijection",
    )
    return {
        "mapping": mapping,
        "n_vectors_voting": vectors,
        "unanimous": True,
        "is_identity": all(
            str(index) == str(batch).lstrip("Uu") for index, batch in mapping.items()
        ),
        "votes": detail,
    }


# --- fold-safe pieces -----------------------------------------------------
def standardise(train, test):
    train = np.asarray(train, dtype=float)
    test = np.asarray(test, dtype=float)
    _check(train.ndim == test.ndim == 2, "blocks must be two-dimensional")
    _check(train.shape[1] == test.shape[1], "blocks must share their columns")
    mean = train.mean(axis=0)
    sd = train.std(axis=0, ddof=0)
    safe = np.where(sd > 0.0, sd, 1.0)
    return (train - mean) / safe, (test - mean) / safe


def proteome_components(train, test, k):
    """Training-fold principal directions of the log proteome block.

    Computed through the n x n Gram matrix, which is exact and cheap for a
    4,000-protein block at n < 16. Returns the training scores, the held-out
    scores and how many components were actually available.
    """
    _check(k >= 0, "a non-negative number of components is required")
    train = np.asarray(train, dtype=float)
    test = np.asarray(test, dtype=float)
    _check(np.isfinite(train).all(), "the training proteome block has a missing cell")
    _check(np.isfinite(test).all(), "the held-out proteome row has a missing cell")
    if k == 0 or train.size == 0:
        return np.zeros((train.shape[0], 0)), np.zeros((test.shape[0], 0)), 0

    centred, held = standardise(train, test)
    shift = centred.mean(axis=0)
    centred = centred - shift
    held = held - shift

    gram = centred @ centred.T
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    top = max(float(values[0]), 0.0)
    available = int((values > max(1e-10 * top, 1e-12)).sum())
    taken = int(min(k, available))
    if taken == 0:
        return np.zeros((train.shape[0], 0)), np.zeros((test.shape[0], 0)), 0
    basis = vectors[:, :taken]
    singular = np.sqrt(values[:taken])
    return basis * singular, (held @ centred.T) @ basis / singular, taken


def _ridge_fit_predict(design_train, target_train, design_test, alpha):
    """Ridge with an intercept, on standardised predictors. No library needed."""
    design_train = np.asarray(design_train, dtype=float)
    design_test = np.asarray(design_test, dtype=float)
    target_train = np.asarray(target_train, dtype=float)
    centre = float(target_train.mean())
    if design_train.shape[1] == 0:
        return np.full(design_test.shape[0], centre)
    z_train, z_test = standardise(design_train, design_test)
    centred = target_train - centre
    gram = z_train.T @ z_train + alpha * np.eye(z_train.shape[1])
    coefficients = np.linalg.solve(gram, z_train.T @ centred)
    return z_test @ coefficients + centre


def select_penalty(design, target, alphas):
    """Inner leave-one-out inside the training donors only."""
    design = np.asarray(design, dtype=float)
    target = np.asarray(target, dtype=float)
    n = design.shape[0]
    _check(n >= 3, f"need at least three training donors, got {n}")
    if design.shape[1] == 0:
        return float(alphas[0])
    errors = np.zeros(len(alphas))
    for hold in range(n):
        keep = np.arange(n) != hold
        for index, alpha in enumerate(alphas):
            predicted = _ridge_fit_predict(
                design[keep], target[keep], design[hold:hold + 1], alpha
            )
            errors[index] += (predicted[0] - target[hold]) ** 2
    return float(np.asarray(alphas, dtype=float)[int(np.argmin(errors))])


# --- losses ---------------------------------------------------------------
def losses(absolute_log_errors, native_absolute_errors):
    """Three readings of the same error, as the freeze requires."""
    absolute_log_errors = np.asarray(absolute_log_errors, dtype=float)
    native = np.asarray(native_absolute_errors, dtype=float)
    mae_log = float(absolute_log_errors.mean())
    return {
        "n_donors": int(absolute_log_errors.size),
        "mae_log10": mae_log,
        "geometric_fold_error": float(10.0 ** mae_log),
        "mae_ml_per_min_per_kg": float(native.mean()),
        "how_to_read": (
            "mae_log10 is the primary loss. geometric_fold_error is the same "
            "quantity read as 'off by a factor of x on average'. "
            "mae_ml_per_min_per_kg is the back-transformed error in the "
            "endpoint's own units. The last two carry no separate test."
        ),
    }


def paired_bootstrap(errors_new, errors_baseline, *, n_boot=None, seed=None,
                     alpha=0.05, family_size=1):
    """Paired difference over DONORS, with an optional Bonferroni level.

    A fold is a donor here, so the resampling unit and the fold unit coincide.
    """
    n_boot = FROZEN["n_bootstrap"] if n_boot is None else n_boot
    seed = FROZEN["seed"] if seed is None else seed
    new = np.asarray(errors_new, dtype=float)
    baseline = np.asarray(errors_baseline, dtype=float)
    _check(new.shape == baseline.shape, "the two arms must cover the same donors")
    difference = new - baseline
    point = float(difference.mean())
    rng = np.random.default_rng(seed)
    n = difference.size
    draws = difference[rng.integers(0, n, size=(n_boot, n))].mean(axis=1)
    half = 100.0 * alpha / (2.0 * family_size)
    low, high = (float(v) for v in np.percentile(draws, [half, 100.0 - half]))
    reference = float(baseline.mean())
    return {
        "n_independent_donors": int(n),
        "mean_difference_log10": point,
        "relative_change_vs_baseline": point / reference if reference else None,
        "level": f"{100.0 * (1.0 - alpha / family_size):.2f} %",
        "family_size": int(family_size),
        "interval": [low, high],
        "interval_excludes_zero": bool(low * high > 0.0),
        "direction": "improvement" if point < 0 else "degradation",
        "n_bootstrap": int(n_boot),
        "resampling_unit": "donor",
        "exploratory": True,
        "how_to_read": (
            "An interval containing zero means no difference was DETECTED. It "
            "is never a demonstration that there is none, and at n=15 the "
            "design's power is low by construction. The estimate and its "
            "uncertainty are the deliverable, not a significance verdict."
        ),
    }


def permutation_rank(observed, null_values):
    """Rank p-value. The observed difference is an IMPROVEMENT when negative,
    so the null tail of interest is the lower one."""
    null_values = [float(v) for v in null_values]
    at_most = sum(1 for value in null_values if value <= float(observed))
    return {
        "observed": float(observed),
        "n_replicates": len(null_values),
        "null_min": min(null_values) if null_values else None,
        "null_median": float(np.median(null_values)) if null_values else None,
        "null_max": max(null_values) if null_values else None,
        "n_null_at_least_as_good_as_observed": at_most,
        "rank_p_value": (at_most + 1) / (len(null_values) + 1) if null_values else None,
        "better_than_every_replicate": bool(null_values) and at_most == 0,
        "how_to_read": (
            "The null permutes the donor labels of the proteome block only, "
            "keeping the target CYP arm and the endpoint paired, and refits "
            "everything. It will not be centred on zero: the wider arm shrinks "
            "harder under an inner-CV penalty and so scores better with no "
            "information. An improvement inside the null is an artefact."
        ),
    }


# --- the sample-size projection ------------------------------------------
def sample_size_projection(paired_differences, *, effect_scenarios,
                           powers=(0.80, 0.90), alpha=0.05,
                           family_sizes=(1, 5)):
    """Donors required in a future cohort, per scenario, assumptions stated.

    A paired two-sided test on the per-donor difference in log10 absolute
    error. The dispersion is the observed standard deviation of those
    differences, which is itself estimated from the pilot's donors and is
    therefore imprecise; that is reported, not hidden.
    """
    differences = np.asarray(paired_differences, dtype=float)
    _check(differences.size >= 3, "need at least three paired differences")
    spread = float(differences.std(ddof=1))
    normal = NormalDist()
    rows = {}
    for label, effect in effect_scenarios.items():
        for power in powers:
            for family in family_sizes:
                if not effect or spread <= 0:
                    required = None
                else:
                    z_alpha = normal.inv_cdf(1.0 - alpha / (2.0 * family))
                    z_beta = normal.inv_cdf(power)
                    required = int(
                        math.ceil(((z_alpha + z_beta) * spread / abs(effect)) ** 2)
                    )
                rows[f"{label}|power={power:.2f}|family={family}"] = {
                    "assumed_true_mean_difference_log10": effect,
                    "power": power,
                    "alpha": alpha,
                    "family_size_for_bonferroni": family,
                    "donors_required": required,
                }
    return {
        "observed_sd_of_paired_differences_log10": spread,
        "n_donors_the_sd_was_estimated_from": int(differences.size),
        "test_assumed": "paired two-sided test on the per-donor difference in log10 absolute error",
        "scenarios": rows,
        "assumptions": [
            "The dispersion of the paired differences in a future cohort equals "
            "the one observed here, which is estimated from a pilot and is "
            "itself imprecise.",
            "The future cohort uses the same assay, the same preparation, the "
            "same scaling equation and the same arms.",
            "Donors are independent and the effect is constant across them.",
            "The projection is for ONE substrate; the family_size column "
            "carries the Bonferroni correction for five.",
        ],
        "how_to_read": (
            "An order of magnitude with its assumptions attached, never a "
            "single number. A scenario whose assumed effect is near zero "
            "requires an unbounded sample size, which is reported as null "
            "rather than as a large integer."
        ),
    }


def spearman_square_is_not_a_variance_fraction():
    """The sentence this project owes itself after withdrawing one.

    Kept as a function so a test can assert the prohibition is stated in the
    code and not only in a document.
    """
    return (
        "No percentage of explained variance may be derived from a Spearman "
        "coefficient. The square of a rank correlation is not a variance "
        "fraction, because the ranks are not the quantity whose variance is "
        "being apportioned."
    )


def no_clinical_claim():
    return (
        "The endpoint is an in-vitro intrinsic clearance in 15 liver donors. "
        "There is no patient, no dose and no exposure in this source, so no "
        "result here is a clinical-pharmacokinetic proof, and none may be "
        "read as one. The five substrates are not five independent cohorts, "
        "and the microsomal and hepatocyte readings of the same donors are "
        "not two cohorts."
    )
