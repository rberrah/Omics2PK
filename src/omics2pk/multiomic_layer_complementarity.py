"""Fold-safe reconstruction of one masked molecular layer from the others.

Axis A only. This module knows nothing about drugs, doses, concentrations,
enzyme activities or pharmacokinetics, and nothing computed with it may be
reported as evidence about any of them. Its cohort (EATRIS-Plus v4) is blood
from a healthy reference set and carries no pharmacological endpoint at all.

The design is frozen in
``docs/OMICS2PK_AXIS_A_EATRIS_PLUS_LAYER_COMPLEMENTARITY_PRESPECIFICATION_20261110.md``
and this module's :data:`FROZEN` is read back against that document's freeze
table by the runner, which stops before loading any data if the two disagree.

Three properties the primitives here must have, because the prespecification
makes them decision rules rather than niceties:

* **Nothing crosses from held-out to training.** Layer means, layer SDs,
  component loadings, design means, design SDs, target means, target SDs and
  the ridge penalty are all estimated inside the training rows only.
* **The target layer is masked entirely.** It is never a predictor, never a
  standardisation constant and never a component basis.
* **Intra-family re-measurement is kept separable.** The target's same-family
  sibling is excluded from the primary arm and measured on its own, because a
  plasma proteome predicted from a second plasma-proteome assay is assay
  agreement, not inter-omic complementarity.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np

#: The freeze table of the prespecification, §13. The runner re-reads the
#: document and refuses to start if any value here differs from it.
FROZEN = {
    "seed": 20261110,
    "outer_folds": 5,
    "outer_repeats": 2,
    "inner_folds": 3,
    "n_bootstrap": 2000,
    "ridge_alphas": (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0),
    "components_per_layer": 5,
    "component_sensitivities": (2, 10),
    "permutation_replicates": 20,
    "parsimony_tolerance": 0.01,
    "primary_target": "proteins",
    "secondary_target": "acylcarnitines",
    "primary_arm": "A1_other_families",
    "primary_comparator": "A0_covariates",
    "complementarity_comparator": "A2_same_family_sibling",
    "primary_configuration": "C11_primary",
    "log_ratio_threshold": 1000.0,
    "standardised_clip": 4.0,
    "status_ceiling": "exploratory_methodological",
}

#: Twelve deposited layers, six biological families, an exact partition.
FAMILIES = {
    "transcriptome": ("mRNA",),
    "small_rna": ("miRNA-mature", "miRNA-qRT-PCR"),
    "methylome": ("EM-seq",),
    "proteome": ("proteins", "proteomics-ipp"),
    "lipidome": (
        "lipidomics-neg",
        "lipidomics-pos",
        "lipidomics-spline",
        "fattyacids",
    ),
    "metabolome_targeted": ("acylcarnitines", "aminoacids"),
}

#: Families whose molecules overlap chemically with another family's. A
#: positive result across an adjacent pair is chemical overlap, not evidence of
#: cross-layer regulation, so the secondary target's headline arm excludes it.
BIOCHEMICALLY_ADJACENT = {
    "metabolome_targeted": ("lipidome",),
    "lipidome": ("metabolome_targeted",),
}

CONFIGURATIONS = {
    "C11_primary": {
        "drop_layers": ("miRNA-qRT-PCR",),
        "n": 79,
        "individual_set_sha256_prefix": "265ffd52c9ce1f9f25c5e26c",
        "role": "primary",
    },
    "C12_all_layers": {
        "drop_layers": (),
        "n": 55,
        "individual_set_sha256_prefix": "b10acdb92dd6df92cba10585",
        "role": "secondary_what_the_qrtpcr_layer_adds",
    },
    "C9_no_rna": {
        "drop_layers": ("miRNA-qRT-PCR", "miRNA-mature", "mRNA"),
        "n": 114,
        "individual_set_sha256_prefix": "3ba65c29b0c988fd6d9e4a66",
        "role": "secondary_wider_coverage_without_the_transcriptome",
    },
}

#: The prespecified baseline block. Complete on all three individual sets, so
#: the primary contrast involves no imputation whatsoever.
COVARIATES_K0 = ("Age", "Sex", "BMI", "Height", "Weight")

#: Haematology. Incomplete, imputed inside training folds, and admitted only
#: to one declared secondary reference arm.
COVARIATES_H = (
    "Leukocytes",
    "Neutrophils",
    "Lymphocytes",
    "Hematocrit",
    "Erythrocytes",
    "Hemoglobin",
    "MCV",
    "MCH",
    "MCHC",
    "MXD",
    "Platelets",
)

PRIMARY_ARM = FROZEN["primary_arm"]
PRIMARY_COMPARATOR = FROZEN["primary_comparator"]
COMPLEMENTARITY_COMPARATOR = FROZEN["complementarity_comparator"]


class LayerComplementarityError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise LayerComplementarityError(message)


# --- family algebra -------------------------------------------------------
def family_of(layer):
    for family, layers in FAMILIES.items():
        if layer in layers:
            return family
    raise LayerComplementarityError(f"layer {layer!r} belongs to no declared family")


def assert_exact_partition(deposited_layers):
    """The declared grouping must cover the deposited layers exactly."""
    declared = sorted(x for layers in FAMILIES.values() for x in layers)
    _check(
        len(declared) == len(set(declared)),
        "the declared family grouping repeats a layer, so it is not a partition",
    )
    _check(
        declared == sorted(deposited_layers),
        "the declared family grouping is not an exact partition of the "
        f"deposited layers: declared {declared}, deposited {sorted(deposited_layers)}",
    )
    return True


def siblings_of(layer):
    """Layers in the target's own family, excluding the target itself."""
    return tuple(x for x in FAMILIES[family_of(layer)] if x != layer)


def families_present(layers_available):
    out = {}
    for family, layers in FAMILIES.items():
        kept = tuple(x for x in layers if x in layers_available)
        if kept:
            out[family] = kept
    return out


def arms_for(target, layers_available):
    """Return {arm name: tuple of predictor layers}, prespecification §7.

    ``A0_covariates`` and ``A4_covariates_haematology`` carry no layer, which
    is why their value is the empty tuple rather than a missing key. Every arm
    carries the ``K0`` covariate block implicitly; ``A4`` additionally carries
    the haematology block, which the runner adds by name.
    """
    _check(target in layers_available, f"target {target!r} is not in this configuration")
    available = tuple(x for x in layers_available if x != target)
    present = families_present(available)
    target_family = family_of(target)
    sibling = tuple(x for x in siblings_of(target) if x in available)
    other_family_layers = tuple(
        x for x in available if family_of(x) != target_family
    )

    arms = {
        "A0_covariates": (),
        "A1_other_families": other_family_layers,
        "A2_same_family_sibling": sibling,
        "A3_all_but_target": available,
        "A4_covariates_haematology": (),
    }
    non_target_families = tuple(
        family for family in present if family != target_family
    )
    for family in non_target_families:
        arms[f"B_{family}"] = tuple(
            x for x in other_family_layers if family_of(x) == family
        )
        arms[f"L_{family}"] = tuple(
            x for x in other_family_layers if family_of(x) != family
        )
    return arms


def family_subsets(target, layers_available):
    """All subsets of the non-target families, for the parsimony enumeration.

    The empty subset is included and is identical to ``A0_covariates``; it is
    named here so that the enumeration table has a row for "no omic layer".
    """
    available = tuple(x for x in layers_available if x != target)
    target_family = family_of(target)
    non_target = sorted(
        family for family in families_present(available) if family != target_family
    )
    out = {}
    for size in range(len(non_target) + 1):
        for subset in combinations(non_target, size):
            name = "S_" + ("none" if not subset else "+".join(subset))
            out[name] = (
                tuple(subset),
                tuple(x for x in available if family_of(x) in subset),
            )
    return out


def parsimonious_subset(table, tolerance=None):
    """Apply the prespecified parsimony rule, §10.

    ``table`` maps subset name -> {"families": tuple, "n_features": int,
    "median_delta_r2": float}. Returns the winner and the rule, so that a test
    can prove the choice was computed rather than written down.
    """
    tolerance = FROZEN["parsimony_tolerance"] if tolerance is None else tolerance
    _check(bool(table), "parsimony enumeration is empty")
    best = max(float(row["median_delta_r2"]) for row in table.values())
    floor = best - float(tolerance)
    eligible = [
        (len(row["families"]), int(row["n_features"]), tuple(sorted(row["families"])), name)
        for name, row in table.items()
        if float(row["median_delta_r2"]) >= floor
    ]
    _check(bool(eligible), "no subset reached the parsimony floor, which cannot happen")
    eligible.sort()
    n_families, n_features, families, name = eligible[0]
    return {
        "selected": name,
        "families": list(families),
        "n_families": n_families,
        "n_features": n_features,
        "best_median_delta_r2_in_enumeration": best,
        "floor": floor,
        "n_subsets_enumerated": len(table),
        "n_subsets_eligible": len(eligible),
        "rule": (
            "Among subsets whose median delta R-squared is at least the best "
            f"subset's minus {tolerance}, take the fewest families; ties are "
            "broken by fewest features, then lexicographic family order."
        ),
        "performance_of_the_selected_subset_is_optimistically_biased": True,
        "why": (
            "The subset was chosen on out-of-fold performance across "
            f"{len(table)} candidates, so its own number is not an unbiased "
            "estimate and is never reported as one."
        ),
    }


# --- layer preparation ---------------------------------------------------
def log_rule_applies(values, threshold=None):
    """§8.2.1 rule 1, evaluated on the FULL deposited layer.

    A property of the assay's units: strictly positive and spanning more than
    ``threshold``-fold. It uses no individual's relation to the target and is
    identical in every fold and every configuration.
    """
    threshold = FROZEN["log_ratio_threshold"] if threshold is None else threshold
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return False
    if finite.min() <= 0:
        return False
    return bool(finite.max() / finite.min() > threshold)


def complete_columns(values):
    """§8.2.1: columns with any missing cell are dropped, never imputed.

    Reads the pattern of missingness only, so it cannot carry a value or the
    target across the fold boundary.
    """
    array = np.asarray(values, dtype=float)
    return ~np.isnan(array).any(axis=0)


def layer_components(train, test, k, clip=None):
    """Training-fold principal directions of one layer, applied to held-out rows.

    Computed through the n x n Gram matrix, which is exact and stays cheap for
    a 100,000-feature methylome at n < 130. Returns the training scores, the
    held-out scores and the number of components actually available, which can
    fall below ``k`` when the training block is rank-deficient.
    """
    clip = FROZEN["standardised_clip"] if clip is None else clip
    train = np.asarray(train, dtype=float)
    test = np.asarray(test, dtype=float)
    _check(train.ndim == 2 and test.ndim == 2, "layer blocks must be two-dimensional")
    _check(
        train.shape[1] == test.shape[1],
        "training and held-out rows of a layer must share their columns",
    )
    _check(np.isfinite(train).all(), "a training layer block carries a missing cell")
    _check(np.isfinite(test).all(), "a held-out layer block carries a missing cell")
    _check(k >= 1, "at least one component must be requested")

    mean = train.mean(axis=0)
    sd = train.std(axis=0, ddof=0)
    keep = sd > 0.0
    if not keep.any():
        return (
            np.zeros((train.shape[0], 0)),
            np.zeros((test.shape[0], 0)),
            0,
        )
    centred = (train[:, keep] - mean[keep]) / sd[keep]
    held = (test[:, keep] - mean[keep]) / sd[keep]
    centred = np.clip(centred, -clip, clip)
    held = np.clip(held, -clip, clip)
    # Re-centre after clipping so that the Gram matrix is of a centred block.
    shift = centred.mean(axis=0)
    centred = centred - shift
    held = held - shift

    gram = centred @ centred.T
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    top = max(float(values[0]), 0.0)
    positive = values > max(1e-10 * top, 1e-12)
    available = int(positive.sum())
    taken = int(min(k, available))
    if taken == 0:
        return (
            np.zeros((train.shape[0], 0)),
            np.zeros((test.shape[0], 0)),
            0,
        )
    basis = vectors[:, :taken]
    singular = np.sqrt(values[:taken])
    train_scores = basis * singular
    test_scores = (held @ centred.T) @ basis / singular
    return train_scores, test_scores, taken


def standardise_design(train, test):
    """Standardise an assembled design with training constants only."""
    train = np.asarray(train, dtype=float)
    test = np.asarray(test, dtype=float)
    _check(np.isfinite(train).all() and np.isfinite(test).all(),
           "a design matrix carries a non-finite value")
    mean = train.mean(axis=0)
    sd = train.std(axis=0, ddof=0)
    sd = np.where(sd > 0.0, sd, 1.0)
    return (train - mean) / sd, (test - mean) / sd


# --- estimator ------------------------------------------------------------
def _ridge_paths(design, targets, held, alphas):
    """Held-out predictions for every penalty, in standardised target units."""
    left, singular, right = np.linalg.svd(design, full_matrices=False)
    projected = left.T @ targets
    rotated = held @ right.T
    out = np.empty((len(alphas), held.shape[0], targets.shape[1]))
    for index, alpha in enumerate(alphas):
        shrunk = singular / (singular ** 2 + alpha)
        out[index] = rotated @ (shrunk[:, None] * projected)
    return out


def ridge_per_target(design_train, targets_train, design_test, *,
                     alphas=None, inner_folds=None, seed=None):
    """Ridge with the penalty chosen per target feature by inner CV.

    Returns the held-out predictions on the original target scale, the chosen
    penalty per target feature, and the training means of the target features
    (which the prespecified R-squared uses as its reference).
    """
    from sklearn.model_selection import KFold

    alphas = tuple(FROZEN["ridge_alphas"]) if alphas is None else tuple(alphas)
    inner_folds = FROZEN["inner_folds"] if inner_folds is None else inner_folds
    seed = FROZEN["seed"] if seed is None else seed

    design_train = np.asarray(design_train, dtype=float)
    design_test = np.asarray(design_test, dtype=float)
    targets_train = np.asarray(targets_train, dtype=float)
    _check(design_train.shape[0] == targets_train.shape[0],
           "design and target rows disagree")
    _check(np.isfinite(targets_train).all(), "a training target carries a missing cell")

    n_train = design_train.shape[0]
    n_targets = targets_train.shape[1]
    target_mean = targets_train.mean(axis=0)
    target_sd = targets_train.std(axis=0, ddof=0)
    safe_sd = np.where(target_sd > 0.0, target_sd, 1.0)

    # --- inner CV, inside the training rows only -------------------------
    inner = KFold(n_splits=min(inner_folds, n_train), shuffle=True, random_state=seed)
    squared = np.zeros((len(alphas), n_targets))
    for inner_train, inner_val in inner.split(np.arange(n_train)):
        a_design, b_design = standardise_design(
            design_train[inner_train], design_train[inner_val]
        )
        a_targets = targets_train[inner_train]
        inner_mean = a_targets.mean(axis=0)
        inner_sd = a_targets.std(axis=0, ddof=0)
        inner_safe = np.where(inner_sd > 0.0, inner_sd, 1.0)
        standardised = (a_targets - inner_mean) / inner_safe
        paths = _ridge_paths(a_design, standardised, b_design, alphas)
        observed = targets_train[inner_val]
        for index in range(len(alphas)):
            predicted = paths[index] * inner_safe + inner_mean
            squared[index] += ((predicted - observed) ** 2).sum(axis=0)
    chosen = np.asarray(alphas, dtype=float)[np.argmin(squared, axis=0)]

    # --- refit on the whole training fold --------------------------------
    z_train, z_test = standardise_design(design_train, design_test)
    standardised = (targets_train - target_mean) / safe_sd
    paths = _ridge_paths(z_train, standardised, z_test, alphas)
    predictions = np.empty((design_test.shape[0], n_targets))
    for index, alpha in enumerate(alphas):
        columns = np.flatnonzero(chosen == alpha)
        if columns.size:
            predictions[:, columns] = (
                paths[index][:, columns] * safe_sd[columns] + target_mean[columns]
            )
    _check(np.isfinite(predictions).all(), "a held-out prediction is not finite")
    return predictions, chosen, target_mean


# --- scoring --------------------------------------------------------------
def r2_per_feature(squared_error, squared_total):
    """Out-of-fold R-squared per target feature, individual as the unit.

    ``squared_error`` and ``squared_total`` are (features, individuals): the
    per-individual squared prediction error and the per-individual squared
    deviation from that individual's training-fold mean, each already averaged
    over the repeats.
    """
    error = np.asarray(squared_error, dtype=float)
    total = np.asarray(squared_total, dtype=float)
    _check(error.shape == total.shape, "error and total blocks disagree in shape")
    denominator = total.sum(axis=1)
    safe = np.where(denominator > 0.0, denominator, np.nan)
    return 1.0 - error.sum(axis=1) / safe


def paired_bootstrap_median_delta(error_a, error_b, total, *,
                                  n_boot=None, seed=None):
    """Interval on the median across features of R2(a) - R2(b).

    Resamples **individuals**, never folds and never repeats, and recomputes
    both arms' R-squared from stored out-of-fold errors without refitting.
    """
    n_boot = FROZEN["n_bootstrap"] if n_boot is None else n_boot
    seed = FROZEN["seed"] if seed is None else seed
    error_a = np.asarray(error_a, dtype=float)
    error_b = np.asarray(error_b, dtype=float)
    total = np.asarray(total, dtype=float)
    _check(error_a.shape == error_b.shape == total.shape,
           "the two arms and the total block must share their shape")
    n_individuals = error_a.shape[1]

    point = float(np.nanmedian(r2_per_feature(error_a, total)
                               - r2_per_feature(error_b, total)))
    rng = np.random.default_rng(seed)
    picked = rng.integers(0, n_individuals, size=(n_boot, n_individuals))
    # One resample is a multiset of individuals, so its per-feature sums are
    # the error block times the individuals' multiplicities. Writing it as a
    # matrix product is the same estimator as summing the resampled columns,
    # and it keeps 2,000 draws over 1,393 features to one BLAS call per arm.
    weights = np.zeros((n_individuals, n_boot))
    for index in range(n_boot):
        weights[:, index] = np.bincount(picked[index], minlength=n_individuals)
    total_sums = total @ weights
    safe = np.where(total_sums > 0.0, total_sums, np.nan)
    draws = np.nanmedian(((error_b @ weights) - (error_a @ weights)) / safe, axis=0)
    low, high = (float(v) for v in np.percentile(draws, [2.5, 97.5]))
    return {
        "median_delta_r2": point,
        "descriptive_95_interval": [low, high],
        "interval_excludes_zero": bool(low * high > 0.0),
        "n_independent_individuals": int(n_individuals),
        "n_bootstrap": int(n_boot),
        "resampling_unit": "individual",
        "how_to_read": (
            "interval_excludes_zero false means no difference was DETECTED, "
            "never that there is none. A repeat of the cross-validation is "
            "not an additional individual and is never resampled."
        ),
    }


def permutation_rank(observed, null_values):
    """Rank p-value of the observed statistic against the permutation null."""
    null_values = [float(v) for v in null_values]
    at_least = sum(1 for v in null_values if v >= float(observed))
    return {
        "observed": float(observed),
        "n_replicates": len(null_values),
        "null_min": min(null_values) if null_values else None,
        "null_max": max(null_values) if null_values else None,
        "null_median": float(np.median(null_values)) if null_values else None,
        "n_null_at_least_observed": at_least,
        "rank_p_value": (at_least + 1) / (len(null_values) + 1) if null_values else None,
        "exceeds_every_replicate": bool(null_values) and at_least == 0,
    }


def no_pharmacological_claim():
    """The sentence every report built on this module must carry."""
    return (
        "EATRIS-Plus holds no drug, dose, concentration, time, enzyme activity "
        "or clinical endpoint. Nothing in this analysis is evidence about drug "
        "metabolism, enzyme activity, intrinsic clearance or pharmacokinetics, "
        "and the predictors are leukocyte and plasma measurements while the "
        "target is a plasma protein panel, so no hepatic interpretation is "
        "available either."
    )
