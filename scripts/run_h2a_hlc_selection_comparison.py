#!/usr/bin/env python
"""H2a: does a mechanistic ordering select better than a statistical one?

Frozen by protocol version 1.4.1, section 5. The endpoint is measured
microsomal cytochrome activity, which is `functional_in_vitro` and is therefore
the functional analogue of H2, never pooled with a clinical pharmacokinetic
result and never reported in the same effect column.

Every arm sees the same donors, the same folds, the same predictor and the same
number of variables. The only difference between arms is where the ordering of
candidate genes comes from. Selection, preprocessing, imputation and tuning all
happen inside training folds; no arm sees a test fold.

The source archive is read in memory and never copied. Only aggregate
quantities are written out: no donor identifier, no per-donor value, no
expression row.

Run from the repository root, with the archive directory passed in:

    python scripts/run_h2a_hlc_selection_comparison.py \
        --hlc "<directory with curatedExpressionLiver.zip and the phenotype zip>" \
        --knowledge "<directory with the four knowledge sources>" \
        --out reports/data/h2a_hlc_selection.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import date
from pathlib import Path
from zipfile import ZipFile

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from omics2pk.mechanistic_graph_selection import (  # noqa: E402
    FROZEN as SELECTOR_FROZEN,
    load_string_symbol_graph,
    rank_by_mechanistic_proximity,
    seeds_from_ctd_substrate,
    seeds_from_endpoint_label,
    seeds_from_reactome,
)

# --------------------------------------------------------------------------
# Frozen analysis parameters (protocol 1.4.1 section 5)
# --------------------------------------------------------------------------
SEED = 20261108
OUTER_FOLDS = 5
OUTER_REPEATS = 2
INNER_FOLDS = 3
PANEL_SIZES = (1, 2, 3, 5, 10, 20, 50, 100)
RANDOM_REPLICATIONS = 100
N_DECOYS = 200
RIDGE_ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
LASSO_PATH_LIMIT = 150
MAX_PROBE_MISSING_FRACTION = 0.20
TRIM_FRACTION = 0.10
N_BOOTSTRAP = 2000
NON_INFERIORITY_MARGIN = 0.05

ENDPOINTS = ("CYP1A2_activity", "CYP2A6_activity", "CYP2B6_activity",
             "CYP2C8_activity", "CYP2C9_activity", "CYP2C19_activity",
             "CYP2D6_activity", "CYP2E1_activity", "CYP3A4M_activity",
             "CYP3A4T_activity")

# Fixed before running: each endpoint's deliberately wrong mechanism is the
# next endpoint's enzyme. The rotation depends on no result.
WRONG_SEED_ROTATION = {ENDPOINTS[i]: ENDPOINTS[(i + 1) % len(ENDPOINTS)]
                       for i in range(len(ENDPOINTS))}

# The only two endpoints whose probe substrate the source documentation names.
DOCUMENTED_SUBSTRATE = {"CYP3A4M_activity": "midazolam",
                        "CYP3A4T_activity": "testosterone"}

REACTOME_PATHWAY_PATTERN = r"^Phase I - Functionalization of compounds$"

EXPRESSION_MEMBER = "expression.txt"
FEATURES_MEMBER = "features.txt"
PHENOTYPE_MEMBER = "phenotype.txt"
EXPECTED_PROBES = 40638
EXPECTED_EXPRESSION_DONORS = 466
AGE_ROW = "AGE_(YRS)"
SEX_ROW = "GENDER"
SEX_CODES = {"Male": 1.0, "Female": 0.0}


class H2aError(RuntimeError):
    pass


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_hlc(hlc_dir: Path):
    """Read the archives in memory and verify the dimensions the loader assumes."""
    import pandas as pd

    expression_zip = hlc_dir / "curatedExpressionLiver.zip"
    phenotype_zip = hlc_dir / "curatedPhenotypesynapse.zip"
    for path in (expression_zip, phenotype_zip):
        if not path.is_file():
            raise H2aError(f"expected archive absent: {path}")

    with ZipFile(expression_zip) as archive:
        with archive.open(FEATURES_MEMBER) as handle:
            features = pd.read_csv(handle, sep="\t", dtype=str)
        with archive.open(EXPRESSION_MEMBER) as handle:
            expression = pd.read_csv(handle, sep="\t", index_col=0)
    with ZipFile(phenotype_zip) as archive:
        with archive.open(PHENOTYPE_MEMBER) as handle:
            phenotype = pd.read_csv(handle, sep="\t", index_col=0, dtype=str)

    if expression.shape != (EXPECTED_PROBES, EXPECTED_EXPRESSION_DONORS):
        raise H2aError(f"expression dimensions changed: {expression.shape}")
    if len(features) != EXPECTED_PROBES:
        raise H2aError(f"annotation rows changed: {len(features)}")
    if not np.array_equal(expression.index.to_numpy(dtype=str),
                          features["feature_id"].to_numpy(dtype=str)):
        raise H2aError("probe identifiers no longer align with the annotation")

    donors = expression.columns.astype(str).to_numpy()
    if len(set(donors)) != EXPECTED_EXPRESSION_DONORS:
        raise H2aError("duplicate donor identifiers in the expression table")
    missing = sorted(set(donors) - set(phenotype.columns.astype(str)))
    if missing:
        raise H2aError(f"{len(missing)} expression donors absent from the phenotype table")

    phenotype = phenotype.loc[:, donors]
    age = pd.to_numeric(phenotype.loc[AGE_ROW], errors="coerce").to_numpy(float)
    sex = phenotype.loc[SEX_ROW].map(SEX_CODES).to_numpy(float)
    # The leading digit of the donor identifier separates the cohort into three
    # groups. The archive readme states three independent tissue collection
    # centres. This is an INFERENCE from the identifier, not a documented field.
    centre = np.array([int(str(d)[0]) for d in donors], dtype=np.int64)

    return {
        "values": expression.to_numpy(dtype=np.float64),
        "symbols": features["genesymbol"].fillna("").to_numpy(dtype=str),
        "donors": donors,
        "age": age,
        "sex": sex,
        "centre": centre,
        "phenotype": phenotype,
    }


def build_universe(data, graph, rng):
    """One candidate feature per gene symbol, plus the frozen decoy block.

    The choice of probe within a gene is made per training fold, so this returns
    every probe position grouped by symbol and leaves the pick to the fold.
    """
    symbols = data["symbols"]
    annotated = np.nonzero(symbols != "")[0]
    by_symbol: dict[str, list[int]] = {}
    for position in annotated:
        by_symbol.setdefault(str(symbols[position]), []).append(int(position))

    on_network = sorted(set(by_symbol) & set(graph.index))
    off_array = sorted(set(graph.index) - set(by_symbol))
    if len(off_array) < N_DECOYS:
        raise H2aError("too few network genes absent from the array to build decoys")

    # A decoy carries the symbol of a network gene that the array does not
    # measure, and the values of a real probe with its donors shuffled. A
    # selector that retains one has made an error: the signal is destroyed but
    # the gene is still a plausible network neighbour.
    decoy_symbols = [off_array[i] for i in
                     rng.choice(len(off_array), size=N_DECOYS, replace=False)]
    donor_count = data["values"].shape[1]
    decoy_sources = rng.choice(annotated, size=N_DECOYS, replace=False)
    decoy_values = np.empty((N_DECOYS, donor_count), dtype=np.float64)
    for row, source in enumerate(decoy_sources):
        decoy_values[row] = data["values"][source][rng.permutation(donor_count)]

    return {
        "by_symbol": by_symbol,
        "genes_on_network": on_network,
        "decoy_symbols": decoy_symbols,
        "decoy_values": decoy_values,
    }


# --------------------------------------------------------------------------
# Folds
# --------------------------------------------------------------------------
def repeated_folds(n, *, n_splits, n_repeats, seed):
    """Donor-disjoint folds, identical for every arm, model and panel size."""
    out = []
    for repeat in range(n_repeats):
        order = np.random.default_rng(seed + repeat).permutation(n)
        for split in np.array_split(order, n_splits):
            test = np.sort(split)
            train = np.setdiff1d(np.arange(n), test)
            out.append((train, test))
    return out


def centre_folds(centre):
    out = []
    for value in np.unique(centre):
        test = np.nonzero(centre == value)[0]
        out.append((np.setdiff1d(np.arange(centre.size), test), test))
    return out


# --------------------------------------------------------------------------
# Per-fold feature matrix
# --------------------------------------------------------------------------
def fold_matrix(data, universe, train, donor_index):
    """Assemble the candidate matrix with every choice made on `train` only.

    Returns the matrix over all donors in `donor_index`, the gene name of each
    column, and a mask marking the decoy columns.
    """
    values = data["values"][:, donor_index]
    names: list[str] = []
    columns: list[np.ndarray] = []
    is_decoy: list[bool] = []

    for symbol, positions in universe["by_symbol"].items():
        best_position, best_mean = None, -np.inf
        for position in positions:
            row = values[position]
            observed = row[train]
            present = np.isfinite(observed)
            if present.mean() < 1.0 - MAX_PROBE_MISSING_FRACTION:
                continue
            mean = float(observed[present].mean())
            if mean > best_mean:
                best_position, best_mean = position, mean
        if best_position is None:
            continue
        row = values[best_position].copy()
        absent = ~np.isfinite(row)
        if absent.any():
            row[absent] = best_mean
        names.append(symbol)
        columns.append(row)
        is_decoy.append(False)

    for symbol, row in zip(universe["decoy_symbols"], universe["decoy_values"]):
        column = row[donor_index].copy()
        absent = ~np.isfinite(column)
        if absent.any():
            column[absent] = float(np.nanmean(column[train]))
        names.append(symbol)
        columns.append(column)
        is_decoy.append(True)

    matrix = np.column_stack(columns)
    return matrix, np.asarray(names, dtype=object), np.asarray(is_decoy, dtype=bool)


def baseline_matrix(data, donor_index, train):
    age = data["age"][donor_index].copy()
    sex = data["sex"][donor_index].copy()
    for vector in (age, sex):
        absent = ~np.isfinite(vector)
        if absent.any():
            observed = vector[train][np.isfinite(vector[train])]
            if observed.size == 0:
                raise H2aError("a baseline covariate is absent throughout a training fold")
            # Age takes the training median; sex, being binary, takes the
            # training majority.
            fill = (np.median(observed) if vector is age
                    else float(np.round(observed.mean())))
            vector[absent] = fill
    return np.column_stack([age, sex])


# --------------------------------------------------------------------------
# Rankings
# --------------------------------------------------------------------------
def _full_order(ranking):
    """Seeds and candidates in one list, ordered by proximity.

    Seeds must not be ordered alphabetically. A multi-seed variant can carry a
    hundred seeds, so its top-k for small k would otherwise be decided by the
    alphabet rather than by the mechanism, which would make the equal-dimension
    comparison meaningless. Restart mass keeps the seeds ahead of the
    candidates; within the seeds, the walk's own score decides.
    """
    scores = ranking.scores
    return sorted(scores, key=lambda gene: (-scores[gene], gene))


def mechanistic_orders(graph, knowledge_dir: Path):
    """Compute every outcome-blind ordering once; they do not vary by fold."""
    orders = {}
    pathway_seeds = seeds_from_reactome(
        knowledge_dir / "ReactomePathways.gmt.zip", REACTOME_PATHWAY_PATTERN, graph)
    shared = rank_by_mechanistic_proximity(graph, pathway_seeds)
    for endpoint in ENDPOINTS:
        own = rank_by_mechanistic_proximity(
            graph, seeds_from_endpoint_label(endpoint, graph))
        wrong_label = WRONG_SEED_ROTATION[endpoint]
        wrong = rank_by_mechanistic_proximity(
            graph, seeds_from_endpoint_label(wrong_label, graph))
        orders[endpoint] = {
            "netpk_target_seeded": _full_order(own),
            "netpk_wrong_mechanism": _full_order(wrong),
            "netpk_pathway_seeded": _full_order(shared),
            "target_gene_only": list(own.seeds),
            "_wrong_seed_label": wrong_label,
            "_pathway_seed_count": len(shared.seeds),
        }
        substrate = DOCUMENTED_SUBSTRATE.get(endpoint)
        if substrate:
            seeds = seeds_from_ctd_substrate(
                knowledge_dir / "CTD_chem_gene_ixns.tsv.gz", substrate, graph)
            chain = rank_by_mechanistic_proximity(graph, seeds)
            orders[endpoint]["netpk_substrate_seeded"] = _full_order(chain)
            orders[endpoint]["_substrate"] = substrate
            orders[endpoint]["_substrate_seed_count"] = len(chain.seeds)
    return orders


def _spearman_strength(matrix, y):
    ranked = np.argsort(np.argsort(matrix, axis=0), axis=0).astype(np.float64)
    ranked -= ranked.mean(axis=0)
    target = np.argsort(np.argsort(y)).astype(np.float64)
    target -= target.mean()
    denominator = np.sqrt((ranked ** 2).sum(axis=0) * (target ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = (ranked * target[:, None]).sum(axis=0) / denominator
    return np.nan_to_num(np.abs(rho))


def statistical_orders(matrix, names, y, *, seed):
    """Orderings that read the endpoint, computed on a training fold only."""
    from sklearn.linear_model import lars_path
    from sklearn.preprocessing import StandardScaler

    orders = {}
    orders["univariate_spearman"] = [
        str(names[i]) for i in np.argsort(-_spearman_strength(matrix, y), kind="stable")]

    scaler = StandardScaler().fit(matrix)
    scaled = scaler.transform(matrix)
    centred = y - y.mean()
    _, _, coefficients = lars_path(scaled, centred, method="lasso",
                                   max_iter=LASSO_PATH_LIMIT)
    entry_order, seen = [], set()
    for step in range(coefficients.shape[1]):
        for index in np.nonzero(coefficients[:, step])[0]:
            if index not in seen:
                seen.add(int(index))
                entry_order.append(int(index))
    remaining = [i for i in np.argsort(-np.abs(coefficients[:, -1]), kind="stable")
                 if i not in seen]
    orders["lasso_entry_order"] = [str(names[i]) for i in entry_order + remaining]

    import xgboost as xgb
    model = xgb.XGBRegressor(n_estimators=200, max_depth=3, learning_rate=0.1,
                             subsample=0.9, colsample_bytree=0.3,
                             tree_method="hist", random_state=seed, n_jobs=0,
                             verbosity=0)
    model.fit(matrix, y)
    booster = model.get_booster()
    booster.feature_names = [f"f{i}" for i in range(matrix.shape[1])]
    gain = booster.get_score(importance_type="total_gain")
    strength = np.zeros(matrix.shape[1], dtype=np.float64)
    for key, value in gain.items():
        strength[int(key[1:])] = float(value)
    orders["xgboost_gain"] = [
        str(names[i]) for i in np.argsort(-strength, kind="stable")]
    return orders


# --------------------------------------------------------------------------
# Prediction
# --------------------------------------------------------------------------
def ridge_prediction(x_train, y_train, x_test, *, inner_folds, seed):
    """Ridge with the penalty chosen by inner cross-validation on the training fold."""
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    n = x_train.shape[0]
    order = np.random.default_rng(seed).permutation(n)
    splits = np.array_split(order, inner_folds)
    errors = np.zeros(len(RIDGE_ALPHAS), dtype=np.float64)
    for held in splits:
        keep = np.setdiff1d(np.arange(n), held)
        scaler = StandardScaler().fit(x_train[keep])
        inner_train = scaler.transform(x_train[keep])
        inner_test = scaler.transform(x_train[held])
        for position, alpha in enumerate(RIDGE_ALPHAS):
            model = Ridge(alpha=alpha).fit(inner_train, y_train[keep])
            errors[position] += np.abs(
                model.predict(inner_test) - y_train[held]).sum()
    alpha = RIDGE_ALPHAS[int(np.argmin(errors))]
    scaler = StandardScaler().fit(x_train)
    model = Ridge(alpha=alpha).fit(scaler.transform(x_train), y_train)
    return model.predict(scaler.transform(x_test)), alpha


# --------------------------------------------------------------------------
# Losses and contrasts
# --------------------------------------------------------------------------
def trimmed_mean(values, fraction):
    values = np.sort(np.asarray(values, dtype=np.float64))
    keep = values.size - int(np.floor(fraction * values.size))
    return float(values[:keep].mean()) if keep > 0 else float("nan")


def loss_summary(errors, y_true, y_pred):
    positive = (y_true > 0) & (y_pred > 0)
    log_errors = (np.abs(np.log(y_pred[positive]) - np.log(y_true[positive]))
                  if positive.any() else np.array([np.nan]))
    return {
        "mae": float(np.mean(errors)),
        "trimmed_mae": trimmed_mean(errors, TRIM_FRACTION),
        "log_mae": float(np.nanmean(log_errors)),
        "log_mae_pairs": int(positive.sum()),
        "n": int(errors.size),
    }


def paired_contrast(errors_a, errors_b, *, seed, label):
    """Paired difference a - b with a donor bootstrap interval."""
    errors_a = np.asarray(errors_a, dtype=np.float64)
    errors_b = np.asarray(errors_b, dtype=np.float64)
    difference = errors_a - errors_b
    observed = float(difference.mean())
    base = float(errors_b.mean())
    rng = np.random.default_rng(seed)
    n = difference.size
    draws = np.empty(N_BOOTSTRAP, dtype=np.float64)
    relative = np.empty(N_BOOTSTRAP, dtype=np.float64)
    for draw in range(N_BOOTSTRAP):
        pick = rng.integers(0, n, size=n)
        draws[draw] = difference[pick].mean()
        denominator = errors_b[pick].mean()
        relative[draw] = (draws[draw] / denominator
                          if denominator > 0 else np.nan)
    low, high = np.percentile(draws, [2.5, 97.5])
    rel_low, rel_high = np.nanpercentile(relative, [2.5, 97.5])
    return {
        "label": label,
        "delta_mae": observed,
        "reference_mae": base,
        "relative": observed / base if base > 0 else None,
        "ci_low": float(low),
        "ci_high": float(high),
        "relative_ci_low": float(rel_low),
        "relative_ci_high": float(rel_high),
        "interval_excludes_zero": bool(low > 0 or high < 0),
        "non_inferior_within_margin": bool(rel_high <= NON_INFERIORITY_MARGIN),
        "how_to_read": (
            "Negative means the first arm has the smaller error. "
            "interval_excludes_zero false means no difference was DETECTED, "
            "never that there is none. non_inferior_within_margin is true only "
            "when the whole interval stays below the frozen 5 per cent margin."),
    }


# --------------------------------------------------------------------------
# One endpoint
# --------------------------------------------------------------------------
def run_endpoint(endpoint, data, universe, orders, *, folds, fold_label, log):
    import pandas as pd

    y_all = pd.to_numeric(data["phenotype"].loc[endpoint],
                          errors="coerce").to_numpy(float)
    keep = np.nonzero(np.isfinite(y_all))[0]
    y = y_all[keep]
    if y.size < 50:
        raise H2aError(f"{endpoint}: only {y.size} donors with a measured value")

    donor_folds = [(np.searchsorted(keep, np.intersect1d(keep, train)),
                    np.searchsorted(keep, np.intersect1d(keep, test)))
                   for train, test in folds]
    donor_folds = [(tr, te) for tr, te in donor_folds if te.size and tr.size]

    arm_names = ["netpk_target_seeded", "netpk_wrong_mechanism",
                 "netpk_pathway_seeded", "univariate_spearman",
                 "lasso_entry_order", "xgboost_gain"]
    if "netpk_substrate_seeded" in orders[endpoint]:
        arm_names.append("netpk_substrate_seeded")

    errors = {arm: {k: np.full(y.size, np.nan) for k in PANEL_SIZES}
              for arm in arm_names}
    errors["random_panel"] = {k: np.full((RANDOM_REPLICATIONS, y.size), np.nan)
                              for k in PANEL_SIZES}
    predictions = {arm: {k: np.full(y.size, np.nan) for k in PANEL_SIZES}
                   for arm in arm_names}
    fixed_arms = {"target_gene_only": np.full(y.size, np.nan),
                  "baseline_age_sex": np.full(y.size, np.nan),
                  "no_selection": np.full(y.size, np.nan)}
    fixed_predictions = {arm: np.full(y.size, np.nan) for arm in fixed_arms}
    decoys_retained = {arm: {k: 0 for k in PANEL_SIZES} for arm in arm_names}
    selected_sets = {arm: {k: [] for k in PANEL_SIZES} for arm in arm_names}

    started = time.time()
    for fold_number, (train, test) in enumerate(donor_folds):
        matrix, names, is_decoy = fold_matrix(data, universe, train, keep)
        position_of = {str(name): index for index, name in enumerate(names)}
        decoy_names = {str(names[i]) for i in np.nonzero(is_decoy)[0]}
        clinical = baseline_matrix(data, keep, train)

        stat = statistical_orders(matrix[train], names, y[train],
                                  seed=SEED + fold_number)
        fold_orders = dict(stat)
        for arm in ("netpk_target_seeded", "netpk_wrong_mechanism",
                    "netpk_pathway_seeded", "netpk_substrate_seeded"):
            if arm in orders[endpoint]:
                fold_orders[arm] = [g for g in orders[endpoint][arm]
                                    if g in position_of]

        def predict(column_names):
            picked = [position_of[name] for name in column_names
                      if name in position_of]
            block = (np.column_stack([clinical, matrix[:, picked]])
                     if picked else clinical)
            return ridge_prediction(block[train], y[train], block[test],
                                    inner_folds=INNER_FOLDS,
                                    seed=SEED + fold_number)[0]

        for arm in arm_names:
            order = fold_orders[arm]
            for k in PANEL_SIZES:
                chosen = order[:k]
                predicted = predict(chosen)
                predictions[arm][k][test] = predicted
                errors[arm][k][test] = np.abs(predicted - y[test])
                decoys_retained[arm][k] += len(set(chosen) & decoy_names)
                if fold_number == 0:
                    selected_sets[arm][k] = list(chosen)

        target = [g for g in orders[endpoint]["target_gene_only"] if g in position_of]
        fixed_predictions["target_gene_only"][test] = predict(target)
        fixed_arms["target_gene_only"][test] = np.abs(
            fixed_predictions["target_gene_only"][test] - y[test])
        fixed_predictions["baseline_age_sex"][test] = predict([])
        fixed_arms["baseline_age_sex"][test] = np.abs(
            fixed_predictions["baseline_age_sex"][test] - y[test])
        every = [str(name) for name in names]
        fixed_predictions["no_selection"][test] = predict(every)
        fixed_arms["no_selection"][test] = np.abs(
            fixed_predictions["no_selection"][test] - y[test])

        rng = np.random.default_rng(SEED + 7919 * fold_number)
        for replication in range(RANDOM_REPLICATIONS):
            draw = rng.choice(len(names), size=max(PANEL_SIZES), replace=False)
            for k in PANEL_SIZES:
                chosen = [str(names[i]) for i in draw[:k]]
                predicted = predict(chosen)
                errors["random_panel"][k][replication, test] = np.abs(
                    predicted - y[test])
        log(f"    fold {fold_number + 1}/{len(donor_folds)} "
            f"({time.time() - started:.0f}s)")

    summary = {
        "endpoint": endpoint,
        "fold_design": fold_label,
        "n_donors": int(y.size),
        "folds": len(donor_folds),
        "wrong_seed_label": orders[endpoint]["_wrong_seed_label"],
        "target_gene_count": len(orders[endpoint]["target_gene_only"]),
        "pathway_seed_count": orders[endpoint]["_pathway_seed_count"],
        "substrate": orders[endpoint].get("_substrate"),
        "substrate_seed_count": orders[endpoint].get("_substrate_seed_count"),
        "losses": {},
        "decoys_retained": decoys_retained,
        "first_fold_selection": {arm: {k: selected_sets[arm][k]
                                       for k in (1, 5, 20)}
                                 for arm in arm_names},
    }
    for arm in arm_names:
        summary["losses"][arm] = {
            str(k): loss_summary(errors[arm][k], y, predictions[arm][k])
            for k in PANEL_SIZES}
    for arm, vector in fixed_arms.items():
        summary["losses"][arm] = {
            "fixed": loss_summary(vector, y, fixed_predictions[arm])}
    summary["losses"]["random_panel"] = {
        str(k): {
            "mae_mean_over_replications": float(
                np.mean([np.mean(errors["random_panel"][k][r])
                         for r in range(RANDOM_REPLICATIONS)])),
            "mae_replication_2_5_pct": float(np.percentile(
                [np.mean(errors["random_panel"][k][r])
                 for r in range(RANDOM_REPLICATIONS)], 2.5)),
            "mae_replication_97_5_pct": float(np.percentile(
                [np.mean(errors["random_panel"][k][r])
                 for r in range(RANDOM_REPLICATIONS)], 97.5)),
            "replications": RANDOM_REPLICATIONS,
        } for k in PANEL_SIZES}

    contrasts = []
    for k in PANEL_SIZES:
        best_statistical = min(
            ("univariate_spearman", "lasso_entry_order", "xgboost_gain"),
            key=lambda arm: np.mean(errors[arm][k]))
        contrasts.append(paired_contrast(
            errors["netpk_target_seeded"][k], errors[best_statistical][k],
            seed=SEED + k, label=f"netpk_target_seeded_minus_best_statistical_k{k}"))
        contrasts[-1]["best_statistical_arm"] = best_statistical
        contrasts.append(paired_contrast(
            errors["netpk_target_seeded"][k], errors["netpk_wrong_mechanism"][k],
            seed=SEED + 100 + k, label=f"netpk_target_minus_wrong_mechanism_k{k}"))
        contrasts.append(paired_contrast(
            errors["netpk_target_seeded"][k],
            errors["random_panel"][k].mean(axis=0),
            seed=SEED + 200 + k, label=f"netpk_target_minus_random_mean_k{k}"))
        contrasts.append(paired_contrast(
            errors["netpk_target_seeded"][k], fixed_arms["target_gene_only"],
            seed=SEED + 300 + k, label=f"netpk_target_minus_target_only_k{k}"))
        contrasts.append(paired_contrast(
            errors["netpk_pathway_seeded"][k], errors["netpk_target_seeded"][k],
            seed=SEED + 400 + k, label=f"netpk_pathway_minus_target_seeded_k{k}"))
    contrasts.append(paired_contrast(
        fixed_arms["target_gene_only"], fixed_arms["baseline_age_sex"],
        seed=SEED + 500, label="target_only_minus_baseline"))
    contrasts.append(paired_contrast(
        fixed_arms["no_selection"], fixed_arms["baseline_age_sex"],
        seed=SEED + 600, label="no_selection_minus_baseline"))
    summary["contrasts"] = contrasts
    return summary


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------
def main(argv=None) -> int:
    global RANDOM_REPLICATIONS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc", required=True)
    parser.add_argument("--knowledge", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--endpoints", nargs="*", default=list(ENDPOINTS))
    parser.add_argument("--design", choices=("repeated", "centre"),
                        default="repeated")
    parser.add_argument("--random-replications", type=int,
                        default=RANDOM_REPLICATIONS,
                        help="Lower only for a smoke test; the protocol "
                             "requires at least 100 for a reported run.")
    args = parser.parse_args(argv)
    RANDOM_REPLICATIONS = int(args.random_replications)

    def log(message):
        print(message, flush=True)

    knowledge = Path(args.knowledge)
    graph = load_string_symbol_graph(
        knowledge / "9606.protein.info.v12.0.txt.gz",
        knowledge / "9606.protein.links.v12.0.txt.gz")
    log(f"network: {graph.n_nodes} nodes, {graph.n_edges} edges")

    data = load_hlc(Path(args.hlc))
    log(f"donors with expression: {data['values'].shape[1]}; "
        f"centres inferred from the identifier: "
        f"{dict(zip(*np.unique(data['centre'], return_counts=True)))}")

    universe = build_universe(data, graph, np.random.default_rng(SEED))
    log(f"candidate genes: {len(universe['by_symbol'])}; "
        f"on the network: {len(universe['genes_on_network'])}; "
        f"decoys: {N_DECOYS}")

    orders = mechanistic_orders(graph, knowledge)
    log("mechanistic orderings computed once, outcome-blind")

    if args.design == "centre":
        folds, fold_label = centre_folds(data["centre"]), "leave_one_centre_out"
    else:
        folds = repeated_folds(data["values"].shape[1], n_splits=OUTER_FOLDS,
                               n_repeats=OUTER_REPEATS, seed=SEED)
        fold_label = f"{OUTER_FOLDS}_fold_x{OUTER_REPEATS}_repeats"

    report = {
        "analysis": "H2a: mechanistic versus statistical selection, equal dimension",
        "protocol_version": "1.4.1",
        "run_on": date.today().isoformat(),
        "endpoint_class": "functional_in_vitro",
        "endpoint_class_note": (
            "Measured microsomal activity, not clinical pharmacokinetics. This is "
            "the functional analogue of H2 and is never pooled with a clinical PK "
            "result nor reported in the same effect column."),
        "frozen": {
            "seed": SEED, "outer_folds": OUTER_FOLDS,
            "outer_repeats": OUTER_REPEATS, "inner_folds": INNER_FOLDS,
            "panel_sizes": list(PANEL_SIZES),
            "random_replications": RANDOM_REPLICATIONS,
            "decoys": N_DECOYS, "bootstrap": N_BOOTSTRAP,
            "non_inferiority_margin": NON_INFERIORITY_MARGIN,
            "trim_fraction": TRIM_FRACTION,
            "selector": dict(SELECTOR_FROZEN),
            "baseline_covariates": ["age", "sex"],
            "fold_design": fold_label,
        },
        "network": dict(graph.provenance),
        "universe": {
            "candidate_genes": len(universe["by_symbol"]),
            "candidate_genes_on_network": len(universe["genes_on_network"]),
            "decoy_block": N_DECOYS,
        },
        "endpoints": {},
    }

    for endpoint in args.endpoints:
        log(f"  {endpoint}")
        report["endpoints"][endpoint] = run_endpoint(
            endpoint, data, universe, orders, folds=folds,
            fold_label=fold_label, log=log)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    log(f"wrote {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
