#!/usr/bin/env python
"""H2a: does a mechanistic ordering select better than a statistical one?

Frozen by protocol version 1.4.3, section 5. The endpoint is measured
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
CONTROL_SYMBOL_PATTERN = __import__("re").compile(
    r"^(NegCont|PosCont|Corner|DarkCorner|SpikeIn|GE_Bright|E1A_r60|ETG)",
    __import__("re").IGNORECASE)
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
        "feature_ids": features["feature_id"].astype(str).to_numpy(dtype=str),
        "donors": donors,
        "age": age,
        "sex": sex,
        "centre": centre,
        "phenotype": phenotype,
    }


# A decoy column carries a network gene's symbol and a real probe's values with
# the donors shuffled. It is measured on nothing, so it has no feature.
DECOY_FEATURE_ID = "decoy_no_feature"


def qualified_pool_features(pool_path, pool_name):
    """The admissible feature identifiers of one named pool, as a set."""
    pool = json.loads(Path(pool_path).read_text(encoding="utf-8"))
    admissible = pool.get("admissible_features", {}).get(pool_name)
    if not admissible:
        raise H2aError(f"the pool file carries no admissible set named {pool_name!r}")
    return set(admissible)


def restrict_universe_to_qualified_pool(universe, data, pool_path, pool_name):
    """Keep only the probes a graded identity qualification admits.

    Applied at the **probe** level and only then grouped, so a gene with one
    admissible and one inadmissible probe keeps the admissible one instead of
    being dropped whole. The pool is produced by
    `scripts/build_qualified_probe_pools.py`, which grades each feature by how
    strongly two independent annotations agree on a human gene, resolves
    accessions to genes before comparing them, and keeps the undetermined and
    conflicting strata separate rather than deleting them.
    """
    pool = json.loads(Path(pool_path).read_text(encoding="utf-8"))
    admissible = pool.get("admissible_features", {}).get(pool_name)
    if not admissible:
        raise H2aError(f"the pool file carries no admissible set named {pool_name!r}")
    allowed = set(admissible)
    position_of = {str(fid): index
                   for index, fid in enumerate(data["feature_ids"])}
    keep_positions = {position_of[fid] for fid in allowed if fid in position_of}
    if not keep_positions:
        raise H2aError("no admissible probe matches this expression table")

    kept, dropped_probes = {}, 0
    for symbol, positions in universe["by_symbol"].items():
        surviving = [p for p in positions if p in keep_positions]
        dropped_probes += len(positions) - len(surviving)
        if surviving:
            kept[symbol] = surviving
    if not kept:
        raise H2aError("restricting to the qualified pool left no candidate")

    # A technical control carries no transcript, so its presence here would mean
    # the pool was built wrongly. Checked rather than assumed.
    controls = [s for s in kept if CONTROL_SYMBOL_PATTERN.match(s)]
    if controls:
        raise H2aError(
            f"{len(controls)} technical control symbols survived the pool, "
            f"starting with {controls[:3]}")

    restricted = dict(universe)
    restricted["by_symbol"] = kept
    restricted["qualified_pool"] = {
        "pool_file": Path(pool_path).name,
        "pool_name": pool_name,
        "rule_version": pool.get("rule_version"),
        "admitted_verdicts": pool.get("pools", {}).get(pool_name, {}).get(
            "admitted_verdicts"),
        "symbols_before": len(universe["by_symbol"]),
        "symbols_after": len(kept),
        "probes_dropped": dropped_probes,
        "technical_controls_present": 0,
    }
    return restricted


def restrict_universe_to_network(universe, graph):
    """Keep only candidates the mechanistic ordering could ever reach.

    Sensitivity analysis, added after the primary run and labelled as such. In
    the primary comparison every arm selects the same *number* of genes but not
    from the same *pool*: the mechanistic ordering can only reach genes present
    in the public network, while a statistical selector ranks the whole array --
    including 1,007 technical control probes and about eleven thousand entries
    carrying an accession rather than a gene symbol, none of which is in the
    network. A statistical arm selected a negative control probe on one
    endpoint, so the asymmetry is not hypothetical.

    That confounds ordering quality with pool size. Restricting every arm to the
    genes the network contains separates the two. It does not replace the
    primary result, which stays as frozen and as run.
    """
    kept = {symbol: positions for symbol, positions in universe["by_symbol"].items()
            if symbol in graph.index}
    if not kept:
        raise H2aError("restricting to the network left no candidate")
    restricted = dict(universe)
    restricted["by_symbol"] = kept
    restricted["restricted_to_network"] = True
    restricted["candidates_before_restriction"] = len(universe["by_symbol"])
    restricted["candidates_after_restriction"] = len(kept)
    return restricted


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


def repeat_assignment(n_folds, *, n_splits, repeated):
    """Which repeat each fold belongs to.

    Out-of-fold results must be stored per repeat. A donor appears in exactly
    one test fold per repeat, so writing every fold into one vector per arm
    lets the last repeat overwrite the others and silently reduces a design
    declared as `n_splits` by `n_repeats` to a single pass.
    """
    if n_folds < 1:
        raise H2aError("a design needs at least one fold")
    if not repeated:
        return 1, [0] * n_folds
    if n_splits < 1 or n_folds % n_splits:
        raise H2aError("fold count does not factor into the declared repeats")
    repeats = n_folds // n_splits
    return repeats, [index // n_splits for index in range(n_folds)]


def centre_folds(centre):
    out = []
    for value in np.unique(centre):
        test = np.nonzero(centre == value)[0]
        out.append((np.setdiff1d(np.arange(centre.size), test), test))
    return out


# --------------------------------------------------------------------------
# Per-fold feature matrix
# --------------------------------------------------------------------------
def fold_matrix(data, universe, train, donor_index, chosen=None):
    """Assemble the candidate matrix with every choice made on `train` only.

    Returns the matrix over all donors in `donor_index`, the gene name of each
    column, and a mask marking the decoy columns.

    When a dict is passed as `chosen`, it is filled with the feature identifier
    actually used for each symbol in this fold, so any selection can be traced
    back to the probe it was measured on and checked against the qualified
    pool. A decoy column has no feature identifier and is recorded as such.
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
        if chosen is not None:
            chosen[symbol] = str(data["feature_ids"][best_position])

    for symbol, row in zip(universe["decoy_symbols"], universe["decoy_values"]):
        column = row[donor_index].copy()
        absent = ~np.isfinite(column)
        if absent.any():
            column[absent] = float(np.nanmean(column[train]))
        names.append(symbol)
        columns.append(column)
        is_decoy.append(True)
        if chosen is not None:
            chosen[symbol] = DECOY_FEATURE_ID

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
    """Absolute rank correlation of every column with the endpoint.

    Ties take the average rank. Breaking them by position instead, as a
    double argsort does, hands a constant column a perfect correlation with any
    monotone endpoint -- it would turn an uninformative probe into the
    univariate arm's first choice. The expression values are clipped to a
    narrow range and carry real ties, so this is not a corner case.
    """
    from scipy.stats import rankdata

    ranked = rankdata(matrix, axis=0).astype(np.float64)
    ranked -= ranked.mean(axis=0)
    target = rankdata(y).astype(np.float64)
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
def ridge_path_predictions(x_train, y_train, x_test, alphas):
    """Ridge predictions for every penalty in `alphas`, sharing one decomposition.

    Statistically identical to fitting each penalty separately: the training
    fold's mean and standard deviation standardise both blocks, the response is
    centred to carry the intercept, and the penalty is applied to the
    standardised coefficients. One decomposition serves the whole path, which
    matters because the comparison needs a hundred random panels per fold and a
    fresh solver object per penalty would dominate the run.

    Two branches, because the arms differ by four orders of magnitude in width:
    the primal form through a singular value decomposition when there are at
    least as many donors as variables, and the dual form through the Gram
    matrix when there are not -- as in the no-selection arm, which carries
    every candidate gene.
    """
    mean = x_train.mean(axis=0)
    spread = x_train.std(axis=0)
    spread = np.where(spread > 0.0, spread, 1.0)
    train = (x_train - mean) / spread
    test = (x_test - mean) / spread
    centre = float(y_train.mean())
    centred = y_train - centre

    n, p = train.shape
    out = []
    if p <= n:
        u, singular, vt = np.linalg.svd(train, full_matrices=False)
        projected = u.T @ centred
        for alpha in alphas:
            shrunk = singular / (singular ** 2 + alpha)
            out.append(test @ (vt.T @ (shrunk * projected)) + centre)
    else:
        gram = train @ train.T
        eigenvalues, vectors = np.linalg.eigh(gram)
        eigenvalues = np.clip(eigenvalues, 0.0, None)
        projected = vectors.T @ centred
        cross = test @ train.T
        for alpha in alphas:
            dual = vectors @ (projected / (eigenvalues + alpha))
            out.append(cross @ dual + centre)
    return out


def ridge_prediction(x_train, y_train, x_test, *, inner_folds, seed):
    """Ridge with the penalty chosen by inner cross-validation on the training fold."""
    n = x_train.shape[0]
    order = np.random.default_rng(seed).permutation(n)
    splits = np.array_split(order, inner_folds)
    errors = np.zeros(len(RIDGE_ALPHAS), dtype=np.float64)
    for held in splits:
        keep = np.setdiff1d(np.arange(n), held)
        predictions = ridge_path_predictions(
            x_train[keep], y_train[keep], x_train[held], RIDGE_ALPHAS)
        for position, predicted in enumerate(predictions):
            errors[position] += np.abs(predicted - y_train[held]).sum()
    best = int(np.argmin(errors))
    alpha = RIDGE_ALPHAS[best]
    final = ridge_path_predictions(x_train, y_train, x_test, (alpha,))[0]
    return final, alpha


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

    # Every out-of-fold quantity is indexed by repeat as well as by donor.
    # Writing into one vector per arm would make the second repeat overwrite the
    # first, so a design declared as five folds by two repeats would silently
    # report the last repeat alone. The donor's error is averaged across repeats
    # before any loss or bootstrap, so a repeat never counts as an extra donor.
    repeats, repeat_of = repeat_assignment(
        len(donor_folds), n_splits=OUTER_FOLDS,
        repeated=fold_label.endswith("repeats"))

    def blank():
        return {k: np.full((repeats, y.size), np.nan) for k in PANEL_SIZES}

    errors = {arm: blank() for arm in arm_names}
    errors["random_panel"] = {
        k: np.full((repeats, RANDOM_REPLICATIONS, y.size), np.nan)
        for k in PANEL_SIZES}
    predictions = {arm: blank() for arm in arm_names}
    fixed_arms = {name: np.full((repeats, y.size), np.nan)
                  for name in ("target_gene_only", "baseline_age_sex",
                               "no_selection")}
    fixed_predictions = {arm: np.full((repeats, y.size), np.nan)
                         for arm in fixed_arms}
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

        repeat = repeat_of[fold_number]
        for arm in arm_names:
            order = fold_orders[arm]
            for k in PANEL_SIZES:
                chosen = order[:k]
                predicted = predict(chosen)
                predictions[arm][k][repeat, test] = predicted
                errors[arm][k][repeat, test] = np.abs(predicted - y[test])
                decoys_retained[arm][k] += len(set(chosen) & decoy_names)
                if fold_number == 0:
                    selected_sets[arm][k] = list(chosen)

        target = [g for g in orders[endpoint]["target_gene_only"] if g in position_of]
        for name, columns in (("target_gene_only", target),
                              ("baseline_age_sex", []),
                              ("no_selection", [str(n) for n in names])):
            predicted = predict(columns)
            fixed_predictions[name][repeat, test] = predicted
            fixed_arms[name][repeat, test] = np.abs(predicted - y[test])

        rng = np.random.default_rng(SEED + 7919 * fold_number)
        for replication in range(RANDOM_REPLICATIONS):
            draw = rng.choice(len(names), size=max(PANEL_SIZES), replace=False)
            for k in PANEL_SIZES:
                chosen = [str(names[i]) for i in draw[:k]]
                predicted = predict(chosen)
                errors["random_panel"][k][repeat, replication, test] = np.abs(
                    predicted - y[test])
        log(f"    fold {fold_number + 1}/{len(donor_folds)} "
            f"({time.time() - started:.0f}s)")

    def per_donor(block):
        """Collapse the repeat dimension: one value per donor, averaged.

        The donor is the independent unit. Averaging across repeats first keeps
        it that way, so the bootstrap cannot treat a repeated measurement of the
        same donor as a second donor.
        """
        return np.nanmean(np.asarray(block, dtype=np.float64), axis=0)

    summary = {
        "endpoint": endpoint,
        "fold_design": fold_label,
        "repeats": repeats,
        "n_donors": int(y.size),
        "folds": len(donor_folds),
        "analysis_unit": (
            "one averaged absolute error per donor; repeats are averaged before "
            "any loss or bootstrap, so a repeat is never a second donor"),
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
    donor_errors = {arm: {k: per_donor(errors[arm][k]) for k in PANEL_SIZES}
                    for arm in arm_names}
    donor_errors["random_panel"] = {
        k: per_donor(errors["random_panel"][k].mean(axis=1))
        for k in PANEL_SIZES}
    fixed_errors = {arm: per_donor(block) for arm, block in fixed_arms.items()}

    for arm in arm_names:
        summary["losses"][arm] = {
            str(k): loss_summary(donor_errors[arm][k], y,
                                 per_donor(predictions[arm][k]))
            for k in PANEL_SIZES}
    for arm, block in fixed_arms.items():
        summary["losses"][arm] = {
            "fixed": loss_summary(fixed_errors[arm], y,
                                  per_donor(fixed_predictions[arm]))}
    summary["losses"]["random_panel"] = {}
    for k in PANEL_SIZES:
        # Each replication keeps its own out-of-fold error, averaged over
        # repeats, so the spread below is across panels and not across folds.
        by_replication = [float(np.nanmean(per_donor(
            errors["random_panel"][k][:, replication, :])))
            for replication in range(RANDOM_REPLICATIONS)]
        summary["losses"]["random_panel"][str(k)] = {
            "mae_mean_over_replications": float(np.mean(by_replication)),
            "mae_replication_2_5_pct": float(np.percentile(by_replication, 2.5)),
            "mae_replication_97_5_pct": float(np.percentile(by_replication, 97.5)),
            "replications": RANDOM_REPLICATIONS,
        }

    # ------------------------------------------------------------------
    # Destroyed-signal diagnostic.
    #
    # The random-symbol decoy block turns out to measure susceptibility to a
    # spurious correlation, which only an arm that reads the endpoint can have.
    # A mechanistic arm scores zero on it for an uninteresting reason: a decoy
    # carries the symbol of a randomly drawn network gene, which is almost
    # always far from the seed, so it was never a candidate for the top of the
    # ordering. That is not discrimination.
    #
    # This is the test that does expose the weakness. The target gene's own
    # column has its donor order permuted, so the gene keeps its position in
    # the mechanistic ordering while carrying no information at all. An arm
    # that reads the endpoint must notice and choose something else; an arm
    # that does not read it cannot.
    # ------------------------------------------------------------------
    target_gene = orders[endpoint]["target_gene_only"][0]
    destroyed = {"gene": target_gene,
                 "blind_errors": np.full((repeats, y.size), np.nan),
                 "univariate_errors": np.full((repeats, y.size), np.nan),
                 "univariate_still_picked_the_dead_gene": 0,
                 "folds_examined": 0}
    for fold_number, (train, test) in enumerate(donor_folds):
        matrix, names, _ = fold_matrix(data, universe, train, keep)
        position_of = {str(name): index for index, name in enumerate(names)}
        if target_gene not in position_of:
            break
        clinical = baseline_matrix(data, keep, train)
        wrecked = matrix.copy()
        shuffle = np.random.default_rng(SEED + 31337 + fold_number).permutation(y.size)
        wrecked[:, position_of[target_gene]] = matrix[shuffle,
                                                      position_of[target_gene]]

        def predict_on(block_columns):
            picked = [position_of[name] for name in block_columns
                      if name in position_of]
            block = (np.column_stack([clinical, wrecked[:, picked]])
                     if picked else clinical)
            return ridge_prediction(block[train], y[train], block[test],
                                    inner_folds=INNER_FOLDS,
                                    seed=SEED + fold_number)[0]

        destroyed["blind_errors"][repeat_of[fold_number], test] = np.abs(
            predict_on([target_gene]) - y[test])
        strength = _spearman_strength(wrecked[train], y[train])
        chosen = str(names[int(np.argmax(strength))])
        destroyed["univariate_errors"][repeat_of[fold_number], test] = np.abs(
            predict_on([chosen]) - y[test])
        destroyed["univariate_still_picked_the_dead_gene"] += int(
            chosen == target_gene)
        destroyed["folds_examined"] += 1
    examined = destroyed["folds_examined"]
    summary["destroyed_signal_diagnostic"] = {
        "gene": target_gene,
        "folds_examined": examined,
        "mechanistic_arm_k1_mae": (
            float(np.nanmean(per_donor(destroyed["blind_errors"])))
            if examined else None),
        "univariate_arm_k1_mae": (
            float(np.nanmean(per_donor(destroyed["univariate_errors"])))
            if examined else None),
        # With no fold examined these counters are not zero findings, they are
        # no findings. Reporting "0 of 10 folds fooled" for a test that never
        # ran would read as a test the arm passed.
        "folds_in_which_the_univariate_arm_could_have_been_fooled":
            (len(donor_folds) if examined else None),
        "univariate_still_picked_the_dead_gene":
            (destroyed["univariate_still_picked_the_dead_gene"]
             if examined else None),
        "how_to_read": (
            "The target gene's donor order is permuted, so it carries no "
            "information while keeping its place in the mechanistic ordering. "
            "The mechanistic arm must keep choosing it and must therefore lose "
            "its advantage; an arm that reads the endpoint should choose "
            "something else."),
        "what_the_number_alone_does_not_tell_you": (
            "A lower error for the arm that reads the endpoint is compatible "
            "with three incompatible situations, and this cohort contains all "
            "three: a non-human duplicate of the destroyed transcript, an "
            "admissible paralogue whose probe plausibly measures the same "
            "transcript, and a different gene whose association may be real. "
            "The replacement must be named before the number is read, which "
            "scripts/name_destroyed_signal_replacement.py does. See protocol "
            "section 5.5 as amended in 1.4.7."),
    }
    if not examined:
        summary["destroyed_signal_diagnostic"]["not_assessable"] = (
            f"{target_gene} has no feature in this universe, so its donor "
            "order cannot be permuted and this diagnostic did not run for this "
            "endpoint. Its counters are null rather than zero: nothing was "
            "examined, so nothing was resisted.")

    # ------------------------------------------------------------------
    # Contrasts.
    #
    # There is deliberately no "best statistical comparator" arm. Picking the
    # statistical method with the lowest out-of-fold error and then comparing
    # against it uses the test folds twice: the winner of three draws beats any
    # single method by construction, so the comparator looks better than any of
    # them really is. Each method is therefore contrasted separately, and the
    # frozen non-inferiority margin must hold against EVERY one of them -- a
    # strictly harder requirement than beating a test-chosen winner, and one
    # that needs no extra computation.
    # ------------------------------------------------------------------
    STATISTICAL = ("univariate_spearman", "lasso_entry_order", "xgboost_gain")
    contrasts = []
    for k in PANEL_SIZES:
        for offset, comparator in enumerate(STATISTICAL):
            contrasts.append(paired_contrast(
                donor_errors["netpk_target_seeded"][k], donor_errors[comparator][k],
                seed=SEED + 1000 * offset + k,
                label=f"netpk_target_minus_{comparator}_k{k}"))
        contrasts.append(paired_contrast(
            donor_errors["netpk_target_seeded"][k],
            donor_errors["netpk_wrong_mechanism"][k],
            seed=SEED + 4000 + k, label=f"netpk_target_minus_wrong_mechanism_k{k}"))
        contrasts.append(paired_contrast(
            donor_errors["netpk_target_seeded"][k], donor_errors["random_panel"][k],
            seed=SEED + 5000 + k, label=f"netpk_target_minus_random_mean_k{k}"))
        contrasts.append(paired_contrast(
            donor_errors["netpk_target_seeded"][k], fixed_errors["target_gene_only"],
            seed=SEED + 6000 + k, label=f"netpk_target_minus_target_only_k{k}"))
        contrasts.append(paired_contrast(
            donor_errors["netpk_pathway_seeded"][k],
            donor_errors["netpk_target_seeded"][k],
            seed=SEED + 7000 + k, label=f"netpk_pathway_minus_target_seeded_k{k}"))
    contrasts.append(paired_contrast(
        fixed_errors["target_gene_only"], fixed_errors["baseline_age_sex"],
        seed=SEED + 8000, label="target_only_minus_baseline"))
    contrasts.append(paired_contrast(
        fixed_errors["no_selection"], fixed_errors["baseline_age_sex"],
        seed=SEED + 8001, label="no_selection_minus_baseline"))
    summary["contrasts"] = contrasts

    summary["non_inferiority_against_every_statistical_arm"] = {
        str(k): {
            "comparators": list(STATISTICAL),
            "within_margin_against_all": all(
                next(c for c in contrasts
                     if c["label"] == f"netpk_target_minus_{comparator}_k{k}"
                     )["non_inferior_within_margin"]
                for comparator in STATISTICAL),
            "worst_relative_upper_bound": max(
                next(c for c in contrasts
                     if c["label"] == f"netpk_target_minus_{comparator}_k{k}"
                     )["relative_ci_high"]
                for comparator in STATISTICAL),
        } for k in PANEL_SIZES}
    summary["why_there_is_no_best_statistical_arm"] = (
        "Choosing the statistical comparator by its out-of-fold error would use "
        "the test folds twice and credit the winner of three draws. Each method "
        "is contrasted separately and the margin must hold against all three.")
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
    parser.add_argument("--probe-pool", default=None,
                        help="A qualified-probe-pool file; restricts every "
                             "arm to the same analytically admissible "
                             "probes.")
    parser.add_argument("--pool-name", default="A_qualified_human",
                        help="Which admissible set inside the pool file "
                             "(default: %(default)s).")
    parser.add_argument("--restrict-to-network", action="store_true",
                        help="Sensitivity analysis: let every arm select "
                             "only from genes the public network contains, "
                             "so ordering quality is not confounded with "
                             "pool size. Does not replace the primary run.")
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
    if args.probe_pool:
        universe = restrict_universe_to_qualified_pool(
            universe, data, args.probe_pool, args.pool_name)
        log(f"QUALIFIED POOL {args.pool_name}: "
            f"{universe['qualified_pool']['symbols_after']} symbols "
            f"(was {universe['qualified_pool']['symbols_before']}), "
            f"{universe['qualified_pool']['probes_dropped']} probes dropped, "
            "no technical control present")
    if args.restrict_to_network:
        universe = restrict_universe_to_network(universe, graph)
        log(f"SENSITIVITY: restricted to the network -> "
            f"{universe['candidates_after_restriction']} candidates "
            f"(was {universe['candidates_before_restriction']})")

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
        "protocol_version": "1.4.7",
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
            # Derived from the kept candidates, not from the unrestricted
            # universe: after a pool restriction some genes on the network are
            # no longer candidates, and this field is exactly the mechanistic
            # ordering's reach.
            "candidate_genes_on_network": len(
                set(universe["by_symbol"]) & set(universe["genes_on_network"])),
            "decoy_block": N_DECOYS,
            "restricted_to_network": bool(args.restrict_to_network),
            "qualified_pool": universe.get("qualified_pool"),
            "pool_asymmetry_note": (
                "In the primary run every arm selects the same NUMBER of genes "
                "but not from the same POOL: the mechanistic ordering reaches "
                "only genes present in the public network, while a statistical "
                "selector ranks the whole array, technical control probes "
                "included. The restricted run separates ordering quality from "
                "pool size."),
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
