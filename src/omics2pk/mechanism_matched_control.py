"""Measured-gene function-prior falsification against matched source-feature controls.

This is *cell-line* doxorubicin response (GI50), NOT clinical human PK,
specific drug transport, metabolic enzyme activity, or an external cohort.
Original gene labels and independent literature motivate the fixed 7-gene
panel; randomized matched RNA expression panels are NEGATIVE CONTROLS.
All matching and model standardization occur on outer training lines only.
"""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer

DRUG_GENES = ("TOP2A", "TOP2B", "ABCB1", "ABCC1", "ABCG2", "CBR1", "CBR3")
SEED = 20261104
ALPHA = 10.0


class FunctionalControlError(ValueError):
    pass


def matched_source_control_panels(train, fixed_indices, *, n_controls=64,
                                  pool_size=40, seed=SEED):
    """Train-only detection and dispersion matched controls (not STRING degree).

    Each target gene receives a matched RNA candidate from the same *measured*
    feature universe. No heldout outcome, test feature or test lineage enters
    this choice. Matched assay properties are not evidence of biological
    equivalence, and transcript levels are not functional enzyme activity.
    """
    train = np.asarray(train, dtype=float)
    if train.ndim != 2 or len(train) < 15 or train.shape[1] < 20:
        raise FunctionalControlError("Measured original training expression required")
    if np.isinf(train).any():
        raise FunctionalControlError("Infinite original RNA value")
    fixed_indices = tuple(fixed_indices)
    if (not fixed_indices or len(set(fixed_indices)) != len(fixed_indices)
            or any(not isinstance(j, (int, np.integer)) or isinstance(j, (bool, np.bool_))
                   or not 0 <= j < train.shape[1] for j in fixed_indices)):
        raise FunctionalControlError("Distinct source genes with valid source indexes required")
    if (not isinstance(n_controls, int) or isinstance(n_controls, bool)
            or n_controls < 1 or not isinstance(pool_size, int)
            or isinstance(pool_size, bool) or pool_size < 2):
        raise FunctionalControlError("Control number and candidate pool must be fixed positive integers")
    observed = np.isfinite(train)
    coverage = observed.mean(axis=0)
    good = np.flatnonzero(coverage >= 0.65)
    if not set(fixed_indices).issubset(set(good)):
        raise FunctionalControlError("Prior gene missing in >=35% of TRAIN samples")
    vals = train[:, good]
    med = np.nanmedian(vals, axis=0)
    if not np.isfinite(med).all():
        raise FunctionalControlError("Unqualified train-only molecular medians")
    imputed = np.where(np.isfinite(vals), vals, med[None, :])
    variance = np.var(imputed, axis=0)
    viable = good[variance > 1e-10]
    if len(viable) < len(fixed_indices) + pool_size:
        raise FunctionalControlError("Insufficient independently measured gene controls")
    if not set(fixed_indices).issubset(set(viable)):
        raise FunctionalControlError("Prior source gene lacks measurable train-only variation")
    all_var = dict(zip(good.tolist(), variance.tolist()))
    candidates = np.array([i for i in viable if i not in fixed_indices], dtype=int)
    logvars = np.log(np.maximum([all_var[i] for i in candidates], 1e-10))
    windows = []
    for j in fixed_indices:
        distance = np.abs(logvars - np.log(all_var[j])) + 3.0 * np.abs(coverage[candidates] - coverage[j])
        windows.append(candidates[np.argsort(distance, kind="stable")[:pool_size]])
    rng = np.random.default_rng(seed)
    panels = []
    for _ in range(n_controls):
        chosen = set()
        panel = []
        for window in windows:
            available = [int(j) for j in window if int(j) not in chosen]
            if not available:
                raise FunctionalControlError("Source-matched random control pool exhausted")
            idx = int(rng.choice(available))
            chosen.add(idx)
            panel.append(idx)
        panels.append(tuple(panel))
    return panels


def _ridge_predict(train, test, y, lineage_train, lineage_test,
                   molecular, indices):
    encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    class_train = encoder.fit_transform(np.asarray(lineage_train)[:, None])
    class_test = encoder.transform(np.asarray(lineage_test)[:, None])
    if indices:
        x = molecular[np.asarray(train)][:, list(indices)]
        xt = molecular[np.asarray(test)][:, list(indices)]
        imp = SimpleImputer(strategy="median", keep_empty_features=True)
        a = imp.fit_transform(x)
        b = imp.transform(xt)
        if a.shape[1] != len(indices) or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise FunctionalControlError("Incomplete train-only gene normalization")
        sc = StandardScaler().fit(a)
        class_train = np.column_stack((class_train, sc.transform(a)))
        class_test = np.column_stack((class_test, sc.transform(b)))
    fitted = Ridge(alpha=ALPHA).fit(class_train, y[np.asarray(train)])
    pred = np.asarray(fitted.predict(class_test), dtype=float)
    if not np.isfinite(pred).all():
        raise FunctionalControlError("Invalid heldout measured GI50 prediction")
    return pred


def benchmark_prior_against_controls(source_ids, lineage, response, gene_names,
                                     original_rna, *, repeats=3, n_controls=64):
    """Fixed source-lineage heldout predictions on n biologically distinct lines."""
    ids = tuple(source_ids)
    names = tuple(gene_names)
    x = np.asarray(original_rna, dtype=float)
    y = np.asarray(response, dtype=float)
    cls = np.asarray(lineage)
    n = len(ids)
    if (n < 25 or len(set(ids)) != n or len(set(names)) != len(names)
            or x.shape != (n, len(names)) or len(cls) != n
            or y.shape != (n,) or not np.isfinite(y).all()
            or np.isinf(x).any() or any(not c for c in cls)):
        raise FunctionalControlError("Matched original same cell-line RNA/response source required")
    if not isinstance(repeats, int) or isinstance(repeats, bool) or not 1 <= repeats <= 5:
        raise FunctionalControlError("Prespecified 1-5 repeat count required")
    if not isinstance(n_controls, int) or isinstance(n_controls, bool) or not 2 <= n_controls <= 300:
        raise FunctionalControlError("Prespecified control count required")
    if not all(names.count(g) == 1 for g in DRUG_GENES):
        raise FunctionalControlError("ABSTAIN original source missing a prior drug gene")
    fixed = tuple(names.index(g) for g in DRUG_GENES)
    baseline = np.full((repeats, n), np.nan)
    drug_prior = np.full((repeats, n), np.nan)
    controls = np.full((n_controls, repeats, n), np.nan)
    for repeat in range(repeats):
        folds = KFold(n_splits=5, shuffle=True, random_state=SEED + repeat)
        for fold, (tr, te) in enumerate(folds.split(np.arange(n))):
            panels = matched_source_control_panels(
                x[tr], fixed, n_controls=n_controls,
                seed=SEED + repeat * 100 + fold)
            baseline[repeat, te] = _ridge_predict(
                tr, te, y, cls[tr], cls[te], x, ())
            drug_prior[repeat, te] = _ridge_predict(
                tr, te, y, cls[tr], cls[te], x, fixed)
            for j, panel in enumerate(panels):
                controls[j, repeat, te] = _ridge_predict(
                    tr, te, y, cls[tr], cls[te], x, panel)
    if (not np.isfinite(baseline).all() or not np.isfinite(drug_prior).all()
            or not np.isfinite(controls).all()):
        raise FunctionalControlError("Incomplete source-cell-line heldout coverage")
    mae_baseline = np.abs(baseline-y[None, :]).mean()
    mae_prior = np.abs(drug_prior-y[None, :]).mean()
    ctrl_mae = np.abs(controls-y[None, None, :]).mean(axis=(1, 2))
    # Controls use matched random panels in every fold: null *procedure*,
    # not an exact permutation test of fixed gene labels.
    q05,q50,q95 = np.quantile(ctrl_mae,[.05,.5,.95])
    return {
        "status": "SOURCE_MATCHED_RNA_PRIOR_VS_TRAIN_MATCHED_RANDOM_GENE_CONTROLS_PD_ONLY",
        "n_source_cell_lines": n,
        "n_correlated_repeat_test_records": n * repeats,
        "n_control_pipelines": n_controls,
        "predeclared_doxorubicin_gene_symbols": list(DRUG_GENES),
        "n_predeclared_genes_measured": len(fixed),
        "assay_matched_on": "training-fold-only RNA detection rate and log variance; NOT external STRING network degree",
        "observed_lineage_only_MAE": round(float(mae_baseline),6),
        "observed_lineage_plus_prior_RNA_MAE": round(float(mae_prior),6),
        "delta_prior_minus_lineage_MAE_negative_better": round(float(mae_prior-mae_baseline),6),
        "random_gene_control_MAE_5_50_95_percentiles": [round(float(v),6) for v in (q05,q50,q95)],
        "fraction_matched_random_control_pipelines_with_lower_MAE_than_prior":
            round(float(np.mean(ctrl_mae < mae_prior)),6),
        "same_original_cell_line_heldout_test_for_each_arm": True,
        "omic_features_and_preprocessing_selected_training_folds_only": True,
        "original_gene_expression_is_not_drug_transporter_or_enzyme_activity": True,
        "original_cell_line_GI50_is_not_human_PK": True,
        "no_original_cell_line_omics_or_predictions_returned": True,
        "mechanism_function_confirmed": False,
        "external_biological_replication": False,
        "confirmatory_significance_claimed": False,
    }
