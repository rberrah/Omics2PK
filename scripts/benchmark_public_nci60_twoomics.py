#!/usr/bin/env python3
"""Prespecified exploratory NCI-60 cross-study GI50 and RNA/metabolome comparison.

DO NOT interpret cell-line-name matched measurements as contemporaneous or
same-specimen multiomics. Measured in-vitro GI50 is not clinical PK/PD.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from xgboost import XGBRegressor

SEED = 20261016
RIDGE_ALPHA = 10.0
RNA_CAP = 300
METAB_CAP = 120
COMPONENTS = 3
XGB_OPTIONS = dict(n_estimators=90, max_depth=2, learning_rate=0.05,
                   min_child_weight=4, subsample=0.8, colsample_bytree=0.8,
                   reg_lambda=10.0, n_jobs=1, random_state=SEED,
                   objective="reg:squarederror")

MODEL_NAMES = (
    "fold_mean", "tumor_only",
    "rna_pca3", "metab_pca3", "rna_metab_pca3",
    "tumor_rna_pca3", "tumor_metab_pca3", "tumor_rna_metab_pca3",
    "tumor_rna_raw300", "tumor_metab_raw120", "tumor_rna_metab_raw",
    "tumor_rna_shuffled_metab_negative",
    "rna_xgboost", "rna_metab_xgboost",
    "tumor_rna_xgboost", "tumor_rna_metab_xgboost",
)


def read_table(path: Path, numeric=True):
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = csv.reader(f)
        header = next(rows)
        if not header or header[0] != "subject_id":
            raise ValueError(f"Unexpected subject_id field in {path.name}")
        fields = header[1:]
        if len(set(fields)) != len(fields) or not fields:
            raise ValueError(f"Duplicate or empty feature schema: {path.name}")
        ids, values = [], []
        for row in rows:
            if len(row) != len(header):
                raise ValueError(f"Wrong field count in {path.name}")
            ids.append(row[0])
            if numeric:
                values.append([float(v) if v not in ("NA", "", "NaN") else np.nan
                               for v in row[1:]])
            else:
                values.append(row[1:])
    if not ids or len(ids) != len(set(ids)) or any(not x for x in ids):
        raise ValueError(f"Invalid cell-line identity in {path.name}")
    return tuple(ids), tuple(fields), np.asarray(values, dtype=float if numeric else str)


def load_dataset(folder: Path):
    with (folder / "response.csv").open(newline="", encoding="utf-8") as f:
        rows = csv.DictReader(f)
        if rows.fieldnames != ["subject_id", "cancer_class", "y"]:
            raise ValueError("Unexpected response/cancer-class schema")
        record = list(rows)
    ids = tuple(r["subject_id"] for r in record)
    classes = tuple(r["cancer_class"] for r in record)
    y = np.array([float(r["y"]) for r in record], dtype=float)
    if (len(ids) != 48 or len(set(ids)) != 48 or
            not np.isfinite(y).all() or any(not c for c in classes)):
        raise ValueError("Expected 48 unique, measured, lineage-qualified cell lines")
    rid, rn, rna = read_table(folder / "rna.csv")
    mid, mn, met = read_table(folder / "metabolite.csv")
    if rid != ids or mid != ids:
        raise ValueError("RNA, metabolite, and response rows are misaligned")
    if (rna.shape[1] < RNA_CAP or met.shape[1] < METAB_CAP or
            np.isinf(rna).any() or np.isinf(met).any()):
        raise ValueError("Unexpected source modality feature number or infinity")
    return ids, classes, y, rn, mn, rna, met


def fit_view(training, testing, max_features, seed=SEED):
    """Complete-case threshold, median imputation, variance ranking and PCA on TRAIN only.

    Return low-dimensional PCs and matched preprocessed raw features.
    """
    if training.ndim != 2 or testing.ndim != 2 or training.shape[1] != testing.shape[1]:
        raise ValueError("Incompatible modality matrices")
    if np.isinf(training).any() or np.isinf(testing).any():
        raise ValueError("Infinite feature value")
    available = np.mean(np.isfinite(training), axis=0) >= 0.5
    if not np.any(available):
        raise ValueError("No train-observed features")
    indices = np.flatnonzero(available)
    a = training[:, indices]
    b = testing[:, indices]
    med = np.nanmedian(a, axis=0)
    if not np.isfinite(med).all():
        raise ValueError("Nonfinite TRAIN medians")
    a = np.where(np.isfinite(a), a, med[None, :])
    b = np.where(np.isfinite(b), b, med[None, :])
    variances = np.var(a, axis=0)
    order = np.argsort(-variances, kind="stable")
    selected = order[variances[order] > 1e-12][:max_features]
    if len(selected) < COMPONENTS:
        raise ValueError("Too few nonconstant training features")
    picked_indices = indices[selected]
    a, b = a[:, selected], b[:, selected]
    scale = StandardScaler().fit(a)
    at, bt = scale.transform(a), scale.transform(b)
    pca = PCA(n_components=COMPONENTS, svd_solver="randomized", random_state=seed).fit(at)
    pc_train, pc_test = pca.transform(at), pca.transform(bt)
    if not all(np.isfinite(z).all() for z in (pc_train, pc_test, at, bt)):
        raise ValueError("Nonfinite final fold features")
    return (pc_train, pc_test, at, bt, picked_indices)


def predict(a, b, y, method="ridge"):
    norm = StandardScaler().fit(a)
    at, bt = norm.transform(a), norm.transform(b)
    if method == "xgb":
        model = XGBRegressor(**XGB_OPTIONS)
    else:
        model = Ridge(alpha=RIDGE_ALPHA)
    model.fit(at, y)
    guess = np.asarray(model.predict(bt), dtype=float)
    if not np.isfinite(guess).all():
        raise ValueError("Nonfinite model predictions")
    return guess


def evaluate(ids, classes, y, rna, met, repeats=3):
    if (len(ids) != 48 or len(set(ids)) != len(ids) or len(classes) != 48
            or y.shape != (48,) or rna.shape[0] != 48 or met.shape[0] != 48):
        raise ValueError("Cross-study cohort must have 48 unique cell lines")
    if repeats not in (1, 3):
        raise ValueError("Only pre-fixed 1 or 3 repeats allowed")
    # Every held-out subject appears exactly once in each 5-fold repetition.
    preds = {name: np.full((repeats, 48), np.nan) for name in MODEL_NAMES}
    n_selected = {"rna": [], "metabolite": []}
    for rep in range(repeats):
        cv = KFold(n_splits=5, shuffle=True, random_state=SEED + rep)
        for fold, (train, test) in enumerate(cv.split(np.arange(48))):
            yt = y[train]
            pr, qr, xr, xrt, select_r = fit_view(rna[train], rna[test], RNA_CAP)
            pm, qm, xm, xmt, select_m = fit_view(met[train], met[test], METAB_CAP)
            n_selected["rna"].append(len(select_r))
            n_selected["metabolite"].append(len(select_m))
            encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
            ct = encoder.fit_transform(np.asarray(classes)[train, None])
            cq = encoder.transform(np.asarray(classes)[test, None])
            stack = lambda *blocks: np.concatenate(blocks, axis=1)
            rng = np.random.default_rng(SEED + rep * 10 + fold)
            # Shuffle only TRAIN metabolomics, keeping the test fixed.
            shuffled = pm[rng.permutation(len(train))]
            options = {
                "fold_mean": np.full(len(test), yt.mean()),
                "tumor_only": predict(ct, cq, yt),
                "rna_pca3": predict(pr, qr, yt),
                "metab_pca3": predict(pm, qm, yt),
                "rna_metab_pca3": predict(stack(pr, pm), stack(qr, qm), yt),
                "tumor_rna_pca3": predict(stack(ct, pr), stack(cq, qr), yt),
                "tumor_metab_pca3": predict(stack(ct, pm), stack(cq, qm), yt),
                "tumor_rna_metab_pca3": predict(stack(ct, pr, pm), stack(cq, qr, qm), yt),
                "tumor_rna_raw300": predict(stack(ct, xr), stack(cq, xrt), yt),
                "tumor_metab_raw120": predict(stack(ct, xm), stack(cq, xmt), yt),
                "tumor_rna_metab_raw": predict(stack(ct, xr, xm), stack(cq, xrt, xmt), yt),
                "tumor_rna_shuffled_metab_negative": predict(
                    stack(ct, pr, shuffled), stack(cq, qr, qm), yt),
                "rna_xgboost": predict(pr, qr, yt, "xgb"),
                "rna_metab_xgboost": predict(stack(pr, pm), stack(qr, qm), yt, "xgb"),
                "tumor_rna_xgboost": predict(stack(ct, pr), stack(cq, qr), yt, "xgb"),
                "tumor_rna_metab_xgboost": predict(
                    stack(ct, pr, pm), stack(cq, qr, qm), yt, "xgb"),
            }
            for name in MODEL_NAMES:
                if len(options[name]) != len(test):
                    raise ValueError(f"Prediction count mismatch: {name}")
                preds[name][rep, test] = options[name]
        for name, matrix in preds.items():
            if not np.isfinite(matrix[rep]).all():
                raise RuntimeError(f"Missing OOF values for {name}, repetition {rep}")
    stats = {}
    tumor_error = np.abs(preds["tumor_only"] - y[None, :]).mean(axis=0)
    rng = np.random.default_rng(SEED)
    draws = rng.integers(0, 48, size=(1000, 48))
    for name, p in preds.items():
        errors = np.abs(p - y[None, :])
        per_line = errors.mean(axis=0)  # correlated repetitions averaged per biological cell line
        diff = per_line - tumor_error
        lo, hi = np.quantile(diff[draws].mean(axis=1), [0.025, 0.975])
        stats[name] = {
            "mae": round(float(errors.mean()), 6),
            "rmse": round(float(math.sqrt(np.mean((p - y[None, :])**2))), 6),
            "paired_delta_mae_vs_tumor_only": round(float(diff.mean()), 6),
            "descriptive_cell_line_bootstrap_95pct_delta": [round(float(lo), 6), round(float(hi), 6)]
        }
    # Outputs contain no row-level assay measurements, labels, or observed responses.
    return {
        "analysis": "exploratory cross-study cell-line identity-linked data; NOT contemporaneous paired omics",
        "target": "pre-registered doxorubicin measured NCI-60 GI50 (NOT human PK/clearance/AUC)",
        "sources": {
            "mixOmics": "mixOmics-org/mixOmics@9f9fd66b69105fdb618b54588b638b567e85c0d6",
            "IntLIM": "Mathelab/NCI60_GeneMetabolite_Data@4c140708c18c41910458414a626b1a51d3dd2a85",
        },
        "join": "48 exact punctuation/case-only unique keys, independently agreeing cancer-lineage labels; 12 excluded",
        "n_biological_cell_lines": 48, "n_repetitions_correlated": repeats, "n_splits_per_repetition": 5,
        "random_seed": SEED,
        "fixed_design": {"ridge_alpha": RIDGE_ALPHA, "rna_max_train_variance_features": RNA_CAP,
                         "metabolite_max_train_variance_features": METAB_CAP,
                         "feature_min_train_observed_fraction": 0.5,
                         "pca_components_per_view": COMPONENTS, "xgb": XGB_OPTIONS,
                         "all_preprocessing_is_fold_train_only": True},
        "selected_feature_count_per_outer_fold": n_selected,
        "uncertainty_warning": "Bootstrap resamples 48 cell lines, not repeat-fold observations; descriptive only.",
        "major_limitations": [
            "Public assays from different studies and culture batches are NOT same-sample multiomics",
            "Independent source genome-wide normalization/probe processing predates outer folds",
            "NCI-60 drug sensitivity is in-vitro PD, not PopPK or clinical precision dosing",
            "Only one prespecified compound, small n and no external cohort",
            "Technical replicates collapsed upstream; IntLIM metabolite preprocessing is upstream",
            "Nonlinear and full-omics results are exploratory, not a mechanism or clinical claim",
        ],
        "metrics": stats
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--repeats", type=int, choices=[1, 3], default=3)
    args = parser.parse_args()
    ids, classes, y, _, _, rna, met = load_dataset(args.data)
    result = evaluate(ids, classes, y, rna, met, args.repeats)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"analysis":result["analysis"],"n":48,"repeats":args.repeats,
                      "models":{k: {"mae":v["mae"],"delta_vs_tumor":v["paired_delta_mae_vs_tumor_only"]}
                                for k,v in result["metrics"].items()}}, indent=2))


if __name__ == "__main__":
    main()
