#!/usr/bin/env python3
"""Research-only paired validation on actual NCI-60 doxorubicin response.

No individual patient data or simulated pharmacological endpoints. No
clinical PK/PD or causal mechanistic assertions. All modeling and preprocessing
use ONLY each outer training set. Every model uses the same splits.
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
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from xgboost import XGBRegressor

from omics2pk.inductive_sweet import (
    InductiveNetworkError,
    build_outer_split_features,
)

COMPOUND = "Doxorubicin"
SEED = 20261016
ALPHA = 10.0
MAX_PCA = 3
XGB_KW = dict(n_estimators=90, max_depth=2, learning_rate=0.05,
              min_child_weight=4, subsample=0.8, colsample_bytree=0.8,
              reg_lambda=10.0, n_jobs=1, random_state=SEED, objective="reg:squarederror")


def load_public(source_dir: Path):
    """Fail rather than join the wrong cell lines, endpoints or fields."""
    with (source_dir/"expression.csv").open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        if fields[:2] != ["subject_id", "cancer_class"] or len(fields) != 50:
            raise ValueError("Unexpected public ABC transporter expression schema")
        samples = list(reader)
    with (source_dir/"measured_doxorubicin_response.csv").open(
            newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != ["subject_id", "y"]:
            raise ValueError("Unexpected pre-registered measured GI50 endpoint")
        outcomes = list(reader)
    genes = tuple(fields[2:])
    ids = tuple(r["subject_id"] for r in samples)
    if (len(ids)!=60 or len(set(ids))!=60
            or len(genes)!=48 or len(set(genes))!=48
            or len(outcomes)!=60
            or [r["subject_id"] for r in outcomes]!=list(ids)):
        raise ValueError("Missing, duplicated or misaligned public cell lines")
    def number(v):
        return np.nan if v == "NA" else float(v)
    x = np.asarray([[number(r[g]) for g in genes] for r in samples])
    y = np.asarray([float(r["y"]) for r in outcomes])
    if not np.isfinite(y).all() or np.isinf(x).any():
        raise ValueError("Non-finite original observed outcomes or expressions")
    if len(set(r["cancer_class"] for r in samples)) < 2:
        raise ValueError("Cancer-class metadata unexpectedly degenerate")
    classes = tuple(r["cancer_class"] for r in samples)
    return ids, genes, x, y, classes


def model_predict(x_train, y_train, x_test, kind: str):
    """Fixed a priori regression settings; no outcome-driven hyperparameter search."""
    if kind == "xgb":
        model = make_pipeline(StandardScaler(), XGBRegressor(**XGB_KW))
    elif kind == "pca_ridge":
        model = make_pipeline(StandardScaler(),
                              PCA(n_components=min(MAX_PCA, x_train.shape[0]-1, x_train.shape[1])),
                              Ridge(alpha=ALPHA))
    else:
        model = make_pipeline(StandardScaler(), Ridge(alpha=ALPHA))
    model.fit(x_train, y_train)
    return model.predict(x_test)


def network_split(x, ids, genes, train, test):
    """Fit network transforms with training-only imputation and gene selection."""
    xi = SimpleImputer(strategy="median", keep_empty_features=True)
    observed_train = xi.fit_transform(x[train])
    observed_test = xi.transform(x[test])
    if not np.isfinite(observed_train).all() or not np.isfinite(observed_test).all():
        raise ValueError("Entirely missing source gene: cannot build network")
    # Variant selection by TRAIN-ONLY gene variance, independent of drug target.
    select = np.std(observed_train, axis=0) > 1e-12
    if select.sum() < 3:
        raise ValueError("Too few nonconstant ABC genes in training")
    use = tuple(g for g, chosen in zip(genes, select) if chosen)
    a = np.vstack([observed_train[:, select], observed_test[:, select]]).T
    order = [ids[i] for i in list(train) + list(test)]
    split = build_outer_split_features(
        a, genes=use, subject_ids=order,
        train_ids=[ids[i] for i in train],
        test_ids=[ids[i] for i in test],
        tissue="NCI60_cell_lines", folds=0, k=0.1)
    ntr = np.asarray([[r["features"][g] for g in use] for r in split["training"]])
    nte = np.asarray([[r["features"][g] for g in use] for r in split["test"]])
    # The exact same selected expression features enter the raw comparator.
    return observed_train[:, select], observed_test[:, select], ntr, nte, split


def evaluate(x, y, ids, genes, repeats=1, classes=None):
    """Fixed, pre-registered grouped-by-cell-line outer CV. No PK data."""
    if x.shape != (60, 48) or len(y)!=60 or len(ids)!=60:
        raise ValueError("Require published 60-by-48 NCI-60 dataset")
    if repeats < 1 or repeats > 5:
        raise ValueError("Allowed repeat count 1..5")
    if classes is not None and len(classes) != len(ids):
        raise ValueError("Cancer-class length must match biological subjects")
    model_names = ("mean", "raw_ridge", "pca_ridge", "network_ridge",
                   "raw_plus_network_ridge", "negative_shuffled_network",
                   "M1_raw_xgboost", "M1bis_network_xgboost", "M2_combined_xgboost")
    if classes is not None:
        model_names += ("tumor_class_ridge", "tumor_class_raw_ridge",
                        "tumor_class_network_ridge", "tumor_class_raw_network_ridge")
    predictions = {name: [] for name in model_names}
    observed = []
    clipping = []
    n_genes_by_fold = []
    split_fingerprints = []
    for rep in range(repeats):
        cv = KFold(n_splits=5, shuffle=True, random_state=SEED + rep)
        for fold, (train, test) in enumerate(cv.split(x)):
            raw_train, raw_test, net_train, net_test, split = network_split(
                x, ids, genes, train, test)
            if len({r["reference_fingerprint"] for r in split["test"]}) != 1:
                raise RuntimeError("Unfrozen reference across evaluation subjects")
            n_genes_by_fold.append(len(raw_train[0]))
            split_fingerprints.append(split["outer_reference_fingerprint"])
            clipping.append(sum(v["weight_clamped_outside_reference"] for v in split["test"]))
            ytr, yte = y[train], y[test]
            if not np.isfinite(ytr).all() or not np.isfinite(yte).all():
                raise ValueError("Missing measured response")
            rng = np.random.default_rng(SEED + rep*10+fold)
            net_shuffled_train = net_train[rng.permutation(len(train))]
            raw_plus = np.hstack([raw_train, net_train])
            test_plus = np.hstack([raw_test, net_test])
            pred = {
                "mean": np.full(len(test), ytr.mean()),
                "raw_ridge": model_predict(raw_train, ytr, raw_test, "ridge"),
                "pca_ridge": model_predict(raw_train, ytr, raw_test, "pca_ridge"),
                "network_ridge": model_predict(net_train, ytr, net_test, "ridge"),
                "raw_plus_network_ridge": model_predict(raw_plus, ytr, test_plus, "ridge"),
                "negative_shuffled_network": model_predict(
                    net_shuffled_train, ytr, net_test, "ridge"),
                "M1_raw_xgboost": model_predict(raw_train, ytr, raw_test, "xgb"),
                "M1bis_network_xgboost": model_predict(net_train, ytr, net_test, "xgb"),
                "M2_combined_xgboost": model_predict(raw_plus, ytr, test_plus, "xgb"),
            }
            if classes is not None:
                labels = np.array(classes, dtype=str)
                encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
                cancer_train = encoder.fit_transform(labels[train, None])
                cancer_test = encoder.transform(labels[test, None])
                pred["tumor_class_ridge"] = model_predict(
                    cancer_train, ytr, cancer_test, "ridge")
                pred["tumor_class_raw_ridge"] = model_predict(
                    np.hstack([cancer_train, raw_train]), ytr,
                    np.hstack([cancer_test, raw_test]), "ridge")
                pred["tumor_class_network_ridge"] = model_predict(
                    np.hstack([cancer_train, net_train]), ytr,
                    np.hstack([cancer_test, net_test]), "ridge")
                pred["tumor_class_raw_network_ridge"] = model_predict(
                    np.hstack([cancer_train, raw_train, net_train]), ytr,
                    np.hstack([cancer_test, raw_test, net_test]), "ridge")
            observed.extend(yte.tolist())
            for name in model_names:
                if len(pred[name])!=len(yte) or not np.isfinite(pred[name]).all():
                    raise RuntimeError(f"Invalid predictions: {name}")
                predictions[name].extend(np.asarray(pred[name]).tolist())
    actual = np.asarray(observed)
    scores = {}
    base_mae = mean_absolute_error(actual, predictions["mean"])
    for name in model_names:
        guessed = np.asarray(predictions[name])
        scores[name] = {
            "mae": round(float(mean_absolute_error(actual, guessed)), 6),
            "rmse": round(float(math.sqrt(mean_squared_error(actual, guessed))), 6),
            "mae_relative_to_mean_baseline": round(
                float(mean_absolute_error(actual, guessed) / base_mae), 4),
        }
    data_fingerprint = hashlib.sha256(
        np.nan_to_num(x, nan=-99999).tobytes() + y.tobytes()).hexdigest()
    return {
        "dataset": "mixOmics multidrug: NCI-60, ABC RT-qPCR and observed drug response",
        "source_commit": "9f9fd66b69105fdb618b54588b638b567e85c0d6",
        "drug_pre_registered": COMPOUND,
        "endpoint": "original processed drug GI50 activity (not PK exposure)",
        "biological_unit": "independent cancer cell lines",
        "tumor_class_adjusted_secondary_analysis": classes is not None,
        "repeats_are_correlated_within_cell_line": repeats > 1,
        "n_cell_lines": 60, "n_transporter_genes": 48,
        "outer_splits": 5, "repeats": repeats, "seed": SEED,
        "outer_prediction_count": len(actual),
        "ridge_alpha_pre_fixed": ALPHA,
        "xgboost_pre_fixed": XGB_KW,
        "target_access_in_network": False,
        "no_pk_pd_or_clinical_validation": True,
        "negative_control": "training network rows shuffled, test network unchanged",
        "genes_per_fold": n_genes_by_fold,
        "test_weight_clippings_by_fold": clipping,
        "reference_fingerprints": split_fingerprints,
        "source_matrix_sha256": data_fingerprint,
        "metrics": scores,
        "warning": "POC only, 60 cell lines; correlations, cancer-line confounding, and "
                   "nested tuning/external validation remain to be qualified. "
                   "Do not claim causal, mechanistic or PopPK gains."
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--repeats", type=int, default=1)
    args = p.parse_args()
    ids, genes, x, y, classes = load_public(args.data)
    results = evaluate(x, y, ids, genes, repeats=args.repeats, classes=classes)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    # Explicitly do not print or upload any individual expression or response.
    print(json.dumps({"drug":COMPOUND,"n":len(y),"models":results["metrics"],
                      "clipped_weights":sum(results["test_weight_clippings_by_fold"]),
                      "result":"RESEARCH_ONLY_OBSERVED_GI50_PILOT"}, indent=2))


if __name__ == "__main__":
    main()
