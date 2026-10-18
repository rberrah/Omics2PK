"""Leakage-resistant comparison of transcriptomic representations for measured PK.

Research-only EXPLORATORY benchmark. Inputs are verified one-row-per-biological-
person transcriptomes and original predose-to-dose timestamps. The original
source join/units/permissions must be audited BEFORE calling this function.
No original source records, predictions or identifiers are returned.
The SWEET-inspired network is NOT official SWEET nor a validated PK covariate.
"""
from __future__ import annotations

from datetime import datetime
from hashlib import sha256
from math import isfinite

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from .inductive_sweet import build_outer_split_features, InductiveNetworkError


class RepresentationBenchmarkError(ValueError):
    pass


def _matrix(value, shape=None):
    try:
        arr = np.asarray(value, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise RepresentationBenchmarkError("Numerical source matrix is invalid") from exc
    if arr.ndim != 2 or not np.isfinite(arr).all() or (shape is not None and arr.shape != shape):
        raise RepresentationBenchmarkError("Incorrect shape or nonfinite source measurement")
    return arr


def _date(value):
    if not isinstance(value, str):
        raise RepresentationBenchmarkError("Original absolute timestamp required")
    try:
        dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise RepresentationBenchmarkError("Malformed absolute timestamp") from exc
    if dt.tzinfo is None or dt.utcoffset() is None:
        raise RepresentationBenchmarkError("Source timestamps must specify timezones")
    return dt


def _predict(train, test, y_train):
    scaler = StandardScaler().fit(train)
    reg = Ridge(alpha=1.0).fit(scaler.transform(train), np.log(y_train))
    pred = float(np.exp(reg.predict(scaler.transform(test))[0]))
    if not isfinite(pred) or pred <= 0:
        raise RepresentationBenchmarkError("Non-finite or invalid PK prediction")
    return pred


def compare_predose_rna_pk(
    *, subject_ids, expression, gene_ids, clinical, measured_pk,
    collected_at, first_dose_at, tissue, target_name, target_unit,
    source_kind, source_sha256, genotype_dosages=None,
    prespecified_gene_weights=None, max_lookback_days=30,
    outer_limit=60,
):
    """Same-person outer LOPO ablation, transformations always TRAIN-only.

    expression: gene-by-subject normalized individual RNA values, not DE group
    coefficients. genotype_dosages: optional subject-by-variant observed 0/1/2.
    prespecified_gene_weights: externally fixed gene-to-weight map, not
    weights selected after inspecting this cohort's PK or heldout scores.
    No cohort-level source deconvolution, normalization or imputation is done;
    the caller must demonstrate those did not include heldout subjects.
    """
    ids = tuple(subject_ids)
    genes = tuple(gene_ids)
    n = len(ids)
    if not 12 <= n <= outer_limit or len(set(ids)) != n or any(
        not isinstance(v, str) or not v.strip() for v in ids
    ):
        raise RepresentationBenchmarkError("Expected 12-60 unique biological persons")
    if len(genes) < 3 or len(genes) > 128 or len(set(genes)) != len(genes) or any(
        not isinstance(g, str) or not g.strip() for g in genes
    ):
        raise RepresentationBenchmarkError("Prespecified unique 3-128 transcript genes required")
    if not isinstance(tissue, str) or not tissue.strip() or not target_name or not target_unit:
        raise RepresentationBenchmarkError("Measured PK target, physical units and tissue required")
    if source_kind not in ("synthetic", "original_source_declared"):
        raise RepresentationBenchmarkError("Source classification must be explicit")
    if (not isinstance(source_sha256, str) or len(source_sha256) != 64
            or any(c not in "0123456789abcdef" for c in source_sha256)):
        raise RepresentationBenchmarkError("Original-source SHA256 declaration required")
    if not isinstance(max_lookback_days, (int, float)) or not 0 < max_lookback_days <= 365:
        raise RepresentationBenchmarkError("Invalid declared predose window")
    if len(collected_at) != n or len(first_dose_at) != n:
        raise RepresentationBenchmarkError("Timestamps not matched to biological subjects")
    for collected, dosed in zip(collected_at, first_dose_at):
        delta = (_date(dosed) - _date(collected)).total_seconds() / 86400
        if not 0 <= delta <= max_lookback_days:
            raise RepresentationBenchmarkError("Postdose or out-of-window omics is forbidden")
    X = _matrix(expression, (len(genes), n)).T
    C = _matrix(clinical)
    if C.shape[0] != n or C.shape[1] < 1:
        raise RepresentationBenchmarkError("Clinical covariates must be matched person-by-feature")
    y = np.asarray(measured_pk, dtype=float)
    if y.shape != (n,) or not np.isfinite(y).all() or np.any(y <= 0):
        raise RepresentationBenchmarkError("Original measured PK endpoint must be positive per person")
    G = None
    if genotype_dosages is not None:
        G = _matrix(genotype_dosages)
        if G.shape[0] != n or G.shape[1] < 1 or not np.all(G == np.floor(G)) or np.any((G < 0) | (G > 2)):
            raise RepresentationBenchmarkError("Genotypes must be authentic coded 0/1/2 per person")
    W = None
    if prespecified_gene_weights is not None:
        if (not isinstance(prespecified_gene_weights, dict)
                or not prespecified_gene_weights
                or not set(prespecified_gene_weights).issubset(genes)):
            raise RepresentationBenchmarkError("Biological score weights must map to prespecified genes")
        W = np.array([float(prespecified_gene_weights.get(g, 0.0)) for g in genes])
        if not np.isfinite(W).all() or not np.any(W):
            raise RepresentationBenchmarkError("Finite nonzero fixed gene-score weights required")
        W /= np.sum(np.abs(W))

    names = ["clinical"]
    if G is not None:
        names.append("clinical_pgx")
    names.extend(("reference_rna_pca2", "reference_direct_rna", "reference_sweet_inspired_graph",
                  "reference_rna_plus_graph"))
    if W is not None:
        names.append("reference_prespecified_gene_score")
    errors = {m: [] for m in names}
    biases = {m: [] for m in names}
    reference = "clinical_pgx" if G is not None else "clinical"
    for held in range(n):
        train_idx = np.array([i for i in range(n) if i != held])
        test_idx = np.array([held])
        base = C if G is None else np.column_stack((C, G))
        c_train, c_test = C[train_idx], C[test_idx]
        xtrain, xtest = X[train_idx], X[test_idx]
        btrain, btest = base[train_idx], base[test_idx]

        predictions = {"clinical": _predict(c_train, c_test, y[train_idx])}
        if G is not None:
            predictions["clinical_pgx"] = _predict(btrain, btest, y[train_idx])
        scale = StandardScaler().fit(xtrain)
        ztrain, ztest = scale.transform(xtrain), scale.transform(xtest)
        pca = PCA(n_components=2, svd_solver="full").fit(ztrain)
        predictions["reference_rna_pca2"] = _predict(
            np.column_stack((btrain, pca.transform(ztrain))),
            np.column_stack((btest, pca.transform(ztest))), y[train_idx])
        predictions["reference_direct_rna"] = _predict(
            np.column_stack((btrain, ztrain)),
            np.column_stack((btest, ztest)), y[train_idx])
        if W is not None:
            predictions["reference_prespecified_gene_score"] = _predict(
                np.column_stack((btrain, ztrain @ W)),
                np.column_stack((btest, ztest @ W)), y[train_idx])
        # Each outer training individual's graph is projected out-of-fold.
        # The heldout individual's graph sees only the frozen outer training
        # transcriptomes, NEVER other heldout subjects or PK labels.
        try:
            graphs = build_outer_split_features(
                X.T, genes=genes, subject_ids=ids,
                train_ids=[ids[i] for i in train_idx], test_ids=[ids[held]],
                tissue=tissue, folds=5)
        except InductiveNetworkError as exc:
            raise RepresentationBenchmarkError("Out-of-fold graph could not be qualified") from exc
        train_map = {r["subject_id"]: r["features"] for r in graphs["training"]}
        test_map = graphs["test"][0]["features"]
        if len(train_map) != len(train_idx) or set(test_map) != set(genes):
            raise RepresentationBenchmarkError("Graph features not aligned across biological subjects")
        gt = np.array([[train_map[ids[i]][g] for g in genes] for i in train_idx])
        gv = np.array([[test_map[g] for g in genes]])
        predictions["reference_sweet_inspired_graph"] = _predict(
            np.column_stack((btrain, gt)), np.column_stack((btest, gv)), y[train_idx])
        predictions["reference_rna_plus_graph"] = _predict(
            np.column_stack((btrain, ztrain, gt)),
            np.column_stack((btest, ztest, gv)), y[train_idx])
        for name in names:
            errors[name].append(abs(predictions[name] - y[held]))
            biases[name].append(predictions[name] - y[held])
    baseline_mae = float(np.mean(errors[reference]))
    results = {}
    for name in names:
        err, bias = np.asarray(errors[name]), np.asarray(biases[name])
        results[name] = {
            "mae": round(float(np.mean(err)), 8),
            "rmse": round(float(np.sqrt(np.mean(bias ** 2))), 8),
            "bias": round(float(np.mean(bias)), 8),
            "delta_mae_vs_reference": round(float(np.mean(err) - baseline_mae), 8),
            "n_biological_test_people": n,
        }
    return {
        "status": "SYNTHETIC_SOFTWARE_CHECK_NOT_PK_EVIDENCE" if source_kind == "synthetic"
                  else "ORIGINAL_SOURCE_DECLARED_EXPLORATORY_LOPO_NOT_EXTERNAL_VALIDATION",
        "method": "fixed_ridge_log_pk_leave_one_person_out",
        "n_independent_subjects": n,
        "pk_target_name": str(target_name),
        "pk_target_unit_as_declared": str(target_unit),
        "n_prespecified_genes": len(genes),
        "source_sha256": source_sha256,
        "genotype_available": G is not None,
        "biological_score_weights_preregistered_by_caller": W is not None,
        "primary_reference": reference,
        "models": results,
        "graph_qualifier": "SWEET-inspired, not official SWEET, not causal or PK-validated",
        "upstream_omic_normalization_out_of_sample_not_independently_verified": True,
        "data_source_identity_and_measurement_units_not_independently_verified": True,
        "no_individual_original_source_values_or_predictions_returned": True,
    }
