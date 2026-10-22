"""Exploratory *SWEET-inspired* inductive graph, NOT the original SWEET method.

Only transcriptomic measurements and a training-frozen reference are accepted.
This is NOT fit for clinical use, not yet numerically calibrated against SWEET,
and not a gene-mechanism or CYP enzyme activity estimate.

Requires optional NumPy dependency; no PK/PD outcomes may enter this module.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Any

import numpy as np

MAX_GENES = 128
MIN_REFERENCE_SAMPLES = 5
EPSILON = 0.01
CONSTRUCTION_V0 = "reference_add_query_v0"
CONSTRUCTION_V1 = "query_anchored_duplicate_v1"
CONSTRUCTIONS = {CONSTRUCTION_V0, CONSTRUCTION_V1}


class InductiveNetworkError(ValueError):
    pass


@dataclass(frozen=True)
class FrozenReference:
    genes: tuple[str, ...]
    subject_ids: tuple[str, ...]
    tissue: str
    expression: np.ndarray
    gene_correlation: np.ndarray
    sample_similarity_min: float
    sample_similarity_max: float
    k: float
    fingerprint: str
    weight_expression: np.ndarray
    weight_genes: tuple[str, ...]


def _matrix(values: Any) -> np.ndarray:
    try:
        a = np.array(values, dtype=np.float64, copy=True)
    except (TypeError, ValueError) as err:
        raise InductiveNetworkError("Expression must be a numeric matrix") from err
    if a.ndim != 2 or not np.isfinite(a).all():
        raise InductiveNetworkError("Expression must be a 2D finite matrix")
    return a


def fit_frozen_reference(
    expression: Any,
    *,
    genes: list[str] | tuple[str, ...],
    subject_ids: list[str] | tuple[str, ...],
    tissue: str,
    k: float = 0.1,
    weight_expression: Any | None = None,
    weight_genes: list[str] | tuple[str, ...] | None = None,
) -> FrozenReference:
    """Fit only on training samples, with a prespecified feature/gene universe.

    expression: genes x training-subjects. No outcome is accepted.
    This independently calibrated candidate is not identical to upstream SWEET.
    """
    x = _matrix(expression)
    genes, subjects = tuple(genes), tuple(subject_ids)
    if not (3 <= len(genes) <= MAX_GENES):
        raise InductiveNetworkError(f"Research POC supports 3–{MAX_GENES} prespecified genes")
    if len(subjects) < MIN_REFERENCE_SAMPLES:
        raise InductiveNetworkError("Need at least five independent reference subjects")
    if (x.shape != (len(genes), len(subjects)) or
            len(set(genes)) != len(genes) or len(set(subjects)) != len(subjects) or
            any(not isinstance(z, str) or not z.strip() for z in genes + subjects) or
            not tissue or not tissue.strip()):
        raise InductiveNetworkError("Invalid genes, subjects, tissue or expression shape")
    if not isinstance(k, (int, float)) or not math.isfinite(k) or not (0 < k <= 1):
        raise InductiveNetworkError("k must be finite and in (0, 1]")
    # Each gene varies across reference subjects, each subject across selected genes.
    if np.any(np.std(x, axis=1) <= 1e-12) or np.any(np.std(x, axis=0) <= 1e-12):
        raise InductiveNetworkError("Constant gene or sample in reference")
    if (weight_expression is None) != (weight_genes is None):
        raise InductiveNetworkError("Full weight expression and gene IDs must be supplied together")
    wx = x.copy() if weight_expression is None else _matrix(weight_expression)
    wgenes = genes if weight_genes is None else tuple(weight_genes)
    if (wx.shape != (len(wgenes), len(subjects))
            or not wgenes or len(set(wgenes)) != len(wgenes)
            or not set(genes).issubset(set(wgenes))
            or any(not isinstance(g, str) or not g.strip() for g in wgenes)):
        raise InductiveNetworkError("Weight genes must include edge genes and align with training subjects")
    if np.any(np.std(wx, axis=0) <= 1e-12):
        raise InductiveNetworkError("Constant sample in weight-expression reference")
    gene_corr = np.corrcoef(x)
    subject_corr = np.corrcoef(wx.T)
    if not np.isfinite(gene_corr).all() or not np.isfinite(subject_corr).all():
        raise InductiveNetworkError("Undefined Pearson correlation")
    n = len(subjects)
    means = (np.sum(subject_corr, axis=1) - 1.0) / (n - 1)
    lo, hi = float(means.min()), float(means.max())
    if hi - lo <= 1e-12:
        raise InductiveNetworkError("Reference similarity range is degenerate")
    digest = hashlib.sha256()
    digest.update("\n".join(genes).encode())
    digest.update("\n".join(subjects).encode())
    digest.update(tissue.encode())
    digest.update(x.tobytes(order="C"))
    digest.update("\n".join(wgenes).encode())
    digest.update(wx.tobytes(order="C"))
    digest.update(str(float(k)).encode())
    x.setflags(write=False)
    wx.setflags(write=False)
    gene_corr.setflags(write=False)
    return FrozenReference(genes, subjects, tissue, x, gene_corr,
                           lo, hi, float(k), digest.hexdigest(), wx, wgenes)


def predict_network(
    reference: FrozenReference,
    query: Any,
    *,
    subject_id: str,
    genes: list[str] | tuple[str, ...],
    tissue: str,
    weight_query: Any | None = None,
    weight_genes: list[str] | tuple[str, ...] | None = None,
    construction: str = CONSTRUCTION_V1,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Compute ONE unseen individual's network from a training-frozen reference.

    v1 mirrors the published SWEET perturbation more closely: the query is
    first added once to the frozen training cohort to form its query-specific
    anchor, then duplicated once for the perturbation. No other held-out
    subject and no outcome enters either matrix.

    v0 is retained only for reproducibility of earlier Omics2PK experiments.
    """
    if not isinstance(subject_id, str) or not subject_id or subject_id in reference.subject_ids:
        raise InductiveNetworkError("Query subject must be new and have a valid identifier")
    if construction not in CONSTRUCTIONS:
        raise InductiveNetworkError("Unknown inductive network construction")
    if tuple(genes) != reference.genes or tissue != reference.tissue:
        raise InductiveNetworkError("Gene order and tissue must match the frozen reference")
    try:
        x = np.asarray(query, dtype=float)
    except (TypeError, ValueError) as err:
        raise InductiveNetworkError("Invalid query expression values") from err
    if x.shape != (len(reference.genes),) or not np.isfinite(x).all():
        raise InductiveNetworkError("Invalid query gene values or shape")
    if np.std(x) <= 1e-12:
        raise InductiveNetworkError("Query profile is constant across genes")
    if reference.weight_genes != reference.genes:
        if weight_query is None or tuple(weight_genes or ()) != reference.weight_genes:
            raise InductiveNetworkError("Unseen subject needs full aligned weight-expression genes")
    elif weight_query is None:
        weight_query = x
    elif weight_genes is not None and tuple(weight_genes) != reference.weight_genes:
        raise InductiveNetworkError("Weight-gene order mismatch")
    try:
        wquery = np.asarray(weight_query, dtype=float)
    except (TypeError, ValueError) as err:
        raise InductiveNetworkError("Invalid weight-query expression") from err
    if (wquery.shape != (len(reference.weight_genes),)
            or not np.isfinite(wquery).all() or np.std(wquery) <= 1e-12):
        raise InductiveNetworkError("Invalid or constant weight-query profile")
    # Weight correlations use the complete, training-frozen gene universe.
    query_weight_matrix = np.column_stack((reference.weight_expression, wquery))
    query_similarities = np.corrcoef(query_weight_matrix.T)[-1, :-1]
    # Edge correlations use only the prespecified mechanism/panel genes.
    query_anchor = np.corrcoef(np.column_stack((reference.expression, x)))
    if construction == CONSTRUCTION_V1:
        # Published SWEET compares a cohort containing the target once against
        # the same cohort with that target duplicated. For an unseen query, the
        # frozen training cohort replaces the published transductive cohort.
        perturbed = np.corrcoef(np.column_stack((reference.expression, x, x)))
        anchor = query_anchor
        effective_cohort_n = len(reference.subject_ids) + 1
    else:
        # Historical Omics2PK v0: compare training-only correlation against
        # training plus the query once.
        perturbed = query_anchor
        anchor = reference.gene_correlation
        effective_cohort_n = len(reference.subject_ids)
    if (not np.isfinite(query_similarities).all()
            or not np.isfinite(anchor).all() or not np.isfinite(perturbed).all()):
        raise InductiveNetworkError("Undefined augmented Pearson correlations")
    similarity = float(query_similarities.mean())
    position = (similarity - reference.sample_similarity_min + EPSILON) / (
        reference.sample_similarity_max - reference.sample_similarity_min + EPSILON)
    clamped = not (0 <= position <= 1)
    position = float(np.clip(position, 0.0, 1.0))
    weight = reference.k * effective_cohort_n * position
    graph = anchor + weight * (perturbed - anchor)
    if not np.isfinite(graph).all():
        raise InductiveNetworkError("Non-finite graph")
    graph = (graph + graph.T) * 0.5
    np.fill_diagonal(graph, 0.0)
    return graph, {
        "reference_fingerprint": reference.fingerprint,
        "reference_n": len(reference.subject_ids),
        "subject_id": subject_id,
        "tissue": tissue,
        "genes": len(reference.genes),
        "weight_genes": len(reference.weight_genes),
        "uses_full_weight_expression": reference.weight_genes != reference.genes,
        "sample_similarity": similarity,
        "weight": weight,
        "weight_clamped_outside_reference": clamped,
        "construction": construction,
        "query_specific_anchor_includes_query_once": construction == CONSTRUCTION_V1,
        "perturbation_duplicates_query_once": construction == CONSTRUCTION_V1,
        "effective_cohort_n_for_weight": effective_cohort_n,
        "method": (
            "experimental_frozen_reference_sweet_inspired_v1"
            if construction == CONSTRUCTION_V1
            else "experimental_frozen_reference_sweet_inspired_v0"
        ),
        "equivalent_to_official_sweet": False,
        "validated_for_pk": False,
    }


def gene_strength_features(graph: Any, *, genes: list[str] | tuple[str, ...]) -> dict[str, float]:
    """Gene-wise unsigned strength; descriptive only, never enzyme activity."""
    a = _matrix(graph)
    n = len(genes)
    if a.shape != (n, n) or n < 3 or not np.allclose(a, a.T, rtol=0, atol=1e-12):
        raise InductiveNetworkError("Expected symmetric graph with matching genes")
    if not np.allclose(np.diag(a), 0, atol=1e-12):
        raise InductiveNetworkError("Expected zero graph diagonal")
    return {g: float(np.sum(np.abs(a[i, :])) / (n - 1))
            for i, g in enumerate(genes)}


def crossfit_training_features(
    expression: Any,
    *,
    genes: list[str] | tuple[str, ...],
    subject_ids: list[str] | tuple[str, ...],
    tissue: str,
    folds: int = 0,
    k: float = 0.1,
    weight_expression: Any | None = None,
    weight_genes: list[str] | tuple[str, ...] | None = None,
    construction: str = CONSTRUCTION_V1,
) -> list[dict[str, Any]]:
    """Build honest, out-of-reference features for TRAINING subjects only.

    The reference for every subject EXCLUDES the subject (and its whole fold).
    folds=0: strict leave-one-out; folds>=2: blocked K-fold cross-fitting.
    No PK/PD target, test subject or metadata outside the supplied training
    cohort is accepted. Preprocessing and gene selection must be frozen by
    the caller before this stage.
    """
    x = _matrix(expression)
    names, subjects = tuple(genes), tuple(subject_ids)
    n = len(subjects)
    if x.shape != (len(names), n):
        raise InductiveNetworkError("Mismatched expression and training identifiers")
    if (weight_expression is None) != (weight_genes is None):
        raise InductiveNetworkError("Supply weight-expression matrix and genes together")
    wx = x if weight_expression is None else _matrix(weight_expression)
    wnames = names if weight_genes is None else tuple(weight_genes)
    if wx.shape != (len(wnames), n):
        raise InductiveNetworkError("Weight-expression sample count mismatched")
    if n < MIN_REFERENCE_SAMPLES + 1:
        raise InductiveNetworkError("At least six training samples required for leave-one-out")
    if len(set(subjects)) != n or any(not isinstance(s, str) or not s for s in subjects):
        raise InductiveNetworkError("Duplicate or invalid training subject IDs")
    if isinstance(folds, bool) or not isinstance(folds, int) or (
        folds != 0 and not 2 <= folds <= n
    ):
        raise InductiveNetworkError("folds must be 0 (leave-one-out) or 2..n")
    if folds != 0 and n - math.ceil(n / folds) < MIN_REFERENCE_SAMPLES:
        raise InductiveNetworkError("Too few training samples after withholding a fold")
    # IDs determine deterministic folds, independent of input column order.
    assignment = {sid: j % folds for j, sid in enumerate(sorted(subjects))} if folds else {}
    records = []
    for pos, sid in enumerate(subjects):
        test_group = assignment.get(sid, pos) if folds else pos
        included = [j for j, other in enumerate(subjects)
                    if (assignment[other] != test_group if folds else j != pos)]
        try:
            ref = fit_frozen_reference(
                x[:, included], genes=names,
                subject_ids=[subjects[j] for j in included], tissue=tissue, k=k,
                weight_expression=wx[:, included], weight_genes=wnames)
            graph, info = predict_network(
                ref, x[:, pos], subject_id=sid, genes=names, tissue=tissue,
                weight_query=wx[:, pos], weight_genes=wnames,
                construction=construction)
        except InductiveNetworkError as exc:
            raise InductiveNetworkError(f"Training subject {sid}: {exc}") from exc
        records.append({
            "subject_id": sid,
            "set": "training",
            "reference_n": info["reference_n"],
            "reference_fingerprint": info["reference_fingerprint"],
            "excluded_fold": test_group,
            "weight": info["weight"],
            "weight_clamped_outside_reference": info["weight_clamped_outside_reference"],
            "features": gene_strength_features(graph, genes=names),
            "construction": info["construction"],
            "effective_cohort_n_for_weight": info["effective_cohort_n_for_weight"],
            "method": info["method"],
            "equivalent_to_official_sweet": False,
        })
    return records


def build_outer_split_features(
    expression: Any,
    *,
    genes: list[str] | tuple[str, ...],
    subject_ids: list[str] | tuple[str, ...],
    train_ids: list[str] | tuple[str, ...],
    test_ids: list[str] | tuple[str, ...],
    tissue: str,
    folds: int = 0,
    k: float = 0.1,
    weight_expression: Any | None = None,
    weight_genes: list[str] | tuple[str, ...] | None = None,
    construction: str = CONSTRUCTION_V1,
) -> dict[str, Any]:
    """Use only the outer TRAINING cohort to fit the held-out reference.

    Training features are cross-fitted. Test features never use any other
    held-out subject. No endpoint or dose/concentration data can enter here.
    A subject is the biological unit; repeated assays must be aggregated or
    handled upstream before invoking this API.
    """
    x = _matrix(expression)
    names, subjects = tuple(genes), tuple(subject_ids)
    train, test = tuple(train_ids), tuple(test_ids)
    if x.shape != (len(names), len(subjects)):
        raise InductiveNetworkError("Inconsistent expression matrix / subject IDs")
    if (weight_expression is None) != (weight_genes is None):
        raise InductiveNetworkError("Supply weight-expression matrix and genes together")
    wx = x if weight_expression is None else _matrix(weight_expression)
    wnames = names if weight_genes is None else tuple(weight_genes)
    if wx.shape != (len(wnames), len(subjects)):
        raise InductiveNetworkError("Weight-expression matrix must match sample count")
    if (len(set(subjects)) != len(subjects)
            or len(set(train)) != len(train) or len(set(test)) != len(test)
            or not train or not test
            or set(train) & set(test)
            or set(train) | set(test) != set(subjects)):
        raise InductiveNetworkError("Outer split must partition unique subjects exactly")
    positions = {sid: i for i, sid in enumerate(subjects)}
    train_idx = [positions[sid] for sid in train]
    test_idx = [positions[sid] for sid in test]
    train_values = x[:, train_idx]
    training = crossfit_training_features(
        train_values, genes=names, subject_ids=train, tissue=tissue, folds=folds, k=k,
        weight_expression=wx[:, train_idx], weight_genes=wnames,
        construction=construction)
    reference = fit_frozen_reference(
        train_values, genes=names, subject_ids=train, tissue=tissue, k=k,
        weight_expression=wx[:, train_idx], weight_genes=wnames)
    held_out = []
    for sid, col in zip(test, test_idx):
        graph, info = predict_network(
            reference, x[:, col], subject_id=sid, genes=names, tissue=tissue,
            weight_query=wx[:, col], weight_genes=wnames,
            construction=construction)
        held_out.append({
            "subject_id": sid,
            "set": "test",
            "reference_n": info["reference_n"],
            "reference_fingerprint": info["reference_fingerprint"],
            "excluded_fold": None,
            "weight": info["weight"],
            "weight_clamped_outside_reference": info["weight_clamped_outside_reference"],
            "features": gene_strength_features(graph, genes=names),
            "construction": info["construction"],
            "effective_cohort_n_for_weight": info["effective_cohort_n_for_weight"],
            "method": info["method"],
            "equivalent_to_official_sweet": False,
        })
    return {
        "training": training,
        "test": held_out,
        "outer_reference_fingerprint": reference.fingerprint,
        "n_train": len(train),
        "n_test": len(test),
        "n_weight_genes": len(wnames),
        "folds": folds,
        "construction": construction,
        "method": (
            "experimental_frozen_reference_sweet_inspired_v1"
            if construction == CONSTRUCTION_V1
            else "experimental_frozen_reference_sweet_inspired_v0"
        ),
        "equivalent_to_official_sweet": False,
        "trained_with_pk_pd_target": False,
    }
