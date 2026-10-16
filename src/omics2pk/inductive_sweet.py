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
    gene_corr = np.corrcoef(x)
    subject_corr = np.corrcoef(x.T)
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
    digest.update(str(float(k)).encode())
    x.setflags(write=False)
    gene_corr.setflags(write=False)
    return FrozenReference(genes, subjects, tissue, x, gene_corr,
                           lo, hi, float(k), digest.hexdigest())


def predict_network(
    reference: FrozenReference,
    query: Any,
    *,
    subject_id: str,
    genes: list[str] | tuple[str, ...],
    tissue: str,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Compute ONE unseen individual's network without updating the reference."""
    if not isinstance(subject_id, str) or not subject_id or subject_id in reference.subject_ids:
        raise InductiveNetworkError("Query subject must be new and have a valid identifier")
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
    augmented = np.column_stack((reference.expression, x))
    # Frozen reference remains unchanged. Query enters only its own contrast.
    query_similarities = np.corrcoef(augmented.T)[-1, :-1]
    new_corr = np.corrcoef(augmented)
    if not np.isfinite(query_similarities).all() or not np.isfinite(new_corr).all():
        raise InductiveNetworkError("Undefined augmented Pearson correlations")
    similarity = float(query_similarities.mean())
    position = (similarity - reference.sample_similarity_min + EPSILON) / (
        reference.sample_similarity_max - reference.sample_similarity_min + EPSILON)
    clamped = not (0 <= position <= 1)
    position = float(np.clip(position, 0.0, 1.0))
    weight = reference.k * len(reference.subject_ids) * position
    graph = reference.gene_correlation + weight * (
        new_corr - reference.gene_correlation)
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
        "sample_similarity": similarity,
        "weight": weight,
        "weight_clamped_outside_reference": clamped,
        "method": "experimental_frozen_reference_sweet_inspired_v0",
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
