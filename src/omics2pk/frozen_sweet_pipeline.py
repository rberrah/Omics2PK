"""Prospective train-frozen SWEET-inspired pipeline for Omics2PK research.

This module bundles every cohort-fitted component needed to transform one new
transcriptomic profile into an inductive SWEET-inspired feature vector. It is
an Omics2PK extension, not the original transductive SWEET algorithm, and it
accepts no PK/PD outcome.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from omics2pk.inductive_sweet import (
    CONSTRUCTION_V1,
    CONSTRUCTIONS,
    FrozenReference,
    fit_frozen_reference,
    gene_strength_features,
    predict_network,
)
from omics2pk.transcriptome_preprocessing import (
    FrozenTranscriptomePreprocessor,
    fit_transcriptome_preprocessor,
    transform_transcriptome,
)


class FrozenSweetPipelineError(ValueError):
    """Invalid fit or prediction request for the frozen inductive pipeline."""


@dataclass(frozen=True)
class FrozenInductiveSweetPipeline:
    """Immutable training artifact used for genuinely unseen RNA profiles."""

    genes: tuple[str, ...]
    edge_genes: tuple[str, ...]
    edge_indices: tuple[int, ...]
    tissue: str
    preprocessing_mode: str
    construction: str
    k: float
    preprocessor: FrozenTranscriptomePreprocessor
    reference: FrozenReference
    fingerprint: str


def _fingerprint_pipeline(
    *,
    genes: tuple[str, ...],
    edge_genes: tuple[str, ...],
    tissue: str,
    preprocessing_mode: str,
    construction: str,
    k: float,
    preprocessor_fingerprint: str,
    reference_fingerprint: str,
) -> str:
    digest = hashlib.sha256()
    digest.update(b"omics2pk.frozen-inductive-sweet.v1\n")
    for value in (
        "\n".join(genes),
        "\n".join(edge_genes),
        tissue,
        preprocessing_mode,
        construction,
        repr(float(k)),
        preprocessor_fingerprint,
        reference_fingerprint,
    ):
        digest.update(value.encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def fit_frozen_inductive_sweet_pipeline(
    expression: Any,
    *,
    genes: list[str] | tuple[str, ...],
    subject_ids: list[str] | tuple[str, ...],
    edge_genes: list[str] | tuple[str, ...],
    tissue: str,
    mode: str,
    k: float = 0.1,
    construction: str = CONSTRUCTION_V1,
) -> FrozenInductiveSweetPipeline:
    """Fit one complete RNA-to-network pipeline on training subjects only.

    expression is full-gene expression for TRAINING subjects only
    (genes x subjects). edge_genes is a prespecified ordered subset.
    The same frozen transformed full-gene matrix is used for sample-weight
    calibration and the subset is used for network edges.

    No endpoint, test cohort or hyperparameter search is accepted here.
    """
    full_genes = tuple(genes)
    edges = tuple(edge_genes)
    subjects = tuple(subject_ids)
    if construction not in CONSTRUCTIONS:
        raise FrozenSweetPipelineError("Unknown inductive SWEET construction")
    if (
        not isinstance(k, (int, float))
        or isinstance(k, bool)
        or not math.isfinite(k)
        or not (0 < k <= 1)
    ):
        raise FrozenSweetPipelineError("k must be finite and in (0, 1]")
    if not isinstance(tissue, str) or not tissue.strip():
        raise FrozenSweetPipelineError("Tissue must be a nonempty prespecified label")
    if len(edges) < 3 or len(set(edges)) != len(edges):
        raise FrozenSweetPipelineError("Need at least three unique prespecified edge genes")
    positions = {gene: i for i, gene in enumerate(full_genes)}
    if len(positions) != len(full_genes) or any(gene not in positions for gene in edges):
        raise FrozenSweetPipelineError(
            "Edge genes must be a unique subset of the full gene universe"
        )
    edge_indices = tuple(positions[gene] for gene in edges)

    try:
        preprocessor = fit_transcriptome_preprocessor(
            expression,
            genes=full_genes,
            subject_ids=subjects,
            mode=mode,
        )
        transformed = transform_transcriptome(
            preprocessor,
            expression,
            genes=full_genes,
        )
        edge_expression = transformed[np.asarray(edge_indices, dtype=int), :]
        reference = fit_frozen_reference(
            edge_expression,
            genes=edges,
            subject_ids=subjects,
            tissue=tissue,
            k=float(k),
            weight_expression=transformed,
            weight_genes=full_genes,
        )
    except (ValueError, TypeError) as exc:
        raise FrozenSweetPipelineError(str(exc)) from exc

    fingerprint = _fingerprint_pipeline(
        genes=full_genes,
        edge_genes=edges,
        tissue=tissue,
        preprocessing_mode=mode,
        construction=construction,
        k=float(k),
        preprocessor_fingerprint=preprocessor.fingerprint,
        reference_fingerprint=reference.fingerprint,
    )
    return FrozenInductiveSweetPipeline(
        genes=full_genes,
        edge_genes=edges,
        edge_indices=edge_indices,
        tissue=tissue,
        preprocessing_mode=mode,
        construction=construction,
        k=float(k),
        preprocessor=preprocessor,
        reference=reference,
        fingerprint=fingerprint,
    )


def predict_frozen_inductive_sweet_pipeline(
    pipeline: FrozenInductiveSweetPipeline,
    query: Any,
    *,
    subject_id: str,
    genes: list[str] | tuple[str, ...],
) -> dict[str, Any]:
    """Transform and predict one new subject without refitting any component."""
    if tuple(genes) != pipeline.genes:
        raise FrozenSweetPipelineError(
            "Gene identity/order must match the frozen pipeline"
        )
    try:
        values = np.asarray(query, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise FrozenSweetPipelineError("Query transcriptome must be numeric") from exc
    if values.ndim == 1:
        values = values[:, None]
    if values.shape != (len(pipeline.genes), 1):
        raise FrozenSweetPipelineError(
            "Prediction accepts exactly one full-gene subject profile"
        )

    try:
        transformed = transform_transcriptome(
            pipeline.preprocessor,
            values,
            genes=pipeline.genes,
        )
        weight_query = transformed[:, 0]
        edge_query = weight_query[np.asarray(pipeline.edge_indices, dtype=int)]
        graph, network_info = predict_network(
            pipeline.reference,
            edge_query,
            subject_id=subject_id,
            genes=pipeline.edge_genes,
            tissue=pipeline.tissue,
            weight_query=weight_query,
            weight_genes=pipeline.genes,
            construction=pipeline.construction,
        )
    except (ValueError, TypeError) as exc:
        raise FrozenSweetPipelineError(str(exc)) from exc

    return {
        "subject_id": subject_id,
        "pipeline_fingerprint": pipeline.fingerprint,
        "preprocessing_fingerprint": pipeline.preprocessor.fingerprint,
        "reference_fingerprint": pipeline.reference.fingerprint,
        "training_subject_count": len(pipeline.reference.subject_ids),
        "preprocessing_mode": pipeline.preprocessing_mode,
        "construction": pipeline.construction,
        "k": pipeline.k,
        "tissue": pipeline.tissue,
        "n_weight_genes": len(pipeline.genes),
        "n_edge_genes": len(pipeline.edge_genes),
        "features": gene_strength_features(graph, genes=pipeline.edge_genes),
        "network_qc": {
            "sample_similarity": network_info["sample_similarity"],
            "sample_similarity_z_vs_reference": network_info[
                "sample_similarity_z_vs_reference"
            ],
            "weight_position_unclipped": network_info["weight_position_unclipped"],
            "weight_position_clipped": network_info["weight_position_clipped"],
            "weight_clamped_outside_reference": network_info[
                "weight_clamped_outside_reference"
            ],
            "weight_clamp_direction": network_info["weight_clamp_direction"],
            "weight": network_info["weight"],
        },
        "preprocessor_refit_for_query": False,
        "reference_refit_for_query": False,
        "other_query_subjects_used": False,
        "pk_pd_target_used": False,
        "equivalent_to_official_sweet": False,
        "validated_for_pk": False,
    }


def frozen_pipeline_manifest(
    pipeline: FrozenInductiveSweetPipeline,
) -> dict[str, Any]:
    """Return non-expression metadata suitable for provenance records."""
    return {
        "schema": "omics2pk.frozen-inductive-sweet.v1",
        "pipeline_fingerprint": pipeline.fingerprint,
        "preprocessing_fingerprint": pipeline.preprocessor.fingerprint,
        "reference_fingerprint": pipeline.reference.fingerprint,
        "training_subject_count": len(pipeline.reference.subject_ids),
        "genes": list(pipeline.genes),
        "edge_genes": list(pipeline.edge_genes),
        "tissue": pipeline.tissue,
        "preprocessing_mode": pipeline.preprocessing_mode,
        "construction": pipeline.construction,
        "k": pipeline.k,
        "training_only_fit": True,
        "pk_pd_target_used": False,
        "equivalent_to_official_sweet": False,
        "validated_for_pk": False,
    }
