"""Prospective train-frozen SWEET-inspired pipeline for Omics2PK research.

This module bundles every cohort-fitted component needed to transform one new
transcriptomic profile into an inductive SWEET-inspired feature vector. It is
an Omics2PK extension, not the original transductive SWEET algorithm, and it
accepts no PK/PD outcome.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
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
    EDGE_SELECTION_MODES,
    EDGE_SELECTION_PRESPECIFIED,
    MODE_PASSTHROUGH,
    FrozenTranscriptomePreprocessor,
    fit_edge_gene_selector,
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
    edge_selection_mode: str
    edge_selection_top_k: int
    edge_selection_fingerprint: str
    preprocessor: FrozenTranscriptomePreprocessor
    reference: FrozenReference
    fingerprint: str
    fingerprint_version: int


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
    edge_selection_mode: str | None = None,
    edge_selection_top_k: int | None = None,
) -> str:
    """Fingerprint prediction state plus the frozen edge-selection protocol.

    Supplying no selection metadata reproduces the historical v1 fingerprint,
    which keeps old persisted artifacts loadable. New fits always use v2.
    """
    digest = hashlib.sha256()
    if edge_selection_mode is None:
        digest.update(b"omics2pk.frozen-inductive-sweet.v1\n")
        values = (
            "\n".join(genes),
            "\n".join(edge_genes),
            tissue,
            preprocessing_mode,
            construction,
            repr(float(k)),
            preprocessor_fingerprint,
            reference_fingerprint,
        )
    else:
        digest.update(b"omics2pk.frozen-inductive-sweet.v2\n")
        values = (
            "\n".join(genes),
            "\n".join(edge_genes),
            tissue,
            preprocessing_mode,
            construction,
            repr(float(k)),
            edge_selection_mode,
            str(edge_selection_top_k),
            preprocessor_fingerprint,
            reference_fingerprint,
        )
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def _prespecified_selector_fingerprint(
    *,
    preprocessing_mode: str,
    subject_ids: tuple[str, ...],
    edge_genes: tuple[str, ...],
) -> str:
    """Reproduce the selector fingerprint for a prespecified panel."""
    digest = hashlib.sha256()
    digest.update(EDGE_SELECTION_PRESPECIFIED.encode())
    digest.update(preprocessing_mode.encode())
    digest.update("\n".join(subject_ids).encode())
    digest.update("\n".join(edge_genes).encode())
    digest.update(str(len(edge_genes)).encode())
    return digest.hexdigest()


def fit_frozen_inductive_sweet_pipeline(
    expression: Any,
    *,
    genes: list[str] | tuple[str, ...],
    subject_ids: list[str] | tuple[str, ...],
    edge_genes: list[str] | tuple[str, ...] | None = None,
    tissue: str,
    mode: str,
    k: float = 0.1,
    construction: str = CONSTRUCTION_V1,
    edge_selection_mode: str = EDGE_SELECTION_PRESPECIFIED,
    edge_top_k: int | None = None,
) -> FrozenInductiveSweetPipeline:
    """Fit one complete RNA-to-network pipeline on training subjects only.

    expression contains TRAINING subjects only. Edge genes are either a
    prospectively prespecified panel or an outcome-free variance panel learned
    exclusively from those training molecular measurements. The selected genes,
    selector fingerprint, preprocessor and SWEET reference are then frozen
    together before any new subject is transformed.

    No endpoint, test cohort or hyperparameter search is accepted here.
    """
    full_genes = tuple(genes)
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
    if edge_selection_mode not in EDGE_SELECTION_MODES:
        raise FrozenSweetPipelineError("Unknown edge-gene selection mode")
    positions = {gene: i for i, gene in enumerate(full_genes)}
    if len(positions) != len(full_genes):
        raise FrozenSweetPipelineError("Full gene universe must be unique")

    try:
        selector = fit_edge_gene_selector(
            expression,
            genes=full_genes,
            subject_ids=subjects,
            mode=mode,
            selection_mode=edge_selection_mode,
            prespecified_genes=edge_genes,
            top_k=edge_top_k,
        )
        edges = selector.selected_genes
        edge_indices = tuple(positions[gene] for gene in edges)
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
        edge_selection_mode=selector.mode,
        edge_selection_top_k=selector.top_k,
    )
    return FrozenInductiveSweetPipeline(
        genes=full_genes,
        edge_genes=edges,
        edge_indices=edge_indices,
        tissue=tissue,
        preprocessing_mode=mode,
        construction=construction,
        k=float(k),
        edge_selection_mode=selector.mode,
        edge_selection_top_k=selector.top_k,
        edge_selection_fingerprint=selector.fingerprint,
        preprocessor=preprocessor,
        reference=reference,
        fingerprint=fingerprint,
        fingerprint_version=2,
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
        "pipeline_fingerprint_version": pipeline.fingerprint_version,
        "preprocessing_fingerprint": pipeline.preprocessor.fingerprint,
        "reference_fingerprint": pipeline.reference.fingerprint,
        "training_subject_count": len(pipeline.reference.subject_ids),
        "preprocessing_mode": pipeline.preprocessing_mode,
        "construction": pipeline.construction,
        "k": pipeline.k,
        "tissue": pipeline.tissue,
        "n_weight_genes": len(pipeline.genes),
        "n_edge_genes": len(pipeline.edge_genes),
        "edge_selection_mode": pipeline.edge_selection_mode,
        "edge_selection_top_k": pipeline.edge_selection_top_k,
        "edge_selection_fingerprint": pipeline.edge_selection_fingerprint,
        "edge_selector_refit_for_query": False,
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
    base = {
        "schema": (
            "omics2pk.frozen-inductive-sweet.v1"
            if pipeline.fingerprint_version == 1
            else "omics2pk.frozen-inductive-sweet.v2"
        ),
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
    if pipeline.fingerprint_version == 2:
        base.update({
            "edge_selection_mode": pipeline.edge_selection_mode,
            "edge_selection_top_k": pipeline.edge_selection_top_k,
            "edge_selection_fingerprint": pipeline.edge_selection_fingerprint,
            "edge_selection_fitted_on_training_only": True,
        })
    elif pipeline.fingerprint_version != 1:
        raise FrozenSweetPipelineError("Unknown frozen pipeline fingerprint version")
    return base


ARTIFACT_SCHEMA_V1 = "omics2pk.frozen-inductive-sweet.artifact.v1"
ARTIFACT_SCHEMA = "omics2pk.frozen-inductive-sweet.artifact.v2"
ARTIFACT_ARRAYS = "arrays.npz"
ARTIFACT_MANIFEST = "manifest.json"


def _fingerprint_preprocessor_state(
    *,
    mode: str,
    genes: tuple[str, ...],
    subject_ids: tuple[str, ...],
    medians: np.ndarray | None,
    means: np.ndarray | None,
    scales: np.ndarray | None,
) -> str:
    """Recompute the exact training-preprocessor fingerprint from frozen state."""
    digest = hashlib.sha256()
    digest.update(mode.encode())
    digest.update("\n".join(genes).encode())
    digest.update("\n".join(subject_ids).encode())
    if medians is not None:
        if means is None or scales is None:
            raise FrozenSweetPipelineError("Incomplete frozen preprocessing state")
        digest.update(medians.tobytes(order="C"))
        digest.update(means.tobytes(order="C"))
        digest.update(scales.tobytes(order="C"))
    return digest.hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_frozen_inductive_sweet_pipeline(
    pipeline: FrozenInductiveSweetPipeline,
    directory: str | Path,
) -> dict[str, Any]:
    """Persist a fitted artifact without pickle or implicit refitting.

    The artifact contains transformed training-reference expression and training
    identifiers. It therefore inherits the source dataset's access restrictions
    and must not be published merely because Omics2PK source code is public.
    """
    root = Path(directory).expanduser().resolve()
    if root.exists():
        if not root.is_dir():
            raise FrozenSweetPipelineError("Frozen artifact destination is not a directory")
        if any(root.iterdir()):
            raise FrozenSweetPipelineError("Frozen artifact destination must be empty")
    else:
        root.mkdir(parents=True, exist_ok=False)

    if (
        pipeline.preprocessor.train_subject_ids != pipeline.reference.subject_ids
        or pipeline.preprocessor.genes != pipeline.genes
        or pipeline.reference.weight_genes != pipeline.genes
        or pipeline.reference.genes != pipeline.edge_genes
        or pipeline.fingerprint_version not in (1, 2)
        or pipeline.edge_selection_mode not in EDGE_SELECTION_MODES
        or pipeline.edge_selection_top_k != len(pipeline.edge_genes)
        or len(pipeline.edge_selection_fingerprint) != 64
        or any(ch not in "0123456789abcdef" for ch in pipeline.edge_selection_fingerprint)
    ):
        raise FrozenSweetPipelineError("Inconsistent in-memory frozen pipeline state")

    medians = pipeline.preprocessor.medians
    means = pipeline.preprocessor.means
    scales = pipeline.preprocessor.scales
    has_scaling = medians is not None
    if has_scaling != (means is not None and scales is not None):
        raise FrozenSweetPipelineError("Incomplete frozen preprocessing state")

    arrays_path = root / ARTIFACT_ARRAYS
    np.savez_compressed(
        arrays_path,
        preprocessor_medians=(
            np.asarray(medians, dtype=np.float64)
            if medians is not None else np.empty(0, dtype=np.float64)
        ),
        preprocessor_means=(
            np.asarray(means, dtype=np.float64)
            if means is not None else np.empty(0, dtype=np.float64)
        ),
        preprocessor_scales=(
            np.asarray(scales, dtype=np.float64)
            if scales is not None else np.empty(0, dtype=np.float64)
        ),
        reference_expression=np.asarray(
            pipeline.reference.expression, dtype=np.float64
        ),
        reference_weight_expression=np.asarray(
            pipeline.reference.weight_expression, dtype=np.float64
        ),
    )
    arrays_sha256 = _sha256_file(arrays_path)
    pipeline_manifest = frozen_pipeline_manifest(pipeline)
    manifest = {
        "schema": (
            ARTIFACT_SCHEMA_V1
            if pipeline.fingerprint_version == 1
            else ARTIFACT_SCHEMA
        ),
        "pipeline": pipeline_manifest,
        "training_subject_ids": list(pipeline.reference.subject_ids),
        "has_scaling_state": has_scaling,
        "arrays": {
            "file": ARTIFACT_ARRAYS,
            "sha256": arrays_sha256,
            "allow_pickle": False,
        },
        "contains_training_derived_molecular_data": True,
        "contains_training_identifiers": True,
        "redistribution_requires_source_permission": True,
    }
    manifest_path = root / ARTIFACT_MANIFEST
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "directory": str(root),
        "manifest": str(manifest_path),
        "arrays": str(arrays_path),
        "arrays_sha256": arrays_sha256,
        "pipeline_fingerprint": pipeline.fingerprint,
        "contains_training_derived_molecular_data": True,
    }


def load_frozen_inductive_sweet_pipeline(
    directory: str | Path,
) -> FrozenInductiveSweetPipeline:
    """Load and integrity-check a persisted artifact without pickle."""
    root = Path(directory).expanduser().resolve()
    manifest_path = root / ARTIFACT_MANIFEST
    arrays_path = root / ARTIFACT_ARRAYS
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, TypeError) as exc:
        raise FrozenSweetPipelineError("Invalid frozen artifact manifest") from exc
    if manifest.get("schema") != ARTIFACT_SCHEMA:
        raise FrozenSweetPipelineError("Unknown frozen artifact schema")
    arrays_meta = manifest.get("arrays")
    if (
        not isinstance(arrays_meta, dict)
        or arrays_meta.get("file") != ARTIFACT_ARRAYS
        or arrays_meta.get("allow_pickle") is not False
        or not isinstance(arrays_meta.get("sha256"), str)
    ):
        raise FrozenSweetPipelineError("Invalid frozen artifact array metadata")
    if not arrays_path.is_file():
        raise FrozenSweetPipelineError("Frozen artifact arrays are missing")
    if _sha256_file(arrays_path) != arrays_meta["sha256"]:
        raise FrozenSweetPipelineError("Frozen artifact array checksum mismatch")

    pipeline_meta = manifest.get("pipeline")
    subject_values = manifest.get("training_subject_ids")
    if not isinstance(pipeline_meta, dict) or not isinstance(subject_values, list):
        raise FrozenSweetPipelineError("Invalid frozen artifact metadata")
    try:
        genes = tuple(str(x) for x in pipeline_meta["genes"])
        edge_genes = tuple(str(x) for x in pipeline_meta["edge_genes"])
        subjects = tuple(str(x) for x in subject_values)
        tissue = str(pipeline_meta["tissue"])
        mode = str(pipeline_meta["preprocessing_mode"])
        construction = str(pipeline_meta["construction"])
        k = float(pipeline_meta["k"])
        expected_preprocessor = str(pipeline_meta["preprocessing_fingerprint"])
        expected_reference = str(pipeline_meta["reference_fingerprint"])
        expected_pipeline = str(pipeline_meta["pipeline_fingerprint"])
    except (KeyError, TypeError, ValueError) as exc:
        raise FrozenSweetPipelineError("Incomplete frozen artifact metadata") from exc
    if (
        not genes
        or len(set(genes)) != len(genes)
        or len(edge_genes) < 3
        or len(set(edge_genes)) != len(edge_genes)
        or not set(edge_genes).issubset(set(genes))
        or not subjects
        or len(set(subjects)) != len(subjects)
        or any(not x.strip() for x in genes + edge_genes + subjects)
    ):
        raise FrozenSweetPipelineError("Invalid frozen artifact identities")

    try:
        with np.load(arrays_path, allow_pickle=False) as arrays:
            expected_keys = {
                "preprocessor_medians",
                "preprocessor_means",
                "preprocessor_scales",
                "reference_expression",
                "reference_weight_expression",
            }
            if set(arrays.files) != expected_keys:
                raise FrozenSweetPipelineError("Unexpected frozen artifact arrays")
            medians = np.array(arrays["preprocessor_medians"], dtype=np.float64, copy=True)
            means = np.array(arrays["preprocessor_means"], dtype=np.float64, copy=True)
            scales = np.array(arrays["preprocessor_scales"], dtype=np.float64, copy=True)
            reference_expression = np.array(
                arrays["reference_expression"], dtype=np.float64, copy=True
            )
            weight_expression = np.array(
                arrays["reference_weight_expression"], dtype=np.float64, copy=True
            )
    except (OSError, ValueError, KeyError) as exc:
        if isinstance(exc, FrozenSweetPipelineError):
            raise
        raise FrozenSweetPipelineError("Invalid frozen artifact arrays") from exc

    has_scaling = manifest.get("has_scaling_state")
    if not isinstance(has_scaling, bool):
        raise FrozenSweetPipelineError("Invalid preprocessing-state flag")
    if has_scaling:
        for name, value in (
            ("medians", medians), ("means", means), ("scales", scales)
        ):
            if value.shape != (len(genes),) or not np.isfinite(value).all():
                raise FrozenSweetPipelineError(
                    f"Invalid frozen preprocessing {name}"
                )
        if np.any(scales <= 1e-12):
            raise FrozenSweetPipelineError("Invalid frozen preprocessing scales")
        for value in (medians, means, scales):
            value.setflags(write=False)
        pre_medians, pre_means, pre_scales = medians, means, scales
    else:
        if mode != MODE_PASSTHROUGH or medians.size or means.size or scales.size:
            raise FrozenSweetPipelineError("Unexpected preprocessing state")
        pre_medians = pre_means = pre_scales = None

    preprocessor_fingerprint = _fingerprint_preprocessor_state(
        mode=mode,
        genes=genes,
        subject_ids=subjects,
        medians=pre_medians,
        means=pre_means,
        scales=pre_scales,
    )
    if preprocessor_fingerprint != expected_preprocessor:
        raise FrozenSweetPipelineError("Frozen preprocessor fingerprint mismatch")
    preprocessor = FrozenTranscriptomePreprocessor(
        genes=genes,
        mode=mode,
        train_subject_ids=subjects,
        medians=pre_medians,
        means=pre_means,
        scales=pre_scales,
        fingerprint=preprocessor_fingerprint,
    )

    positions = {gene: i for i, gene in enumerate(genes)}
    edge_indices = tuple(positions[gene] for gene in edge_genes)
    try:
        reference = fit_frozen_reference(
            reference_expression,
            genes=edge_genes,
            subject_ids=subjects,
            tissue=tissue,
            k=k,
            weight_expression=weight_expression,
            weight_genes=genes,
        )
    except (ValueError, TypeError) as exc:
        raise FrozenSweetPipelineError("Invalid frozen reference state") from exc
    if reference.fingerprint != expected_reference:
        raise FrozenSweetPipelineError("Frozen reference fingerprint mismatch")

    pipeline_fingerprint = _fingerprint_pipeline(
        genes=genes,
        edge_genes=edge_genes,
        tissue=tissue,
        preprocessing_mode=mode,
        construction=construction,
        k=k,
        preprocessor_fingerprint=preprocessor.fingerprint,
        reference_fingerprint=reference.fingerprint,
    )
    if pipeline_fingerprint != expected_pipeline:
        raise FrozenSweetPipelineError("Frozen pipeline fingerprint mismatch")
    if construction not in CONSTRUCTIONS:
        raise FrozenSweetPipelineError("Unknown frozen pipeline construction")

    return FrozenInductiveSweetPipeline(
        genes=genes,
        edge_genes=edge_genes,
        edge_indices=edge_indices,
        tissue=tissue,
        preprocessing_mode=mode,
        construction=construction,
        k=k,
        preprocessor=preprocessor,
        reference=reference,
        fingerprint=pipeline_fingerprint,
    )
