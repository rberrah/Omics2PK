"""Train-only transcriptome preprocessing for inductive Omics2PK research.

No PK/PD target is accepted. Cohort-learned statistics are fitted on outer
training subjects only and frozen before any held-out subject is transformed.
This module does not claim to replace assay-specific QC or batch correction.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Any

import numpy as np


MODE_PASSTHROUGH = "source_normalized_passthrough"
MODE_TRAIN_ZSCORE = "source_normalized_train_zscore"
MODE_COUNTS_CPM_LOG1P_ZSCORE = "raw_counts_cpm_log1p_train_zscore"
MODES={MODE_PASSTHROUGH,MODE_TRAIN_ZSCORE,MODE_COUNTS_CPM_LOG1P_ZSCORE}


class TranscriptomePreprocessingError(ValueError):
    pass


@dataclass(frozen=True)
class FrozenTranscriptomePreprocessor:
    genes: tuple[str, ...]
    mode: str
    train_subject_ids: tuple[str, ...]
    medians: np.ndarray | None
    means: np.ndarray | None
    scales: np.ndarray | None
    fingerprint: str


def _input(values: Any, *, allow_nan: bool) -> np.ndarray:
    try:
        x=np.array(values,dtype=np.float64,copy=True)
    except (TypeError,ValueError) as err:
        raise TranscriptomePreprocessingError("Transcriptome values must be numeric") from err
    if x.ndim!=2 or np.isinf(x).any() or (not allow_nan and np.isnan(x).any()):
        raise TranscriptomePreprocessingError("Transcriptome must be a finite genes-by-subjects matrix")
    return x


def _validate_ids(genes,subject_ids,shape):
    genes=tuple(genes); subjects=tuple(subject_ids)
    if (shape!=(len(genes),len(subjects)) or len(genes)<3 or len(subjects)<5
            or len(set(genes))!=len(genes) or len(set(subjects))!=len(subjects)
            or any(not isinstance(x,str) or not x.strip() for x in genes+subjects)):
        raise TranscriptomePreprocessingError("Invalid genes, subjects or matrix alignment")
    return genes,subjects


def _per_sample_cpm_log1p(x: np.ndarray) -> np.ndarray:
    if np.isnan(x).any() or np.any(x<0):
        raise TranscriptomePreprocessingError("Raw counts must be finite and nonnegative")
    totals=x.sum(axis=0,keepdims=True)
    if np.any(totals<=0):
        raise TranscriptomePreprocessingError("Each raw-count sample needs positive library size")
    return np.log1p(x/totals*1_000_000.0)


def fit_transcriptome_preprocessor(
    expression: Any, *, genes, subject_ids, mode: str
) -> FrozenTranscriptomePreprocessor:
    """Fit preprocessing statistics on training subjects only."""
    if mode not in MODES:
        raise TranscriptomePreprocessingError("Unknown transcriptome preprocessing mode")
    allow_nan=mode==MODE_TRAIN_ZSCORE
    x=_input(expression,allow_nan=allow_nan)
    genes,subjects=_validate_ids(genes,subject_ids,x.shape)
    if mode==MODE_COUNTS_CPM_LOG1P_ZSCORE:
        x=_per_sample_cpm_log1p(x)
    if mode==MODE_PASSTHROUGH:
        medians=means=scales=None
    else:
        medians=np.nanmedian(x,axis=1)
        if not np.isfinite(medians).all():
            raise TranscriptomePreprocessingError("Training gene entirely missing")
        missing=np.isnan(x)
        if missing.any():
            x=x.copy()
            x[missing]=np.take(medians,np.nonzero(missing)[0])
        means=x.mean(axis=1)
        scales=x.std(axis=1)
        if np.any(scales<=1e-12) or not np.isfinite(scales).all():
            raise TranscriptomePreprocessingError("Constant training gene cannot be standardized")
        for a in (medians,means,scales):
            a.setflags(write=False)
    digest=hashlib.sha256()
    digest.update(mode.encode())
    digest.update("\n".join(genes).encode())
    digest.update("\n".join(subjects).encode())
    if medians is not None:
        digest.update(medians.tobytes(order="C"))
        digest.update(means.tobytes(order="C"))
        digest.update(scales.tobytes(order="C"))
    return FrozenTranscriptomePreprocessor(
        genes,mode,subjects,medians,means,scales,digest.hexdigest())


def transform_transcriptome(
    frozen: FrozenTranscriptomePreprocessor, expression: Any, *, genes
) -> np.ndarray:
    """Transform any number of subjects using only frozen training statistics."""
    allow_nan=frozen.mode==MODE_TRAIN_ZSCORE
    x=_input(expression,allow_nan=allow_nan)
    if tuple(genes)!=frozen.genes or x.shape[0]!=len(frozen.genes):
        raise TranscriptomePreprocessingError("Gene identity/order must match frozen preprocessing")
    if frozen.mode==MODE_COUNTS_CPM_LOG1P_ZSCORE:
        x=_per_sample_cpm_log1p(x)
    if frozen.mode==MODE_PASSTHROUGH:
        return x
    missing=np.isnan(x)
    if missing.any():
        x=x.copy()
        x[missing]=np.take(frozen.medians,np.nonzero(missing)[0])
    z=(x-frozen.means[:,None])/frozen.scales[:,None]
    if not np.isfinite(z).all():
        raise TranscriptomePreprocessingError("Nonfinite frozen transcriptome transform")
    return z


def preprocess_outer_split(
    expression: Any, *, genes, subject_ids, train_ids, test_ids, mode: str
) -> dict[str, object]:
    """Fit on outer train only and transform train/test without outcome access."""
    allow_nan=mode==MODE_TRAIN_ZSCORE
    x=_input(expression,allow_nan=allow_nan)
    genes,subjects=_validate_ids(genes,subject_ids,x.shape)
    train=tuple(train_ids); test=tuple(test_ids)
    if (not train or not test or len(set(train))!=len(train) or len(set(test))!=len(test)
            or set(train)&set(test) or set(train)|set(test)!=set(subjects)):
        raise TranscriptomePreprocessingError("Outer train/test IDs must partition subjects exactly")
    pos={sid:i for i,sid in enumerate(subjects)}
    ti=[pos[sid] for sid in train]; qi=[pos[sid] for sid in test]
    frozen=fit_transcriptome_preprocessor(
        x[:,ti],genes=genes,subject_ids=train,mode=mode)
    return {
        "training":transform_transcriptome(frozen,x[:,ti],genes=genes),
        "test":transform_transcriptome(frozen,x[:,qi],genes=genes),
        "train_ids":train,
        "test_ids":test,
        "genes":genes,
        "mode":mode,
        "preprocessing_fingerprint":frozen.fingerprint,
        "uses_test_subjects_to_fit_preprocessing":False,
        "batch_correction_executed":False,
        "pk_pd_target_used":False,
        "frozen":frozen,
    }

def build_preprocessed_inductive_outer_split(
    expression: Any, *, genes, edge_genes, subject_ids, train_ids, test_ids,
    tissue: str, mode: str, folds: int = 0, k: float = 0.1,
    construction: str = "query_anchored_duplicate_v1",
) -> dict[str, object]:
    """Preprocess on outer train, then build leakage-controlled inductive features.

    The full gene universe is used for sample-weight similarity after the same
    frozen preprocessing; ``edge_genes`` is a prespecified subset used for the
    network edges. No endpoint is accepted anywhere in this API.
    """
    prep=preprocess_outer_split(
        expression,genes=genes,subject_ids=subject_ids,
        train_ids=train_ids,test_ids=test_ids,mode=mode)
    all_genes=tuple(genes)
    edges=tuple(edge_genes)
    if (len(edges)<3 or len(set(edges))!=len(edges)
            or not set(edges).issubset(set(all_genes))):
        raise TranscriptomePreprocessingError(
            "Edge genes must be a unique prespecified subset of the frozen gene universe")
    combined=np.column_stack((prep["training"],prep["test"]))
    combined_ids=tuple(prep["train_ids"])+tuple(prep["test_ids"])
    pos={g:i for i,g in enumerate(all_genes)}
    edge_values=combined[[pos[g] for g in edges],:]
    from omics2pk.inductive_sweet import build_outer_split_features
    out=build_outer_split_features(
        edge_values,genes=edges,subject_ids=combined_ids,
        train_ids=prep["train_ids"],test_ids=prep["test_ids"],
        tissue=tissue,folds=folds,k=k,
        weight_expression=combined,weight_genes=all_genes,
        construction=construction)
    out.update({
        "preprocessing_mode":mode,
        "preprocessing_fingerprint":prep["preprocessing_fingerprint"],
        "preprocessing_fitted_on_outer_train_only":True,
        "batch_correction_executed":False,
        "preprocessing_pk_pd_target_used":False,
        "n_preprocessing_genes":len(all_genes),
        "n_edge_genes":len(edges),
    })
    return out
