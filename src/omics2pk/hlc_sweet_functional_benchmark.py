"""Site-disjoint HLC transcriptome -> measured CYP activity benchmark.

This module is deliberately source-agnostic: it accepts an already authenticated
same-donor Human Liver Cohort matrix and never downloads or redistributes donor
rows. The primary design is leave-one-study-site-out so site/batch structure
cannot be mistaken for biological generalization.

Research scope: E2 human-liver molecular phenotype -> in-vitro CYP activity.
This is not clinical PK, PopPK, or validation of original transductive SWEET.
"""
from __future__ import annotations

from math import expm1, isfinite, log1p
from statistics import mean
from typing import Any

import numpy as np

from omics2pk.epa_liver_cyp_validation import _ridge
from omics2pk.frozen_sweet_pipeline import (
    fit_frozen_inductive_sweet_pipeline,
    predict_frozen_inductive_sweet_pipeline,
)
from omics2pk.inductive_sweet import CONSTRUCTION_V1
from omics2pk.transcriptome_preprocessing import MODE_TRAIN_ZSCORE


class HLCFunctionalBenchmarkError(ValueError):
    pass


def _validate(
    expression: Any, *, genes, subject_ids, sites, activity,
    target_gene: str, edge_genes
):
    try:
        x=np.asarray(expression,dtype=np.float64)
        y=np.asarray(activity,dtype=np.float64)
    except (TypeError,ValueError) as exc:
        raise HLCFunctionalBenchmarkError("Expression and activity must be numeric") from exc
    genes=tuple(str(g) for g in genes)
    subjects=tuple(str(s) for s in subject_ids)
    sites=tuple(str(s) for s in sites)
    edges=tuple(str(g) for g in edge_genes)
    if x.ndim!=2 or x.shape!=(len(genes),len(subjects)):
        raise HLCFunctionalBenchmarkError("Expression must be genes x subjects")
    if len(y)!=len(subjects):
        raise HLCFunctionalBenchmarkError("Activity length must match subjects")
    if len(set(genes))!=len(genes) or len(set(subjects))!=len(subjects):
        raise HLCFunctionalBenchmarkError("Gene and subject IDs must be unique")
    if any(not s.strip() for s in subjects+sites):
        raise HLCFunctionalBenchmarkError("Subject/site labels cannot be blank")
    if len(sites)!=len(subjects) or len(set(sites))<3:
        raise HLCFunctionalBenchmarkError("Need at least three source study sites")
    if not np.isfinite(y).all() or np.any(y<0):
        raise HLCFunctionalBenchmarkError("Activity must be finite and nonnegative")
    if target_gene not in genes:
        raise HLCFunctionalBenchmarkError("Target CYP transcript absent from expression")
    if len(edges)<3 or len(set(edges))!=len(edges) or not set(edges).issubset(set(genes)):
        raise HLCFunctionalBenchmarkError("Edge genes must be a unique expression subset")
    counts={site:sites.count(site) for site in set(sites)}
    if min(counts.values())<3:
        raise HLCFunctionalBenchmarkError("Each source site needs at least three donors")
    return x,y,genes,subjects,sites,edges


def _pca_fit_transform(train_x: np.ndarray, test_x: np.ndarray, n_components: int):
    """Training-only PCA using rows=features, columns=subjects."""
    mu=np.nanmedian(train_x,axis=1)
    if not np.isfinite(mu).all():
        raise HLCFunctionalBenchmarkError("Training gene entirely missing")
    a=train_x.copy()
    b=test_x.copy()
    for matrix in (a,b):
        miss=np.isnan(matrix)
        if miss.any():
            matrix[miss]=np.take(mu,np.nonzero(miss)[0])
    center=np.mean(a,axis=1)
    a=a-center[:,None]
    b=b-center[:,None]
    # SVD over training subjects only. Components are frozen before test projection.
    _,_,vt=np.linalg.svd(a.T,full_matrices=False)
    k=min(n_components,vt.shape[0])
    basis=vt[:k,:]
    return a.T@basis.T, b.T@basis.T


def _mae(values):
    return mean(values) if values else float("nan")


def evaluate_site_disjoint_hlc(
    expression: Any, *, genes, subject_ids, sites, activity,
    target_gene: str, edge_genes, tissue: str="liver",
    preprocessing_mode: str=MODE_TRAIN_ZSCORE,
    pca_components: int=3, k: float=0.1,
):
    """Evaluate fixed comparators with one entire source site held out at a time.

    Arms:
      - training_logmean: outcome-only reference
      - target_RNA: prespecified CYP transcript
      - PCA3: broad transcriptome, PCA fitted on training sites only
      - SWEET_v1: train-frozen inductive SWEET-inspired gene-strength features

    No endpoint is used to fit preprocessing, PCA, SWEET reference, or edge panel.
    Regression is fixed Ridge(alpha=1) with train-only standardization.
    """
    x,y,genes,subjects,sites,edges=_validate(
        expression,genes=genes,subject_ids=subject_ids,sites=sites,activity=activity,
        target_gene=target_gene,edge_genes=edge_genes)
    if isinstance(pca_components,bool) or not isinstance(pca_components,int) or not 1<=pca_components<=10:
        raise HLCFunctionalBenchmarkError("pca_components must be 1..10")
    target_ix=genes.index(target_gene)
    errors={name:[] for name in ("training_logmean","target_RNA","PCA3","SWEET_v1")}
    site_results={}
    all_clamps=0
    all_queries=0

    for held_site in sorted(set(sites)):
        train_idx=[i for i,s in enumerate(sites) if s!=held_site]
        test_idx=[i for i,s in enumerate(sites) if s==held_site]
        train_ids=tuple(subjects[i] for i in train_idx)
        train_y=[log1p(float(y[i])) for i in train_idx]

        # Target transcript: train-only median imputation and train-only ridge scaling.
        target_train=np.asarray(x[target_ix,train_idx],dtype=float)
        med=float(np.nanmedian(target_train))
        if not isfinite(med):
            raise HLCFunctionalBenchmarkError("Target CYP transcript entirely missing in training sites")
        target_train=np.where(np.isnan(target_train),med,target_train)

        # PCA broad RNA: all learned on training sites.
        pca_train,pca_test=_pca_fit_transform(
            x[:,train_idx],x[:,test_idx],pca_components)

        # SWEET v1: fit only on training sites.
        pipe=fit_frozen_inductive_sweet_pipeline(
            x[:,train_idx],genes=genes,subject_ids=train_ids,
            edge_genes=edges,tissue=tissue,mode=preprocessing_mode,
            k=k,construction=CONSTRUCTION_V1)

        local={name:[] for name in errors}
        for pos,idx in enumerate(test_idx):
            truth=float(y[idx])
            pred_ref=max(0.0,expm1(min(30.0,mean(train_y))))
            target_test=float(x[target_ix,idx])
            if not isfinite(target_test):
                target_test=med
            pred_target=max(0.0,expm1(min(30.0,_ridge(
                [[float(v)] for v in target_train],train_y,[target_test]))))
            pred_pca=max(0.0,expm1(min(30.0,_ridge(
                pca_train.tolist(),train_y,pca_test[pos].tolist()))))
            sweet=predict_frozen_inductive_sweet_pipeline(
                pipe,x[:,idx],subject_id=subjects[idx],genes=genes)
            sweet_train=[]
            # Cross-fitted training SWEET features are required for the outcome model:
            # each training subject is represented without being in its own reference.
            for j,train_subject_idx in enumerate(train_idx):
                inner=[q for q in train_idx if q!=train_subject_idx]
                inner_ids=tuple(subjects[q] for q in inner)
                inner_pipe=fit_frozen_inductive_sweet_pipeline(
                    x[:,inner],genes=genes,subject_ids=inner_ids,
                    edge_genes=edges,tissue=tissue,mode=preprocessing_mode,
                    k=k,construction=CONSTRUCTION_V1)
                feat=predict_frozen_inductive_sweet_pipeline(
                    inner_pipe,x[:,train_subject_idx],
                    subject_id=subjects[train_subject_idx],genes=genes)["features"]
                sweet_train.append([feat[g] for g in edges])
            pred_sweet=max(0.0,expm1(min(30.0,_ridge(
                sweet_train,train_y,[sweet["features"][g] for g in edges]))))
            all_queries+=1
            all_clamps+=bool(sweet["network_qc"]["weight_clamped_outside_reference"])

            preds={
                "training_logmean":pred_ref,
                "target_RNA":pred_target,
                "PCA3":pred_pca,
                "SWEET_v1":pred_sweet,
            }
            for arm,pred in preds.items():
                e=abs(pred-truth)
                errors[arm].append(e)
                local[arm].append(e)

        site_results[held_site]={
            "n_heldout_donors":len(test_idx),
            "mae":{arm:round(_mae(vals),6) for arm,vals in local.items()},
            "train_subject_count":len(train_idx),
            "sweET_pipeline_fingerprint":pipe.fingerprint,
        }

    mae={arm:round(_mae(vals),6) for arm,vals in errors.items()}
    return {
        "status":"SITE_DISJOINT_HLC_FUNCTIONAL_BENCHMARK_E2",
        "n_subjects":len(subjects),
        "n_sites":len(set(sites)),
        "target_gene":target_gene,
        "n_edge_genes":len(edges),
        "models":"fixed Ridge alpha=1; activity log1p; no target-driven feature selection",
        "outer_validation":"leave-one-source-site-out",
        "mae_original_activity_units":mae,
        "delta_MAE_SWEET_minus_target_RNA":round(mae["SWEET_v1"]-mae["target_RNA"],6),
        "delta_MAE_PCA_minus_target_RNA":round(mae["PCA3"]-mae["target_RNA"],6),
        "site_results":site_results,
        "sweET_query_clamp_fraction":round(all_clamps/all_queries,6) if all_queries else None,
        "preprocessing_and_network_fit_on_training_sites_only":True,
        "training_SWEET_features_cross_fitted":True,
        "pk_pd_target_used_for_representation_learning":False,
        "not_original_transductive_SWEET":True,
        "not_clinical_PK_or_PopPK":True,
        "individual_donor_predictions_exported":False,
    }
