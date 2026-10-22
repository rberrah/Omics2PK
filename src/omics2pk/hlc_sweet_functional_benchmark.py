"""Site-disjoint HLC RNA -> measured CYP activity benchmark.

Consumes an already authenticated same-donor Human Liver Cohort matrix.
No source rows are bundled or redistributed. Primary validation is
leave-one-study-site-out to prevent center/batch structure from masquerading
as biological generalization.

Scope: E2 human-liver molecular phenotype -> in-vitro CYP activity.
Not clinical PK/PopPK and not original transductive SWEET.
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


def _validate(expression: Any, *, genes, subject_ids, sites, activity,
              target_gene: str, edge_genes):
    x=np.asarray(expression,dtype=np.float64)
    y=np.asarray(activity,dtype=np.float64)
    genes=tuple(map(str,genes)); subjects=tuple(map(str,subject_ids))
    sites=tuple(map(str,sites)); edges=tuple(map(str,edge_genes))
    if x.ndim!=2 or x.shape!=(len(genes),len(subjects)):
        raise HLCFunctionalBenchmarkError("Expression must be genes x subjects")
    if len(y)!=len(subjects) or len(sites)!=len(subjects):
        raise HLCFunctionalBenchmarkError("Activity/site length must match subjects")
    if len(set(genes))!=len(genes) or len(set(subjects))!=len(subjects):
        raise HLCFunctionalBenchmarkError("Gene and subject IDs must be unique")
    if len(set(sites))<3 or min(sites.count(s) for s in set(sites))<3:
        raise HLCFunctionalBenchmarkError("Need >=3 sites and >=3 donors per site")
    if not np.isfinite(y).all() or np.any(y<0):
        raise HLCFunctionalBenchmarkError("Activity must be finite and nonnegative")
    if target_gene not in genes:
        raise HLCFunctionalBenchmarkError("Target CYP transcript absent")
    if len(edges)<3 or len(set(edges))!=len(edges) or not set(edges).issubset(set(genes)):
        raise HLCFunctionalBenchmarkError("Edge genes must be a unique expression subset")
    return x,y,genes,subjects,sites,edges


def _impute_train_median(train: np.ndarray, test: np.ndarray):
    med=np.nanmedian(train,axis=1)
    if not np.isfinite(med).all():
        raise HLCFunctionalBenchmarkError("Training gene entirely missing")
    a=train.copy(); b=test.copy()
    for z in (a,b):
        miss=np.isnan(z)
        if miss.any():
            z[miss]=np.take(med,np.nonzero(miss)[0])
    return a,b


def _pca(train: np.ndarray, test: np.ndarray, n_components: int):
    a,b=_impute_train_median(train,test)
    center=np.mean(a,axis=1)
    a=a-center[:,None]; b=b-center[:,None]
    _,_,vt=np.linalg.svd(a.T,full_matrices=False)
    basis=vt[:min(n_components,vt.shape[0]),:]
    return a.T@basis.T,b.T@basis.T


def _crossfit_sweet_train(x, genes, subjects, train_idx, edges, tissue, mode, k):
    """Outcome-free leave-one-training-subject-out SWEET features."""
    rows=[]
    for held in train_idx:
        inner=[j for j in train_idx if j!=held]
        pipe=fit_frozen_inductive_sweet_pipeline(
            x[:,inner],genes=genes,subject_ids=tuple(subjects[j] for j in inner),
            edge_genes=edges,tissue=tissue,mode=mode,k=k,
            construction=CONSTRUCTION_V1)
        pred=predict_frozen_inductive_sweet_pipeline(
            pipe,x[:,held],subject_id=subjects[held],genes=genes)
        rows.append([pred["features"][g] for g in edges])
    return rows


def evaluate_site_disjoint_hlc(
    expression: Any, *, genes, subject_ids, sites, activity,
    target_gene: str, edge_genes, tissue: str="liver",
    preprocessing_mode: str=MODE_TRAIN_ZSCORE,
    pca_components: int=3, k: float=0.1,
):
    """Fixed comparison under leave-one-site-out external-style validation.

    Arms are training mean, prespecified target CYP RNA, train-only PCA,
    and inductive train-frozen SWEET v1. Representation learning never sees
    activity values.
    """
    x,y,genes,subjects,sites,edges=_validate(
        expression,genes=genes,subject_ids=subject_ids,sites=sites,
        activity=activity,target_gene=target_gene,edge_genes=edge_genes)
    if isinstance(pca_components,bool) or not isinstance(pca_components,int) or not 1<=pca_components<=10:
        raise HLCFunctionalBenchmarkError("pca_components must be 1..10")

    target_ix=genes.index(target_gene)
    arms=("training_logmean","target_RNA","PCA","SWEET_v1")
    errors={a:[] for a in arms}; site_results={}
    clamp_count=query_count=0

    for held_site in sorted(set(sites)):
        train_idx=[i for i,s in enumerate(sites) if s!=held_site]
        test_idx=[i for i,s in enumerate(sites) if s==held_site]
        train_ids=tuple(subjects[i] for i in train_idx)
        train_y=[log1p(float(y[i])) for i in train_idx]

        # Target-RNA comparator: median learned on outer training sites only.
        target_train=np.asarray(x[target_ix,train_idx],float)
        target_med=float(np.nanmedian(target_train))
        if not isfinite(target_med):
            raise HLCFunctionalBenchmarkError("Target RNA missing in all training donors")
        target_train=np.where(np.isnan(target_train),target_med,target_train)

        pca_train,pca_test=_pca(x[:,train_idx],x[:,test_idx],pca_components)

        outer_pipe=fit_frozen_inductive_sweet_pipeline(
            x[:,train_idx],genes=genes,subject_ids=train_ids,edge_genes=edges,
            tissue=tissue,mode=preprocessing_mode,k=k,construction=CONSTRUCTION_V1)
        sweet_train=_crossfit_sweet_train(
            x,genes,subjects,train_idx,edges,tissue,preprocessing_mode,k)

        local={a:[] for a in arms}
        for p,idx in enumerate(test_idx):
            truth=float(y[idx])
            ref=max(0.0,expm1(min(30.0,mean(train_y))))
            tv=float(x[target_ix,idx])
            if not isfinite(tv): tv=target_med
            tr=max(0.0,expm1(min(30.0,_ridge(
                [[float(v)] for v in target_train],train_y,[tv]))))
            pc=max(0.0,expm1(min(30.0,_ridge(
                pca_train.tolist(),train_y,pca_test[p].tolist()))))
            sp=predict_frozen_inductive_sweet_pipeline(
                outer_pipe,x[:,idx],subject_id=subjects[idx],genes=genes)
            sw=max(0.0,expm1(min(30.0,_ridge(
                sweet_train,train_y,[sp["features"][g] for g in edges]))))
            query_count+=1
            clamp_count+=bool(sp["network_qc"]["weight_clamped_outside_reference"])
            for arm,pred in zip(arms,(ref,tr,pc,sw)):
                e=abs(pred-truth); errors[arm].append(e); local[arm].append(e)

        site_results[held_site]={
            "n_heldout_donors":len(test_idx),
            "n_training_donors":len(train_idx),
            "mae":{a:round(mean(v),6) for a,v in local.items()},
            "sweet_pipeline_fingerprint":outer_pipe.fingerprint,
        }

    mae={a:round(mean(v),6) for a,v in errors.items()}
    return {
        "status":"SITE_DISJOINT_HLC_FUNCTIONAL_BENCHMARK_E2",
        "n_subjects":len(subjects),
        "n_sites":len(set(sites)),
        "target_gene":target_gene,
        "n_edge_genes":len(edges),
        "outer_validation":"leave-one-source-site-out",
        "mae_original_activity_units":mae,
        "delta_MAE_SWEET_minus_target_RNA":round(mae["SWEET_v1"]-mae["target_RNA"],6),
        "delta_MAE_PCA_minus_target_RNA":round(mae["PCA"]-mae["target_RNA"],6),
        "site_results":site_results,
        "sweet_query_clamp_fraction":round(clamp_count/query_count,6),
        "preprocessing_and_network_fit_on_training_sites_only":True,
        "training_SWEET_features_cross_fitted":True,
        "pk_pd_target_used_for_representation_learning":False,
        "not_original_transductive_SWEET":True,
        "not_clinical_PK_or_PopPK":True,
        "individual_donor_predictions_exported":False,
    }
