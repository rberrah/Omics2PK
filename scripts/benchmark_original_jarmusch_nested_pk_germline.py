"""Nested, source-verified real human CYP3A PK: CYP3A4/5 PGx vs pre-dose MS1.

Jarmusch et al (2020), DOI 10.1111/cts.12785; official 2019 GNPS task
b3966445abad4658a3cdc63c8198a666. Baseline TRUE observed sex/weight +
original categorical CYP3A4 and CYP3A5 genotype; outcome original measured
day1/day9 oral midazolam CL/F (units not authoritatively qualified).
Predose DAY1 plasma only. Source data and volunteer identifiers never logged.

Strict LOO outer person holdout + 4-fold inner hyperparameter and feature
number selection, repeated independently in every outer model's train.
Target pathway candidates are TENTATIVE GNPS spectral hits (not analytical
standard-validated endogenous metabolites); remain ABSTAIN when absent.
These are only preliminary n13/n12, no tacrolimus evidence.
"""
from __future__ import annotations

from collections import Counter
import csv
from io import BytesIO,TextIOWrapper
import json
from math import isfinite
from statistics import mean
from zipfile import ZipFile

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder,StandardScaler

from scripts.validate_local_jarmusch14_gnps_biom_pk import (
    load_original_source,GNPSOriginalSourceError,
)
from scripts.inspect_public_jarmusch_gnps_annotation_schema import STEROID_SIGNATURES

ALPHAS=(1.,20.,100.)
OMIC_K=(0,1,3)
INNER_FOLDS=4
SEED=20261106


class NestedPKValidationError(ValueError):
    pass


def _original_spectral_features(archive):
    """Source-confirmed original BIOM IDs and GNPS putative MS2 library matches."""
    import h5py
    with ZipFile(archive) as z:
        db=[s for s in z.namelist() if s.startswith("DB_result/") and s.endswith(".tsv")]
        if len(db)!=1:
            raise NestedPKValidationError("Not exactly one source DB spectral hit table")
        with TextIOWrapper(z.open(db[0]),encoding="utf-8-sig") as f:
            hits=list(csv.DictReader(f,delimiter="\t"))
        qza=z.read("qiime2_output/qiime2_table.qza")
    with ZipFile(BytesIO(qza)) as q:
        members=[s for s in q.namelist() if s.endswith("/data/feature-table.biom")]
        if len(members)!=1:
            raise NestedPKValidationError("No unique source BIOM table in original QZA")
        biom=q.read(members[0])
    with h5py.File(BytesIO(biom),"r") as h:
        source_ids=[x.decode() for x in h["observation/ids"][:]]
    if len(source_ids)!=6349 or len(set(source_ids))!=len(source_ids):
        raise NestedPKValidationError("Original MS1 feature identifier universe changed")
    id_to_col={s:i for i,s in enumerate(source_ids)}
    eligible=set()
    steroid=set()
    uncertain=[]
    best_library_by_source_idx={}
    for row in hits:
        scan=str(row.get("#Scan#","")).strip()
        if scan not in id_to_col:
            continue
        try:
            mq=float(row.get("MQScore",""))
            ppm=float(row.get("MZErrorPPM",""))
            shared=float(row.get("SharedPeaks",""))
        except (TypeError,ValueError):
            continue
        if not all(isfinite(v) for v in (mq,ppm,shared)):
            continue
        if mq<0.7 or abs(ppm)>20 or shared<6:
            continue
        feature=id_to_col[scan]
        name=str(row.get("Compound_Name","")).casefold()
        eligible.add(feature)
        # The GNPS MS2 library annotation is *putative*. Never assume that
        # library membership proves exact metabolite identity or CYP3A flux.
        current=best_library_by_source_idx.get(int(feature))
        if current is None or mq>current["library_MQScore"]:
            best_library_by_source_idx[int(feature)]={
                "putative_GNPS_compound_name":str(row.get("Compound_Name",""))[:180],
                "library_MQScore":round(float(mq),4),
            }
        if any(t in name for t in STEROID_SIGNATURES):
            steroid.add(feature)
            uncertain.append(name)
    return {
        "spectrally_annotated_qc_feature_indices":np.array(sorted(eligible),int),
        "putative_steroid_bile_feature_indices":np.array(sorted(steroid),int),
        "source_biom_n_features":len(source_ids),
        "source_public_MS1_scan_ids":source_ids,
        "putative_MS2_annotation_by_source_index":best_library_by_source_idx,
        "n_spectral_library_rows":len(hits),
        "annotation_quality_gate":"MQScore>=0.7; shared peaks>=6; |ppm|<=20, putative MS2, not reference-validated",
        "n_putative_steroid_or_bile_features":len(steroid),
        "no_claim_of_unique_identity_or_CYP3A_specificity":True,
    }


def _rank_source_features_inside_train(x,log_target,train,allowed,k):
    """Actual source public ion index; only original TRAIN outcomes enter ranking."""
    if k not in OMIC_K or k<1:
        raise NestedPKValidationError("No nonzero source feature rank requested")
    allowed_idx=np.asarray(allowed,dtype=int)
    if len(allowed_idx)<k:
        raise NestedPKValidationError("No sufficient independent assigned spectral feature")
    valid=allowed_idx[np.mean(x[train][:,allowed_idx]>0,axis=0)>=0.5]
    if len(valid)<k:
        raise NestedPKValidationError("No enough training-only measurable molecular features")
    ztr=x[train][:,valid]
    yd=log_target[train]-log_target[train].mean()
    centered=ztr-ztr.mean(axis=0)
    norms=np.sqrt(np.square(centered).sum(axis=0))
    val=np.abs(np.dot(yd,centered))/np.maximum(norms,1e-12)
    # Stable source feature-index tie break; no heldout person outcome.
    return valid[np.argsort(-val,kind="stable")[:k]]


def _nested_predict_once(x,clinical,geno,log_target,train,test,*,alpha,k,allowed):
    if alpha not in ALPHAS or not 0<=k<=3 or k not in OMIC_K:
        raise NestedPKValidationError("Exploratory unqualified model parameters")
    if len(test)!=1 and len(test)>len(train):
        raise NestedPKValidationError("Invalid validation folds")
    if set(train)&set(test):
        raise NestedPKValidationError("Train/test person overlap")
    encoder=OneHotEncoder(handle_unknown="ignore")
    gt=encoder.fit_transform(geno[train]).toarray()
    gs=encoder.transform(geno[test]).toarray()
    clinical_scale=StandardScaler()
    ct=clinical_scale.fit_transform(clinical[train])
    cs=clinical_scale.transform(clinical[test])
    train_base=np.column_stack((ct,gt))
    test_base=np.column_stack((cs,gs))
    if k>0:
        select=_rank_source_features_inside_train(x,log_target,train,allowed,k)
        sc=StandardScaler()
        train_base=np.column_stack((train_base,sc.fit_transform(x[train][:,select])))
        test_base=np.column_stack((test_base,sc.transform(x[test][:,select])))
    model=Ridge(alpha=alpha)
    model.fit(train_base,log_target[train])
    return np.exp(np.clip(model.predict(test_base),-10,10))


def _nested_cv(original_x,clinical,genotypes,observed,
               *,candidates,source_label):
    x=np.asarray(original_x,dtype=float)
    c=np.asarray(clinical,dtype=float)
    g=np.asarray(genotypes,dtype=object)
    y=np.asarray(observed,dtype=float)
    if (not 10<=len(y)<=30 or x.ndim!=2 or x.shape[0]!=len(y)
         or c.shape!=(len(y),2) or g.shape!=(len(y),2)
         or np.any(x<0) or not np.isfinite(x).all()
         or not np.isfinite(c).all() or not np.isfinite(y).all()
         or np.any(y<=0)):
        raise NestedPKValidationError("Original source person-level inputs invalid")
    # This normalization is per-individual, never using other people or PK.
    xx=np.log1p(x/np.maximum(x.sum(axis=1,keepdims=True),1.)*1_000_000)
    logy=np.log(y)
    models={
        "weight_sex_PGx_only":np.array([],dtype=int),
        "PGx_plus_putative_CYP3A_steroid_bile":np.asarray(
            candidates["putative_steroid_bile_feature_indices"],dtype=int),
        "PGx_plus_all_GNPS_MS2_annotated":np.asarray(
            candidates["spectrally_annotated_qc_feature_indices"],dtype=int),
        "PGx_plus_untargeted_MS1":np.arange(x.shape[1],dtype=int),
    }
    results={}
    for name,allowed in models.items():
        if name!="weight_sex_PGx_only" and len(allowed)<1:
            results[name]={"status":"ABSTAIN_NO_SOURCE_GROUNDED_FEATURES",
                           "n_candidate_molecular_features":0}
            continue
        predictions=np.empty(len(y))
        selected=[]
        ranked_feature_counts=Counter()
        inner_margins_nonzero_minus_zero=[]
        for held in range(len(y)):
            outer_train=np.array([i for i in range(len(y)) if i!=held],int)
            outer_test=np.array([held],int)
            params=[]
            # No molecular candidate extraction, ranking, imputation or
            # regularization selection outside this outer TRAIN.
            k_set=(0,) if name=="weight_sex_PGx_only" else tuple(
                k for k in OMIC_K if k<=len(allowed))
            for k in k_set:
                for alpha in ALPHAS:
                    errors=[]
                    for inner_train_id,inner_val_id in KFold(
                        n_splits=INNER_FOLDS,shuffle=True,random_state=SEED
                    ).split(outer_train):
                        inner_train=outer_train[inner_train_id]
                        inner_val=outer_train[inner_val_id]
                        # One-by-one prediction supports multi-subject inner folds
                        # and ensures identical complete-case basis.
                        for v in inner_val:
                            try:
                                pred=_nested_predict_once(
                                    xx,c,g,logy,inner_train,np.array([v]),
                                    alpha=alpha,k=k,allowed=allowed)[0]
                                errors.append(abs(np.log(pred)-logy[v]))
                            except NestedPKValidationError:
                                errors.append(float("inf"))
                    params.append((mean(errors),k,alpha))
            # Prefer parsimonious k and larger alpha for ties.
            valid=[z for z in params if np.isfinite(z[0])]
            if not valid:
                raise NestedPKValidationError("All nested fits fail original source validation")
            # Explain each k-selection with the *actual inner training-only*
            # log-MAE gain versus k=0, not with the held-out patient's label.
            # Positive margin means the best nonzero-omic configuration was
            # worse than k=0 on inner validation; negative means it was better.
            zeros=[v[0] for v in valid if v[1]==0]
            nonzeros=[v[0] for v in valid if v[1]>0]
            if not zeros:
                raise NestedPKValidationError("Model must always permit zero molecular features")
            if nonzeros:
                inner_margins_nonzero_minus_zero.append(min(nonzeros)-min(zeros))
            chosen=sorted(valid,key=lambda z:(z[0],z[1],-z[2]))[0]
            _,k,alpha=chosen
            pred=_nested_predict_once(
                xx,c,g,logy,outer_train,outer_test,
                alpha=alpha,k=k,allowed=allowed)[0]
            predictions[held]=pred
            selected.append((k,alpha))
            if k>0:
                # Report which publicly identified ion *features* were
                # selected across independently fitted original outer folds.
                # Do not export anything from the heldout person's sample.
                ranked_feature_counts.update(
                    int(v) for v in _rank_source_features_inside_train(
                        xx,logy,outer_train,allowed,k))
        losses=np.abs(predictions-y)
        results[name]={
            "status":"MEASURED_ORIGINAL_PATIENT_HELDOUT_NESTED",
            "n_original_independent_people":len(y),
            "n_candidate_molecular_features":len(allowed),
            "n_nested_outer_LOPO":len(y),
            "inner_model_choice":"4fold train-only; 0/1/3 spectral features; Ridge alpha1/20/100; log CL/F error",
            "MAE_original_source_CLF_numeric_scale":round(float(losses.mean()),6),
            "mean_absolute_log_error":round(float(np.mean(abs(np.log(predictions)-logy))),6),
            "p20":round(float(np.mean(losses/y<=.2)),6),
            "k_selected_counts":{str(k):int(sum(q[0]==k for q in selected)) for k in OMIC_K},
            "selected_original_MS1_ion_feature_frequency_top10":[
                {
                    "original_public_scan_id":str(
                        candidates.get("source_public_MS1_scan_ids",[])[int(idx)]
                        if int(idx)<len(candidates.get("source_public_MS1_scan_ids",[]))
                        else idx),
                    "n_outer_training_models_that_selected_ion":int(frequency),
                    "n_total_outer_training_models":len(y),
                    "putative_MS2_library_hit":candidates.get(
                        "putative_MS2_annotation_by_source_index",{}
                    ).get(int(idx)),
                    "identity_standard_verified":False,
                    "biological_pathway_fully_qualified":False,
                }
                for idx,frequency in sorted(
                    ranked_feature_counts.items(),key=lambda z:(-z[1],z[0])
                )[:10]
            ],
            "selection_frequency_is_model_training_stability_not_biological_validity":True,
            "alpha_selected_counts":{str(a):int(sum(q[1]==a for q in selected)) for a in ALPHAS},
            "selection_rule":{
                "objective":"minimum inner fourfold mean absolute log(CL/F) prediction error",
                "models_compared":"Ridge alpha1/20/100 x k0/1/3, training-only metabolite prevalence and ranking",
                "zero_features_is_eligible_at_every_outer_fold":True,
                "tie_break":"smaller k; then stronger ridge alpha",
                "selection_is_statistical_not_proof_of_absent_biology":True,
                "PK_MAE_in_original_units_is_outer_evaluation_not_inner_selection_metric":True,
            },
            "inner_cv_nonzero_minus_zero_logMAE_summary":(
                {
                    "n_outer_folds":len(inner_margins_nonzero_minus_zero),
                    "n_folds_nonzero_beats_zero":int(sum(v<0 for v in inner_margins_nonzero_minus_zero)),
                    "n_folds_zero_better_or_equal":int(sum(v>=0 for v in inner_margins_nonzero_minus_zero)),
                    "mean_margin":round(float(mean(inner_margins_nonzero_minus_zero)),6),
                    "median_margin":round(float(np.median(inner_margins_nonzero_minus_zero)),6),
                    "min_margin":round(float(min(inner_margins_nonzero_minus_zero)),6),
                    "max_margin":round(float(max(inner_margins_nonzero_minus_zero)),6),
                    "positive_means_unselected_omics_increases_inner_validation_error":True,
                    "no_subject_level_margins_exported":True,
                }
                if inner_margins_nonzero_minus_zero else None),
            "no_source_individual_outcomes_or_predictions_printed":True,
        }
        # Internal only, used for paired model comparisons, never logged.
        results[name]["_losses_private"]=losses
    ref=results["weight_sex_PGx_only"]
    if ref["status"]!="MEASURED_ORIGINAL_PATIENT_HELDOUT_NESTED":
        raise NestedPKValidationError("Invalid baseline genomic prediction")
    r=ref.pop("_losses_private")
    rng=np.random.default_rng(SEED+len(y))
    for name,entry in results.items():
        if entry.get("status")!="MEASURED_ORIGINAL_PATIENT_HELDOUT_NESTED" or name=="weight_sex_PGx_only":
            continue
        d=entry.pop("_losses_private")-r
        bs=np.mean(d[rng.integers(0,len(y),size=(1999,len(y)))],axis=1)
        entry["paired_delta_MAE_vs_clinical_PGx_negative_better"]=round(float(d.mean()),6)
        entry["descriptive_95pct_paired_patient_residual_bootstrap"]=[
            round(float(v),6) for v in np.quantile(bs,[.025,.975])]
        entry["bootstrap_is_descriptive_no_nested_refits_or_multiplicity_control"]=True
    return {"endpoint":source_label,"n_source_people":len(y),"models":results,
            "original_subject_covariates_only_no_future_dose_leakage":True}


def actual_public_source_benchmark(archive):
    x,c,g,outcomes,qc=load_original_source(archive)
    candidates=_original_spectral_features(archive)
    source={
        "source_abstract":"original GNPS Jarmusch et al 2020 research 13 predose plasma biological people",
        "public_archive_is_official_repack_source_verified_by_internal_SHA":True,
        "original_source_QC":qc,
        "n_qc_MS2_library_annotated_features":len(
            candidates["spectrally_annotated_qc_feature_indices"]),
        "n_putative_CYP3A_steroid_bile_MS2_features":candidates["n_putative_steroid_or_bile_features"],
        "n_putative_steroid_bile_features_detected_in_predose_plasma":int(
            np.any(x[:,candidates["putative_steroid_bile_feature_indices"]]>0,axis=0).sum()),
        "n_putative_steroid_bile_features_present_in_at_least_half_predose_plasma":int(
            (np.mean(x[:,candidates["putative_steroid_bile_feature_indices"]]>0,axis=0)>=.5).sum()),
        "n_library_annotated_features_present_in_at_least_half_predose_plasma":int(
            (np.mean(x[:,candidates["spectrally_annotated_qc_feature_indices"]]>0,axis=0)>=.5).sum()),
        "library_annotation_quality_assumption":candidates["annotation_quality_gate"],
        "qualifier":"CYP3A4/5 source PGx limited categorical calls; not comprehensive CYP star allele sequencing",
        "results":{},
    }
    for day in ("day1","day9"):
        mask=np.isfinite(outcomes[day])&(outcomes[day]>0)
        source["results"][day]=_nested_cv(x[mask],c[mask],g[mask],outcomes[day][mask],
                                          candidates=candidates,
                                          source_label="observed oral midazolam CL/F "+day)
    source.update({
        "status":"ACTUAL_ORIGINAL_JARMUSCH_HUMAN_PK_NESTED_PGX_OMICS_EXPLORATORY",
        "not_tacrolimus_PK_or_complete_genotype":True,
        "n_treated_humans_not_independent_cross_period_replicates":13,
        "omic_source_is_pretreatment_plasma_day1_not_liver":True,
        "source_GNPS_MS2_hits_are_putative_not_MS1_chemically_confirmed":True,
        "no_volunteer_rows_or_PK_predictions_redistributed":True,
    })
    return source


if __name__=="__main__":
    import sys
    if len(sys.argv)!=2:
        raise SystemExit("Usage: python -m scripts.benchmark_original_jarmusch_nested_pk_germline /private/source.zip")
    try:
        print(json.dumps(actual_public_source_benchmark(sys.argv[1]),sort_keys=True,indent=2))
    except (NestedPKValidationError,GNPSOriginalSourceError,OSError,ValueError) as exc:
        print(json.dumps({"status":"ABSTAIN_SOURCE_OR_NESTED_BIOLOGICAL_COHORT",
                          "reason_class":type(exc).__name__}))
        raise SystemExit(2)
