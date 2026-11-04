"""Nested source-matched metformin PK: minimal molecular budget preserving 80/90% inner gain.

Outcome is observed individual 0-24h AUC or renal CLr, not PD.
Original Rotroff 2016, 32 source-verified matched human predose GC-TOF.
The SHAP ranking, k choice and ALL preprocessing are training-only.
k=0 means no incremental verified omic value, not a software failure.
Inner gain targets are not claims of outer test improvement.
"""
from __future__ import annotations

from collections import Counter
from random import Random
import numpy as np
from sklearn.model_selection import KFold
from sklearn.impute import SimpleImputer

from omics2pk.information_budget import (
    DEFAULT_BUDGETS,TARGETS,cumulative_shap_coverage,
    pareto_frontier, select_gain_retention_budget,validate_budgets,
)
from omics2pk.rotroff33_clinical_e3 import (
    E3SourceError,SEED,actual_original_tables,
)
from omics2pk.rotroff33_xgb_shap_whitebox import _fit_xgb,_model_pred_shap
from omics2pk.rotroff33_pathway_e3 import qualify_original_feature_labels


def _fold_models(clinical_train,molecular_train,y_train,clinical_test,
                 molecular_test,budgets):
    """One fitting fold, original 384 fields ranked by TRAIN-only XGB TreeSHAP.

    Returns heldout predictions for every prespecified k and full-model
    training SHAP importance (not outer-test SHAP).
    """
    median=SimpleImputer(strategy="median",keep_empty_features=True)
    a=median.fit_transform(molecular_train)
    b=median.transform(molecular_test)
    nfields=molecular_train.shape[1]
    if a.shape[1]!=nfields or b.shape[1]!=nfields:
        raise E3SourceError("Source molecular field count drift after imputation")
    full_a=np.column_stack((clinical_train,a))
    full_b=np.column_stack((clinical_test,b))
    model=_fit_xgb(full_a,full_b,y_train)
    _,contrib=_model_pred_shap(model,full_a)
    values=np.abs(contrib[:,clinical_train.shape[1]:]).mean(axis=0)
    ranking=np.argsort(-values,kind="stable")
    predictions={}
    for k in budgets:
        if k==nfields:
            use=model
            features=full_b
        else:
            idx=sorted(int(j) for j in ranking[:k])
            aa=np.column_stack((clinical_train,a[:,idx])) if idx else clinical_train
            bb=np.column_stack((clinical_test,b[:,idx])) if idx else clinical_test
            use=_fit_xgb(aa,bb,y_train)
            features=bb
        pred,_=_model_pred_shap(use,features)
        if not np.isfinite(pred).all():
            raise E3SourceError("Non-finite observed clinical PK predictions")
        predictions[k]=pred
    return predictions,ranking,values


def _heldout_metrics(predictions,y):
    if predictions.shape!=y.shape or not np.isfinite(predictions).all():
        raise E3SourceError("Incomplete original patient-heldout model")
    err=predictions-y
    return {
        "MAE_original_PK_table_units":round(float(np.mean(np.abs(err))),6),
        "RMSE_original_PK_table_units":round(float(np.sqrt(np.mean(err*err))),6),
        "P20_percent":round(float(np.mean(np.abs(err)/y<=.2)*100),4),
        "bias_original_PK_table_units":round(float(np.mean(err)),6),
    }


def _paired_delta(pred,baseline,y):
    d=np.abs(pred-y)-np.abs(baseline-y)
    rnd=Random(SEED)
    boots=sorted(float(np.mean([d[rnd.randrange(len(d))] for _ in d]))
                 for _ in range(1999))
    return {"delta_MAE_negative_is_better":round(float(np.mean(d)),6),
            "descriptive_paired_bootstrap95_no_model_refit":[round(boots[49],6),round(boots[1949],6)],
            "n_original_people_lower_error":int(np.sum(d<0))}


def run_nested_budget(clinical,molecular,y,source_headers,*,candidate_budgets=DEFAULT_BUDGETS,
                      inner_n_splits=4,seed=SEED,source_status=None):
    """Tiny-sample proof-of-principle; ALL threshold choices heldout nested.

    Research outcomes are numeric and *not* normalised drug PK unit claims.
    Subject disjointness must be qualified upstream, before calling this.
    """
    clinical=np.asarray(clinical,dtype=float)
    molecular=np.asarray(molecular,dtype=float)
    y=np.asarray(y,dtype=float)
    if clinical.ndim!=2 or molecular.ndim!=2 or y.ndim!=1:
        raise E3SourceError("Clinical / omics / PK source input must be 2D/2D/1D")
    n,p=molecular.shape
    if (n<12 or clinical.shape[0]!=n or len(y)!=n or len(source_headers)!=p or
        len(set(source_headers))!=p or not np.isfinite(clinical).all() or
        np.isinf(molecular).any() or not np.isfinite(y).all() or np.any(y<=0)):
        raise E3SourceError("Missing original same-subject PK or molecular feature values")
    budgets=validate_budgets(candidate_budgets,p)
    if not isinstance(inner_n_splits,int) or not 3<=inner_n_splits<=min(n-2,7):
        raise E3SourceError("Inner CV must be group-disjoint and sufficiently populated")
    results={k:np.full(n,np.nan) for k in budgets}
    by_fraction={q:np.full(n,np.nan) for q in TARGETS}
    retained={q:Counter() for q in TARGETS}
    proportions={q:[] for q in TARGETS}
    rank_frequency=Counter()
    full_shap_acc=np.zeros(p,dtype=float)
    per_k_shap_training={k:[] for k in budgets}
    n_replicates=0
    for outer_index in range(n):
        train_outer=np.array([i for i in range(n) if i!=outer_index],dtype=int)
        inner_pred={k:np.full(len(train_outer),np.nan) for k in budgets}
        inner_splitter=KFold(n_splits=inner_n_splits,shuffle=True,
                             random_state=seed+outer_index)
        # Within each OUTER TRAINING set, every candidate includes independent
        # ranking + model refit on each INNER train fold.
        for inner_train,inner_valid in inner_splitter.split(train_outer):
            it=train_outer[inner_train];iv=train_outer[inner_valid]
            pred,_,_=_fold_models(clinical[it],molecular[it],y[it],
                                  clinical[iv],molecular[iv],budgets)
            for k in budgets:inner_pred[k][inner_valid]=pred[k]
        if any(not np.isfinite(vec).all() for vec in inner_pred.values()):
            raise E3SourceError("Leaked or incomplete inner patient predictions")
        inner_losses={k:float(np.mean(np.abs(pred-y[train_outer])))
                      for k,pred in inner_pred.items()}
        selected={fraction:select_gain_retention_budget(inner_losses,fraction)
                  for fraction in TARGETS}
        # ONE honest outer prediction per source human, with whole nested
        # training pipeline repeated without ever consulting outer test y.
        outer_preds,rank,importances=_fold_models(
            clinical[train_outer],molecular[train_outer],y[train_outer],
            clinical[outer_index:outer_index+1],
            molecular[outer_index:outer_index+1],budgets)
        for k in budgets:
            results[k][outer_index]=float(outer_preds[k][0])
            coverage=cumulative_shap_coverage(importances,[int(x) for x in rank[:k]])
            if coverage["fraction_of_absolute_SHAP"] is not None:
                per_k_shap_training[k].append(coverage["fraction_of_absolute_SHAP"])
        full_shap_acc+=importances
        n_replicates+=1
        for j in rank[:min(20,p)]:
            if importances[j]>0:rank_frequency[int(j)]+=1
        for fraction,choice in selected.items():
            k=choice["k"]
            by_fraction[fraction][outer_index]=float(outer_preds[k][0])
            retained[fraction][k]+=1
            if choice["best_positive_inner_gain"]>0:
                proportions[fraction].append(
                    choice["selected_inner_gain"]/choice["best_positive_inner_gain"])
    if n_replicates!=n:raise E3SourceError("Wrong original donor-disjoint heldout coverage")
    baseline=results[0]
    fixed_curves={}
    pareto_inputs={}
    for k in budgets:
        info=_heldout_metrics(results[k],y)
        fixed_curves[str(k)]={
            "n_source_omics_features":k,
            **info,
            "paired_vs_clinical":_paired_delta(results[k],baseline,y),
            "mean_fraction_of_training_fold_absolute_SHAP_coverage":
                round(float(np.mean(per_k_shap_training[k])),6)
                if per_k_shap_training[k] else None,
            "SHAP_coverage_does_not_measure_PK_gain":True
        }
        pareto_inputs[str(k)]={"k":k,"loss":info["MAE_original_PK_table_units"]}
    auto={}
    for fraction,predictions in by_fraction.items():
        info=_heldout_metrics(predictions,y)
        auto[str(int(100*fraction))]={
            "gain_retention_inner_cv_target_percent":int(100*fraction),
            "observed_outer_patient_heldout_performance":info,
            "paired_outer_vs_clinical":_paired_delta(predictions,baseline,y),
            "outer_fold_selected_k_frequency":{str(k):retained[fraction][k]
                                               for k in budgets if retained[fraction][k]>0},
            "n_outer_folds_using_zero_omic_features":retained[fraction].get(0,0),
            "mean_selected_omic_features":round(
                sum(k*v for k,v in retained[fraction].items())/n,5),
            "mean_inner_gain_fraction_retained_conditionally_on_positive_inner_gain":
                round(float(np.mean(proportions[fraction])),5) if proportions[fraction] else None,
            "fraction_of_gain_is_training_nested_not_outer_test_claim":True,
        }
    top=sorted(range(p),key=lambda i:(-full_shap_acc[i],i))[:min(25,p)]
    return {
       "status":"NESTED_ORIGINAL_PATIENT_HELDOUT_INFORMATION_BUDGET_EXPLORATORY",
       "n_unique_people":n,"n_source_original_omic_fields":p,
       "n_source_true_pk_outcomes":n,
       "number_outer_leave_one_person_out_folds":n,
       "number_inner_folds_per_outer":inner_n_splits,
       "source_annotation_provenance":source_status,
       "candidate_budgets":list(budgets),
       "inner_selection_target":"smallest k with >=80/90 percent of BEST POSITIVE INNER-CV gain vs clinical-only",
       "outer_patient_biomarker_used_in_feature_ranking":False,
       "outer_patient_biomarker_used_in_choice_of_k":False,
       "fixed_k_full_nested_outer_heldout_curve":fixed_curves,
       "adaptive_gain_retention_budgets":auto,
       "posthoc_descriptive_heldout_pareto_do_not_use_for_selection":pareto_frontier(pareto_inputs),
       "top25_training_fitted_SHAP_ranking_frequency":{
           str(source_headers[i]):{"n_outer_training_models_in_top20":int(rank_frequency.get(i,0)),
                                   "average_train_SHAP_abs_log1p_PK":round(float(full_shap_acc[i]/n),7)}
           for i in top if full_shap_acc[i]>0
       },
       "SHAP_fraction_is_not_fraction_of_model_improvement":True,
       "network_module_biological_link_requires_independent_source_annotations":True,
       "individual_patient_concentrations_intensities_or_predictions_exported":False,
    }


def run_original_source(raw):
    people,source,annotations=qualify_original_feature_labels(raw)
    _,omics=actual_original_tables(raw)
    names=[str(omics[0][0].get(j,"")).strip() for j in range(1,385)]
    if len(names)!=384 or len(set(names))!=384 or len(annotations["unknown"])!=229:
        raise E3SourceError("Wrong source original molecule headers")
    endpoint_results={}
    for endpoint in sorted(people[0]["targets"]):
        sampled=[p for p in people if p["targets"][endpoint] is not None]
        if len(sampled)!=32 or len({p["internal_patient_id"] for p in sampled})!=32:
            raise E3SourceError("Expected 32 original unique human PK subjects")
        baseline=np.asarray([p["clinical_plus_original_serum_creatinine"] for p in sampled])
        omic=np.asarray([[np.nan if v is None else v for v in p["predose_omics"]]
                         for p in sampled])
        y=np.asarray([p["targets"][endpoint] for p in sampled])
        endpoint_results[endpoint]=run_nested_budget(
             baseline,omic,y,names,
             source_status="Rotroff2016 plasma predose authentic original published bin headers")
    return {
        "status":"ORIGINAL_ROTR0FF_METFORMIN_NESTED_80_90_GAIN_BUDGET_SOURCE_EXECUTED",
        "original_DOI":"10.3389/fphar.2016.00135",
        "source_original_clinical_SHA256":source["original_clinical_workbook_SHA256"],
        "source_original_omics_SHA256":source["original_untargeted_metabolomic_workbook_SHA256"],
        "source_tissue":"circulating_human_plasma_predose_GC_TOF",
        "n_actual_genotypes_or_PopPK_reference_models":0,
        "pk_physical_units_confirmed":False,
        "source_original_non_patient_descriptive_metadata_only":True,
        "endpoints":endpoint_results,
        "clinical_dosing_value_demonstrated":False,
    }
