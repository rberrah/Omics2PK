"""Exploratory *original human* metformin PK: XGBoost Tree SHAP + white-box ablation.

Uses ONLY same published subject-matched 32 original Rotroff people and
pre-dose A plasma GC-TOF signals; clinical source age/sex/weight/serum creatinine.
All fit/imputation/preselection are TRAIN-ONLY. Tree SHAP is the exact
XGBoost pred_contribs=False? No: pred_contribs=True, approx_contribs=False.
SHAP explains the fitted log1p(PK) prediction, not causality.
Not PGx, PopPK, external validation, transporter activity, or PTM.
"""
from __future__ import annotations
from collections import Counter
from random import Random

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor, DMatrix

from omics2pk.rotroff33_clinical_e3 import actual_original_tables, E3SourceError, SEED
from omics2pk.rotroff33_pathway_e3 import qualify_original_feature_labels

# Parameters set BEFORE seeing any tree results; no outer-fold tuning.
PARAMS=dict(n_estimators=100, max_depth=2, learning_rate=.045,
            min_child_weight=3, subsample=1.0,colsample_bytree=.8,
            reg_lambda=20.0,reg_alpha=1.0,objective="reg:squarederror",
            tree_method="hist",random_state=SEED,n_jobs=1,verbosity=0)
ARMS=("clinical","clinical_all_384","clinical_named_155",
      "clinical_unknown_229","clinical_urea_purines_5",
      "train_SHAP_top5_Ridge_whitebox","fixed_biochemical5_Ridge_whitebox")

def _fit_xgb(train,test,target):
    model=XGBRegressor(**PARAMS)
    model.fit(train,np.log1p(target))
    return model

def _model_pred_shap(model, x):
    ylog=np.asarray(model.predict(x),dtype=float).reshape(-1)
    contributions=np.asarray(model.get_booster().predict(
        DMatrix(x),pred_contribs=True,approx_contribs=False),dtype=float)
    if contributions.shape!=(len(x),x.shape[1]+1):
        raise E3SourceError("Exact Tree SHAP contribution dimension mismatch")
    if not np.allclose(contributions.sum(axis=1),ylog,rtol=1e-5,atol=1e-4):
        raise E3SourceError("Exact SHAP contributions not additive to model raw log PK output")
    return np.maximum(0.0,np.expm1(np.clip(ylog,-100.0,25.0))),contributions[:,:-1]

def _loocv(records,annotation,headers,target,clinical_key):
    source=[r for r in records if r["targets"][target] is not None]
    if len(source)!=32:raise E3SourceError("Need 32 observed original source PK donors")
    ids=[r["internal_patient_id"] for r in source]
    if len(set(ids))!=32:raise E3SourceError("Identity leak")
    y=np.asarray([r["targets"][target] for r in source],dtype=float)
    clinical=np.asarray([r[clinical_key] for r in source],dtype=float)
    X=np.asarray([[np.nan if v is None else v for v in r["predose_omics"]] for r in source],dtype=float)
    if X.shape!=(32,384) or not np.isfinite(clinical).all() or np.any(y<=0):
        raise E3SourceError("Original unqualified source shape/positive PK")
    groups={
      "clinical":[],
      "clinical_all_384":list(range(384)),
      "clinical_named_155":annotation["named"],
      "clinical_unknown_229":annotation["unknown"],
      "clinical_urea_purines_5":sorted(set(annotation["paths"]["urea_cycle_partial"]
                                               +annotation["paths"]["purine_degradation_partial"])),
    }
    result={arm:np.full(32,np.nan,dtype=float) for arm in ARMS}
    heldout_shap=np.zeros((32,384),dtype=float)
    fold_rank=Counter()
    fold_selected=Counter()
    for i in range(32):
        tr=np.arange(32)!=i
        c_train=clinical[tr];c_test=clinical[i:i+1]
        # Imputation in one given outer training fold only, no precomputed
        # feature selection or model-dependent patient sample QC.
        median=SimpleImputer(strategy="median",keep_empty_features=True)
        train=median.fit_transform(X[tr]);test=median.transform(X[i:i+1])
        models={}
        for arm in groups:
            chosen=groups[arm]
            if chosen:
                a=np.column_stack([c_train,train[:,chosen]])
                b=np.column_stack([c_test,test[:,chosen]])
            else:
                a,b=c_train,c_test
            model=_fit_xgb(a,b,y[tr])
            pred,_=_model_pred_shap(model,b)
            result[arm][i]=pred[0]
            models[arm]=model
        full=models["clinical_all_384"]
        big_t=np.column_stack([c_test,test])
        _,held=_model_pred_shap(full,big_t)
        heldout_shap[i]=held[0,len(c_train[0]):]
        # Train-only SHAP representation to select five markers for a
        # fresh interpretable Ridge, in this same outer held-out iteration.
        _,train_contrib=_model_pred_shap(full,np.column_stack([c_train,train]))
        train_import=np.abs(train_contrib[:,len(c_train[0]):]).mean(axis=0)
        ranked=np.argsort(-train_import,kind="stable")
        selected=ranked[:5]
        if not np.any(train_import>0):
            raise E3SourceError("XGBoost never used molecular features: cannot run SHAP-guided white-box stage")
        for j in ranked[:10]:
            if train_import[j]>0:fold_rank[int(j)]+=1
        for j in selected:fold_selected[int(j)]+=1
        scaler=StandardScaler().fit(np.column_stack([c_train,train[:,selected]]))
        a=scaler.transform(np.column_stack([c_train,train[:,selected]]))
        b=scaler.transform(np.column_stack([c_test,test[:,selected]]))
        # Low-dimensional fixed alpha penalized regression, no optimization
        # of coefficient stability on heldout patients.
        white=Ridge(alpha=10.0).fit(a,np.log1p(y[tr]))
        result["train_SHAP_top5_Ridge_whitebox"][i]=max(0.0,float(np.expm1(np.clip(white.predict(b)[0],-100,25))))
        selected_fixed=groups["clinical_urea_purines_5"]
        scalar=StandardScaler().fit(np.column_stack([c_train,train[:,selected_fixed]]))
        white_fixed=Ridge(alpha=10.0).fit(
            scalar.transform(np.column_stack([c_train,train[:,selected_fixed]])),np.log1p(y[tr]))
        result["fixed_biochemical5_Ridge_whitebox"][i]=max(
            0.0,float(np.expm1(np.clip(white_fixed.predict(
                scalar.transform(np.column_stack([c_test,test[:,selected_fixed]])))[0],-100,25))))
    scores={}
    for arm,yp in result.items():
        if not np.all(np.isfinite(yp)):raise E3SourceError("Incomplete heldout prediction")
        err=np.abs(yp-y)
        scores[arm]={"mae_original_PK_numeric_units":round(float(err.mean()),5),
                     "rmse_original_PK_numeric_units":round(float(np.sqrt(np.mean((yp-y)**2))),5),
                     "p20_percent":round(float(np.mean(err/y<=.2)*100),2),
                     "bias_original_PK_numeric_units":round(float(np.mean(yp-y)),5)}
    errors={key:np.abs(val-y) for key,val in result.items()}
    control=errors["clinical"]
    differences={}
    for arm,err in errors.items():
        if arm=="clinical":continue
        d=err-control
        rng=Random(SEED)
        boot=sorted(float(np.mean([d[rng.randrange(32)] for _ in range(32)])) for _ in range(1999))
        differences[arm]={"delta_MAE_vs_clinical_negative_is_better":round(float(d.mean()),5),
                          "descriptive_95pct_paired_bootstrap_no_refit":[round(boot[49],5),round(boot[1949],5)],
                          "n_people_with_lower_absolute_error":int(np.sum(d<0))}
    imp=np.mean(np.abs(heldout_shap),axis=0)
    ordered=np.argsort(-imp,kind="stable")
    top=[]
    for j in ordered[:20]:
        if imp[j]<=0:continue
        top.append({
          "source_header_public_chemical_or_unknown_bin":headers[j],
          "annotation_status":"UNRESOLVED_BIN_NOT_MASS" if int(j) in set(annotation["unknown"]) else
                              "SOURCE_CHEMICAL_NAME_NOT_CONFIRMED_BY_STANDARD",
          "mean_abs_heldout_SHAP_on_log1p_PK_model":round(float(imp[j]),7),
          "n_outer_train_folds_selected_in_train_SHAP_top5":int(fold_selected.get(int(j),0)),
          "n_outer_train_folds_in_train_SHAP_top10":int(fold_rank.get(int(j),0)),
          "model_specific_SHAP_not_causal_activity":True,
        })
    total=float(imp.sum())
    group_sums={
        "clinical_4_vars_SHAP_not_in_metabolite_total":True,
        "named_source_chemicals_fraction_omics_abs_SHAP":round(float(imp[annotation["named"]].sum()/total),5) if total>0 else None,
        "numbered_unknown_fraction_omics_abs_SHAP":round(float(imp[annotation["unknown"]].sum()/total),5) if total>0 else None,
        "sum_heldout_omics_SHAP_log_output":round(total,6),
    }
    return {"n_unique_original_patients":32,
            "n_exact_original_LOOCV_test_folds":32,
            "source_PK_outcome_original_label":target,
            "clinical_reference_original_source":"age, sex, weight"+(", serum creatinine" if "creatinine" in clinical_key else ""),
            "fixed_XGB_parameters":PARAMS,
            "models_same_subject_same_outer_folds":scores,
            "paired_deltas_vs_XGB_clinical_reference":differences,
            "top20_original_source_molecules_or_unknown_bin_SHAP":top,
            "heldout_all_384_omics_SHAP_aggregate":group_sums,
            "training_fold_SHAP_top5_for_Ridge_only":True,
            "outer_test_SHAP_never_used_to_select_scored_features":True,
            "clinical_GPx_or_PopPK_available":False,
            "no_patient_values_or_fits_exported":True}

def run_real_source(raw):
    records,source,annotation=qualify_original_feature_labels(raw)
    _,omics=actual_original_tables(raw)
    headers=[str(omics[0][0].get(j,"")).strip() for j in range(1,385)]
    if (len(headers)!=384 or len(set(headers))!=384 or
        len(annotation["named"])!=155 or len(annotation["unknown"])!=229):
        raise E3SourceError("Original spectrometry labels not qualified")
    targets=sorted(records[0]["targets"])
    return {
      "status":"EXECUTED_ORIGINAL_HUMAN_METFORMIN_BLACK_BOX_XGB_TREE_SHAP_AND_FOLDSAFE_WHITE_BOX_EXPLORATORY",
      "original_source_DOI":"10.3389/fphar.2016.00135",
      "original_clinical_SHA256":source["original_clinical_workbook_SHA256"],
      "original_metabolome_SHA256":source["original_untargeted_metabolomic_workbook_SHA256"],
      "source_specimen":"human_plasma_predose_GC_TOF",
      "n_source_omics_features":384,"n_source_named_chemical_headers":155,
      "n_unknown_numeric_bins_unidentified":229,
      "endpoints":{
        target:_loocv(records,annotation,headers,target,"clinical_plus_original_serum_creatinine")
        for target in targets
      },
      "endpoint_table_physical_units":"UNRESOLVED_original_numeric_PK_units",
      "source_actual_proteome_RNA_available":False,
      "actual_genotype_and_PopPK_baseline_available":False,
      "external_clinical_validation":False,
      "SHAP_is_predictor_attribution_not_enzyme_transporter_molecule_identity":True,
      "white_box_features_are_hypotheses_selected_in_train_folds_not_independent_replications":True,
      "individual_PK_or_molecular_data_exported":False,
    }
