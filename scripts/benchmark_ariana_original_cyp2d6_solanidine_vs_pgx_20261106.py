"""Real ORIGINAL 47-patient ARIANA CYP2D6 PGx vs identified plasma solanidine ratios.

Paper Sarömba et al., Front Pharmacol 2026 DOI 10.3389/fphar.2026.1812007.
Original author depositor pinned git tree / XLSX SHA; 47 genuine rows contain
CYP2D6 diplotypes, genotype-derived AS and CNV, identified solanidine-derived
OH-solanidine ratio and semi-quantitative SSDA ratio, plasma metoprolol,
daily dose and clinical covariates. Only aggregate results are emitted.

**Crucial**: Plasma measured during ongoing stable metoprolol therapy, not
predose treatment-naive metabolomics; clinically observed single plasma
metoprolol concentration is NOT a full individual AUC or PK clearance.
Drug-derived OH-metoprolol/metoprolol measured on the SAME sample must NEVER
enter predictors as its denominator would leak the outcome.
Source ARIANA data already published/used by authors -> this is an
INDEPENDENT MODEL REPLAY, not independent external population replication.
"""
from __future__ import annotations
from collections import Counter
from statistics import mean
import json

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler,OneHotEncoder

from scripts.audit_ariana_original_source import (
    _download,_read_shared,_get_workbook_sheets,_read_sheet
)
from hashlib import sha256
from io import BytesIO
from zipfile import ZipFile

EXPECTED_SHA256="7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"
ALPHAS=(.1,1.,10.,100.)
SEED=20261106
ARMS={
    "clinical_only":("age","weight","log_daily_dose"),
    "clinical_PGx_CYP2D6_AS":("age","weight","log_daily_dose","activity_score"),
    "clinical_PGx_plus_OH_solanidine_ratio":(
        "age","weight","log_daily_dose","activity_score","ln_oh_solanidine_ratio"),
    "clinical_PGx_plus_SSDA_solanidine_ratio":(
        "age","weight","log_daily_dose","activity_score","ln_ssda_solanidine_ratio"),
    "clinical_PGx_plus_both_identified_metabolic_ratios":(
        "age","weight","log_daily_dose","activity_score",
        "ln_oh_solanidine_ratio","ln_ssda_solanidine_ratio"),
    "clinical_plus_OH_solanidine_without_PGx":(
        "age","weight","log_daily_dose","ln_oh_solanidine_ratio"),
}
CATEGORICAL=("sex","CYP2D6_inhibitor","copy_number_variation")

class OriginalARIANADataError(ValueError):
    pass

def _float(s):
    s=str(s or "").strip().replace("\u00a0","").replace(",",".")
    if s in ("","-","NA","n/a","NaN","None","missing"):
        return np.nan
    try:return float(s)
    except ValueError:return np.nan

def _load_ari_47():
    raw=_download()
    if sha256(raw).hexdigest()!=EXPECTED_SHA256:
        raise OriginalARIANADataError("Original author deposited XLSX source bytes changed")
    with ZipFile(BytesIO(raw)) as z:
        shared=_read_shared(z)
        sheets=_get_workbook_sheets(z)
        if len(sheets)!=1 or sheets[0][0]!="Tabelle1":
            raise OriginalARIANADataError("Original study author workbook sheet layout changed")
        rows=_read_sheet(z,sheets[0][1],shared)
    if len(rows)!=48:
        raise OriginalARIANADataError("Original author study 47 patient source rows changed")
    headings=["ID","Allel 1","Allel 2","AS","CNV",
        "c(Metoprolol)[ng/ml]","Daily Dose Metoprolol",
        "c(metoprolol)weighted [ng/ml]",
        "lnMR OH-metoprolol/metoprolol",
        "lnMR OH-solanidine/solanidine",
        "lnMR SSDA/solanidine",
        "Age (years)","male sex","BMI","weight","CYP2D6 inhibitor intake"]
    actual=[" ".join(str(rows[0].get(i,"")).split()) for i in range(16)]
    if actual!=headings:
        raise OriginalARIANADataError("Author clinical+omics+genotype XLSX heading source contract changed")
    records=[]
    ids=[]
    for row in rows[1:]:
        ids.append(str(row.get(0,"")).strip())
        y=_float(row.get(5))
        dose=_float(row.get(6))
        if not np.isfinite(y) or y<=0 or not np.isfinite(dose) or dose<=0:
            raise OriginalARIANADataError("Missing metoprolol concentration/dose in reported 47")
        r={
          "plasma_metoprolol_ng_ml":y,
          "age":_float(row.get(11)),
          "weight":_float(row.get(14)),
          "log_daily_dose":float(np.log(dose)),
          "activity_score":_float(row.get(3)),
          "ln_oh_solanidine_ratio":_float(row.get(9)),
          "ln_ssda_solanidine_ratio":_float(row.get(10)),
          "sex":str(row.get(12,"")).strip(),
          "CYP2D6_inhibitor":str(row.get(15,"")).strip(),
          "copy_number_variation":str(row.get(4,"")).strip(),
        }
        # Original free-text CNV labels are categorical, never numerically
        # imputed or guessed as diplotypes; AS is original source PGx consensus.
        if not r["sex"] or not r["CYP2D6_inhibitor"]:
            raise OriginalARIANADataError("Unqualified original clinical categorical values")
        records.append(r)
    if len(set(ids))!=47 or any(not x for x in ids):
        raise OriginalARIANADataError("Original author study person IDs not exactly 47 unique")
    # We don't train an imputed all-n47 model if an original ratio is absent
    # (different arms would have different heldout biological populations).
    numeric=("age","weight","log_daily_dose","activity_score",
             "ln_oh_solanidine_ratio","ln_ssda_solanidine_ratio")
    qualified=[r for r in records if all(np.isfinite(r[k]) for k in numeric)]
    if not 30<=len(qualified)<=47:
        raise OriginalARIANADataError("Too few fully observed original patient PGx/ratio rows")
    # No ID persisted in qualified records or model output.
    return qualified,{"n_original_source_rows":47,
        "n_complete_case_human_patients":len(qualified),
        "n_excluded_missing_any_clinical_AS_or_omics":47-len(qualified),
        "source_CNV_n_categories":len(set(r["copy_number_variation"] for r in qualified)),
        "source_genotype_CYP2D6_activity_score_distinct_values":len(set(r["activity_score"] for r in qualified)),
        "source_plasma_omics":{"OH_solanidine":"identified and LC-MS targeted","SSDA":"semi-quantitative, no reference standard"},
        "drug_derived_OH_metoprolol_ratio_excluded_to_avoid_outcome_leakage":True,
    }

def _fit_predict(rows,cols,train,test,alpha):
    if alpha not in ALPHAS:
        raise OriginalARIANADataError("Tuning outside predefined Ridge grid")
    if set(train)&set(test):
        raise OriginalARIANADataError("Source patient leakage")
    y=np.log(np.asarray([r["plasma_metoprolol_ng_ml"] for r in rows]))
    x=np.array([[r[col] for col in cols] for r in rows],float)
    scaler=StandardScaler()
    xt=scaler.fit_transform(x[train])
    xv=scaler.transform(x[test])
    # Copy-number variation, sex and inhibitor as source categorical effects,
    # not biological surrogate metabolites, and one-hot fit training-only.
    encoder=OneHotEncoder(handle_unknown="ignore")
    cat=np.array([[r[col] for col in CATEGORICAL] for r in rows],object)
    a=encoder.fit_transform(cat[train]).toarray()
    b=encoder.transform(cat[test]).toarray()
    Xtr=np.column_stack((xt,a))
    Xte=np.column_stack((xv,b))
    model=Ridge(alpha=alpha)
    model.fit(Xtr,y[train])
    return np.exp(np.clip(model.predict(Xte),-8,14))

def _cv(rows):
    n=len(rows)
    if not 30<=n<=47:
        raise OriginalARIANADataError("Unqualified real original patient sample")
    y=np.asarray([r["plasma_metoprolol_ng_ml"] for r in rows])
    outcomes={}
    private_losses={}
    for name,features in ARMS.items():
        yhat=np.empty(n)
        tuned=Counter()
        for held in range(n):
            train=np.array([i for i in range(n) if i!=held],int)
            test=np.array([held],int)
            candidates=[]
            for alpha in ALPHAS:
                errs=[]
                kf=KFold(n_splits=4,shuffle=True,random_state=SEED)
                for it,iv in kf.split(train):
                    prediction=_fit_predict(rows,features,train[it],train[iv],alpha)
                    truth=np.log(y[train[iv]])
                    errs.extend(np.abs(np.log(prediction)-truth))
                candidates.append((float(mean(errs)),-alpha,alpha))
            _,_,chosen=min(candidates)
            tuned[str(chosen)]+=1
            yhat[held]=_fit_predict(rows,features,train,test,chosen)[0]
        errors=np.abs(yhat-y)
        private_losses[name]=errors
        outcomes[name]={
            "n_real_independent_original_humans":n,
            "outer_evaluation":"leave_one_original_patient_out",
            "inner_Ridge_regularization":"train-only 4fold, alpha .1 1 10 100, log-concentration MAE",
            "MAE_observed_plasma_metoprolol_ng_ml":round(float(errors.mean()),6),
            "mean_absolute_log_concentration_error":round(float(np.mean(np.abs(np.log(yhat)-np.log(y)))),6),
            "p20_relative_error":round(float(np.mean(errors/y<=.2)),6),
            "inner_alpha_counts":dict(sorted(tuned.items())),
            "source_ratio_screening_is_mechanism_prespecified_not_supervised":True,
            "all_variants_same_complete_cases":True,
            "original_individual_patient_values_never_exported":True,
        }
    ref=private_losses["clinical_PGx_CYP2D6_AS"]
    rng=np.random.default_rng(SEED+len(y))
    for name,e in private_losses.items():
        if name=="clinical_PGx_CYP2D6_AS":continue
        delta=e-ref
        draws=np.mean(delta[rng.integers(0,len(y),size=(1999,len(y)))],axis=1)
        outcomes[name]["delta_MAE_vs_clinical_and_genotype_negative_better_ng_ml"]=round(float(delta.mean()),6)
        outcomes[name]["paired_original_patient_loss_bootstrap_95pct_descriptive"]=[
            round(float(q),6) for q in np.quantile(draws,[.025,.975])]
        outcomes[name]["bootstrap_not_nested_model_refit_or_multiple_testing_corrected"]=True
    return outcomes

def original_benchmark():
    rows,qualification=_load_ari_47()
    return {
      "status":"ACTUAL_PUBLIC_ARIANA_ORIGINAL_CYP2D6_GENOTYPE_VS_IDENTIFIED_SOLANIDINE_PATIENT_HELDOUT",
      "original_study":"Sarömba et al Front Pharmacol 2026 10.3389/fphar.2026.1812007",
      "original_xlsx_SHA256":EXPECTED_SHA256,
      "source_original_author_repository_committed":True,
      "outcome":"plasma metoprolol ng/mL on stable metoprolol BID, not individual PK AUC or CL",
      "sampling":"plasma during stable metoprolol treatment; contemporaneous endogenous/exogenous metabolic proxy, not treatment-naive pre-exposure",
      "solanidine_is_plant_derived_dietary_alkaloid_not_endogenous":True,
      "SSDA_lacks_authentic_standard_and_is_semiquantitative":True,
      "original_source_PGx":"CYP2D6 activity score and source measured CNV category, not full resequencing",
      "source_qualification":qualification,
      "nested_independent_patient_validation_results":_cv(rows),
      "exploratory_secondary_reanalysis_not_new_external_population":True,
      "original_patient_identifiers_and_records_not_printed_or_committed":True,
    }

if __name__=="__main__":
    try:print(json.dumps(original_benchmark(),sort_keys=True,indent=2))
    except (OriginalARIANADataError,ValueError,OSError) as exc:
        print(json.dumps({"status":"ABSTAIN_SOURCE_SCHEMA_OR_NUMERIC_OR_QUALIFICATION",
                          "reason_class":type(exc).__name__}))
        raise SystemExit(2)
