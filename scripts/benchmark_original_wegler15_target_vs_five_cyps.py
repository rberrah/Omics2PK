#!/usr/bin/env python3
"""Real Wegler 2021 original 15-donor CYP proteome -> 5-drug intrinsic clearance.

Compare a mechanistically chosen corresponding CYP vs all five probe CYP
proteins vs prespecified strongly regularized tree ensemble. Molecular
proteins and *independent in-vitro substrate depletion* are observed on
matched human liver donors. This is NOT clinical PK, nor multiview omics.
Official primary PMC AWS supplements are SHA256-pinned by the existing
strict source reader. Donor values/identifiers never enter logs or git.
"""
from __future__ import annotations

import io
import json
from math import log
from statistics import mean
from collections import defaultdict

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

from scripts.audit_wegler_source_structure import (
    DRUG_TERMS, RAW_CLINT_MAPPING, _parse_wegler_kinetics_records,
    _validate_raw_clint_column_mapping, _paired_loss_bootstrap_ci,
    get_verified_supplements, EXPECTED_HASHES,
)

SEED=20261105
ALPHAS=(1.0, 20.0)
MODELS=("training_mean","target_CYP_ridge1","five_CYP_ridge1",
        "five_CYP_ridge20","non_target_four_CYP_ridge1",
        "five_CYP_regularized_XGB")
SYSTEMS=("HH","HLM")
# Probe-specific CYP protein label in the authenticated Table S3 drug block.
CYP_PROBES={
    "midazolam":"CYP3A4",
    "omeprazole":"CYP2C19",
    "bupropion":"CYP2B6",
    "diclofenac":"CYP2C9",
    "bufuralol":"CYP2D6",
}


class WeglerPanelAbstention(ValueError):
    pass


def _source_panel(records,system,drug):
    if system not in SYSTEMS or drug not in DRUG_TERMS:
        raise WeglerPanelAbstention("Unknown source assay system/substrate")
    pkey,ckey=RAW_CLINT_MAPPING[system]
    by_drug={d:{} for d in DRUG_TERMS}
    for row in records:
        d=row["drug"]
        if d not in by_drug:
            raise WeglerPanelAbstention("Unexpected original substrate")
        donor=row["donor_id"]
        if donor in by_drug[d]:
            raise WeglerPanelAbstention("Duplicate source donor-drug protein")
        by_drug[d][donor]=row
    shared=set.intersection(*[set(v) for v in by_drug.values()])
    if len(shared)!=15:
        raise WeglerPanelAbstention("Expected 15 identical original source donors")
    rows=[]
    for donor in sorted(shared):
        yrow=by_drug[drug][donor]
        y=yrow[ckey]
        if y is None:
            continue  # Original ND: exclude corresponding donor for ALL arms, never zero.
        if not np.isfinite(y) or float(y)<=0:
            raise WeglerPanelAbstention("Invalid observed intrinsic clearance")
        five=[]
        for other in DRUG_TERMS:
            p=by_drug[other][donor][pkey]
            # No unknown LLOQ protein replacement, and original protein >0.
            if p is None or not np.isfinite(p) or float(p)<=0:
                five=[]
                break
            five.append(log(float(p)))
        if len(five)==len(DRUG_TERMS):
            rows.append((donor,np.array(five,dtype=float),log(float(y))))
    if not 8<=len(rows)<=15:
        raise WeglerPanelAbstention("Original complete-case panel too small")
    ids=[r[0] for r in rows]
    if len(set(ids))!=len(ids):
        raise WeglerPanelAbstention("Duplicate original donor")
    X=np.stack([r[1] for r in rows])
    y=np.array([r[2] for r in rows],dtype=float)
    idx=DRUG_TERMS.index(drug)
    return ids,X,y,idx


def _fit_heldout(X,y,idx, *, seed=SEED):
    n=X.shape[0]
    if X.ndim!=2 or X.shape[1]!=5 or y.shape!=(n,) or not 8<=n<=15:
        raise WeglerPanelAbstention("Unqualified per-substrate matrix dimensions")
    if not (np.isfinite(X).all() and np.isfinite(y).all()) or idx not in range(5):
        raise WeglerPanelAbstention("Nonfinite molecular input/outcome")
    errors={key:np.empty(n) for key in MODELS}
    for held in range(n):
        train=np.array([i for i in range(n) if i!=held],dtype=int)
        test=np.array([held])
        mu=float(y[train].mean())
        errors["training_mean"][held]=abs(float(y[held])-mu)
        for name,cols,alpha in (
            ("target_CYP_ridge1",[idx],1.),
            ("five_CYP_ridge1",list(range(5)),1.),
            ("five_CYP_ridge20",list(range(5)),20.),
            ("non_target_four_CYP_ridge1",[i for i in range(5) if i!=idx],1.)):
            scaler=StandardScaler()
            tr=scaler.fit_transform(X[train][:,cols])
            te=scaler.transform(X[test][:,cols])
            m=Ridge(alpha=alpha).fit(tr,y[train])
            errors[name][held]=abs(float(m.predict(te)[0])-float(y[held]))
        m=XGBRegressor(n_estimators=60,learning_rate=.05,max_depth=1,
                       min_child_weight=3,subsample=1.,colsample_bytree=1.,
                       reg_lambda=12.,reg_alpha=0.,objective="reg:squarederror",
                       n_jobs=1,random_state=seed,tree_method="hist")
        # Tree models need no scaling, and no feature selection on test.
        m.fit(X[train],y[train])
        errors["five_CYP_regularized_XGB"][held]=abs(
            float(m.predict(X[test])[0])-float(y[held]))
    if any(not np.isfinite(v).all() for v in errors.values()):
        raise WeglerPanelAbstention("Held-out error is not finite")
    results={name:round(float(err.mean()),6) for name,err in errors.items()}
    ci={}
    for name in MODELS[1:]:
        delta=(errors[name]-errors["target_CYP_ridge1"]).tolist()
        ci[name]={
            "difference_logMAE_vs_target_CYP_negative_better":
                round(mean(delta),6),
            "paired_donor_residual_bootstrap95":
                _paired_loss_bootstrap_ci(delta,seed=seed+idx+len(name),draws=1999),
            "bootstrap_does_not_refit_models_or_correct_multiplicity":True,
        }
    return {
        "n_independent_original_donors":n,
        "heldout_evaluation":"leave-one-biological-donor-out",
        "MAE_log_original_intrinsic_clearance":results,
        "paired_within_drug_delta_vs_target_CYP":ci,
        "number_target_only_CYP_features":1,
        "number_broad_CYP_features":5,
        "same_complete_case_donor_set_all_models":True,
        "drug_clint_not_computed_from_CYP_protein_predictor":True,
        "protein_quantification_cannot_substitute_active_enzyme_fraction":True,
        "no_individual_original_donor_data_emitted":True,
    }


def benchmark_records(records):
    _validate_raw_clint_column_mapping(records)
    out={}
    for system in SYSTEMS:
        out[system]={}
        for drug in DRUG_TERMS:
            ids,X,y,idx=_source_panel(records,system,drug)
            entry=_fit_heldout(X,y,idx,seed=SEED+sum(ord(c) for c in (system+drug)))
            entry["source_target_probed_CYP"]=CYP_PROBES[drug]
            entry["source_assay"]="human_hepatocytes" if system=="HH" else "liver_microsomes"
            entry["in_vitro_clearance_source_column"]=RAW_CLINT_MAPPING[system][1]
            out[system][drug]=entry
    return out


def run_original():
    from pypdf import PdfReader
    original=get_verified_supplements()
    source=PdfReader(io.BytesIO(original["wegler_si_methods_pdf"]))
    pages=[page.extract_text() or "" for page in source.pages]
    records=_parse_wegler_kinetics_records(pages)
    return {
        "status":"REAL_WEGLER_15_DONOR_TARGET_VS_5_PROTEINS_VS_XGB_NOT_HUMAN_PK",
        "source_publication_DOI":"10.1021/acs.molpharmaceut.1c00053",
        "verified_source_PDF_SHA256":EXPECTED_HASHES["wegler_si_methods_pdf"],
        "original_study_source_matched_donors":15,
        "study_contrast_is_posthoc_exploratory":True,
        "mechanistic_one_gene_CYP_probe_is_simplification_CYP3A4_not_specific":True,
        "no_original_patient_PK_or_full_PGx_or_multilayer_omics":True,
        "data_processing_fold_frozen":True,
        "no_source_donor_IDs_or_measurements_printed":True,
        "results":benchmark_records(records)
    }


if __name__=="__main__":
    try:
        print(json.dumps(run_original(),indent=2,sort_keys=True))
    except (WeglerPanelAbstention,ValueError,OSError) as err:
        print(json.dumps({"status":"ABSTAIN_UNQUALIFIED_SOURCE_OR_PROTEIN_PANEL",
                          "reason_class":type(err).__name__}))
        raise SystemExit(2)
