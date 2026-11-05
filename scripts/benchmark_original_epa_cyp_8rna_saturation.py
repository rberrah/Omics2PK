#!/usr/bin/env python3
"""Real US EPA n100 human liver CYP source: target RNA vs eight-RNA saturation.

The measured endpoint is *in vitro* CYP enzyme probe activity, not drug PK.
This follows examination of prior original EPA source runs and is exploratory.
Original EPA xlsx is downloaded into memory, SHA256 verified, ID-joined using
both original donor keys, and never copied to repo or Actions artifacts.
Output ONLY aggregate CV errors, count and descriptive uncertainty.
"""
from __future__ import annotations

import hashlib
import json
from math import isfinite
from statistics import mean

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

from scripts.audit_epa_liver_rna_activity import _get, SOURCE_DOI, PAPER_DOI
from scripts.inspect_epa_liver_crosswalk import parse_epa_workbook
from scripts.validate_epa_liver_rna_cyp import EXPECTED_SHA
from omics2pk.epa_liver_cyp_validation import (
    EXPECTED_GENES, EPALinkageError, _panel, _numeric, read_paired_epa,
)

GENES = tuple(EXPECTED_GENES)
SEED = 20261105
ALPHAS = (1.0, 100.0)
ARMS = ("clinical", "clinical_plus_target_RNA",
        "clinical_plus_non_target_RNA", "clinical_plus_all_eight_RNAs")


def _target_rows(original, measured, target):
    if target not in GENES:
        raise EPALinkageError("Unexpected CYP endpoint")
    rna, cols, _, _ = _panel(*original["Table 3"])
    if not set(GENES).issubset(cols):
        raise EPALinkageError("Eight original CYP assays absent")
    rows = []
    for key, target_rna, activity, age, male in measured[target]:
        original_rna_row = rna[key]
        values = [_numeric(original_rna_row[cols[gene]])
                  if cols[gene] < len(original_rna_row) else None
                  for gene in GENES]
        if any(v is None or v <= 0 for v in values):
            continue
        if not all(isfinite(z) for z in (age, male, target_rna, activity)):
            raise EPALinkageError("Invalid source measured phenotype")
        if abs(values[GENES.index(target)] - target_rna) > 1e-9:
            raise EPALinkageError("Target transcript diverged from original paired source")
        rows.append((key, float(age), float(male), float(activity),
                     *[float(np.log1p(v)) for v in values]))
    if len(rows) < 60 or len(rows) > 100 or len(set(r[0] for r in rows)) != len(rows):
        raise EPALinkageError("Too few complete original donor RNA panels")
    return sorted(rows, key=lambda x: x[0])


def _metrics_cv(rows, target, repeats=3):
    if target not in GENES or repeats not in (1, 2, 3):
        raise EPALinkageError("Unqualified endpoint/CV specification")
    names=[f"{arm}|alpha_{int(a)}" for arm in ARMS for a in ALPHAS]
    n=len(rows)
    X=np.asarray([r[1:] for r in rows],dtype=float)
    if X.shape!=(n,3+len(GENES)) or not np.isfinite(X).all():
        raise EPALinkageError("Source X dimensions or values invalid")
    age_sex=X[:,:2]
    measured_y=X[:,2]
    if np.any(measured_y<0):
        raise EPALinkageError("Observed CYP activity below zero")
    rnas=X[:,3:]
    tindex=GENES.index(target)
    y=np.log1p(measured_y)
    inputs={
        "clinical":age_sex,
        "clinical_plus_target_RNA":np.column_stack([age_sex,rnas[:,tindex]]),
        "clinical_plus_non_target_RNA":
            np.column_stack([age_sex,np.delete(rnas,tindex,axis=1)]),
        "clinical_plus_all_eight_RNAs":np.column_stack([age_sex,rnas]),
    }
    errors={k:np.full((repeats,n),np.nan) for k in names}
    for rep in range(repeats):
        kfold=KFold(n_splits=5,shuffle=True,random_state=SEED+rep)
        for tr,te in kfold.split(X):
            for arm,x in inputs.items():
                scale=StandardScaler()
                a=scale.fit_transform(x[tr])
                b=scale.transform(x[te])
                for alpha in ALPHAS:
                    fitted=Ridge(alpha=alpha).fit(a,y[tr])
                    pred=np.expm1(np.minimum(fitted.predict(b),30))
                    pred=np.maximum(pred,0.0)
                    errors[f"{arm}|alpha_{int(alpha)}"][rep,te]=(
                        np.abs(pred-measured_y[te]))
    if any(not np.isfinite(v).all() for v in errors.values()):
        raise EPALinkageError("Missing/nonfinite held-out original donor errors")
    person={name:np.mean(mat,axis=0) for name,mat in errors.items()}
    scores={k:round(float(v.mean()),6) for k,v in person.items()}
    def compare(a,b,seed):
        delta=person[a]-person[b]
        rng=np.random.default_rng(seed)
        draws=np.mean(delta[rng.integers(0,n,size=(1999,n))],axis=1)
        return {
            "delta_MAE_tested_minus_reference_negative_better":
                round(float(np.mean(delta)),6),
            "descriptive_95pct_paired_person_residual_bootstrap":
                [round(float(v),6) for v in np.quantile(draws,[.025,.975])],
            "n_unique_independent_persons_bootstrapped":n,
            "bootstrap_does_not_refit_models_or_adjust_multiplicity":True,
        }
    contrasts={}
    for alpha in ALPHAS:
        pref=f"|alpha_{int(alpha)}"
        contrasts[f"target_RNA_vs_clinical{pref}"]=compare(
            "clinical_plus_target_RNA"+pref,"clinical"+pref,seed=SEED+int(alpha))
        contrasts[f"all_eight_RNAs_vs_target_RNA{pref}"]=compare(
            "clinical_plus_all_eight_RNAs"+pref,"clinical_plus_target_RNA"+pref,
            seed=SEED+int(alpha)+900)
        contrasts[f"non_target_RNAs_vs_clinical{pref}"]=compare(
            "clinical_plus_non_target_RNA"+pref,"clinical"+pref,
            seed=SEED+int(alpha)+1800)
    return {
        "n_original_independent_donors":n,
        "heldout_folds":f"{repeats}x5_original_person_disjoint",
        "target":"original measured CYP microsomal probe enzyme activity",
        "source_protein_tissue":"human liver microsomes",
        "predictor_tissue":"same human liver bulk qPCR",
        "matched_CYP":target,
        "eight_original_CYP_RNA_assay_genes":list(GENES),
        "fixed_ridge_alphas":list(ALPHAS),
        "MAE_in_original_activity_units":scores,
        "exploratory_paired_person_contrasts":contrasts,
        "per_donor_values_and_IDs_emitted":False,
    }


def evaluate_original(original,*,repeats=3):
    measurements,_=read_paired_epa(original)
    results={}
    for gene in GENES:
        results[gene]=_metrics_cv(
            _target_rows(original,measurements,gene),gene,repeats=repeats)
    return results


def live_original():
    original_bytes=_get()
    if hashlib.sha256(original_bytes).hexdigest()!=EXPECTED_SHA:
        raise EPALinkageError("EPA official original source SHA changed")
    original=parse_epa_workbook(original_bytes)
    return {
        "status":"REAL_EPA_ORIGINAL_LIVER_TARGET_VS_EIGHT_RNAs_SATURATION_EXPLORATORY_NOT_PK",
        "official_source":SOURCE_DOI,
        "original_publication":PAPER_DOI,
        "verified_source_SHA256":EXPECTED_SHA,
        "n_original_donors_in_source":100,
        "n_original_cyp_probe_endpoints":len(GENES),
        "data_driven_tuning_or_target_outcome_feature_selection":False,
        "design_selected_after_previous_source_results":True,
        "8_endpoints_related_same_liver_bank_no_multitest_adjustment":True,
        "not_human_AUC_CL_or_genotype_proof":True,
        "source_individual_measurements_never_written":True,
        "results":evaluate_original(original),
    }


if __name__=="__main__":
    try:
        print(json.dumps(live_original(),indent=2,sort_keys=True))
    except (EPALinkageError,ValueError,OSError) as e:
        print(json.dumps({"status":"ABSTAIN_SOURCE_OR_EIGHT_RNA_PANEL",
                          "reason_category":type(e).__name__,
                          "original_EPA_model_completed":False}))
        raise SystemExit(2)
