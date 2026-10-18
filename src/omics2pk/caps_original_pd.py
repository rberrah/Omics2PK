"""Original-public CaPS/RACE baseline pharmacometabolomics to 90-day mortality.

Source: Puskarich et al. Clin Transl Sci 2021 DOI 10.1111/cts.13088,
original author release UMichNMR-Metabolomics/CaPS_TimeZero,
commit 2dcacb5d82fe07cf625ae0e4edab34e51546f3e6.

**This is a PD / treatment-outcome RESEARCH proof**, not pharmacokinetics:
RACE randomised placebo vs L-carnitine at doses recorded in original
demographics, pretreatment NMR metabolites and LC-MS acylcarnitines,
90-day observed mortality. Source input is loaded transiently. Original
person IDs, values, predictions and original files NEVER emitted to logs
or committed to GitHub.

Both metabolomics assays are from the METABOLOME (different analytical
platforms); combining them is multi-view metabolomics, NOT demonstrated
multi-omic transcriptome+proteome fusion.

Five held-PERSON folds with fold-local missing value handling, scaling and
PCA; all model arms use EXACT SAME complete-case original 2-platform
person subset. Results describe discrimination, NOT patient-level causal
individual treatment benefit or externally validated clinical prediction.
"""
from __future__ import annotations

import csv
from io import StringIO
from math import isfinite
from hashlib import sha256
from urllib.request import Request,urlopen

SOURCE_COMMIT="2dcacb5d82fe07cf625ae0e4edab34e51546f3e6"
SOURCE_BASE="https://raw.githubusercontent.com/UMichNMR-Metabolomics/CaPS_TimeZero/"+SOURCE_COMMIT+"/"
SOURCE_FILES=("RACE_Demographics.csv","raw_nmr.csv","raw_ac.csv")
SEED=20261018
CV=5
PCS_PER_ASSAY=3


class OriginalPDSourceError(ValueError):
    pass


def _numerical(value):
    if value is None or str(value).strip()=="":
        return float("nan")
    try:
        result=float(value)
    except (TypeError,ValueError):
        return float("nan")
    return result if isfinite(result) else float("nan")


def _unique_rows(csv_content,key):
    rows=list(csv.DictReader(StringIO(csv_content)))
    if not rows or key not in rows[0]:
        raise OriginalPDSourceError("Missing original study person-key field")
    selected={}
    for row in rows:
        person=(row.get(key) or "").strip()
        if not person or person in selected:
            raise OriginalPDSourceError("Duplicate/missing original independent donor ID")
        selected[person]=row
    return selected


def _download_author_original(filename):
    if filename not in SOURCE_FILES:
        raise OriginalPDSourceError("Unexpected original-source file")
    try:
        with urlopen(Request(SOURCE_BASE+filename,headers={
            "User-Agent":"Omics2PK-CaPS-RACE-public-aggregate-only-PD/1.0",
        }),timeout=40) as response:
            original=response.read(3_500_001)
    except OSError as err:
        raise OriginalPDSourceError("Original public RACE study unavailable") from err
    if len(original)>3_500_000 or not original:
        raise OriginalPDSourceError("Source size outside original study limits")
    return original.decode("utf-8-sig"),sha256(original).hexdigest()


def load_original_caps(*,files=None):
    """Public original rows are kept only in RAM; callable with fictional tests."""
    import numpy as np
    if files is None:
        files={f:_download_author_original(f)[0] for f in SOURCE_FILES}
    if set(files)!=set(SOURCE_FILES):
        raise OriginalPDSourceError("All three original assay/outcome tables required")
    demographics=_unique_rows(files["RACE_Demographics.csv"],"studyid")
    nmr=_unique_rows(files["raw_nmr.csv"],"studyid")
    ac=_unique_rows(files["raw_ac.csv"],"studyid")
    # Original source declares raw_nmr first two columns as studyid and
    # t0nmronly; original AC first two as timepoint, studyid. The latter
    # must be PRETREATMENT to avoid treatment-effect/postdose leakage.
    nmr_header=list(next(iter(nmr.values())))
    ac_header=list(next(iter(ac.values())))
    if nmr_header[:2]!=["studyid","t0nmronly"] or ac_header[:2]!=["timepoint","studyid"]:
        raise OriginalPDSourceError("Original metabolite feature schema changed")
    nmr_features=nmr_header[2:]
    ac_features=ac_header[2:]
    if len(nmr_features)!=27 or len(ac_features)!=24:
        raise OriginalPDSourceError("Original assay feature dimensions changed")
    if len(nmr)<200 or len(ac)<200 or len(demographics)<200:
        raise OriginalPDSourceError("Original study source too small")
    common=sorted(set(demographics)&set(nmr)&set(ac))
    original=[]
    for subject in common:
        meta=demographics[subject]
        # RACE dose 0 for placebo; nonzero signifies randomized carnitine.
        dose=_numerical(meta.get("dose"))
        arm=(meta.get("carntx") or "").strip()
        if (dose==0 and arm!="Placebo") or (dose>0 and arm!="CarnTx"):
            raise OriginalPDSourceError("Original dose / randomized arm mismatch")
        if dose not in (0,6,12,18):
            raise OriginalPDSourceError("Unexpected original study dose")
        if (ac[subject].get("timepoint") or "").strip().casefold() not in (
            "t0","time0","baseline","0"):
            raise OriginalPDSourceError("LC-MS sample not explicitly prespecified predose T0")
        age=_numerical(meta.get("Age"))
        sofa=_numerical(meta.get("Calculated.T0.SOFA.Value"))
        survival=_numerical(meta.get("survivaldays"))
        sex={"Male":1.,"Female":0.}.get(meta.get("Sex"))
        if (sex is None or not 18<=age<=110 or not 0<=sofa<=24
            or not 0<=survival<=366):
            continue
        nmr_values=[_numerical(nmr[subject][f]) for f in nmr_features]
        ac_values=[_numerical(ac[subject][f]) for f in ac_features]
        for values in (nmr_values,ac_values):
            if any(isfinite(x) and x<0 for x in values):
                raise OriginalPDSourceError("Negative original metabolite concentration")
            if sum(isfinite(x) for x in values)<len(values)//2:
                raise OriginalPDSourceError("Original assay has too many missing values")
        original.append((age,sex,sofa,dose,1. if dose>0 else 0.,
                         float(survival<=90),nmr_values,ac_values))
    if len(original)<185:
        raise OriginalPDSourceError("Original same-person NMR/LCMS/PD intersection too small")
    clinical=np.asarray([v[:5] for v in original],dtype=float)
    outcome=np.asarray([v[5] for v in original],dtype=int)
    nmr_matrix=np.asarray([v[6] for v in original],dtype=float)
    ac_matrix=np.asarray([v[7] for v in original],dtype=float)
    return {
        "clinical":clinical,"outcome":outcome,
        "nmr":nmr_matrix,"acylcarnitines":ac_matrix,
        "source_counts":{
            "demographics_source_people":len(demographics),
            "NMR_source_people":len(nmr),
            "LCMS_acylcarnitine_source_people":len(ac),
            "n_true_same_original_people_both_omics_and_outcome":len(original),
            "n_NMR_metabolites":len(nmr_features),
            "n_LCMS_acylcarnitines":len(ac_features),
        },
    }


def _safe_pca_fit_transform(train,test):
    import numpy as np
    from sklearn.impute import SimpleImputer
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    # raw pre-treatment metabolite concentration with missing and possible 0
    transform=lambda a:np.log1p(a)
    imp=SimpleImputer(strategy="median",keep_empty_features=True)
    a=imp.fit_transform(transform(train))
    b=imp.transform(transform(test))
    scale=StandardScaler().fit(a)
    a=scale.transform(a);b=scale.transform(b)
    pca=PCA(n_components=PCS_PER_ASSAY,random_state=SEED,
            svd_solver="full")
    return pca.fit_transform(a),pca.transform(b)


def evaluate_original_caps(data,*,folds=CV,seed=SEED):
    import numpy as np
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score,brier_score_loss,log_loss
    c=data["clinical"]
    y=data["outcome"]
    nmr=data["nmr"]
    ac=data["acylcarnitines"]
    if (c.ndim!=2 or c.shape[1]!=5 or
        not len(c)==len(y)==len(nmr)==len(ac) or len(y)<185 or
        not set(np.unique(y))=={0,1} or
        not np.isfinite(c).all()):
        raise OriginalPDSourceError("Not the original two-platform matched public PD study")
    strat=y*2+c[:,4].astype(int)
    cv=StratifiedKFold(n_splits=folds,shuffle=True,random_state=seed)
    arms=("clinical_only","clinical_plus_NMR","clinical_plus_LCMS",
          "clinical_plus_both","clinical_plus_both_treatment_interactions")
    predicted={a:np.full(len(y),np.nan) for a in arms}
    for train,test in cv.split(c,strat):
        standard=StandardScaler().fit(c[train])
        clinical_train=standard.transform(c[train])
        clinical_test=standard.transform(c[test])
        nmr_train,nmr_test=_safe_pca_fit_transform(nmr[train],nmr[test])
        ac_train,ac_test=_safe_pca_fit_transform(ac[train],ac[test])
        train_treat=c[train,4:5]
        test_treat=c[test,4:5]
        specs={
            "clinical_only":(clinical_train,clinical_test),
            "clinical_plus_NMR":(
                np.column_stack([clinical_train,nmr_train]),
                np.column_stack([clinical_test,nmr_test])),
            "clinical_plus_LCMS":(
                np.column_stack([clinical_train,ac_train]),
                np.column_stack([clinical_test,ac_test])),
            "clinical_plus_both":(
                np.column_stack([clinical_train,nmr_train,ac_train]),
                np.column_stack([clinical_test,nmr_test,ac_test])),
            "clinical_plus_both_treatment_interactions":(
                np.column_stack([clinical_train,nmr_train,ac_train,
                                 (nmr_train*train_treat),(ac_train*train_treat)]),
                np.column_stack([clinical_test,nmr_test,ac_test,
                                 (nmr_test*test_treat),(ac_test*test_treat)])),
        }
        for arm,(Xtr,Xte) in specs.items():
            m=LogisticRegression(C=1.0,max_iter=2000,solver="lbfgs")
            m.fit(Xtr,y[train])
            predicted[arm][test]=m.predict_proba(Xte)[:,1]
    if any(not np.isfinite(v).all() for v in predicted.values()):
        raise OriginalPDSourceError("Unobserved heldout mortality probabilities")
    metric={}
    for arm,p in predicted.items():
        metric[arm]={
            "AUROC_original_person_heldout":round(float(roc_auc_score(y,p)),6),
            "Brier_90day_mortality_probability":round(float(brier_score_loss(y,p)),6),
            "logloss_90day_mortality_probability":round(float(log_loss(y,p)),6),
        }
    base=(y-predicted["clinical_only"])**2
    rng=np.random.default_rng(seed)
    comparisons={}
    for arm in arms[1:]:
        losses=(y-predicted[arm])**2
        delta=losses-base
        boot=rng.choice(delta,size=(1000,len(delta)),replace=True).mean(axis=1)
        comparisons[arm]={
            "delta_Brier_vs_clinical_only":round(float(delta.mean()),6),
            "descriptive_paired_95pct_bootstrap_interval_delta_Brier":[
                round(float(v),6) for v in np.quantile(boot,[.025,.975])
            ],
        }
    return {
        "status":"ACTUAL_ORIGINAL_RANDOMIZED_RACE_PRETREATMENT_TWO_METABOLOMIC_PLATFORMS_PD_RESEARCH",
        "original_study":"RACE CaPS L-carnitine randomized dose vs placebo, 90d mortality",
        "n_original_independent_people":len(y),
        "n_original_mortality_by_day90":int(y.sum()),
        "n_placebo_original_participants":int(np.sum(c[:,4]==0)),
        "n_treated_original_participants":int(np.sum(c[:,4]==1)),
        "data_types":"pretreatment NMR 27 metabolites and LCMS 24 acylcarnitines, same metabolome",
        "not_transcriptome_proteome_multiomics":True,
        "n_omic_PCs_trainonly_each_platform":PCS_PER_ASSAY,
        "same_original_people_and_outer_split_all_models":True,
        "no_outcome_guided_pca_or_tuning_using_outer_test":True,
        "predictive_PD_outcome_not_individual_causal_treatment_benefit":True,
        "source_original_donor_assay_values_predictions_not_exported":True,
        "metrics":metric,
        "paired_brier_comparisons":comparisons,
        "original_source_cardinality":data.get("source_counts",{}),
    }


def run_original_public_caps():
    contents={}
    hashes={}
    for filename in SOURCE_FILES:
        contents[filename],hashes[filename]=_download_author_original(filename)
    data=load_original_caps(files=contents)
    results=evaluate_original_caps(data)
    results["pinned_original_source_Git_commit"]=SOURCE_COMMIT
    results["individual_original_source_file_SHA256"]=hashes
    return results
