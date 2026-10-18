"""True ORIGINAL Rotroff 2016 n33 human clinical metformin E3 benchmark.

Reconstruct *only* genuine publisher-distributed 4-digit matched n33 subjects:
molecular source sample A predose (never B/C) and authentic individual
metformin AUClast0–24 or CLr rows. 32 participants have each quantitative
endpoint. Published person identifiers, molecular values, predictions and
concentration records never printed, saved or committed.

Two fixed 32-donor leave-one-person-out regression estimators:
  clinical: age, sex, weight
  experimental: age, sex, weight + 3 TRAIN-ONLY PCA components from up to
                384 untargeted predose original GC-TOF molecular features.
No response-guided feature selection, no hyperparameter search. Measured
original source units not inferred if publisher omits physical AUC units.

This is E3 INTERNAL exploratory CV on a genuine clinical PK cohort, NOT
independent external cohort prediction, clinical dosing support, SWEET or
an omics causal effect.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
from io import BytesIO
from math import isfinite
from statistics import mean
from random import Random
from zipfile import ZipFile
import json

from scripts.audit_rotroff33_metformin_original_pk_omics import (
    LIMIT, MEMBER_LIMIT, _as_float
)
from scripts.audit_rotroff33_crosssheet_pk_omics import decode_workbook_sheets
from scripts.qualify_rotroff33_original_metformin_patient_keys import (
    _sample_index, _clinical_ids,
)

DOI="10.3389/fphar.2016.00135"
SOURCE_CLINICAL_SHA="76ef6d62e1809f37b7fcf0d9a07e38eb55a67331b23bd058890911547fa9681a"
SOURCE_METAB_SHA="1debb35eac7c22383031d4aaeaa47c0362d57609956f16fa4420626efd4a8fe6"
SEED=20261018
COMPONENTS=3
ALPHA=1.0

class E3SourceError(ValueError):
    pass


def actual_original_tables(zipdata):
    """Source SHA binds the exact scientific workbook bytes despite volatile ZIP."""
    if len(zipdata)>LIMIT or not zipdata.startswith(b"PK"):
        raise E3SourceError("Original published source was not downloaded intact")
    found={}
    with ZipFile(BytesIO(zipdata)) as z:
        for member in z.infolist():
            if member.is_dir() or not member.filename.lower().endswith(".xlsx"):
                continue
            if member.file_size>MEMBER_LIMIT or ".." in member.filename.split("/"):
                raise E3SourceError("Original article workbook too large or unsafe")
            blob=z.read(member)
            h=sha256(blob).hexdigest()
            if h in (SOURCE_CLINICAL_SHA,SOURCE_METAB_SHA):
                if h in found:raise E3SourceError("Duplicate original source workbook")
                found[h]=decode_workbook_sheets(blob)
    if set(found)!={SOURCE_CLINICAL_SHA,SOURCE_METAB_SHA}:
        raise E3SourceError("Original source clinical or metabolomic workbook checksum changed")
    return found[SOURCE_CLINICAL_SHA], found[SOURCE_METAB_SHA]


def measured_original_data(zipdata):
    """Returns original patient dataset INTERNALLY, never serialized by script."""
    clinical,omics=actual_original_tables(zipdata)
    if len(clinical)!=6:
        raise E3SourceError("Expected original six clinical/stats worksheets")
    pk=clinical[2]  # original publication clinical Table: AUClast and CLr
    if _clinical_ids(pk) is None:
        raise E3SourceError("Not exactly 33 unique original clinical 4-digit person IDs")
    expected_headers={
        0:"Subject_ID",1:"Age",2:"Gender",3:"Weight",
        4:"AUClast (0-24)",10:"CLr_renal clearance",
    }
    for col,label in expected_headers.items():
        if pk[0].get(col,"").strip()!=label:
            raise E3SourceError("Clinical PK source header identity changed")
    if len(omics)!=2:
        raise E3SourceError("Expected original two metabolomic source worksheets")
    mol=omics[0]
    groups=_sample_index(mol)
    if not groups or groups["ambiguous"] or groups["n_complete_patient_triplets"]!=33:
        raise E3SourceError("Original 33 paired A/B/C mass spectrometry sample IDs not verified")
    # The source codebook itself is four-letter + four-digit + A/B/C.
    # Reconstruct original SUBJECT ID from this embedded original code,
    # never from row position or sorted order.
    from scripts.qualify_rotroff33_original_metformin_patient_keys import SAMPLE
    baseline={}
    matched_n=Counter()
    for sample in mol[1:]:
        name=sample.get(0,"").strip()
        m=SAMPLE.fullmatch(name)
        if m:
            _,patient_id,time=m.groups()
            matched_n[patient_id]+=1
            if time=="A":
                if patient_id in baseline:
                    raise E3SourceError("Duplicate original predose A sample")
                baseline[patient_id]=sample
    pkrows={}
    for row in pk[1:]:
        key=row.get(0,"").strip()
        if key in pkrows:
            raise E3SourceError("Duplicate original PK patient ID")
        pkrows[key]=row
    if not(set(pkrows)==set(baseline)) or len(baseline)!=33:
        raise E3SourceError("Original 33-person clinical and predose omic keys mismatch")
    if any(matched_n[s]!=3 for s in baseline):
        raise E3SourceError("Missing or extra original within-patient A/B/C samples")
    # Require >=300 empirically numeric GC-TOF features among 33 A samples;
    # exclude original sample-ID column. NaN imputation is ALWAYS fold local.
    feature_cols=sorted(set().union(*(set(r) for r in baseline.values()))-{0})
    source_quantified=[
        j for j in feature_cols if sum(_as_float(r.get(j)) is not None for r in baseline.values())>=25
    ]
    if not 300<=len(source_quantified)<=500:
        raise E3SourceError("Original predose quantitative molecular feature count is not supported")
    records=[]
    for key in sorted(pkrows):
        cl=pkrows[key]
        age=_as_float(cl.get(1));weight=_as_float(cl.get(3))
        sex=cl.get(2,"").strip().casefold()
        sexval={"m":1.0,"male":1.0,"f":0.0,"female":0.0}.get(sex)
        if age is None or not (18<=age<=100) or weight is None or not (25<=weight<=250) or sexval is None:
            raise E3SourceError("Ineligible or unrecognized clinical original age, sex, weight")
        features=[_as_float(baseline[key].get(j)) for j in source_quantified]
        targets={
            "metformin_AUClast_0_24_original_units":_as_float(cl.get(4)),
            "metformin_CLr_renal_clearance_original_units":_as_float(cl.get(10)),
        }
        for val in targets.values():
            if val is not None and val<=0:
                raise E3SourceError("Nonpositive original measured PK value")
        records.append({
            "internal_patient_id":key,
            "clinical":[float(age),float(sexval),float(weight)],
            "predose_omics":[float(v) if v is not None else None for v in features],
            "targets":targets,
        })
    complete={key:sum(r["targets"][key] is not None for r in records)
              for key in records[0]["targets"]}
    if any(n != 32 for n in complete.values()):
        raise E3SourceError("Original original individual drug AUC/renal clearance count changed")
    return records, {
        "n_original_unique_source_patients":len(records),
        "n_original_triplet_source_samples":sum(matched_n.values()),
        "original_predose_sample_category":"A, collected BEFORE first drug dose",
        "n_source_predose_molecular_features":len(source_quantified),
        "n_subjects_complete_by_original_measured_PK_target":complete,
        "original_clinical_workbook_SHA256":SOURCE_CLINICAL_SHA,
        "original_untargeted_metabolomic_workbook_SHA256":SOURCE_METAB_SHA,
        "original_subject_rows_or_numeric_measurements_exported":False,
    }


def _report_one_target(records,key):
    try:
        import numpy as np
        from sklearn.decomposition import PCA
        from sklearn.impute import SimpleImputer
        from sklearn.linear_model import Ridge
        from sklearn.preprocessing import StandardScaler
    except ImportError as e:
        raise E3SourceError("Original E3 optional sklearn/numpy dependencies required") from e
    selected=[x for x in records if x["targets"][key] is not None]
    if len(selected)!=32:
        raise E3SourceError("Expected exactly 32 independently measured original PK outcome subjects")
    y=np.array([x["targets"][key] for x in selected],dtype=float)
    if not np.all(np.isfinite(y)) or np.min(y)<=0:
        raise E3SourceError("Original quantitatively measured patient PK target invalid")
    log_y=np.log1p(y)
    clinical=np.asarray([x["clinical"] for x in selected],dtype=float)
    molecular=np.asarray([[np.nan if v is None else v for v in x["predose_omics"]]
                          for x in selected],dtype=float)
    observed=np.isfinite(molecular)
    if min(np.sum(observed,axis=1))<250:
        raise E3SourceError("Some original predose human metabolome has insufficient quantified peaks")

    errs={name:[] for name in ("clinical_age_sex_weight","clinical_plus_3_predose_omics_PCs",
                              "predose_omics_3_PCs_only")}
    preds={name:[] for name in errs}
    for h in range(len(selected)):
        train=np.arange(len(selected))!=h
        train_clin=StandardScaler().fit(clinical[train])
        xx=train_clin.transform(clinical[train])
        testc=train_clin.transform(clinical[h:h+1])
        # CRITICAL: each donor's heldout molecular profile has NO say in
        # training-stage imputation, scaling, PCA directions, or ridge fit.
        imp=SimpleImputer(strategy="median")
        trainm=imp.fit_transform(molecular[train])
        testm=imp.transform(molecular[h:h+1])
        scale=StandardScaler().fit(trainm)
        zm=scale.transform(trainm)
        ztest=scale.transform(testm)
        pca=PCA(n_components=COMPONENTS,svd_solver="randomized",random_state=SEED)
        pcs=pca.fit_transform(zm)
        pctest=pca.transform(ztest)
        plans={
            "clinical_age_sex_weight":(xx,testc),
            "clinical_plus_3_predose_omics_PCs":(
                np.column_stack((xx,pcs)),np.column_stack((testc,pctest))),
            "predose_omics_3_PCs_only":(pcs,pctest),
        }
        for name,(xtr,xt) in plans.items():
            model=Ridge(alpha=ALPHA,fit_intercept=True)
            model.fit(xtr,log_y[train])
            p=float(np.expm1(min(float(model.predict(xt)[0]),30.0)))
            p=max(0.,p)
            if not isfinite(p):
                raise E3SourceError("Nonfinite patient heldout actual-source prediction")
            preds[name].append(p)
            errs[name].append(abs(float(y[h])-p))
    metrics={}
    for name in errs:
        diffs=np.asarray(preds[name])-y
        metrics[name]={
            "MAE_original_source_PK_table_units":round(float(np.mean(np.abs(diffs))),6),
            "RMSE_original_source_PK_table_units":round(float(np.sqrt(np.mean(diffs*diffs))),6),
            "mean_bias_original_source_PK_table_units":round(float(np.mean(diffs)),6),
            "percent_within_20pct_of_actual_patient_PK":round(float(np.mean(np.abs(diffs)/y<=.20)*100),3),
        }
    differences=[a-b for a,b in zip(
        errs["clinical_plus_3_predose_omics_PCs"],errs["clinical_age_sex_weight"]
    )]
    bootstrap=Random(SEED)
    samples=sorted(mean(differences[bootstrap.randrange(len(differences))] for _ in differences)
                   for _ in range(1999))
    return {
        "n_original_PK_measured_subjects":len(selected),
        "n_original_subject_heldout_folds":len(selected),
        "n_original_source_omics_features_considered":molecular.shape[1],
        "models_fixed_before_experiment":"ridge alpha1, trainonly imputer/scaler/PCA3, log1p response",
        "clinical_baseline_features":"age, sex, weight",
        "original_omic_input":"one predose GC-TOF human sample A for each source person",
        "performance":metrics,
        "delta_MAE_plus_omics_minus_clinical_baseline":round(mean(differences),6),
        "exploratory_paired_95pct_bootstrap_delta_MAE":[
            round(samples[49],6),round(samples[1949],6)],
        "n_clinical_subjects_omics_model_has_lower_absolute_error":sum(x<0 for x in differences),
        "n_clinical_subjects_omics_model_has_higher_absolute_error":sum(x>0 for x in differences),
        "no_subject_IDs_numeric_or_predictions_exported":True,
        "not_external_validation":True,
        "not_clinical_dosing_advice":True,
    }


def evaluate_e3_original(records,metadata):
    endpoint_results={}
    for key in sorted(records[0]["targets"]):
        endpoint_results[key]=_report_one_target(records,key)
    return {
        "status":"AUTHENTIC_ORIGINAL_HUMAN_METFORMIN33_OMICS_PREDOSE_TO_MEASURED_PK_E3_INTERNAL_EXPLORATORY",
        "source_original_study_DOI":DOI,
        "original_source_schema_provenance_aggregates":metadata,
        "n_original_related_measured_pk_endpoints":len(endpoint_results),
        "endpoint_results":endpoint_results,
        "original_pk_physical_unit_source_column_metadata_missing":True,
        "limitations":[
            "Original 32 complete measured PK outcomes; one original human lacks PK observation",
            "Internal leave-one-person-out exploratory CV on same study, no independent external cohort",
            "384 untargeted GC-TOF predose features, 3 unsupervised TRAIN-ONLY PC components; no SWEET or XGB fit",
            "Baseline age+sex+weight, no actual source PGx genotype; no claim beyond clinical+PGx",
            "No supervised feature search or response-based PCA/tuning; high dimensional/small cohort uncertainty",
            "Source original physical PK outcome unit not inferred where header omits units",
            "Postdrug sample B/C and original glucose AUC excluded from biomarker predictors and original drug endpoint",
            "Study article describes A as predose; source patient IDs match original metadata exactly",
            "This is observational association, no causal proof of metabolite regulation",
        ],
        "original_subject_IDs_and_medical_measurements_exported":False,
    }
