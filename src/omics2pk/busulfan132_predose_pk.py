"""Subject-heldout time-transfer model for ORIGINAL busulfan plasma predose omics.

Eligibility is CONDITIONAL on authorized numeric SAME-person original n132 / n96
Metabolon molecular signal plus real measured Bayesian NONMEM busulfan CL
(10.1002/cpt.2794). This file is a methodological engine / privacy gate;
NO authentic 132×841 donor data has been obtained by this software.

Comparison:
  C0 age/sex/size clinical vs +omics (13 TRAIN-ONLY selected signals)
  C1 full author-validated GSTA1/GSTM1 PGx vs +omics IF available;
  NO claim of genotype superiority without real source PGx rows.

Train original 'immediate_predose' metabolome on people EXCLUDING original
future test humans. Test ONLY the 'two_week_predose' on the 96 originally
paired people, grouped into five folds by donor, ensuring test humans have
NO matching immediate-predose study rows within that model's training data.

This corrects the published same-person temporal transfer's information
overlap; no inference that model will generalize to an independent hospital.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime
from math import isfinite
from typing import Sequence,Mapping

SOURCE_DOI="10.1002/cpt.2794"
FEATURES=841
N_SUBJECTS=132
N_PAIRED=96
WINDOWS=frozenset({"immediate_predose","two_week_predose"})


class BusulfanSourceNotQualified(ValueError):
    pass


def _s(v,what):
    if not isinstance(v,str) or not v.strip() or len(v)>160:
        raise BusulfanSourceNotQualified("Missing or invalid "+what)
    return v.strip()


def _date(v,what):
    try: x=datetime.fromisoformat(_s(v,what).replace("Z","+00:00"))
    except ValueError:raise BusulfanSourceNotQualified("Original source timestamp invalid") from None
    if x.tzinfo is None or x.utcoffset() is None:
        raise BusulfanSourceNotQualified("Original source sampling timestamps need timezone")
    return x


def _finite(v,what,*,strict=False):
    try:
        n=float(v)
    except (ValueError,TypeError,OverflowError):
        raise BusulfanSourceNotQualified("Original "+what+" nonnumeric") from None
    if not isfinite(n) or (n<=0 if strict else n<0):
        raise BusulfanSourceNotQualified("Original "+what+" invalid/censored")
    return n


def qualify_authorized_original_busulfan(rows:Sequence[Mapping]):
    """Return PRIVATE validated source row objects; NEVER print this return.

    The caller may supply these rows only from legitimately obtained original
    clinical assay and Metabolon linked files. The source provenance fields
    are asserted but NOT cryptographically authenticated by this validator.
    """
    if len(rows)!=N_SUBJECTS+N_PAIRED:
        raise BusulfanSourceNotQualified("Expected exactly 228 original assays in 132 humans")
    seen_sample=set()
    by_human={}
    pgx_mode=None
    for row in rows:
        if row.get("source_doi")!=SOURCE_DOI:
            raise BusulfanSourceNotQualified("Source DOI mismatch")
        person=_s(row.get("subject_id"),"original human source id")
        assay=_s(row.get("sample_id"),"original biological assay id")
        if assay in seen_sample:
            raise BusulfanSourceNotQualified("Repeated original biological assay ID")
        seen_sample.add(assay)
        window=_s(row.get("time_window"),"documented dose-relative window")
        if window not in WINDOWS:
            raise BusulfanSourceNotQualified("Original sample not predose window")
        t=_date(row.get("collected_at"),"sample collection")
        dose=_date(row.get("first_busulfan_dose_at"),"first dose")
        if t>=dose or (dose-t).total_seconds()>45*86400:
            raise BusulfanSourceNotQualified("Omics not verified within prospective PRE-dose interval")
        cl=_finite(row.get("busulfan_CL_per_NFM"),"NONMEM estimated busulfan clearance",strict=True)
        if row.get("CL_unit")!="mL/min/kg_NFM":
            raise BusulfanSourceNotQualified("Incorrect source pharmacokinetic clearance units")
        features=row.get("metabolite_abundances")
        if not isinstance(features,(list,tuple)) or len(features)!=FEATURES:
            raise BusulfanSourceNotQualified("Original 841 metabolite panel unavailable")
        # Prohibit omic leakage from post-drug samples, BLQ and negative MS.
        # No source molecular values are emitted in the QC report.
        for v in features:
            if v is None or (isinstance(v,str) and v.strip().upper() in ("NA","")):
                continue
            _finite(v,"molecular feature",strict=False)
        age=_finite(row.get("age_years"),"age",strict=True)
        size=_finite(row.get("normal_fat_mass_kg"),"NFM",strict=True)
        sex=_s(row.get("sex"),"sex").lower()
        if sex not in ("male","female"):
            raise BusulfanSourceNotQualified("Unsupported sex coding")
        # No imputed fictitious GSTA1/GSTM1 genotype: all subjects need
        # authenticated labels and source assay provenance, or none.
        haspgx=any(row.get(k) is not None for k in (
            "GSTA1_original_function_class","GSTM1_original_copy_number"))
        if pgx_mode is None:pgx_mode=haspgx
        if pgx_mode!=haspgx:
            raise BusulfanSourceNotQualified("Incomplete PGx source coverage cannot be assumed")
        pgx=None
        if haspgx:
            if row.get("PGx_original_assay_verified") is not True:
                raise BusulfanSourceNotQualified("PGx original assay provenance unverified")
            pgx=(
                _finite(row.get("GSTA1_original_function_class"),"GSTA1 class"),
                _finite(row.get("GSTM1_original_copy_number"),"GSTM1 copies")
            )
            if int(pgx[0])!=pgx[0] or int(pgx[1])!=pgx[1] or pgx[0]>2 or pgx[1]>2:
                raise BusulfanSourceNotQualified("Original verified PGx category unavailable")
        entries=by_human.setdefault(person,{})
        if window in entries:
            raise BusulfanSourceNotQualified("Duplicated biological donor in one window")
        if entries:
            other=next(iter(entries.values()))
            if (abs(other["CL"]-cl)>1e-7 or other["age"]!=age
                or other["sex"]!=sex or other["nfm"]!=size or other["dose"]!=dose
                or other["pgx"]!=pgx):
                raise BusulfanSourceNotQualified("Inconsistent original repeated-person outcome or metadata")
        entries[window]={
            "features":features,"CL":cl,"age":age,"sex":sex,
            "nfm":size,"dose":dose,"sample_time":t,"pgx":pgx
        }
    if len(by_human)!=N_SUBJECTS:
        raise BusulfanSourceNotQualified("132 distinct original humans required")
    pair=0
    for windows in by_human.values():
        if "immediate_predose" not in windows:
            raise BusulfanSourceNotQualified("Each original person requires source-immediate predose measurement")
        if "two_week_predose" in windows:
            pair+=1
            if windows["two_week_predose"]["sample_time"]>=windows["immediate_predose"]["sample_time"]:
                raise BusulfanSourceNotQualified("Original 2-week predose assay after immediate predose")
    if pair!=N_PAIRED:
        raise BusulfanSourceNotQualified("Actual same-person two-week predose n96 expected")
    return by_human


def audit_original_busulfan_study(rows):
    source=qualify_authorized_original_busulfan(rows)
    return {
        "status":"PROVISIONAL_AUTHORIZED_SOURCE_SCHEMA_QC_NOT_CERTIFIED_PROVENANCE",
        "original_source_DOI":SOURCE_DOI,
        "n_original_unique_biological_humans":len(source),
        "n_original_pre_dose_biological_samples":sum(len(x) for x in source.values()),
        "n_original_two_week_samples":sum("two_week_predose" in x for x in source.values()),
        "n_reported_metabolites_per_sample":FEATURES,
        "all_samples_documented_before_first_busulfan_dose":True,
        "genuine_source_doc_authorized_and_independently_reviewed":False,
        "original_genotype_available":all(
            entries["immediate_predose"]["pgx"] is not None for entries in source.values()),
        "n_original_individual_PK_outcomes":len(source),
        "n_external_independent_clinical_humans":0,
        "no_original_subjects_omics_or_PK_rows_exported":True,
        "study_level_source_QC_not_clinical_PK_model":True,
    }


def benchmark_authorized_original_busulfan(rows,*,seed=20261018):
    """Return ONLY aggregate PK prediction errors. Requires actual source.

    Never call with simulated records and report as an original n132 result.
    """
    import numpy as np
    from sklearn.model_selection import KFold
    from sklearn.pipeline import make_pipeline
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    from sklearn.feature_selection import SelectKBest,f_regression
    from sklearn.linear_model import Ridge

    source=qualify_authorized_original_busulfan(rows)
    people=list(sorted(source))
    test_people=[p for p in people if "two_week_predose" in source[p]]
    fold=KFold(n_splits=5,shuffle=True,random_state=seed)
    pgx_ready=all(source[p]["immediate_predose"]["pgx"] is not None for p in people)
    errors={name:[] for name in (
        "clinical_baseline","clinical_plus_omics_13_trainonly",
        "clinical_plus_verified_PGx","clinical_plus_verified_PGx_plus_omics_13_trainonly"
    ) if pgx_ready or "PGx" not in name}
    actuals=[]
    predicted={name:[] for name in errors}
    # NO participant train features of test people from another timepoint.
    for train_idx,test_idx in fold.split(test_people):
        held={test_people[i] for i in test_idx}
        train=[p for p in people if p not in held]
        test=[test_people[i] for i in test_idx]
        assert not held.intersection(train)

        def build(ppl,window):
            cc=[];genes=[];yy=[]
            for p in ppl:
                src=source[p][window]
                cc.append([src["age"],1. if src["sex"]=="male" else 0.,src["nfm"]])
                arr=[float(v) if v not in (None,"","NA") else np.nan
                     for v in src["features"]]
                genes.append(np.log1p(arr))
                yy.append(src["CL"])
            return np.asarray(cc),np.asarray(genes),np.asarray(yy)

        ctrain,mtrain,ytrain=build(train,"immediate_predose")
        ctest,mtest,ytest=build(test,"two_week_predose")
        if (ytrain<=0).any() or (ytest<=0).any():
            raise BusulfanSourceNotQualified("Invalid NONMEM original clearance")
        ylog=np.log(ytrain)
        train_sets={"clinical_baseline":ctrain}
        test_sets={"clinical_baseline":ctest}
        features=make_pipeline(
            SimpleImputer(strategy="median",keep_empty_features=True),
            StandardScaler(),SelectKBest(f_regression,k=13))
        # Selection is on TRAIN human source values ONLY.
        mt=features.fit_transform(mtrain,ylog)
        mq=features.transform(mtest)
        train_sets["clinical_plus_omics_13_trainonly"]=np.column_stack((ctrain,mt))
        test_sets["clinical_plus_omics_13_trainonly"]=np.column_stack((ctest,mq))
        if pgx_ready:
            gp_train=np.asarray([source[p]["immediate_predose"]["pgx"] for p in train])
            gp_test=np.asarray([source[p]["two_week_predose"]["pgx"] for p in test])
            train_sets["clinical_plus_verified_PGx"]=np.column_stack((ctrain,gp_train))
            test_sets["clinical_plus_verified_PGx"]=np.column_stack((ctest,gp_test))
            train_sets["clinical_plus_verified_PGx_plus_omics_13_trainonly"]=np.column_stack((ctrain,gp_train,mt))
            test_sets["clinical_plus_verified_PGx_plus_omics_13_trainonly"]=np.column_stack((ctest,gp_test,mq))
        for key in errors:
            mdl=make_pipeline(SimpleImputer(strategy="median",keep_empty_features=True),
                              StandardScaler(),Ridge(alpha=1.0))
            mdl.fit(train_sets[key],ylog)
            pred=np.exp(np.clip(mdl.predict(test_sets[key]),-20.,20.))
            errors[key].extend(np.abs(pred-ytest))
            predicted[key].extend(pred)
        actuals.extend(ytest)
    truth=np.asarray(actuals)
    metric={}
    for key,p in predicted.items():
        pp=np.asarray(p)
        metric[key]={
            "original_busulfan_CL_MAE_mL_min_kg_NFM":round(float(np.mean(np.abs(pp-truth))),6),
            "mean_absolute_percentage_error_pct":round(float(100*np.mean(np.abs(pp/truth-1))),4),
            "n_true_heldout_people":len(truth),
        }
    return {
        "status":"RESEARCH_SOURCE_ONLY_BUSULFAN_132_PAIRED_PERSON_EXCLUDED_TIME_TRANSFER",
        "n_unique_training_source_people":N_SUBJECTS,
        "n_true_person_heldout_two_week_assays":len(truth),
        "n_original_metabolomic_compounds":FEATURES,
        "molecular_feature_selectK13_fold_local":True,
        "each_two_week_predose_heldout_person_excluded_from_all_training_timepoints":True,
        "verified_full_clinical_genotype_available":False,
        "source_original_specific_two_loci_PGx_present":pgx_ready,
        "source_independent_external_clinical_cohort":False,
        "genuine_source_documents_and_crosswalk_independently_reviewed":False,
        "original_subject_genotypes_concentrations_omicro_results_never_exported":True,
        "results":metric,
    }
