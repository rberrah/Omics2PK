"""Descriptive NCA on ORIGINAL 120-person, 198-occasion human IV mavoglurant PK.

Source original nlmixr2data::mavoglurant pinned SHA. n=198 single-IV
infusion occasions; SAME biological subject can contribute two occasions.
This deliberately is NOT nonlinear mixed effects PopPK or omics validation.
AUC over actually observed first-to-last points only, not AUClast from 0
unless t=0 original observed; no intravenous CL/dose derived without AUCinf.
"""
from __future__ import annotations

from collections import defaultdict,Counter
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
from urllib.request import Request,urlopen

COMMIT="f2cfb01ade88d4d30e66f56b19efce2a2c6e26bd"
SOURCE=f"https://raw.githubusercontent.com/nlmixr2/nlmixr2data/{COMMIT}/data/mavoglurant.rda"
EXPECTED_SHA="ba1d75e20bcc09f20a55a05a9e4424244efa15b2dfcef6079cd27385d5876174"
EXPECTED_N_ROWS=2678
EXPECTED_N_PEOPLE=120
EXPECTED_N_OCCASIONS=198
EXPECTED_N_OBS=2480

class OriginalPKCurveError(ValueError):
    pass

def _float(v,name):
    if isinstance(v,bool):
        raise OriginalPKCurveError("Invalid boolean "+name)
    try:f=float(v)
    except (ValueError,TypeError,OverflowError):
        raise OriginalPKCurveError("Invalid original numerical "+name) from None
    if not isfinite(f):raise OriginalPKCurveError("Non-finite original numerical "+name)
    return f

def observed_curve(rows):
    """Already grouped single (person, OCC), native columns never reindexed by TIME."""
    doses=[r for r in rows if int(_float(r["EVID"],"EVID"))==1]
    raw_samples=[r for r in rows if int(_float(r["EVID"],"EVID"))==0]
    samples=[]
    n_mdv_excluded=0
    n_nonpositive_excluded=0
    for row in raw_samples:
        mdv=int(_float(row["MDV"],"sample MDV"))
        if mdv not in (0,1):
            raise OriginalPKCurveError("Unknown original MDV code")
        if mdv==1:
            n_mdv_excluded+=1
            continue
        dv=_float(row["DV"],"original DV")
        if dv<=0:
            n_nonpositive_excluded+=1
            continue
        samples.append(row)
    if len(doses)!=1 or len(samples)<3:
        raise OriginalPKCurveError("Need exactly one actual IV dose and 3 real concentration samples per occasion")
    dose=doses[0]
    if int(_float(dose["MDV"],"dose MDV"))!=1:
        raise OriginalPKCurveError("Original dose must have MDV 1")
    dose_amt=_float(dose["AMT"],"original AMT")
    rate=_float(dose["RATE"],"original RATE")
    origin=_float(dose["TIME"],"dose TIME")
    if dose_amt<=0 or rate<=0:raise OriginalPKCurveError("IV dose and RATE must be positive")
    if dose_amt/rate<0.04 or dose_amt/rate>0.5:
        raise OriginalPKCurveError("Original IV infusion duration outside expected 10-minute-range")
    points=[]
    for row in samples:
        t=_float(row["TIME"],"sample TIME")-origin
        dv=_float(row["DV"],"original DV")
        if t<0 or dv<0:raise OriginalPKCurveError("Negative source concentration or pre-dose sample")
        points.append((t,dv))
    points.sort()
    if any(b[0]<=a[0] for a,b in zip(points,points[1:])):
        raise OriginalPKCurveError("Replicate time within occasion needs assay-specific protocol, cannot blindly average")
    cmax=max(v for _,v in points)
    tmax=next(t for t,v in points if v==cmax)
    auc=sum((b[0]-a[0])*(a[1]+b[1])/2 for a,b in zip(points,points[1:]))
    return {
        "n_original_samples":len(points),
        "n_original_EVID0_MDV1_excluded":n_mdv_excluded,
        "n_original_EVID0_nonpositive_DV_excluded":n_nonpositive_excluded,
        "source_dose_mg":dose_amt,
        "source_infusion_duration_h":dose_amt/rate,
        "first_original_observation_time_h":points[0][0],
        "last_original_observation_time_h":points[-1][0],
        "Cmax_observed_original_DV":cmax,
        "Tmax_observed_h_from_dose":tmax,
        "AUC_first_to_last_original_DV_h":auc,
        "AUC_first_to_last_per_mg_dose_original_DV_h_per_mg":auc/dose_amt,
        "AUC_0_infinity_calculated":False,
    }

def summarize(rows):
    if len(rows)!=EXPECTED_N_ROWS:raise OriginalPKCurveError("Actual original file n2678 changed")
    grouped=defaultdict(list)
    subjects=set()
    n_obs=0
    for r in rows:
        person=r["ID"]
        occasion=r["OCC"]
        if person is None or occasion is None:
            raise OriginalPKCurveError("Original biological ID or occasion unavailable")
        subjects.add(person)
        grouped[(person,occasion)].append(r)
        if int(_float(r["EVID"],"EVID"))==0:n_obs+=1
    if len(subjects)!=EXPECTED_N_PEOPLE or len(grouped)!=EXPECTED_N_OCCASIONS or n_obs!=EXPECTED_N_OBS:
        raise OriginalPKCurveError("Original unique patients, occasions or measured concentrations have changed")
    original_people_by_occ=defaultdict(set)
    for person,occasion in grouped:
        original_people_by_occ[person].add(occasion)
    original_repeated=sum(len(z)==2 for z in original_people_by_occ.values())
    if original_repeated!=78 or any(len(z) not in (1,2) for z in original_people_by_occ.values()):
        raise OriginalPKCurveError("Original 78 paired subjects and 42 single-occasion subjects source changed")
    reports={}
    rejected=Counter()
    ACCEPTABLE_SOURCE_EXCLUSION={
        "Replicate time within occasion needs assay-specific protocol, cannot blindly average",
        "Need exactly one actual IV dose and 3 real concentration samples per occasion",
    }
    for key,values in grouped.items():
        try:
            reports[key]=observed_curve(values)
        except OriginalPKCurveError as e:
            if str(e) not in ACCEPTABLE_SOURCE_EXCLUSION:
                raise
            rejected[str(e)]+=1
    if not reports:
        raise OriginalPKCurveError("No actual individual PK curves eligible for NCA")
    people_with_occ=defaultdict(list)
    for (person,occasion),r in reports.items():
        people_with_occ[person].append(r)
    repeated_qualified=sum(len(z)==2 for z in people_with_occ.values())
    return {
        "status":"REAL_120_PERSON_198_SINGLE_INFUSION_OCCASION_OBSERVED_NCA_NOT_POPPK",
        "n_original_human_subjects":len(subjects),
        "n_original_independent_subject_occasion_groups":len(grouped),
        "n_original_occasion_groups_qualified_for_descriptive_NCA":len(reports),
        "n_original_occasion_groups_abstained_due_to_QC":sum(rejected.values()),
        "abstention_reasons_original_source_QC":dict(rejected),
        "n_original_EVID0_candidate_PK_records":n_obs,
        "n_original_MDV0_positive_PK_concentrations_used":sum(r["n_original_samples"] for r in reports.values()),
        "n_original_MDV1_excluded":sum(r["n_original_EVID0_MDV1_excluded"] for r in reports.values()),
        "n_original_MDV0_nonpositive_excluded":sum(r["n_original_EVID0_nonpositive_DV_excluded"] for r in reports.values()),
        "n_original_subjects_with_two_original_occasion_profiles":original_repeated,
        "n_original_subjects_with_one_original_occasion_profile":len(original_people_by_occ)-original_repeated,
        "n_original_subjects_with_two_NCA_qualified_occasion_profiles":repeated_qualified,
        "n_original_occ_curves_with_measured_time_zero_sample":sum(r["first_original_observation_time_h"]==0 for r in reports.values()),
        "n_original_profiles_with_three_or_more_eligible_observations":len(reports),
        "duplicate_clock_time_occasion_not_silently_averaged":True,
        "source_MDVs_and_nonpositive_DVs_do_not_count_as_true_quantified_concentrations":True,
        "median_original_samples_per_occasion":median(r["n_original_samples"] for r in reports.values()),
        "median_observed_Cmax_ng_ml":round(median(r["Cmax_observed_original_DV"] for r in reports.values()),6),
        "median_observed_Tmax_h":round(median(r["Tmax_observed_h_from_dose"] for r in reports.values()),6),
        "median_auc_first_to_last_ng_h_ml":round(median(r["AUC_first_to_last_original_DV_h"] for r in reports.values()),6),
        "median_first_source_time_h":round(median(r["first_original_observation_time_h"] for r in reports.values()),6),
        "median_last_source_time_h":round(median(r["last_original_observation_time_h"] for r in reports.values()),6),
        "dose_and_infusion_rates_from_original_EVID1_not_generated":True,
        "source_units":"Original mavoglurant source plasma DV ng/mL; time hours; AMT mg",
        "AUC_first_to_last_intervals_may_differ_and_not_directly_be_comparable":True,
        "never_confuse_198_occasions_with_198_independent_patients":True,
        "individual_PK_and_patient_identifiers_NOT_exported":True,
        "individual_omics_measurements":0,
        "actual_PopPK_SAEM_or_PBPK_fit":False,
    }

def run_original():
    import pyreadr
    with urlopen(Request(SOURCE,headers={"User-Agent":"Omics2PK-original-mavoglurant-human-NCA-198-occasions/1.0"}),timeout=40) as res:
        original=res.read(3_000_001)
    if sha256(original).hexdigest()!=EXPECTED_SHA:
        raise OriginalPKCurveError("Source RDA has changed; abstain")
    with TemporaryDirectory() as td:
        p=Path(td)/"mavoglurant.rda"
        p.write_bytes(original)
        df=pyreadr.read_r(str(p))["mavoglurant"]
    report=summarize(df.to_dict("records"))
    report["source_sha256"]=EXPECTED_SHA
    report["source_commit"]=COMMIT
    return report

if __name__=="__main__":print(json.dumps(run_original(),indent=2,sort_keys=True))
