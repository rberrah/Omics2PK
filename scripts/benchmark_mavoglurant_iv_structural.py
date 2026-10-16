"""Original mavoglurant real-human *pooled* IV 1C/2C structural benchmark.

Genuine source RDA is pinned by SHA, processed ephemerally. Fit is fixed
effects (no ETA/IOV, no shrinkage, no EBEs and NOT PopPK). Hold out whole
biological people, including both clinical occasions, for honest prediction.
Omics unavailable in source. Never print individual doses/DVs/IDs/parameters.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import Request, urlopen

from omics2pk.iv_pk_structural import subject_disjoint_folds
from scripts.validate_mavoglurant_120_real_occasion_nca import (
    SOURCE, EXPECTED_SHA, OriginalPKCurveError, observed_curve, summarize
)


def qualified_original_profiles(rows):
    """Accept exactly the 196 QC-qualified single IV occasions, no refit imputation."""
    qc=summarize(rows)
    grouped=defaultdict(list)
    for row in rows:
        grouped[(row["ID"],row["OCC"])].append(row)
    profiles=[]
    exclusions=0
    for (patient,occ), events in grouped.items():
        try:
            check=observed_curve(events)
        except OriginalPKCurveError as exc:
            if str(exc) not in {
                "Replicate time within occasion needs assay-specific protocol, cannot blindly average",
                "Need exactly one actual IV dose and 3 real concentration samples per occasion",
            }:
                raise
            exclusions+=1
            continue
        dose=next(r for r in events if int(r["EVID"])==1)
        zero=float(dose["TIME"])
        points=tuple(sorted(
            (float(r["TIME"])-zero,float(r["DV"]))
            for r in events if int(r["EVID"])==0 and int(r["MDV"])==0 and
            r["DV"] is not None and float(r["DV"])>0
        ))
        if len(points)!=check["n_original_samples"] or any(t<=0 for t, c in points):
            raise OriginalPKCurveError("Measured time-zero or NCA/Pk input inconsistency; abstain")
        profiles.append({
            "subject":str(patient), "occasion":str(occ),
            "dose_mg":check["source_dose_mg"],
            "duration_h":check["source_infusion_duration_h"],
            "points":points,
        })
    if (len(profiles),exclusions)!=(196,2):
        raise OriginalPKCurveError("Expected exact 196 accepted IV occasions and 2 abstentions")
    if sum(len(p["points"]) for p in profiles)!=2403:
        raise OriginalPKCurveError("Original observed DV count changed: abstain")
    if len({p["subject"] for p in profiles})!=120:
        raise OriginalPKCurveError("Original biological subjects changed: abstain")
    if qc["n_original_occasion_groups_qualified_for_descriptive_NCA"]!=196:
        raise OriginalPKCurveError("Original source report disagreement")
    return profiles


def _arrays(profiles):
    import numpy as np
    counts=defaultdict(int)
    for p in profiles:
        counts[p["subject"]]+=len(p["points"])
    subject_weights=[];times=[];dose=[];duration=[];values=[];ids=[]
    for p in profiles:
        for t,y in p["points"]:
            times.append(t);values.append(y);dose.append(p["dose_mg"])
            duration.append(p["duration_h"])
            ids.append(p["subject"])
            subject_weights.append(1.0 / counts[p["subject"]]**0.5)
    return (
        np.asarray(times,float),np.asarray(values,float),
        np.asarray(dose,float),np.asarray(duration,float),
        np.asarray(subject_weights,float),ids
    )


def _predict(which, parameters, packed):
    import numpy as np
    t,_,dose,T,_,_=packed
    CL,V=parameters[:2]
    if which=="1C":
        terms=[(1.0,CL/V)]
    elif which=="2C":
        Q,Vp=parameters[2:]
        k10,k12,k21=CL/V,Q/V,Q/Vp
        trace=k10+k12+k21
        diff=np.sqrt(max(trace*trace-4*k10*k21,0.0))
        alpha,beta=(trace+diff)/2,(trace-diff)/2
        if alpha<=beta or beta<=0:
            return np.zeros_like(t)
        terms=[((alpha-k21)/(alpha-beta),alpha),
               ((k21-beta)/(alpha-beta),beta)]
    else:
        raise ValueError("Unknown structural comparator")
    response=np.zeros_like(t)
    pre=np.minimum(t,T)
    after=np.maximum(t-T,0.0)
    for weight,k in terms:
        response+=weight*(-np.expm1(-k*pre)/k)*np.exp(-k*after)
    return 1000.0*dose/T/V*response


def _fit(which, train):
    import numpy as np
    from scipy.optimize import least_squares
    pack=_arrays(train)
    logdv=np.log(pack[1])
    starts=([20,40], [8,20]) if which=="1C" else (
        [20,30,20,100], [8,20,5,100], [50,100,40,300]
    )
    lo=[0.1,1] + ([0.05,1] if which=="2C" else [])
    hi=[300,5000] + ([500,5000] if which=="2C" else [])
    def residual(logpars):
        estimate=np.maximum(_predict(which,np.exp(logpars),pack),1e-12)
        return pack[4]*(np.log(estimate)-logdv)
    best=None
    for start in starts:
        sol=least_squares(
            residual,np.log(start),bounds=(np.log(lo),np.log(hi)),
            max_nfev=250,method="trf",ftol=1e-8,xtol=1e-8,gtol=1e-8
        )
        score=float(np.dot(sol.fun,sol.fun))
        if not np.isfinite(score):
            raise ValueError("Nonfinite structural optimizer")
        if best is None or score<best[0]:
            best=(score,np.exp(sol.x))
    return best[1]


def _validate(which, parameters, heldout):
    import numpy as np
    pack=_arrays(heldout)
    pred=_predict(which,parameters,pack)
    if np.any(~np.isfinite(pred)) or np.any(pred<=0):
        raise ValueError("Nonpositive or nonfinite held-out structural prediction")
    diff=np.log(pred)-np.log(pack[1])
    by_human=defaultdict(list)
    for subject,r in zip(pack[5],diff):
        by_human[subject].append(float(r))
    # Equal person importance; two occasions cannot count as two independent people.
    loss={person:sum(r*r for r in records)/len(records)
          for person,records in by_human.items()}
    signed={person:sum(records)/len(records) for person,records in by_human.items()}
    return {
        "n_independent_heldout_humans":len(by_human),
        "n_heldout_occasions":len(heldout),
        "n_measured_heldout_concentrations":len(diff),
        "person_mse_log":sum(loss.values())/len(loss),
        "person_mean_bias_log":sum(signed.values())/len(signed),
    }


def heldout_comparison(profiles, folds=5):
    """Training-free validation on each test human: fixed predictions only."""
    import numpy as np
    import scipy
    people={p["subject"] for p in profiles}
    partitions=subject_disjoint_folds(people,n_folds=folds)
    scores={"1C":[],"2C":[]}
    for train_people,test_people in partitions:
        train=[p for p in profiles if p["subject"] in train_people]
        test=[p for p in profiles if p["subject"] in test_people]
        if {p["subject"] for p in train}!=train_people:
            raise ValueError("Missing training biological person")
        if {p["subject"] for p in test}!=test_people:
            raise ValueError("Missing testing biological person")
        for kind in scores:
            fit=_fit(kind,train)
            scores[kind].append(_validate(kind,fit,test))
    n_people=sum(x["n_independent_heldout_humans"] for x in scores["1C"])
    n_obs=sum(x["n_measured_heldout_concentrations"] for x in scores["1C"])
    if n_people != len(people) or n_obs != sum(len(p["points"]) for p in profiles):
        raise ValueError("Whole-human heldout coverage changed")
    report={
        "status":"MEASURED_PK_ONLY_POOLED_STRUCTURAL_BENCHMARK_NOT_POPPK_NOT_OMICS_VALIDATION",
        "original_biological_humans":len(people),
        "original_qualified_occasions":len(profiles),
        "original_measured_concentrations":n_obs,
        "folds":folds,
        "split":"SHA256-seeded whole human disjoint; both occasions in same fold",
        "fit_scope":"pooled fixed effects only; no patient-specific posthoc fit, EBEs, IIV or IOV",
        "model_units":"mg / h / L / ng mL-1; original RATE establishes each infusion duration",
        "scipy_version":scipy.__version__,
        "no_individual_identifiers_or_source_concentrations_exported":True,
        "no_omics_in_original_source":True,
        "no_nlme_or_popPK_claim":True,
    }
    for kind,parts in scores.items():
        # Held-out human weighted, not sample- or occasion-weighted.
        weighted=sum(x["person_mse_log"]*x["n_independent_heldout_humans"] for x in parts)/n_people
        bias=sum(x["person_mean_bias_log"]*x["n_independent_heldout_humans"] for x in parts)/n_people
        report[kind]={
            "root_mean_individual_mean_squared_log_error":round(float(np.sqrt(weighted)),6),
            "mean_heldout_person_log_bias":round(float(bias),6),
            "n_folds_with_measured_test_concentrations":len(parts),
            "n_biological_test_subjects_all_folds":n_people,
        }
    report["difference_2C_minus_1C_RMSLE"]=round(
        report["2C"]["root_mean_individual_mean_squared_log_error"] -
        report["1C"]["root_mean_individual_mean_squared_log_error"],6
    )
    return report


def run_original():
    import pyreadr
    with urlopen(Request(SOURCE,headers={"User-Agent":"Omics2PK-IV-pk-structural-research-only/1.0"}),timeout=45) as resp:
        payload=resp.read(3_000_001)
    if sha256(payload).hexdigest()!=EXPECTED_SHA:
        raise OriginalPKCurveError("Actual original RDA hash mismatch")
    with TemporaryDirectory() as directory:
        path=Path(directory)/"mavoglurant.rda"
        path.write_bytes(payload)
        rows=pyreadr.read_r(str(path))["mavoglurant"].to_dict("records")
    report=heldout_comparison(qualified_original_profiles(rows))
    report["source_sha256_verified"]=EXPECTED_SHA
    return report


if __name__=="__main__":
    print(json.dumps(run_original(),sort_keys=True,indent=2))
