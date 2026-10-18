"""Original Jarmusch healthy-human clinical *measured* CYP probe PK vs genotype.

Original author PUBLIC source:
https://github.com/alan-jarmusch/UntargetedMS-DrugMetabolism
pinned SHA 95586620ead4ac5241dc1e99a98e047c1f99b17b
Materials/Metadata.txt is 2.87MB of original study MS sample metadata, with
REAL matched fourteen-person pharmacogenotypes and targeted NCA PK outcome
metrics for days 1 and 9, replicated across multiple *analytical samples*.
Do NOT treat ~1423 MS rows as subjects/independent response observations.

This benchmark compares the utility of a **PRIOR OBSERVED DRUG-PROBE PK
PHENOTYPE** on day1 versus genotype for prediction of the day9 PK response
after seven days of cefprozil. Unlike the pending Jarmusch metabolomic MS
feature matrix, this is NOT yet an omic-abundance multiomic improvement test.
It is NOT a prediction of the same day's drug outcome from an after-drug assay.
Small n13/14 and internal LOPO only; no clinical deployment recommendation.
No original sample IDs, PK measurements, genotypes or original rows exported.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from io import StringIO
import csv
import json
from math import isfinite
from urllib.request import Request,urlopen

COMMIT="95586620ead4ac5241dc1e99a98e047c1f99b17b"
SOURCE="https://raw.githubusercontent.com/alan-jarmusch/UntargetedMS-DrugMetabolism/"+COMMIT+"/Materials/Metadata.txt"
DOI="10.1111/cts.12785"
ENDPOINTS={
    "CYP3A_midazolam_oral_CLF":{
        "source_column":"PharmKinetics_CYP3A4",
        "PGx_columns":("PharmGeno_CYP3A5","PharmGeno_CYP3A4"),
        "units":"Original NCA metadata scale: apparent oral CL/F (verify physical unit from assay paper)",
    },
    "CYP2C19_omeprazole_metabolite_AUC_ratio":{
        "source_column":"PharmKinetics_CYP2C19",
        "PGx_columns":("PharmGeno_CYP2C19",),
        "units":"original 5OH-OME/OME AUC source ratio, direction pending source codebook check",
    },
    "CYP1A2_caffeine_oral_CLF":{
        "source_column":"PharmKinetics_CYP1A2",
        "PGx_columns":(),
        "units":"Original NCA metadata scale: apparent oral CL/F (verify physical unit from assay paper)",
    }
}


class JarmuschOriginalPKError(ValueError):
    pass


def _original_text():
    try:
        with urlopen(Request(SOURCE,headers={
            "User-Agent":"Omics2PK-research-original-Jarmusch-pinned-metadata/1.0"
        }),timeout=45) as response:
            content=response.read(3_200_001)
    except OSError as exc:
        raise JarmuschOriginalPKError("Public original author file unavailable") from exc
    if len(content)>3_200_000:
        raise JarmuschOriginalPKError("Original research metadata exceeds expected bound")
    return content.decode("utf-8-sig"),sha256(content).hexdigest()


def _usable(x):
    return str(x or "").strip() not in ("","not applicable","NA","N/A","nan","None")


def _float(x):
    try: y=float(x)
    except (TypeError,ValueError):return None
    return y if isfinite(y) else None


def read_author_metadata(text:str):
    original=csv.DictReader(StringIO(text),delimiter="\t")
    if not original.fieldnames or len(original.fieldnames)<180:
        raise JarmuschOriginalPKError("Published source schema not present")
    needed={"ATTRIBUTE_Subject","ATTRIBUTE_Study_Day","ATTRIBUTE_Sample_Type",
            "MS_SampleType","Enrollment_Subject_Age","ATTRIBUTE_sex",
            "Day0_Vitals_weight_kg",
            *[v["source_column"] for v in ENDPOINTS.values()],
            *[g for v in ENDPOINTS.values() for g in v["PGx_columns"]]}
    if not needed.issubset(set(original.fieldnames)):
        raise JarmuschOriginalPKError("Original source genotype/PK fields absent")
    # DictReader duplicates headings for unused legacy columns; selected
    # original subject/day/genotype/PK headers themselves MUST be unique.
    for name in needed:
        if original.fieldnames.count(name)!=1:
            raise JarmuschOriginalPKError("Ambiguous duplicated original required assay field")
    # source per-subject per-day fields duplicated on many MS sample rows:
    # only blood_plasma sample records are eligible to form PK assay tables.
    byperson=defaultdict(lambda:defaultdict(list))
    all_subjects=set()
    ms_rows=0
    for row in original:
        ms_rows+=1
        person=(row.get("ATTRIBUTE_Subject") or "").strip()
        if not _usable(person):continue
        all_subjects.add(person)
        if (row.get("ATTRIBUTE_Study_Day") not in ("1","9") or
            row.get("ATTRIBUTE_Sample_Type")!="blood_plasma" or
            row.get("MS_SampleType")!="sample"):
            continue
        clinical=(
            _float(row.get("Enrollment_Subject_Age")),
            row.get("ATTRIBUTE_sex"),
            _float(row.get("Day0_Vitals_weight_kg"))
        )
        if not all(v is not None for v in (clinical[0],clinical[2])):
            continue
        if clinical[1] not in ("male","female"):
            continue
        pk={k:_float(row.get(v["source_column"])) for k,v in ENDPOINTS.items()}
        pgx={name:row.get(name) for name in (
            "PharmGeno_CYP2C19","PharmGeno_CYP3A5","PharmGeno_CYP3A4")}
        if any(_usable(v) and _float(v) is not None and _float(v)<=0 for v in (
            row.get("PharmKinetics_CYP3A4"),row.get("PharmKinetics_CYP1A2"),
            row.get("PharmKinetics_CYP2C19"))):
            raise JarmuschOriginalPKError("Non-positive published original PK observation")
        byperson[person][row["ATTRIBUTE_Study_Day"]].append({
            "clinical":clinical,"genotypes":pgx,"PK":pk,
        })
    if not 1350<=ms_rows<=1500:
        raise JarmuschOriginalPKError("Original author metadata source number of MS rows changed")
    if len(all_subjects)!=14:
        raise JarmuschOriginalPKError("Original source does not have exactly 14 real volunteers")

    # Require no contradictions among repeated original MS specimens.
    people={}
    for person,days in byperson.items():
        clean={}
        for day,rows in days.items():
            clinical={x["clinical"] for x in rows}
            if len(clinical)!=1:
                raise JarmuschOriginalPKError("Same original person's clinical covariates inconsistent")
            genotype={tuple(sorted(x["genotypes"].items())) for x in rows}
            if len(genotype)!=1:
                raise JarmuschOriginalPKError("Original drug-sample genotype records contradict each other")
            assays={}
            for end in ENDPOINTS:
                vals={x["PK"][end] for x in rows if x["PK"][end] is not None}
                if len(vals)>1:
                    raise JarmuschOriginalPKError("Sample-pseudoreplication contains contradictory target PK")
                assays[end]=next(iter(vals)) if vals else None
            clean[day]={
                "clinical":next(iter(clinical)),
                "genotypes":rows[0]["genotypes"],"PK":assays
            }
        if clean:people[person]=clean
    if len(people)!=14:
        raise JarmuschOriginalPKError("Not all 14 original humans have source plasma assay records")
    return {
        "people":people,
        "n_original_MS_metadata_rows":ms_rows,
        "n_independent_humans":len(people)
    }


def _source_arm_design(people,endpoint):
    import numpy as np
    description=ENDPOINTS[endpoint]
    clinical=[]
    genotypes=[]
    prev=[]
    response=[]
    for person in sorted(people):
        values=people[person]
        if "1" not in values or "9" not in values:continue
        first,last=values["1"],values["9"]
        y=last["PK"][endpoint]
        early=first["PK"][endpoint]
        if y is None or early is None:
            continue
        if first["clinical"]!=last["clinical"] or first["genotypes"]!=last["genotypes"]:
            raise JarmuschOriginalPKError("Timepoint subject clinical/PGx changed")
        codes=tuple(first["genotypes"].get(g) for g in description["PGx_columns"])
        if description["PGx_columns"] and not all(_usable(v) for v in codes):
            raise JarmuschOriginalPKError("Original targeted CYP germline genotype missing")
        clinical.append([
            first["clinical"][0],
            1 if first["clinical"][1]=="male" else 0,
            first["clinical"][2]
        ])
        genotypes.append(codes)
        prev.append(float(early))
        response.append(float(y))
    if len(response)<10:
        raise JarmuschOriginalPKError("Original source not enough authentic paired PK measurements")
    return (np.asarray(clinical,float),genotypes,np.asarray(prev,float),
            np.asarray(response,float))


def _test_holdout(clinical,genotypes,prior,response,*,include_genotypes):
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler,OneHotEncoder
    from sklearn.model_selection import LeaveOneOut
    from sklearn.impute import SimpleImputer

    variants=("clinical_reference","clinical_plus_CYP_genotype",
              "clinical_plus_PRIOR_drug_probe_PK",
              "clinical_plus_CYP_genotype_plus_PRIOR_drug_probe_PK")
    variants=[x for x in variants if include_genotypes or "genotype" not in x]
    results={m:np.empty(len(response)) for m in variants}
    for train,test in LeaveOneOut().split(clinical):
        scaler=StandardScaler().fit(clinical[train])
        c0=scaler.transform(clinical[train]);ct=scaler.transform(clinical[test])
        p=np.log(prior).reshape(-1,1)
        sd=StandardScaler().fit(p[train])
        p0=sd.transform(p[train]);pt=sd.transform(p[test])
        e0=et=np.empty((len(train),0))
        if include_genotypes:
            enc=OneHotEncoder(handle_unknown="ignore",sparse_output=False)
            g=np.asarray(genotypes,dtype=str)
            e0=enc.fit_transform(g[train])
            et=enc.transform(g[test])
        design={
            "clinical_reference":(c0,ct),
            "clinical_plus_CYP_genotype":(np.column_stack([c0,e0]),
                                         np.column_stack([ct,et])),
            "clinical_plus_PRIOR_drug_probe_PK":(
                np.column_stack([c0,p0]),np.column_stack([ct,pt])),
            "clinical_plus_CYP_genotype_plus_PRIOR_drug_probe_PK":(
                np.column_stack([c0,e0,p0]),np.column_stack([ct,et,pt]))
        }
        for key in variants:
            xx,xt=design[key]
            model=Ridge(alpha=1.0).fit(xx,np.log(response[train]))
            results[key][test]=np.exp(np.clip(model.predict(xt),-20,20))
    metrics={}
    for key,v in results.items():
        metrics[key]={
            "MAE_source_NCA_numeric_units":round(float(np.mean(np.abs(v-response))),5),
            "absolute_percentage_error_pct":round(float(100*np.mean(np.abs(v/response-1))),4)
        }
    return metrics


def benchmark_original_author_Jarmusch14(text):
    base=read_author_metadata(text)
    output={}
    for endpoint,meta in ENDPOINTS.items():
        clinical,genotypes,early,later=_source_arm_design(base["people"],endpoint)
        m=_test_holdout(
            clinical,genotypes,early,later,
            include_genotypes=bool(meta["PGx_columns"]))
        output[endpoint]={
            "n_independent_donors_with_day1_and_day9_actual_NCA_PK":len(later),
            "PK_parameter":meta["units"],
            "source_actual_biological_day_1_prior_PK_value_available":True,
            "assay_day9_after_cefprozil_used_as_target":True,
            "genotype_reference_available":bool(meta["PGx_columns"]),
            "held_out_original_subject_group":"LOPO individual 14 volunteers",
            "same_human_all_repeated_MS_records_grouped_before_CV":True,
            "MAE_comparison":m,
            "original_subject_or_PK_values_exported":False
        }
    return {
        "status":"ORIGINAL_JARMUSCH_CLINICAL_PK_AND_GERMLINE_PGX_LOPO_ACTUAL_NCA",
        "study_doi":DOI,"study_metadata_git_commit":COMMIT,
        "n_original_analytical_MS_source_rows":base["n_original_MS_metadata_rows"],
        "n_unique_original_volunteers":base["n_independent_humans"],
        "endpoint_analysis":output,
        "important_scientific_limitations":[
            "n=13-14; descriptive leave-person-out internal prediction, not external statistical proof",
            "prior day1 CYP probe PK is a measured functional phenotype, NOT omics",
            "no per-person molecular MS intensity matrices joined here",
            "day9 PK is after cefprozil, prior day1 PK is available before antibiotic",
            "CYP2C19 output is original metabolite to parent AUC ratio, not clearance",
            "PGx is small target panel; not comprehensive phased metaboliser alleles",
        ],
        "source_study_subject_ids_PGx_calls_PK_values_never_exported":True,
    }


def run_public_jarmusch_clinical_benchmark():
    text,digest=_original_text()
    result=benchmark_original_author_Jarmusch14(text)
    result["source_original_public_metadata_SHA256"]=digest
    return result


if __name__=="__main__":
    print(json.dumps(run_public_jarmusch_clinical_benchmark(),indent=2,sort_keys=True))
