"""Actual source-profile PBPK multi-analyte NCA from published APAP R reference.

Descriptive AUC between FIRST and LAST observed time only, separately for
acetaminophen, APAP-glucuronide and APAP-sulfate, where measured. The source
offers STUDY aggregates, not individual donor profiles; treat (Study,Dose)
as a literature series and NEVER infer individual clearance or uncertainty.
"""
from __future__ import annotations
from collections import defaultdict
from hashlib import sha256
import json
import math
import os
from statistics import median
from tempfile import TemporaryDirectory
from urllib.request import Request,urlopen

COMMIT="a775203eab90fa09059fde7190396d730ae04351"
URL=f"https://raw.githubusercontent.com/nanhung/pksensi/{COMMIT}/data/APAP.rda"
EXPECTED_SHA="08b1d33b96d41cdd12fd48b2f09ad28c68728bf047ef4c93772a97af3dca00e0"
TARGET=("APAP","AG","AS")


def compute_area(series):
    pairs=sorted((float(t),float(c)) for t,c in series if c is not None and math.isfinite(float(c)))
    if len(pairs)<3 or any(t<0 or c<0 for t,c in pairs):
        return None
    if any(b[0]<=a[0] for a,b in zip(pairs,pairs[1:])):
        # Deliberate ABSTAIN if one original literature series has replicate
        # points per identical time: no arbitrary averaging or pseudo-variance.
        return None
    auc=sum((b[0]-a[0])*(a[1]+b[1])/2.0 for a,b in zip(pairs,pairs[1:]))
    peak=max(pairs,key=lambda p:p[1])
    return {"n_source_time_points":len(pairs),"first_time_h":pairs[0][0],
            "last_time_h":pairs[-1][0],"cmax_mg_l":peak[1],"tmax_h":peak[0],
            "observed_AUC_first_to_last_mg_h_l":auc}


def run():
    import pyreadr
    with urlopen(Request(URL,headers={"User-Agent":"Omics2PK-PBPK-original-APAP-AG-AS-profile-validation/1.0"}),timeout=45) as resp:
        b=resp.read(250_001)
    if sha256(b).hexdigest()!=EXPECTED_SHA:
        raise ValueError("Original APAP source RDA changed: all models must abstain")
    with TemporaryDirectory() as td:
        filename=os.path.join(td,"APAP.rda")
        with open(filename,"wb") as f:f.write(b)
        frame=pyreadr.read_r(filename)["APAP"]
    if len(frame)!=32:raise ValueError("Expect 32 original literature time rows")
    # Preserve original-dose identity, not merely trial name.
    groups=defaultdict(list)
    for _,row in frame.iterrows():
        key=(str(row["Study"]),float(row["Dose"]),float(row["Wt"]))
        groups[key].append(row)
    out=defaultdict(list)
    passed_groups=0
    failed_groups=0
    for key,records in groups.items():
        valid={}
        for analyte in TARGET:
            seq=[(x["Time"],x[analyte]) for x in records if not math.isnan(float(x[analyte]))]
            valid[analyte]=compute_area(seq)
            if valid[analyte] is not None:
                out[analyte].append(valid[analyte])
        if all(valid.values()):passed_groups+=1
        else:failed_groups+=1
    return {
        "status":"ORIGINAL_PBPK_PARENT_AND_CONJUGATES_AGGREGATE_LITERATURE_NCA_NOT_POPPK",
        "source_sha256":EXPECTED_SHA,
        "n_original_time_rows":len(frame),
        "n_original_literature_study_labels":int(frame["Study"].nunique()),
        "n_distinct_original_study_dose_weight_profile_groups":len(groups),
        "n_groups_with_all_three_analytes_NCA":passed_groups,
        "n_groups_with_missing_or_ambiguous_molecular_peak":failed_groups,
        "n_eligible_curve_groups_by_analyte":{a:len(out[a]) for a in TARGET},
        "n_subjects_with_individual_PK_curve":0,
        "fixed_biological_unit":"original published STUDY profile, NOT subject",
        "nca_method":"within measured original first-to-last times, linear trapezoid, no AUC extrapolation",
        "no_independent_omics_features":True,
        "no_original_source_values_or_literature_study_row_ids_exported":True,
        "no_prospective_individual_crossvalidation":True,
        "PBPK_mechanistic_model_fitted":False,
    }


if __name__=="__main__":
    print(json.dumps(run(),indent=2,sort_keys=True))
