"""Original public APAP, APAP-glucuronide, APAP-sulfate PBPK profile inventory.

The pksensi R package (Hsieh et al 2020/2022, source Zurlinden
and Reisfeld 2016) provides 32 Time×dose×study concentration rows
for three actual circulating drug-related analytes. Individual participant
IDs are NOT supplied: values are study-level reported literature profiles,
not 32 independent people nor personalized multimodal longitudinal omics.
This audits data provenance, counts and units without redistributing source.
"""
from __future__ import annotations
import hashlib
import json
import os
from tempfile import TemporaryDirectory
from urllib.request import Request,urlopen

COMMIT="a775203eab90fa09059fde7190396d730ae04351"
REPO="nanhung/pksensi"
SOURCE=f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/data/APAP.rda"
DOI="10.1007/s13318-015-0290-x"
MAX=250_000


def audit():
    import pyreadr
    with urlopen(Request(SOURCE,headers={"User-Agent":"Omics2PK-original-public-APAP-PBPK-three-analyte-source-QA/1.0"}),timeout=40) as r:
        data=r.read(MAX+1)
    if len(data)>MAX:raise ValueError("PBPK reference RDA size unexpectedly high")
    with TemporaryDirectory() as td:
        dest=os.path.join(td,"APAP.rda")
        with open(dest,"wb") as f:f.write(data)
        objects=pyreadr.read_r(dest)
    if len(objects)!=1 or "APAP" not in objects:raise ValueError("Original APAP RDA object/source changed")
    frame=objects["APAP"]
    fields={"Wt","Dose","Time","APAP","AG","AS","Study"}
    if set(frame.columns)!=fields:
        raise ValueError("Unexpected original APAP data frame columns")
    if len(frame)!=32:
        raise ValueError("Expected original 32 dose/time multi-analyte source rows")
    if frame["Study"].isna().any() or frame["Time"].isna().any():
        raise ValueError("Original PBPK time and study labels missing")
    if (frame[["Time","APAP","AG","AS","Dose","Wt"]].min(axis=0)<0).any():
        raise ValueError("Original negative drug/metabolite/time/dose")
    n_3=int(frame[["APAP","AG","AS"]].notna().all(axis=1).sum())
    groups=frame.groupby("Study",dropna=False)
    return {
        "status":"ORIGINAL_PUBLIC_HUMAN_APAP_METABOLITE_PBPK_STUDY_LEVEL_NOT_INDIVIDUAL_PK",
        "source_repo":REPO,
        "source_commit":COMMIT,
        "source_file_sha256":hashlib.sha256(data).hexdigest(),
        "original_literature_source":"Zurlinden & Reisfeld (2016), Bayesian population PBPK human acetaminophen and conjugates",
        "original_pksensi_reference":"R pksensi 1.2.3 APAP dataset",
        "n_original_source_time_rows":len(frame),
        "n_distinct_literature_study_labels":int(frame["Study"].nunique()),
        "n_original_rows_with_three_numeric_analytes":n_3,
        "analytes":["acetaminophen","acetaminophen_glucuronide","acetaminophen_sulfate"],
        "original_concentration_unit":"mg/L",
        "original_time_unit":"hours",
        "original_Weight_definition":"study-averaged weight kg",
        "original_Dose_definition":"reported oral dose mg/kg",
        "explicit_individual_patient_ID_field":False,
        "eligible_for_individual_donor_disjoint_PopPK":False,
        "eligible_for_population_mean_PBPK_curve_comparison":True,
        "human_omics_features_in_source":False,
        "not_a_same_patient_omics_to_PK_dataset":True,
        "original_rows_or_measurements_exported":False,
    }

if __name__=="__main__":
    print(json.dumps(audit(),sort_keys=True,indent=2))
