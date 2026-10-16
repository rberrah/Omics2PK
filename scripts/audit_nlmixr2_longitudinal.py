"""Source-inventory of nlmixr2data PK longitudinal RDA, no raw records exported.

Documents the crucial measured-vs-SIMULATED dosing and PK/PD separation:
- theo_sd genuine measured first-day PK
- theo_md includes simulated repeated-dose study sections, not clinical follow-up
- warfarin includes DIFFERENT DV compartments: PK concentrations AND PCA PD
- mavoglurant includes multiple PK OCC and multiple dosing / concentrations

The actual original .rda binary is fetched from one frozen source commit, parsed
in ephemeral GitHub Actions using pyreadr and only counts are printed.
"""
from __future__ import annotations

import hashlib
import json
import os
from tempfile import TemporaryDirectory
from urllib.request import Request,urlopen

import pyreadr

COMMIT="f2cfb01ade88d4d30e66f56b19efce2a2c6e26bd"
SOURCES={
    "theo_sd": ("data/theo_sd.rda",144),
    "theo_md": ("data/theo_md.rda",348),
    "warfarin": ("data/warfarin.rda",519),
    "mavoglurant": ("data/mavoglurant.rda",2678),
}


def source_report():
    result={}
    with TemporaryDirectory() as td:
        for name,(path,expected_n) in SOURCES.items():
            url=f"https://raw.githubusercontent.com/nlmixr2/nlmixr2data/{COMMIT}/{path}"
            with urlopen(Request(url,headers={"User-Agent":"Omics2PK-longitudinal-original-audit/1.0"}),timeout=30) as res:
                raw=res.read(3_000_000)
            if not raw or len(raw)>2_999_999:raise ValueError("Original PK file too large/empty")
            full=os.path.join(td,name+".rda")
            with open(full,"wb") as f:f.write(raw)
            items=pyreadr.read_r(full)
            if len(items)!=1:raise ValueError("Unexpected RDA object number")
            df=next(iter(items.values()))
            if len(df)!=expected_n:raise ValueError(f"Original source PK row count mismatch for {name}: observed={len(df)}, expected={expected_n}")
            columns=list(map(str,df.columns))
            idkey=next((c for c in columns if c.lower()=="id"),None)
            if idkey is None:raise ValueError("Original subject ID column missing")
            event_name=next((c for c in columns if c.lower()=="evid"),None)
            time_name=next((c for c in columns if c.lower()=="time"),None)
            if time_name is None:raise ValueError("Original PK time column missing")
            event_counts={}
            if event_name:
                event_counts={str(k):int(v) for k,v in df[event_name].value_counts(dropna=False).to_dict().items()}
            extra={}
            if name=="theo_md":
                extra["original_source_truly_repeated_dose_simulated_after_day1"]=True
                extra["cannot_count_simulated_events_as_human_clinical_repeats"]=True
            if name=="warfarin":
                dvid=next((x for x in columns if x.lower()=="dvid"),None)
                extra["DV_contains_both_PK_concentration_and_PCA_PD"]=bool(dvid)
                extra["dvid_category_count"]=len(df[dvid].dropna().unique()) if dvid else None
            if name=="mavoglurant":
                occ=next((x for x in columns if x.lower()=="occ"),None)
                extra["original_occasion_column_present"]=bool(occ)
                if occ:
                    extra["n_distinct_occasion_labels"]=int(df[occ].nunique(dropna=True))
            result[name]={
                "original_rda_sha256":hashlib.sha256(raw).hexdigest(),
                "n_source_rows":len(df),
                "n_independent_subjects":int(df[idkey].nunique()),
                "n_unique_time_points":int(df[time_name].nunique()),
                "time_min_original":float(df[time_name].min()),
                "time_max_original":float(df[time_name].max()),
                "columns":columns,
                "event_code_counts_original_not_converted":event_counts,
                "has_original_omics":False,
                "no_individual_source_data_exported":True,
                **extra,
            }
    return {"status":"PUBLIC_LONGITUDINAL_PK_RDA_FORMAT_AUDIT_ONLY_NOT_POPPK_FIT",
            "upstream_commit":COMMIT,"sources":result}


if __name__=="__main__":
    print(json.dumps(source_report(),sort_keys=True,indent=2))
