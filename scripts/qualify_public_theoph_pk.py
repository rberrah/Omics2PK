"""First original human dose-concentration-time profile source qualification.

Canonical experiment nlme::Theoph: 12 original human individuals x 11 observed
serum concentrations following one oral theophylline administration.
Source: Boeckmann/Sheiner/Beal (1994), R built-in Theoph.
No omics is present; this only verifies our PK longitudinal event layer.
"""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import io
import json
from urllib.request import Request, urlopen

from omics2pk.longitudinal_pk_contract import qualify_longitudinal_pk
from omics2pk.single_dose_nca import aggregate_original_theoph

RD_DATA_COMMIT = "cc03c29690889dbe83089f5206e2422db8c3f71f"
PUBLIC_SOURCE = (
    "https://raw.githubusercontent.com/vincentarelbundock/"
    f"Rdatasets/{RD_DATA_COMMIT}/csv/datasets/Theoph.csv"
)
CAP = 100_000
EXPECTED_ORIGINAL_SOURCE_SHA256 = "f9b77033090644cf38c83056c25c8137bca71cf629d7c0159a72ef18095726b0"


def load_original():
    with urlopen(Request(PUBLIC_SOURCE,headers={"User-Agent":"Omics2PK-original-theoph-PK-event-audit/1.0"}),timeout=30) as r:
        raw=r.read(CAP+1)
    if len(raw)>CAP:raise ValueError("Original single-dose PK source unexpectedly large")
    if hashlib.sha256(raw).hexdigest()!=EXPECTED_ORIGINAL_SOURCE_SHA256:
        raise ValueError("Pinned original Theoph numeric PK source SHA changed")
    records=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(records)!=132:
        raise ValueError("Theoph must have 132 actual original concentration samples")
    if set(records[0])!={"rownames","Subject","Wt","Dose","Time","conc"}:
        raise ValueError("Unexpected original published Theoph source format")
    return raw,records


def observed_pk_events(rows):
    group={}
    for row in rows:
        subject=row["Subject"].strip()
        if not subject:raise ValueError("Original subject missing")
        subject_rows=group.setdefault(subject,[])
        subject_rows.append(row)
    if len(group)!=12 or any(len(v)!=11 for v in group.values()):
        raise ValueError("Original Theoph must retain 12 biological subjects x11 observations")
    events=[]
    for subject, sr in group.items():
        weight=float(sr[0]["Wt"])
        dose_mg_kg=float(sr[0]["Dose"])
        if weight<=0 or dose_mg_kg<=0:
            raise ValueError("Original weight/dose must be strictly positive")
        if any(float(x["Wt"])!=weight or float(x["Dose"])!=dose_mg_kg for x in sr):
            raise ValueError("Original dose/weight changed within stated one-dose occasion")
        events.append({
            "subject_id": subject,
            "drug":"theophylline",
            "time_h": 0.0,
            "event_type":"dose",
            "data_origin":"observed",
            "dose_id":f"theoph-dose-{subject}",
            "amount_mg":dose_mg_kg*weight,
            "route":"oral",
            "formulation":"original_theoph_oral_unspecified",
        })
        for x in sr:
            t=float(x["Time"])
            dv=float(x["conc"])
            if not (0<=t<=26) or dv<0:
                raise ValueError("Invalid original time or observed concentration")
            events.append({
                "subject_id":subject,
                "drug":"theophylline",
                "time_h":t,
                "event_type":"observation",
                "data_origin":"observed",
                "sample_id":f"theoph-obs-{x['rownames']}",
                "analyte":"theophylline",
                "specimen":"serum",
                "sampling_phase":"predose" if t==0 else "postdose",
                "concentration_mg_l":dv,
            })
    return events


def run_public_original():
    raw,rows=load_original()
    qualified=qualify_longitudinal_pk(
        observed_pk_events(rows),
        drug="theophylline",
        time_origin="hours_since_original_single_oral_dose",
        dataset_kind="single_dose",
        dosing_history_status="complete",
    )
    if qualified["n_concentration_observations"]!=132 or qualified["n_independent_subjects"]!=12:
        raise ValueError("Theoph original subjects and concentration count mismatch")
    return {
        "source_data_sha256":hashlib.sha256(raw).hexdigest(),
        "source_commit":RD_DATA_COMMIT,
        "source":"R datasets::Theoph, Rdatasets CSV; originals Boeckmann, Sheiner, Beal (1994)",
        "n_original_pk_observations":len(rows),
        "only_pk_no_real_omics_available":True,
        "all_original_dose_and_samples_not_duplicated":True,
        "baseline_t0_observation_note":"Serum C(t=0) can be >0 in original published dataset; retained, not forced to zero.",
        "original_mg_per_kg_dose_explicitly_converted_using_measured_kg_weight":True,
        "pk_event_contract":qualified,
        "observed_single_dose_NCA":aggregate_original_theoph(rows),
    }


if __name__=="__main__":
    print(json.dumps(run_public_original(),sort_keys=True,indent=2))
