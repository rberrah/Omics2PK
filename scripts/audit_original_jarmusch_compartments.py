"""Source-qualified specimen/time provenance of original public Jarmusch GNPS.

Aggregate-only: Never output original individual source IDs, assays,
genotypes, molecular values, or patient-specific timing. This gate audits
the two *predose plasma* time points and other sampling compartments.
The original GNPS output can be re-zipped by the host; verify canonical
inner metadata/params/QZA source member SHA pins instead of ZIP bytes.
"""
from __future__ import annotations
from collections import defaultdict,Counter
import csv
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile
import json

from scripts.validate_local_jarmusch14_gnps_biom_pk import (
    REPACK_CANONICAL_MEMBER_SHA256, GNPS_TASK,
)

def audit_source_specimens(archive):
    src=Path(archive)
    if not 50_000_000 < src.stat().st_size < 150_000_000:
        raise ValueError("Not a source-qualified GNPS original archive")
    with ZipFile(src) as z:
        for member,expected in REPACK_CANONICAL_MEMBER_SHA256.items():
            if sha256(z.read(member)).hexdigest()!=expected:
                raise ValueError("Canonical GNPS internal source member digest changed")
        if GNPS_TASK not in z.read("params.xml").decode("utf-8"):
            raise ValueError("GNPS original task ID mismatch")
        with z.open("qiime2_output/qiime2_metadata.tsv") as handle:
            rows=list(csv.DictReader(
                (line.decode("utf-8-sig") for line in handle),delimiter="\t"))
    required=("sample_name","MS_SampleType","ATTRIBUTE_Study_Day",
              "ATTRIBUTE_Time_Point_Mins","ATTRIBUTE_Sample_Type",
              "ATTRIBUTE_Subject")
    if not rows or not set(required).issubset(rows[0]):
        raise ValueError("Original GNPS original sample provenance fields changed")
    observed={}
    for row in rows:
        label=row["sample_name"]
        if not label:
            continue
        if label in observed and any(row[k]!=observed[label][k] for k in required[1:]):
            raise ValueError("Original repeated analytical sample identifier disagrees")
        observed.setdefault(label,row)
    if len(rows)!=1440 or len(observed)!=1202:
        raise ValueError("Original GNPS sample duplicate counts changed")
    per_day={}
    for day in ("1","9"):
        samples=defaultdict(set)
        all_types=Counter()
        for row in observed.values():
            if row["MS_SampleType"]!="sample" or row["ATTRIBUTE_Study_Day"]!=day:
                continue
            matrix=row["ATTRIBUTE_Sample_Type"]
            all_types[matrix]+=1
            if row["ATTRIBUTE_Time_Point_Mins"]=="0" and row["ATTRIBUTE_Subject"]:
                samples[matrix].add(row["ATTRIBUTE_Subject"])
        plasma=[
            r for r in observed.values()
            if r["MS_SampleType"]=="sample"
            and r["ATTRIBUTE_Study_Day"]==day
            and r["ATTRIBUTE_Time_Point_Mins"]=="0"
            and r["ATTRIBUTE_Sample_Type"]=="blood_plasma"
            and r["ATTRIBUTE_Subject"]
        ]
        subjects=[r["ATTRIBUTE_Subject"] for r in plasma]
        if len(set(subjects))!=len(subjects):
            raise ValueError("Multiple subject/day/time0 original predose plasma assays")
        per_day[day]={
            "n_original_predose_plasma_biological_people":len(subjects),
            "original_sample_types_and_analytical_sample_counts":
                dict(sorted(all_types.items())),
            "predose_time0_biological_people_by_matrix":
                {k:len(v) for k,v in sorted(samples.items())},
        }
    subjects_day={}
    for day in ("1","9"):
        subjects_day[day]={
            r["ATTRIBUTE_Subject"] for r in observed.values()
            if r["MS_SampleType"]=="sample"
            and r["ATTRIBUTE_Study_Day"]==day
            and r["ATTRIBUTE_Time_Point_Mins"]=="0"
            and r["ATTRIBUTE_Sample_Type"]=="blood_plasma"
            and r["ATTRIBUTE_Subject"]
        }
    if len(subjects_day["1"])!=13:
        raise ValueError("Original verified day1 plasma source does not contain 13 unique people")
    return {
        "status":"ORIGINAL_GNPS_SOURCE_QUALIFIED_COMPARTMENT_AND_TIMEPOINT_AUDIT",
        "research_question":"pre-dose plasma molecular predictor provenance",
        "n_original_metadata_rows":len(rows),
        "n_unique_analytical_source_labels":len(observed),
        "predose_plasma_by_day":per_day,
        "n_original_people_with_day1_and_day9_predose_plasma":
            len(subjects_day["1"] & subjects_day["9"]),
        "sample_material_is_blood_plasma_not_liver_or_hepatic_CYP_protein":True,
        "genotype_biological_source_material_not_ascertained_from_GNPS_metadata":True,
        "day1_plasma_molecules_measured_before_antibiotics":True,
        "day9_plasma_molecules_measured_after_7day_cefprozil":True,
        "original_source_individual_sample_labels_never_exported":True,
    }

if __name__=="__main__":
    import sys
    if len(sys.argv)!=2:
        raise SystemExit("Pass qualified original GNPS ZIP path as sole argument")
    print(json.dumps(audit_source_specimens(sys.argv[1]),sort_keys=True,indent=2))
