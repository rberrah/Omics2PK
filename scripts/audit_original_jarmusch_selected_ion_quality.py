"""Original source GNPS Jarmusch selected MS1 ion: specimen detection and time QC.

Only aggregate detection counts and source metadata diagnostics are returned.
This is NOT a new supervised model, biological identity inference, or
independent validation. It uses the exact source-matched n10 day9 PK cohort,
and no individual donor/assay identifiers, intensities, PK or genotype values
are exported.

Uses the official pinned original GNPS archive source. The arbitrary
zero/nonzero presence reflects peak detection/processing (not concentration).
Any MS sample-type blanks/QC diagnostics are reported separately and do not
prove metabolite identity or instrument source absence.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from hashlib import sha256
from io import BytesIO,TextIOWrapper
import csv
import json
from zipfile import ZipFile

import numpy as np

from scripts.benchmark_original_jarmusch_urine_t0_nested import (
    _strict_original_urine_pair
)
from scripts.benchmark_original_jarmusch_nested_pk_germline import (
    _original_spectral_features
)
from scripts.validate_local_jarmusch14_gnps_biom_pk import (
    REPACK_CANONICAL_MEMBER_SHA256,GNPS_TASK
)

SELECTED={
    "J1_plasma_discovery":("6171","936"),
    "J9_urine_discovery":("369","3018","4321","6135"),
    "J1_urine_putative_bile":("740","1889","2496"),
}
PAIR_KEYS=(
    ("plasma_day1_t0","plasma_day9_t0"),
    ("urine_day1_t0_before_midazolam","urine_day9_t0_before_midazolam")
)

class OriginalIonQCError(ValueError):
    pass

def audited_quality(archive):
    original,predictors,genotype,pk=_strict_original_urine_pair(archive)
    if set(original)!={*["plasma_day1_t0","plasma_day9_t0"],
                      "urine_day1_t0_before_midazolam","urine_day9_t0_before_midazolam"}:
        raise OriginalIonQCError("Unknown original source biological sample matrices")
    if any(x.shape!=(10,6349) for x in original.values()):
        raise OriginalIonQCError("Original human sample count/feature shape differs")
    p=_original_spectral_features(archive)
    index={s:i for i,s in enumerate(p["source_public_MS1_scan_ids"])}
    ions=sorted({ion for batch in SELECTED.values() for ion in batch},key=int)
    if any(ion not in index for ion in ions):
        raise OriginalIonQCError("An original selected source MS1 ion absent in authenticated BIOM")

    matrix_names=("plasma_day1_t0","plasma_day9_t0",
                  "urine_day1_t0_before_midazolam",
                  "urine_day9_t0_before_midazolam")
    bpresence={name:original[name][:,[index[ion] for ion in ions]]>0
               for name in matrix_names}
    summary={}
    for i,ion in enumerate(ions):
        current={}
        for mat in matrix_names:
            current[mat]={
                "original_source_number_of_people":10,
                "n_measured_positive_MS1_peak":int(bpresence[mat][:,i].sum()),
                "n_source_below_software_peak_detection":int((~bpresence[mat][:,i]).sum()),
                "peak_presence_is_not_chemical_identification":True
            }
        pairs={}
        for a,b in PAIR_KEYS:
            aa=bpresence[a][:,i]
            bb=bpresence[b][:,i]
            pairs[a.replace("_day1_t0","")]={
                "n_detected_in_both_J1_J9":int((aa&bb).sum()),
                "n_detected_only_J1":int((aa&(~bb)).sum()),
                "n_detected_only_J9":int(((~aa)&bb).sum()),
                "n_detected_neither":int(((~aa)&(~bb)).sum()),
            }
        summary[ion]={
            "ion_source_identifier":ion,
            "tested_source_specimens":current,
            "same_person_day1_day9_peak_detection_concordance":pairs,
            "candidate_block":[key for key,series in SELECTED.items() if ion in series],
        }

    with ZipFile(archive) as z:
        for member,expected in REPACK_CANONICAL_MEMBER_SHA256.items():
            if sha256(z.read(member)).hexdigest()!=expected:
                raise OriginalIonQCError("Source ZIP contents changed")
        if GNPS_TASK not in z.read("params.xml").decode("utf-8"):
            raise OriginalIonQCError("GNPS task changed")
        with TextIOWrapper(z.open("qiime2_output/qiime2_metadata.tsv"),
                           encoding="utf-8-sig") as fh:
            rows=list(csv.DictReader(fh,delimiter="\t"))
        qza=z.read("qiime2_output/qiime2_table.qza")
    # Count nonclinical technical sample classes without leaking patient keys.
    unique={}
    for row in rows:
        sn=row.get("sample_name")
        if sn and sn not in unique:
            unique[sn]=row
    classification=Counter(x.get("MS_SampleType","") for x in unique.values())
    datecols=[c for c in rows[0] if "dateacquired" in c.casefold() or "batch" in c.casefold()]
    # No original patient/sample-specific dates or individual metadata exported.
    daydate={}
    for matrix in ("blood_plasma","urine"):
        daydate[matrix]={}
        for day in ("1","9"):
            group=[v for v in unique.values() if v.get("MS_SampleType")=="sample"
                   and v.get("ATTRIBUTE_Sample_Type")==matrix
                   and v.get("ATTRIBUTE_Time_Point_Mins")=="0"
                   and v.get("ATTRIBUTE_Study_Day")==day]
            daydate[matrix][day]={
                "n_source_predose_assays":len(group),
                "unique_source_acquisition_date_count":{
                    c:len(set(v.get(c,"") for v in group if v.get(c,""))) for c in datecols
                }
            }
    # Source BIOM n1198 assay rows include technical sample categories; no
    # assumptions that technical class equals actual subtraction blank.
    import h5py
    from scipy.sparse import csr_matrix
    with ZipFile(BytesIO(qza)) as q:
        name=[x for x in q.namelist() if x.endswith("/data/feature-table.biom")]
        if len(name)!=1:raise OriginalIonQCError("Original source BIOM uniqueness gate")
        biom=q.read(name[0])
    technical={}
    with h5py.File(BytesIO(biom),"r") as h:
        ids=[x.decode() for x in h["sample/ids"][:]]
        if len(ids)!=1198:raise OriginalIonQCError("Source assay count changed")
        mat=csr_matrix((h["sample/matrix/data"][:],
                        h["sample/matrix/indices"][:],
                        h["sample/matrix/indptr"][:]),shape=(1198,6349))
        source_cols=[index[x] for x in ions]
        cat=defaultdict(list)
        for idx,assay in enumerate(ids):
            rec=unique.get(assay)
            if rec is None:continue
            if rec.get("MS_SampleType")!="sample":
                cat[rec.get("MS_SampleType","unknown")].append(idx)
        for kind,idx in sorted(cat.items()):
            vals=mat[idx,:][:,source_cols].toarray()
            technical[kind]={
                "n_original_technical_assays":len(idx),
                "n_source_MS1_nonzero_by_selected_ion":{
                    ion:int((vals[:,i]>0).sum()) for i,ion in enumerate(ions)
                },
                "technicals_are_not_necessarily_valid_method_blank":True,
            }
    return {
        "status":"ACTUAL_JARMUSCH_N10_ORIGINAL_SELECTED_MS1_SAMPLING_AND_TECHNICAL_QC",
        "source_public_GNPS_task":GNPS_TASK,
        "n_persons_with_4_precise_source_predose_specimens_and_day9_CLF":10,
        "source_assay_classifications":dict(sorted(classification.items())),
        "technical_MS1_peak_detection_summary":technical,
        "source_batch_metadata_columns_available":datecols,
        "source_acquisition_date_uniqueness_per_matrix_day":daydate,
        "selected_original_source_ions":summary,
        "source_biom_peak_detections_are_not_absolute_concentrations":True,
        "cannot_infer_liver_cyp_activity_from_urine_or_plasma_MS1":True,
        "zero_abundance_does_not_prove_biological_absence":True,
        "no_per_person_original_molecular_values_or_identifiers_exported":True,
        "posthoc_reanalysis_is_not_external_clinical_validation":True,
    }

if __name__=="__main__":
    import sys
    if len(sys.argv)!=2:
        raise SystemExit("Usage: python -m scripts.audit_original_jarmusch_selected_ion_quality source.zip")
    print(json.dumps(audited_quality(sys.argv[1]),sort_keys=True,indent=2))
