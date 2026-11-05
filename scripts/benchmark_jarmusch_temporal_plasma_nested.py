"""Real Jarmusch same-subject temporal comparison: predose plasma J1 versus J9.

Source-confirmed original GNPS task, 13 J1 predose plasma humans, 12 J9,
11 original people with BOTH specimens. Original PK endpoint for
both predictor clocks: measured **day9** oral midazolam CL/F. This fixes
the earlier design's mismatch when interpreting J1 plasma as if measured
after the J2-J8 cefprozil perturbation.

No donor IDs, clinical source rows, genotype values, or original ion
intensities are serialized or printed; only aggregate 11-person CV.
"""
from __future__ import annotations

from collections import defaultdict
import csv
import io
import json
from zipfile import ZipFile

import numpy as np
from scripts.validate_local_jarmusch14_gnps_biom_pk import (
    load_original_source,GNPSOriginalSourceError,
)
from scripts.benchmark_original_jarmusch_nested_pk_germline import (
    _original_spectral_features,_nested_cv,
)
from scripts.audit_original_jarmusch_compartments import audit_source_specimens


class JarmuschTemporalSourceError(ValueError):
    pass


def _temporal_source_join(archive):
    import h5py
    from scipy.sparse import csr_matrix
    x1,c1,g1,pk,qc=load_original_source(archive)
    evidence=audit_source_specimens(archive)
    if (evidence["predose_plasma_by_day"]["1"]["n_original_predose_plasma_biological_people"]!=13
        or evidence["predose_plasma_by_day"]["9"]["n_original_predose_plasma_biological_people"]!=12
        or evidence["n_original_people_with_day1_and_day9_predose_plasma"]!=11):
        raise JarmuschTemporalSourceError("Source predose plasma donor crosswalk changed")
    with ZipFile(archive) as z:
        with io.TextIOWrapper(z.open("qiime2_output/qiime2_metadata.tsv"),
                              encoding="utf-8-sig") as fh:
            meta=list(csv.DictReader(fh,delimiter="\t"))
        qza=z.read("qiime2_output/qiime2_table.qza")
    analytic={}
    for row in meta:
        assay=row["sample_name"]
        if assay in analytic:
            for field in ("ATTRIBUTE_Subject","ATTRIBUTE_Study_Day",
                          "ATTRIBUTE_Time_Point_Mins","ATTRIBUTE_Sample_Type",
                          "MS_SampleType","PharmKinetics_CYP3A4"):
                if analytic[assay][field]!=row[field]:
                    raise JarmuschTemporalSourceError("Repeated source assay metadata disagree")
        else:analytic[assay]=row
    def predose(day):
        mapping={}
        for row in analytic.values():
            if (row["MS_SampleType"]=="sample"
                and row["ATTRIBUTE_Study_Day"]==day
                and row["ATTRIBUTE_Time_Point_Mins"]=="0"
                and row["ATTRIBUTE_Sample_Type"]=="blood_plasma"):
                person=row["ATTRIBUTE_Subject"]
                if not person or person in mapping:
                    raise JarmuschTemporalSourceError("Duplicate original predose biological person")
                mapping[person]=row["sample_name"]
        return mapping
    d1=predose("1");d9=predose("9")
    if len(d1)!=13 or len(d9)!=12:
        raise JarmuschTemporalSourceError("Temporal sample source mismatch")
    orig_persons=sorted(d1)
    intersect=sorted(set(d1)&set(d9))
    if len(intersect)!=11:
        raise JarmuschTemporalSourceError("Nonmatching J1/J9 original biological donor keys")
    d1_idx={person:i for i,person in enumerate(orig_persons)}
    # Biological specimen availability and actual observed PK availability
    # are DISTINCT; qualify the real complete-case intersection before fits.
    # A complete day9 plasma assay does not imply a measured day9 CL/F.
    matched=[(p,d1_idx[p]) for p in intersect]
    complete=[(p,i) for p,i in matched
              if np.isfinite(pk["day9"][i]) and pk["day9"][i]>0]
    if not 10<=len(complete)<=11:
        raise JarmuschTemporalSourceError(
            "Too few original day9 clinically observed PK outcomes in matched source")
    valid_people=[p for p,_ in complete]
    indices=[i for _,i in complete]
    y=pk["day9"][indices]
    with ZipFile(io.BytesIO(qza)) as q:
        names=[s for s in q.namelist() if s.endswith("/data/feature-table.biom")]
        if len(names)!=1:
            raise JarmuschTemporalSourceError("No canonical original GNPS BIOM")
        biom=q.read(names[0])
    with h5py.File(io.BytesIO(biom),"r") as h:
        q_samples=[v.decode() for v in h["sample/ids"][:]]
        if len(q_samples)!=1198 or len(set(q_samples))!=1198:
            raise JarmuschTemporalSourceError("Original measured MS sample IDs changed")
        mapping={v:i for i,v in enumerate(q_samples)}
        assays=[d9[p] for p in valid_people]
        if any(s not in mapping for s in assays):
            raise JarmuschTemporalSourceError("Source day9 plasma sample missing from QZA")
        mat=csr_matrix(
            (h["sample/matrix/data"][:], h["sample/matrix/indices"][:],
             h["sample/matrix/indptr"][:]),
            shape=(1198,6349))
        x9=mat[[mapping[s] for s in assays]].toarray()
    x1matched=x1[indices]
    if x9.shape!=(len(complete),6349) or x1matched.shape!=x9.shape:
        raise JarmuschTemporalSourceError("Source original paired metabolomics shape invalid")
    if np.any(x9<0) or not np.isfinite(x9).all():
        raise JarmuschTemporalSourceError("Source original J9 molecular intensities invalid")
    return x1matched,x9,c1[indices],g1[indices],y,evidence


def evaluate_matched_original_day9(archive):
    x1,x9,c,g,y,audit=_temporal_source_join(archive)
    candidates=_original_spectral_features(archive)
    results={}
    for source_clock,original_x in (
        ("plasma_day1_predose_before_cefprozil",x1),
        ("plasma_day9_predose_after_cefprozil",x9),
    ):
        results[source_clock]=_nested_cv(
            original_x,c,g,y,candidates=candidates,
            source_label="observed oral midazolam day9 CL/F")
    return {
        "status":"ACTUAL_JARMUSCH_ORIGINAL_DAY9_COMPLETE_CASE_PREDOSE_DAY1_VS_DAY9_PLASMA",
        "n_original_unique_patients_with_both_predose_times":11,
        "n_original_patients_with_both_predose_times_and_valid_day9_PK":len(y),
        "source_compartment":"peripheral blood-derived plasma, no liver mRNA/protein",
        "day1_source":"before antibiotics, baseline predicted day9 clinical phenotype",
        "day9_source":"after cefprozil course; immediately before oral drug probe dosing",
        "original_day9_pk_end_point":"observed oral midazolam CL/F (not hepatic CYP3A4-specific CL)",
        "same_complete_case_person_labels_heldout_for_both_sample_clocks":True,
        "primary_provenance_gate_passed":audit["status"]=="ORIGINAL_GNPS_SOURCE_QUALIFIED_COMPARTMENT_AND_TIMEPOINT_AUDIT",
        "source_MS1_identity_is_putative_not_functional_metabolite_confirmation":True,
        "origin_and_timing_matter_for_comparison":True,
        "no_individual_original_source_data_exported":True,
        "results":results,
    }


if __name__=="__main__":
    import sys
    if len(sys.argv)!=2:
        raise SystemExit("usage: python -m scripts.benchmark_jarmusch_temporal_plasma_nested original_gnps.zip")
    try:
        print(json.dumps(evaluate_matched_original_day9(sys.argv[1]),
                         indent=2,sort_keys=True))
    except (JarmuschTemporalSourceError, GNPSOriginalSourceError,ValueError,OSError) as exc:
        # All qualification messages in this script are static schema errors;
        # never print original rows, participant keys or patient values.
        safe=str(exc) if type(exc) is JarmuschTemporalSourceError else type(exc).__name__
        print(json.dumps({"status":"ABSTAIN_SOURCE_TEMPORAL_CROSSWALK",
                          "reason_class":type(exc).__name__,
                          "safe_source_gate":safe}))
        raise SystemExit(2)
