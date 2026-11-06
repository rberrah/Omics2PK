"""Original Jarmusch urine versus plasma predose, true same-person midazolam PK.

Only t=0 **pre-dose** spot urine, never 2-8h post-dose substrate or
metabolite intensities (otherwise severe pharmacokinetic target leakage).
Cohort is exact original n10 with BOTH day1/day9 plasma, both day1/day9
urine, and actual observed day9 midazolam CL/F; fail closed otherwise.
DNA source tissue not inferred, urine renal handling not corrected.
No patient/sample IDs or individual feature values are ever exported.
"""
from __future__ import annotations
import csv
from io import BytesIO,TextIOWrapper
import json
from collections import defaultdict
from zipfile import ZipFile

import numpy as np

from scripts.benchmark_jarmusch_temporal_plasma_nested import (
    _temporal_source_join,JarmuschTemporalSourceError,
)
from scripts.validate_local_jarmusch14_gnps_biom_pk import (
    load_original_source,
)
from scripts.benchmark_original_jarmusch_nested_pk_germline import (
    _nested_cv,_original_spectral_features,
)


def _strict_original_urine_pair(archive):
    from scipy.sparse import csr_matrix
    import h5py

    x_plasma_1,x_plasma_9,clin,pgx,y9,_=_temporal_source_join(archive)
    _,clin_all,pgx_all,pk,_=load_original_source(archive)
    with ZipFile(archive) as z:
        with TextIOWrapper(z.open("qiime2_output/qiime2_metadata.tsv"),
                           encoding="utf-8-sig") as f:
            meta=list(csv.DictReader(f,delimiter="\t"))
        qza=z.read("qiime2_output/qiime2_table.qza")
    uniq={}
    for row in meta:
        key=row["sample_name"]
        if key and key not in uniq:
            uniq[key]=row
    def find(day,matrix):
        chosen={}
        for row in uniq.values():
            if (row.get("MS_SampleType")=="sample" and
                row.get("ATTRIBUTE_Study_Day")==day and
                row.get("ATTRIBUTE_Time_Point_Mins")=="0" and
                row.get("ATTRIBUTE_Sample_Type")==matrix):
                person=row["ATTRIBUTE_Subject"]
                if not person or person in chosen:
                    raise JarmuschTemporalSourceError(
                        "Duplicate original biological person in urine/plasma crosswalk")
                chosen[person]=row["sample_name"]
        return chosen
    plasma1=find("1","blood_plasma")
    plasma9=find("9","blood_plasma")
    urine1=find("1","urine")
    urine9=find("9","urine")
    if (len(plasma1)!=13 or len(plasma9)!=12 or
        len(urine1)!=14 or len(urine9)!=14):
        raise JarmuschTemporalSourceError(
            "Original t0 plasma or urine biological source counts changed")

    # Reconstruct EXACT original person/PK order solely in private memory,
    # from source's verified J1 plasma sorted by biological participant.
    order=sorted(plasma1)
    all_ix={p:i for i,p in enumerate(order)}
    both=sorted(set(plasma1)&set(plasma9))
    qualified=[p for p in both if np.isfinite(pk["day9"][all_ix[p]])
               and pk["day9"][all_ix[p]]>0]
    if len(qualified)!=10 or len(set(qualified))!=10:
        raise JarmuschTemporalSourceError("n10 plasma/PK complete original source changed")
    if not np.array_equal(pk["day9"][[all_ix[p] for p in qualified]],y9):
        raise JarmuschTemporalSourceError("Original PK day9 source order discordance")
    if not np.array_equal(clin_all[[all_ix[p] for p in qualified]],clin):
        raise JarmuschTemporalSourceError("Clinical person source order discordance")
    if not np.array_equal(pgx_all[[all_ix[p] for p in qualified]],pgx):
        raise JarmuschTemporalSourceError("Genotype person source order discordance")

    if any(p not in urine1 or p not in urine9 for p in qualified):
        raise JarmuschTemporalSourceError("Missing true source predose urine in n10 PK cohort")
    with ZipFile(BytesIO(qza)) as q:
        names=[s for s in q.namelist() if s.endswith("/data/feature-table.biom")]
        if len(names)!=1:
            raise JarmuschTemporalSourceError("Original unique mass spectra biom missing")
        biom=q.read(names[0])
    with h5py.File(BytesIO(biom),"r") as h:
        source_ids=[s.decode() for s in h["sample/ids"][:]]
        if len(source_ids)!=1198 or len(set(source_ids))!=1198:
            raise JarmuschTemporalSourceError("Original mass spectra sample set changed")
        idx={name:i for i,name in enumerate(source_ids)}
        matrix=csr_matrix(
            (h["sample/matrix/data"][:],h["sample/matrix/indices"][:],
             h["sample/matrix/indptr"][:]),shape=(1198,6349))
        def take(rows):
            names=[rows[p] for p in qualified]
            if not all(name in idx for name in names):
                raise JarmuschTemporalSourceError("Predose urine original sample missing from BIOM")
            z=matrix[[idx[name] for name in names]].toarray()
            if z.shape!=(10,6349) or not np.isfinite(z).all() or np.any(z<0):
                raise JarmuschTemporalSourceError("Invalid original urine molecular intensities")
            return z
        u1,u9=take(urine1),take(urine9)

    return {
        "urine_day1_t0_before_midazolam":u1,
        "urine_day9_t0_before_midazolam":u9,
        "plasma_day1_t0":x_plasma_1,
        "plasma_day9_t0":x_plasma_9,
    },clin,pgx,y9


def benchmark_original_urine(archive):
    x,clin,pgx,y9=_strict_original_urine_pair(archive)
    panels=_original_spectral_features(archive)
    out={}
    for sample_matrix,matrix in x.items():
        # Per-person compositional normalization in the existing estimator;
        # urine has no creatinine/osmolality adjustment; relative MS1 only.
        out[sample_matrix]=_nested_cv(matrix,clin,pgx,y9,
             candidates=panels,source_label="actual observed midazolam day9 oral CL/F")
    return {
        "status":"ORIGINAL_JARMUSCH_SAME10_URINE_PLASMA_T0_NESTED_DAY9_CLF",
        "n_independent_patients_with_all_t0_specimens_and_observed_PK":10,
        "source":"Original SHA-pinned GNPS BIOM and T0 source metadata",
        "specimen_comparison":"urine vs peripheral blood plasma, predose before midazolam",
        "study_day":"1 before antibiotics versus 9 after cefprozil",
        "outcome":"same genuine observed day9 midazolam CL/F for all four specimens",
        "urine_clinical_renal_covariates_and_urine_creatinine_not_available":True,
        "only_pre_drug_urine_t0_no_postdose_metabolite_target_leakage":True,
        "source_ion_MS2_annotations_are_only_putative":True,
        "clinical_germline_genotyping_available_but_not_full_star_allele_panel":True,
        "not_original_tacrolimus_PK":True,
        "reused_source_cohort_exploratory_not_independent_validation":True,
        "no_source_person_or_assay_identifiers_or_original_omics_exported":True,
        "results":out,
    }


if __name__=="__main__":
    import sys
    if len(sys.argv)!=2:
        raise SystemExit("expected original GNPS ZIP path")
    try:
        print(json.dumps(benchmark_original_urine(sys.argv[1]),indent=2,sort_keys=True))
    except (JarmuschTemporalSourceError,ValueError,OSError,KeyError) as exc:
        print(json.dumps({"status":"ABSTAIN_ORIGINAL_URINE_PK_MATCH",
                          "source_gate_class":type(exc).__name__,
                          "static_safe_reason":(str(exc) if type(exc) is JarmuschTemporalSourceError
                                                else type(exc).__name__)}))
        raise SystemExit(2)
