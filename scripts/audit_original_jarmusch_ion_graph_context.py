"""Research only: chemistry/context audit on *source-original public GNPS network*.

Objective: determine whether MS1 ions selected for tiny n=10 midazolam CL/F
models have qualified source molecular-network neighbors, putative annotated
anchors, precursor m/z redundancy or original RT metadata. No patient rows,
model fitting, metabolite identities invented, or clinical effect claims.

Source: original GNPS task b3966445abad4658a3cdc63c8198a666.
Source task params/metadata/BIOM hashes are checked by the existing
source-qualified precursor audit. The whole-study network is *transductive*,
not eligible for test-fold-selected features or independent confirmation.
"""
from __future__ import annotations

import csv
import io
import json
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile

from scripts.audit_original_jarmusch_selected_ion_mz_rt import (
    _source_selected_spectra, FOCUSED_SOURCE_FEATURE_IDS,
)
from scripts.validate_local_jarmusch14_gnps_biom_pk import (
    GNPS_TASK,REPACK_CANONICAL_MEMBER_SHA256,
)

# Differences below are ONLY compatible mass hypotheses, never identification.
MASS_SHIFT_HYPOTHESES_DA={
    "O_monooxygenation":15.994915,
    "CH2_homologation":14.015650,
    "H2O_hydration_dehydration":18.010565,
    "Na_minus_H_adduct":21.981943,
    "NH4_minus_H_adduct":17.026549,
    "13C_isotope":1.003355,
    "glucuronide_net":176.032088,
    "hexose_net":162.052824,
}
TOLERANCE_DA=0.015

def _read_original_tsv(z, prefix):
    names=[n for n in z.namelist() if n.startswith(prefix) and n.lower().endswith(".tsv")]
    if len(names)!=1:
        raise ValueError(f"Expected exactly one original public GNPS {prefix} TSV")
    with z.open(names[0]) as fh:
        reader=csv.DictReader(io.TextIOWrapper(fh,encoding="utf-8-sig"),delimiter="\t")
        rows=list(reader)
        header=reader.fieldnames or []
    return rows,header

def _maybe_float(s):
    try:
        v=float(s)
        return v if abs(v)<1e20 else None
    except (ValueError,TypeError):
        return None

def _mass_hypotheses(delta):
    return [
        {"delta_mass_compatible_with":key,"nominal_delta_Da":round(target,6)}
        for key,target in MASS_SHIFT_HYPOTHESES_DA.items()
        if abs(abs(float(delta))-target)<=TOLERANCE_DA
    ]

def audit_source_graph(archive):
    source=_source_selected_spectra(archive)
    mass=source["feature_mz_rt"]
    wanted=set(FOCUSED_SOURCE_FEATURE_IDS)
    if len(mass)!=len(wanted):
        raise ValueError("All expected source selected MS1 IDs must map to public MGF")
    with ZipFile(archive) as z:
        edges,edge_header=_read_original_tsv(z,"networking_pairs_results_file_filtered/")
        expected={"CLUSTERID1","CLUSTERID2","ComponentIndex","Cosine","DeltaMZ"}
        if not expected.issubset(edge_header):
            raise ValueError("Original GNPS network source column contract changed")
        lib,lib_header=_read_original_tsv(z,"DB_result/")
        if not {"#Scan#","Compound_Name","MQScore","SharedPeaks","MZErrorPPM"}.issubset(lib_header):
            raise ValueError("Original GNPS source library audit schema changed")
        nodes,source_columns=_read_original_tsv(z,"clusterinfo_summary/")
    verified_lib={}
    for x in lib:
        scan=x["#Scan#"].strip()
        score=_maybe_float(x.get("MQScore"))
        peaks=_maybe_float(x.get("SharedPeaks"))
        ppm=_maybe_float(x.get("MZErrorPPM"))
        if (score is None or peaks is None or ppm is None or
            score<0.7 or peaks<6 or abs(ppm)>20):
            continue
        if scan not in verified_lib or score>verified_lib[scan]["source_MS2_library_cosine"]:
            verified_lib[scan]={
                "putative_name":x["Compound_Name"][:140],
                "source_MS2_library_cosine":round(score,4),
                "not_standard_confirmed":True,
            }
    selected_edges=defaultdict(list)
    all_neighbors=defaultdict(list)
    component_by_scan={}
    for e in edges:
        a,b=e["CLUSTERID1"].strip(),e["CLUSTERID2"].strip()
        cosine=_maybe_float(e["Cosine"])
        diff=_maybe_float(e["DeltaMZ"])
        if not a or not b or cosine is None or not 0<=cosine<=1 or diff is None:
            raise ValueError("Original GNPS molecular network edge validation failed")
        component=e.get("ComponentIndex","")
        for scan,partner in ((a,b),(b,a)):
            all_neighbors[scan].append((partner,cosine,diff,component))
            if scan in wanted:
                selected_edges[scan].append((partner,cosine,diff,component))
                component_by_scan[scan]=component
    # Examine original cluster TSV for usable original source-validated retention
    # descriptors (not person-specific GNPSGROUP columns). Never log a donor row.
    lower={c.strip().lower():c for c in source_columns}
    possible_keys=("cluster index","clusterindex","cluster_id","clusterid",
                   "cluster index (original)", "#scan#","scan")
    id_col=next((lower[x] for x in possible_keys if x in lower),None)
    rt_cols=[c for c in source_columns if (
        ("rt" in c.lower() and "attribute" not in c.lower())
        or "retention time" in c.lower()
        or "retention_time" in c.lower()
    ) and not c.startswith("GNPSGROUP:")]
    # Do NOT assume source RT field units or validity from its name alone.
    feature_meta={}
    if id_col is not None:
        for row in nodes:
            key=(row.get(id_col) or "").strip()
            if key in wanted:
                feature_meta[key]={
                    "source_original_cluster_rt_fields":{
                        c:row.get(c) for c in rt_cols[:12]
                        if row.get(c) and row.get(c).strip() not in ("0","nan","NA","N/A")
                    }
                }
    # Similar masses alone are not a proof of duplicate ions without
    # actual RT and fragmentation comparison.
    close_selected=[]
    for i,a in enumerate(sorted(mass,key=int)):
        for b in sorted(mass,key=int)[i+1:]:
            dd=abs(mass[a]["original_precursor_mz"]-mass[b]["original_precursor_mz"])
            if dd<=0.015:
                close_selected.append({
                    "source_MS1_id1":a,"source_MS1_id2":b,
                    "absolute_precursor_mz_difference":round(dd,6),
                    "interpretation":"same-mass candidates, retention/isomer/adduct unresolved"
                })
    candidates={}
    for scan in sorted(wanted,key=int):
        raw=selected_edges.get(scan,[])
        sorted_neighbors=sorted(raw,key=lambda x:-x[1])
        def row(a):
            neigh,score,delta,component=a
            return {
                "neighbor_source_public_cluster":neigh,
                "original_molecular_network_cosine":round(float(score),4),
                "original_edge_delta_MZ_Da":round(float(delta),5),
                "putative_GNPS_library_match":verified_lib.get(neigh),
                "mass_difference_hypothesis_only":_mass_hypotheses(delta),
                "not_chemical_identity_evidence":True,
            }
        candidates[scan]={
            "source_original_precursor_mz":mass[scan]["original_precursor_mz"],
            "source_original_retention_time_from_MGF":mass[scan]["original_MS2_retention_time_seconds"],
            "source_cluster_retention_fields":feature_meta.get(scan,{}).get(
                "source_original_cluster_rt_fields",{}),
            "source_direct_network_degree":len(raw),
            "source_direct_component_index":component_by_scan.get(scan),
            "direct_neighbors_highest_cosine":[row(v) for v in sorted_neighbors[:5]],
            "direct_neighbors_with_source_library_annotation":[
                row(v) for v in sorted_neighbors if v[0] in verified_lib][:5],
            "source_index_has_tentative_library_match":scan in verified_lib,
            "source_library_match_if_any":verified_lib.get(scan),
        }
    return {
        "status":"SOURCE_ORIGINAL_GNPS_SELECTED_ION_MOLECULAR_GRAPH_CONTEXT_ONLY",
        "original_task":GNPS_TASK,
        "n_original_source_graph_filtered_edges":len(edges),
        "n_original_source_library_qualified_scan_matches":len(verified_lib),
        "n_original_source_cluster_summary_rows":len(nodes),
        "cluster_summary_source_index_column":id_col,
        "cluster_summary_retention_related_header_candidates":rt_cols[:12],
        "original_source_selected_ms1_ion_precursor_duplicates":close_selected,
        "studied_source_selected_ions":candidates,
        "original_network_is_not_independent_validation_or_causal_pathway":True,
        "cluster_graph_built_using_entire_public_multimatrix_study_not_train_frozen":True,
        "no_source_individual_patient_samples_printed":True,
        "mass_deltas_and_spectral_similarity_not_metabolite_identification":True,
    }

if __name__=="__main__":
    import sys
    if len(sys.argv)!=2:
        raise SystemExit("usage: python -m scripts.audit_original_jarmusch_ion_graph_context /tmp/original_gnps.zip")
    d=audit_source_graph(sys.argv[1])
    print(json.dumps(d,indent=2,sort_keys=True))
