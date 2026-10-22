"""Real original GSE9588 GEO sample expression -> Omics2PK frozen networks.

NO pharmacokinetic, CYP activity, donor phenotype, PGx or RNA/protein fusion
target data are available in this imported source. This tests only whether
the Omics2PK SWEET-INSPIRED (NOT actual SWEET) frozen-reference algorithm
accepts authentic 427-column historical microarray signal and projects
unseen GSM rows without refitting reference on query array observations.

Features are 24 deterministic VARIABLE ARRAY PROBES assessed on the fixed
training 75 samples only. Array probes are NOT mapped to gene symbols here:
graphs are algorithmic stress-test examples, NOT real interpretable CYP
regulatory networks. GEO GSM names might be technical biological repeat
labels, so held-out GSM != independently held-out human donors.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import isfinite
from statistics import mean
import gzip
import json
import re
from urllib.request import Request,urlopen

from scripts.qualify_gse9588_real_human_liver_expression import (
    URL,N_EXPECTED,MAX_PROBES,MAX_UNCOMPRESSED
)

TRAIN=75
PANEL=24
WEIGHT_PANEL=512
N150_REPEAT_STRATEGIES=tuple(f"sha256_r{rep}_gsm150" for rep in range(1,5))
def _reference_size(strategy):
    if strategy=="first75":
        return TRAIN
    match=re.fullmatch(r"sha256(?:_r\d+)?_gsm(\d+)",strategy)
    if not match:
        raise ValueError("Unknown GSE9588 reference strategy")
    n=int(match.group(1))
    if n<20 or n>N_EXPECTED-20:
        raise ValueError("Reference size outside prespecified sensitivity range")
    return n


def reference_indices(gsm, strategy="first75"):
    if len(gsm)!=N_EXPECTED or len(set(gsm))!=N_EXPECTED:
        raise ValueError("Reference selection requires unique complete GSM IDs")
    n=_reference_size(strategy)
    if strategy=="first75":
        return list(range(n))
    repeat=re.fullmatch(r"sha256_r(\d+)_gsm\d+",strategy)
    if repeat:
        salt=f"Omics2PK:GSE9588:reference:v1:repeat:{repeat.group(1)}:"
        key=lambda i:(sha256((salt+gsm[i]).encode()).hexdigest(),gsm[i])
    else:
        # Preserve the historical deterministic split exactly.
        key=lambda i:(sha256(gsm[i].encode()).hexdigest(),gsm[i])
    ranked=sorted(range(len(gsm)),key=key)
    return sorted(ranked[:n])


def weight_reference_qc_only(gsm,weight_probe_values,*,reference_strategy):
    """Fast exact sample-weight QC without constructing any gene network.

    This reproduces the sample-similarity and clipping equations used by
    fit_frozen_reference/predict_network, but intentionally creates no graph.
    """
    import numpy as np
    from omics2pk.inductive_sweet import EPSILON
    weight=np.asarray(weight_probe_values,dtype=float)
    if weight.shape!=(WEIGHT_PANEL,N_EXPECTED) or not np.isfinite(weight).all():
        raise ValueError("Invalid GSE9588 weight-expression matrix")
    train_idx=reference_indices(gsm,reference_strategy)
    train_set=set(train_idx)
    query_idx=[i for i in range(N_EXPECTED) if i not in train_set]
    x=weight[:,train_idx]
    q=weight[:,query_idx]

    def unit_center(a):
        centered=a-a.mean(axis=0,keepdims=True)
        norms=np.sqrt(np.sum(centered*centered,axis=0,keepdims=True))
        if np.any(norms<=1e-12):
            raise ValueError("Constant sample in repeated reference QC")
        return centered/norms

    ztrain=unit_center(x)
    zquery=unit_center(q)
    corr=ztrain.T@ztrain
    n=len(train_idx)
    ref_means=(corr.sum(axis=1)-1.0)/(n-1)
    similarity=(ztrain.T@zquery).mean(axis=0)
    lo,hi=float(ref_means.min()),float(ref_means.max())
    mu,sd=float(ref_means.mean()),float(ref_means.std())
    if hi-lo<=1e-12 or sd<=1e-12:
        raise ValueError("Degenerate repeated reference similarity distribution")
    position=(similarity-lo+EPSILON)/(hi-lo+EPSILON)
    low=position<0
    high=position>1
    z=(similarity-mu)/sd
    return {
        "reference_strategy":reference_strategy,
        "reference_n":n,
        "query_n":len(query_idx),
        "v1_clamps":int(np.sum(low|high)),
        "v1_clamp_fraction":round(float(np.mean(low|high)),6),
        "v1_low_clamps":int(np.sum(low)),
        "v1_high_clamps":int(np.sum(high)),
        "v1_zero_weights":int(np.sum(low)),
        "v1_abs_similarity_z_gt3":int(np.sum(np.abs(z)>3)),
        "v1_mean_similarity_z":round(float(np.mean(z)),6),
        "v1_median_weight_position":round(float(np.median(position)),6),
        "reference_selection_uses_expression_values":False,
        "reference_selection_uses_PK_or_outcomes":False,
    }

def extract_train_variable_array_probe_sets(lines, *, reference_strategy="first75"):
    gsm=None
    rows=[]
    n_scan=0
    read_bytes=0
    in_data=False
    for raw in lines:
        read_bytes+=len(raw)
        if read_bytes>MAX_UNCOMPRESSED:
            raise ValueError("Real original matrix exceeded safe decompressed limit")
        if raw.startswith(b"!series_matrix_table_begin"):
            in_data=True
            continue
        if not in_data:continue
        if raw.startswith(b"!series_matrix_table_end"):break
        fields=raw.rstrip(b"\r\n").split(b"\t")
        if gsm is None:
            if fields[0].strip(b'"')!=b"ID_REF":
                raise ValueError("Original GEO microarray probe header incompatible")
            gsm=[x.decode().strip('"') for x in fields[1:]]
            if len(gsm)!=N_EXPECTED or len(set(gsm))!=N_EXPECTED:
                raise ValueError("GEO original 427 accession codes inconsistent")
            if any(not re.fullmatch(r"GSM\d{5,9}",x) for x in gsm):
                raise ValueError("Original expression column not a GEO sample ID")
            train_idx=reference_indices(gsm,reference_strategy)
            continue
        n_scan+=1
        if n_scan>MAX_PROBES:
            raise ValueError("Official GEO matrix too many probes")
        if len(fields)!=N_EXPECTED+1:
            raise ValueError("Original microarray input has unequal sample width")
        if len(rows)==WEIGHT_PANEL:continue
        try:
            vals=[float(x.strip().strip(b'"')) for x in fields[1:]]
        except ValueError:
            continue
        if not all(isfinite(x) for x in vals):
            continue
        # Probe selection uses only the selected reference samples. Query
        # variance is never consulted.
        from statistics import pstdev
        if pstdev([vals[i] for i in train_idx])<1e-8:continue
        rows.append(vals)
    if not gsm or len(rows)!=WEIGHT_PANEL or n_scan<1000:
        raise ValueError("Too few valid original GEO probe measurements")
    return gsm,rows[:PANEL],rows,n_scan


def extract_train_variable_array_probes(lines):
    """Backward-compatible historical 24-probe extractor."""
    gsm,edges,_weights,n_scan=extract_train_variable_array_probe_sets(lines)
    return gsm,edges,n_scan


def real_gse_frozen_network_smoke(
    gsm,probe_values,n_probes,*,weight_probe_values=None,
    train_indices=None,reference_strategy="first75"
):
    import numpy as np
    from omics2pk.inductive_sweet import (
        CONSTRUCTION_V0,CONSTRUCTION_V1,
        fit_frozen_reference,predict_network,gene_strength_features
    )
    mat=np.asarray(probe_values,dtype=float)
    if mat.shape!=(PANEL,N_EXPECTED):
        raise ValueError("Incorrect original GEO panel by sample schema")
    train_idx=(reference_indices(gsm,reference_strategy)
               if train_indices is None else list(train_indices))
    expected_n=_reference_size(reference_strategy)
    if (len(train_idx)!=expected_n or len(set(train_idx))!=expected_n
            or any(i<0 or i>=N_EXPECTED for i in train_idx)):
        raise ValueError("Reference size/IDs inconsistent with declared strategy")
    if train_indices is not None and train_idx!=reference_indices(gsm,reference_strategy):
        raise ValueError("Reference indices disagree with declared strategy")
    train_set=set(train_idx)
    query_idx=[i for i in range(N_EXPECTED) if i not in train_set]
    if len(query_idx)!=N_EXPECTED-expected_n:
        raise ValueError("Reference/query partition inconsistent")
    # Names are syntactic PANEL PROXIES, not genomic annotations.
    probe_names=[f"UNANNOTATED_GEO_PROBE_{i+1:02d}" for i in range(PANEL)]
    if weight_probe_values is None:
        weight_mat=mat
        weight_names=probe_names
        weight_universe="legacy_edge_panel_only"
    else:
        weight_mat=np.asarray(weight_probe_values,dtype=float)
        if weight_mat.shape!=(WEIGHT_PANEL,N_EXPECTED):
            raise ValueError("Incorrect original GEO weight universe by sample schema")
        if not np.array_equal(weight_mat[:PANEL,:],mat):
            raise ValueError("Edge panel must be the leading subset of weight universe")
        weight_names=[f"UNANNOTATED_GEO_WEIGHT_PROBE_{i+1:03d}"
                      for i in range(WEIGHT_PANEL)]
        # The first edge probes need the exact same IDs in the weight universe.
        weight_names[:PANEL]=probe_names
        weight_universe="broader_train_variable_probe_universe"
    trained=fit_frozen_reference(
        mat[:,train_idx],genes=probe_names,
        subject_ids=[gsm[i] for i in train_idx],
        tissue="Human liver (GEO GSE9588 microarray, source context)",
        k=0.1,
        weight_expression=weight_mat[:,train_idx],weight_genes=weight_names
    )
    def project(construction):
        strengths=[]
        saturated=0
        zero_weights=0
        low_clamps=0
        high_clamps=0
        similarity_z=[]
        positions_unclipped=[]
        signatures=set()
        # Exactly held-out GSM sample columns, fixed same reference; no training
        # objective is re-fit while predicting any query graph.
        for i in query_idx:
            graph,details=predict_network(
                trained,mat[:,i],subject_id=gsm[i],
                genes=probe_names,tissue=trained.tissue,
                weight_query=weight_mat[:,i],weight_genes=weight_names,
                construction=construction)
            if graph.shape!=(PANEL,PANEL) or not np.isfinite(graph).all():
                raise ValueError("Nonfinite or invalid original query graph")
            if not np.allclose(graph,graph.T,atol=1e-12) or not np.allclose(np.diag(graph),0):
                raise ValueError("Invalid self-edges or asymmetric source sample graph")
            if details["reference_fingerprint"]!=trained.fingerprint:
                raise ValueError("Original data reference was modified by query")
            s=gene_strength_features(graph,genes=probe_names)
            strengths.append(mean(s.values()))
            saturated+=bool(details["weight_clamped_outside_reference"])
            zero_weights+=details["weight"]==0
            low_clamps+=details["weight_clamp_direction"]=="low"
            high_clamps+=details["weight_clamp_direction"]=="high"
            similarity_z.append(details["sample_similarity_z_vs_reference"])
            positions_unclipped.append(details["weight_position_unclipped"])
            signatures.add(sha256(graph.tobytes()).hexdigest())
        if len(signatures)<50:
            raise ValueError("Near-duplicate or constant projected networks")
        return {
            "construction":construction,
            "unseen_GSM_samples_projected":len(strengths),
            "projected_networks_distinct_hash_count":len(signatures),
            "n_original_queries_clamped_reference_weight":saturated,
            "n_original_queries_zero_graph_update_weight":zero_weights,
            "n_original_queries_low_clamped":low_clamps,
            "n_original_queries_high_clamped":high_clamps,
            "min_weight_position_unclipped":round(float(np.min(positions_unclipped)),6),
            "median_weight_position_unclipped":round(float(np.median(positions_unclipped)),6),
            "max_weight_position_unclipped":round(float(np.max(positions_unclipped)),6),
            "mean_similarity_z_vs_reference":round(float(np.mean(similarity_z)),6),
            "min_similarity_z_vs_reference":round(float(np.min(similarity_z)),6),
            "max_similarity_z_vs_reference":round(float(np.max(similarity_z)),6),
            "n_queries_abs_similarity_z_gt3":int(np.sum(np.abs(similarity_z)>3)),
            "n_projected_individual_graph_hash_collisions":len(strengths)-len(signatures),
            "fraction_unique_heldout_network_hashes":round(len(signatures)/len(strengths),6),
            "all_heldout_graph_hashes_unique":len(signatures)==len(strengths),
            "mean_graph_absolute_node_strength_summary_only":round(mean(strengths),6),
        }

    legacy=project(CONSTRUCTION_V0)
    v1=project(CONSTRUCTION_V1)
    if (v1["n_projected_individual_graph_hash_collisions"]
            > legacy["n_projected_individual_graph_hash_collisions"]):
        raise ValueError("Query-anchored v1 unexpectedly increases graph collisions")
    return {
        "status":"REAL_PUBLIC_GSE9588_MICROARRAY_FROZEN_GRAPH_SOFTWARE_SMOKE_ONLY",
        "authentic_source_accession":"GSE9588",
        "genuine_input_sample_columns":N_EXPECTED,
        "genuine_input_original_probe_rows_scanned":n_probes,
        "unsupervised_array_probes_chosen_by_train_only_variance":PANEL,
        "fixed_reference_training_GSM_samples":expected_n,
        "reference_selection_strategy":reference_strategy,
        "reference_selection_uses_expression_values":False,
        "reference_selection_uses_PK_or_outcomes":False,
        "sample_weight_expression_probe_count":len(weight_names),
        "sample_weight_universe":weight_universe,
        "legacy_v0":legacy,
        "query_anchored_v1":v1,
        # Backward-compatible aliases preserve the historical v0 report schema.
        "unseen_GSM_samples_projected":legacy["unseen_GSM_samples_projected"],
        "projected_networks_distinct_hash_count":legacy["projected_networks_distinct_hash_count"],
        "n_original_queries_clamped_reference_weight":legacy["n_original_queries_clamped_reference_weight"],
        "n_original_queries_zero_graph_update_weight":legacy["n_original_queries_zero_graph_update_weight"],
        "n_projected_individual_graph_hash_collisions":legacy["n_projected_individual_graph_hash_collisions"],
        "fraction_unique_heldout_network_hashes":legacy["fraction_unique_heldout_network_hashes"],
        "mean_graph_absolute_node_strength_summary_only":legacy["mean_graph_absolute_node_strength_summary_only"],
        "v1_collision_reduction":(
            legacy["n_projected_individual_graph_hash_collisions"]
            - v1["n_projected_individual_graph_hash_collisions"]),
        "v1_preserves_or_improves_graph_uniqueness":(
            v1["n_projected_individual_graph_hash_collisions"]
            <= legacy["n_projected_individual_graph_hash_collisions"]),
        "v1_query_anchor_remains_subject_specific_when_weight_zero":True,
        "network_features_sufficiently_person_specific_for_individual_PK":False,
        "warning":(
            "v1 graph uniqueness is a software property only; reference-weight clamping, "
            "unannotated probes, unresolved GSM-to-donor identity and absent matched CYP outcome "
            "still block person-specific PK qualification"
        ),
        "same_original_human_liver_donors_proven_independent":False,
        "does_not_verify_GEO_GSM_to_real_biological_donor_ID":True,
        "gene_symbol_annotation_provenance":False,
        "equivalent_to_official_SWEET":False,
        "predictive_CYP_activity_test_executed":False,
        "clinical_or_invitro_PK_performance_claim":False,
        "n_original_expression_probe_rows_or_GSM_ids_exported":0,
    }


def _load_reference_strategy(strategy):
    with urlopen(Request(URL,headers={
        "User-Agent":"Omics2PK-original-427-geo-sweet-inspired-frozen-reference-matrix/1.0"
    }),timeout=90) as resp:
        with gzip.GzipFile(fileobj=resp) as source:
            ids,edges,weights,n=extract_train_variable_array_probe_sets(
                source,reference_strategy=strategy)
    idx=reference_indices(ids,strategy)
    return ids,edges,weights,n,idx


def main():
    ids,edges,weights,n,first_idx=_load_reference_strategy("first75")
    narrow=real_gse_frozen_network_smoke(
        ids,edges,n,train_indices=first_idx,reference_strategy="first75")
    broad=real_gse_frozen_network_smoke(
        ids,edges,n,weight_probe_values=weights,
        train_indices=first_idx,reference_strategy="first75")

    size_runs={}
    for reference_n in (50,75,100,150):
        strategy=f"sha256_gsm{reference_n}"
        h_ids,h_edges,h_weights,h_n,hash_idx=_load_reference_strategy(strategy)
        size_runs[str(reference_n)]=real_gse_frozen_network_smoke(
            h_ids,h_edges,h_n,weight_probe_values=h_weights,
            train_indices=hash_idx,reference_strategy=strategy)
    hashed=size_runs["75"]

    # Repeated n=150 references use four predeclared salted SHA-256 rankings.
    # Only sample-weight calibration is recomputed here; no graph or outcome is
    # needed for this robustness gate.
    repeated_n150={}
    for strategy in N150_REPEAT_STRATEGIES:
        r_ids,_r_edges,r_weights,_r_n,_r_idx=_load_reference_strategy(strategy)
        repeated_n150[strategy]=weight_reference_qc_only(
            r_ids,r_weights,reference_strategy=strategy)
    b_ids,_b_edges,b_weights,_b_n,_b_idx=_load_reference_strategy("sha256_gsm150")
    baseline150_fast=weight_reference_qc_only(
        b_ids,b_weights,reference_strategy="sha256_gsm150")
    baseline150_full=size_runs["150"]["query_anchored_v1"]
    if (
        baseline150_fast["v1_clamps"]
        !=baseline150_full["n_original_queries_clamped_reference_weight"]
        or baseline150_fast["v1_zero_weights"]
        !=baseline150_full["n_original_queries_zero_graph_update_weight"]
        or baseline150_fast["v1_abs_similarity_z_gt3"]
        !=baseline150_full["n_queries_abs_similarity_z_gt3"]
    ):
        raise ValueError("Fast repeated-reference QC disagrees with full v1 baseline")

    narrow_v1=narrow["query_anchored_v1"]
    broad_v1=broad["query_anchored_v1"]
    hash_v1=hashed["query_anchored_v1"]
    size_summary={
        n:{
            "reference_n":int(n),
            "query_n":run["query_anchored_v1"]["unseen_GSM_samples_projected"],
            "v1_clamps":run["query_anchored_v1"]["n_original_queries_clamped_reference_weight"],
            "v1_clamp_fraction":round(
                run["query_anchored_v1"]["n_original_queries_clamped_reference_weight"]
                /run["query_anchored_v1"]["unseen_GSM_samples_projected"],6),
            "v1_low_clamps":run["query_anchored_v1"]["n_original_queries_low_clamped"],
            "v1_high_clamps":run["query_anchored_v1"]["n_original_queries_high_clamped"],
            "v1_zero_weights":run["query_anchored_v1"]["n_original_queries_zero_graph_update_weight"],
            "v1_abs_similarity_z_gt3":run["query_anchored_v1"]["n_queries_abs_similarity_z_gt3"],
            "v1_mean_similarity_z":run["query_anchored_v1"]["mean_similarity_z_vs_reference"],
            "v1_median_weight_position":run["query_anchored_v1"]["median_weight_position_unclipped"],
            "v1_unique_graph_fraction":run["query_anchored_v1"]["fraction_unique_heldout_network_hashes"],
        } for n,run in size_runs.items()
    }
    return {
        "status":"REAL_PUBLIC_GSE9588_WEIGHT_AND_REFERENCE_SENSITIVITY_ONLY",
        "authentic_source_accession":"GSE9588",
        "legacy_24_probe_weighting":narrow,
        "broader_512_probe_weighting":broad,
        "sha256_reference_75_broader_512_probe_weighting":hashed,
        "sha256_reference_size_sensitivity":size_summary,
        "sha256_n150_repeated_reference_weight_qc":{
            "sha256_gsm150":baseline150_fast,
            **repeated_n150,
        },
        "n150_repeat_count_including_baseline":1+len(repeated_n150),
        "n150_all_references_zero_clamp":all(
            x["v1_clamps"]==0
            for x in [baseline150_fast,*repeated_n150.values()]),
        "n150_max_clamp_fraction":max(
            x["v1_clamp_fraction"]
            for x in [baseline150_fast,*repeated_n150.values()]),
        "v1_clamp_reduction_with_broader_weight_universe":(
            narrow_v1["n_original_queries_clamped_reference_weight"]
            - broad_v1["n_original_queries_clamped_reference_weight"]),
        "v1_low_clamp_reduction_with_broader_weight_universe":(
            narrow_v1["n_original_queries_low_clamped"]
            - broad_v1["n_original_queries_low_clamped"]),
        "v1_zero_weight_reduction_with_broader_weight_universe":(
            narrow_v1["n_original_queries_zero_graph_update_weight"]
            - broad_v1["n_original_queries_zero_graph_update_weight"]),
        "v1_hash_reference_clamp_delta_vs_first75_broad":(
            hash_v1["n_original_queries_clamped_reference_weight"]
            - broad_v1["n_original_queries_clamped_reference_weight"]),
        "v1_hash_reference_mean_similarity_z_delta_vs_first75_broad":round(
            hash_v1["mean_similarity_z_vs_reference"]
            - broad_v1["mean_similarity_z_vs_reference"],6),
        "reference_sensitivity_uses_only_GSM_identifiers_not_expression_or_outcome":True,
        "no_PK_or_CYP_outcome_used":True,
        "no_GSM_or_expression_rows_exported":True,
    }


if __name__=="__main__":
    try:
        print(json.dumps(main(),indent=2,sort_keys=True))
    except (OSError,ValueError,ImportError) as e:
        print(json.dumps({"status":"ORIGINAL_GSE9588_NETWORK_SOFTWARE_SMOKE_ABSTAIN",
                          "error_type":type(e).__name__,
                          "predictive_CYP_activity_test_executed":False}))
        raise SystemExit(2)
