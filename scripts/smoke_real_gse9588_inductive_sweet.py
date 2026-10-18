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


def extract_train_variable_array_probes(lines):
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
            continue
        n_scan+=1
        if n_scan>MAX_PROBES:
            raise ValueError("Official GEO matrix too many probes")
        if len(fields)!=N_EXPECTED+1:
            raise ValueError("Original microarray input has unequal sample width")
        if len(rows)==PANEL:continue
        try:
            vals=[float(x.strip().strip(b'"')) for x in fields[1:]]
        except ValueError:
            continue
        if not all(isfinite(x) for x in vals):
            continue
        # This RNA probe universe is fixed solely by the FIRST ORIGINAL
        # TRAINING samples. Heldout sample variance is not consulted.
        from statistics import pstdev
        if pstdev(vals[:TRAIN])<1e-8:continue
        rows.append(vals)
    if not gsm or len(rows)!=PANEL or n_scan<1000:
        raise ValueError("Too few valid original GEO probe measurements")
    return gsm,rows,n_scan


def real_gse_frozen_network_smoke(gsm,probe_values,n_probes):
    import numpy as np
    from omics2pk.inductive_sweet import (
        fit_frozen_reference,predict_network,gene_strength_features
    )
    mat=np.asarray(probe_values,dtype=float)
    if mat.shape!=(PANEL,N_EXPECTED):
        raise ValueError("Incorrect original GEO panel by sample schema")
    # Names are syntactic PANEL PROXIES, not genomic annotations.
    probe_names=[f"UNANNOTATED_GEO_PROBE_{i+1:02d}" for i in range(PANEL)]
    trained=fit_frozen_reference(
        mat[:,:TRAIN],genes=probe_names,subject_ids=gsm[:TRAIN],
        tissue="Human liver (GEO GSE9588 microarray, source context)",
        k=0.1
    )
    strengths=[]
    saturated=0
    example_signatures=set()
    # Exactly held-out GSM sample columns, fixed same reference; no training
    # objective is re-fit while predicting the 352 sample graphs.
    for i in range(TRAIN,N_EXPECTED):
        graph,details=predict_network(
            trained,mat[:,i],subject_id=gsm[i],
            genes=probe_names,tissue=trained.tissue)
        if graph.shape!=(PANEL,PANEL) or not np.isfinite(graph).all():
            raise ValueError("Nonfinite or invalid original query graph")
        if not np.allclose(graph,graph.T,atol=1e-12) or not np.allclose(np.diag(graph),0):
            raise ValueError("Invalid self-edges or asymmetric source sample graph")
        if details["reference_fingerprint"]!=trained.fingerprint:
            raise ValueError("Original data reference was modified by query")
        s=gene_strength_features(graph,genes=probe_names)
        strengths.append(mean(s.values()))
        saturated+=bool(details["weight_clamped_outside_reference"])
        example_signatures.add(sha256(graph.tobytes()).hexdigest())
    if len(example_signatures)<50:
        raise ValueError("Near-duplicate or constant projected networks")
    return {
        "status":"REAL_PUBLIC_GSE9588_MICROARRAY_FROZEN_GRAPH_SOFTWARE_SMOKE_ONLY",
        "authentic_source_accession":"GSE9588",
        "genuine_input_sample_columns":N_EXPECTED,
        "genuine_input_original_probe_rows_scanned":n_probes,
        "unsupervised_array_probes_chosen_by_train_only_variance":PANEL,
        "fixed_reference_training_GSM_samples":TRAIN,
        "unseen_GSM_samples_projected":len(strengths),
        "projected_networks_distinct_hash_count":len(example_signatures),
        "n_original_queries_clamped_reference_weight":saturated,
        "mean_graph_absolute_node_strength_summary_only":round(mean(strengths),6),
        "same_original_human_liver_donors_proven_independent":False,
        "does_not_verify_GEO_GSM_to_real_biological_donor_ID":True,
        "gene_symbol_annotation_provenance":False,
        "equivalent_to_official_SWEET":False,
        "predictive_CYP_activity_test_executed":False,
        "clinical_or_invitro_PK_performance_claim":False,
        "n_original_expression_probe_rows_or_GSM_ids_exported":0,
    }


def main():
    with urlopen(Request(URL,headers={
        "User-Agent":"Omics2PK-original-427-geo-sweet-inspired-frozen-reference-matrix/1.0"
    }),timeout=90) as resp:
        with gzip.GzipFile(fileobj=resp) as source:
            ids,features,n=extract_train_variable_array_probes(source)
    return real_gse_frozen_network_smoke(ids,features,n)


if __name__=="__main__":
    try:
        print(json.dumps(main(),indent=2,sort_keys=True))
    except (OSError,ValueError,ImportError) as e:
        print(json.dumps({"status":"ORIGINAL_GSE9588_NETWORK_SOFTWARE_SMOKE_ABSTAIN",
                          "error_type":type(e).__name__,
                          "predictive_CYP_activity_test_executed":False}))
        raise SystemExit(2)
