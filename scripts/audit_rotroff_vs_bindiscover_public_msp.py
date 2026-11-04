#!/usr/bin/env python3
"""Original Rotroff 229 GC-TOF unknown-bin ID -> pinned public BinDiscover MSP screen.

Read experimental reference MSP FILES directly, not the unavailable historic API.
Public source pinned at metabolomics-us/bindiscover@045e7af0861ca04b97be902f4f557e8abb2d4940
GCBinbase_unknowns.msp ~39.9 MB, knowns_curated ~1.87 MB.
Original Rotroff A-sample sheet retrieved and SHA-locked by existing source code.
Look for exact ORIGINAL ID in *identity-related* MSP header fields; spectra are
NOT exported and matching bin ID alone NEVER proves same instrument/compound.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.request import Request,urlopen

from omics2pk.rotroff33_clinical_e3 import actual_original_tables,measured_original_data
from scripts.audit_rotroff33_metformin_original_pk_omics import LIMIT,SUPPLEMENT
from scripts.audit_rotroff33_predose_feature_labels import UNRESOLVED

UPSTREAM="metabolomics-us/bindiscover@045e7af0861ca04b97be902f4f557e8abb2d4940"
EXPECTED_SIZES={"unknown":39902584,"known":1870093}
TOP=("222169","221574","218765")
LOOKUP_KEYS={"NAME","DB#","ID","BIN","BINID","BIN_ID","ACCESSION","SYNON","SYNONYM","CAS"}
MAX_FIELDS=30

def original_source_labels(raw):
    obs,manifest=measured_original_data(raw)
    _,omics=actual_original_tables(raw)
    if len(obs)!=33:raise ValueError("Missing genuine 33 matched Rotroff samples")
    names=[str(omics[0][0].get(i,"")).strip() for i in range(1,385)]
    if len(set(names))!=384:raise ValueError("Source 384 named chemical/bins not uniquely declared")
    unknown={name for name in names if UNRESOLVED.fullmatch(name)}
    if len(unknown)!=229 or not set(TOP)<=unknown:
        raise ValueError("Source 229 original unknown bins or top prioritization changed")
    return unknown,manifest

def scan_msp(path,target_bins,expected_size):
    path=Path(path)
    if not path.is_file() or path.stat().st_size!=expected_size:
        raise ValueError("Pinned public upstream MSP file size mismatch")
    target=set(target_bins)
    if len(target)!=len(target_bins):
        raise ValueError("Duplicate original source unknown IDs")
    found={}
    n_records=0
    counts=Counter()
    first_metadata=[]
    source=sha256()
    header={}
    header_bytes=0
    in_spectrum=False
    with path.open("rb") as handle:
        for line in handle:
            source.update(line)
            s=line.decode("utf-8",errors="replace").strip()
            # Some real MSP exports have no blank lines between NAME entries.
            # A new NAME header is therefore a mandatory record boundary,
            # even if the previous record was inside its ion-list section.
            if not s or (s.upper().startswith("NAME:") and header):
                if header:
                    n_records+=1
                    if len(first_metadata)<2:
                        first_metadata.append({
                            "header_keys":sorted(header),
                            "name_example":str(header.get("NAME",""))[:80],
                            "declared_num_peaks":str(header.get("NUM PEAKS",""))[:12],
                        })
                    tokens={}
                    for key,value in header.items():
                        if key in LOOKUP_KEYS:
                            # Full numeric ID tokens only; never match a substring
                            for match in re.finditer(r"(?<![0-9])[0-9]{5,8}(?![0-9])",str(value)):
                                number=match.group()
                                if number in target:
                                    tokens[number]=key
                    for original,key in tokens.items():
                        payload={
                          "source_original_feature_identifier":original,
                          "public_msp_identity_header_field":key,
                          "public_MSP_Name":str(header.get("NAME",""))[:90],
                          "public_msp_declared_number_peaks":header.get("NUM PEAKS"),
                          "public_reference_is_same_original_Rotroff_specimen":False,
                          "same_derivatization_and_retention_index_verified":False,
                          "confirmed_molecular_structure":False,
                        }
                        if original in found:found[original]["n_candidate_records"]+=1
                        else:found[original]={**payload,"n_candidate_records":1}
                header={}
                header_bytes=0
                in_spectrum=False
                if not s:continue
            if in_spectrum:continue
            if ":" not in s:continue
            key,value=s.split(":",1)
            key=key.strip().upper()
            if not re.fullmatch(r"[A-Z0-9_ #/-]{1,32}",key):continue
            if key=="NUM PEAKS":
                # This line belongs to header although following ion list is skipped
                header[key]=value.strip()
                in_spectrum=True
            elif len(header)<MAX_FIELDS and header_bytes<6000:
                header[key]=value.strip()[:500]
                counts[key]+=1
                header_bytes+=len(s)
    if header:n_records+=1 # truncated separator at EOF: counts only, never invent identity
    if n_records<1000:
        raise ValueError("Original public MSP not a multi-record spectrum library")
    return {"n_reference_MSP_records":n_records,
            "source_original_MSP_SHA256":source.hexdigest(),
            "n_exact_identity_field_bin_matches":len(found),
            "matched_source_bins":{id:found[id] for id in sorted(found)},
            "top3_source_bins":[{ "source_bin_id":x,
                                  "status":"EXACT_PUBLIC_MSP_IDENTITY_CANDIDATE_UNCONFIRMED" if x in found else "NO_EXACT_PUBLIC_MSP_IDENTITY_FIELD_MATCH"}
                               for x in TOP],
            "first_two_original_MSP_source_record_header_schemas":first_metadata,
            "source_header_key_counts":{k:counts[k] for k in sorted(counts)},
            "full_spectra_ion_pairs_exported":False}

def main(unknown_path,known_path):
    with urlopen(Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK original Rotroff BinDiscover MSP bin accession crosswalk source-only audit"}),timeout=90) as resp:
        raw=resp.read(LIMIT+1)
    unknown,source=original_source_labels(raw)
    scan1=scan_msp(unknown_path,unknown,EXPECTED_SIZES["unknown"])
    scan2=scan_msp(known_path,unknown,EXPECTED_SIZES["known"])
    bins_in_known=set(scan2["matched_source_bins"])
    bins_in_unknown=set(scan1["matched_source_bins"])
    matched=bins_in_known|bins_in_unknown
    return {
        "status":"ORIGINAL_ROTR0FF_229_BINS_SCREENED_AGAINST_PINNED_PUBLIC_MSP_SOURCES_NOT_IDENTIFICATIONS",
        "upstream":UPSTREAM,
        "original_Rotroff_clinical_source_SHA256":source["original_clinical_workbook_SHA256"],
        "original_Rotroff_molecular_source_SHA256":source["original_untargeted_metabolomic_workbook_SHA256"],
        "n_qualified_original_predose_unknown_bins":229,
        "n_exact_public_MSP_header_accession_candidates":len(matched),
        "n_original_unknown_bins_without_exact_public_MSP_ID":229-len(matched),
        "unknown_spectra_reference":scan1,
        "knowns_curated_reference":scan2,
        "n_assigned_verified_structural_identities":0,
        "n_assigned_verified_enzyme_transporters":0,
        "n_original_patient_measurements_exported":0,
        "no_source_spectral_ion_pair_data_exported":True,
        "binbase_ID_namespace_and_instrument_compatibility_not_automatically_proved":True,
    }
if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("--unknown-msp",required=True)
    p.add_argument("--known-msp",required=True)
    args=p.parse_args()
    try:
        print(json.dumps(main(args.unknown_msp,args.known_msp),indent=2,sort_keys=True))
    except (OSError,ValueError) as e:
        print(json.dumps({"status":"ABSTAIN_ORIGINAL_BIN_ID_CROSSWALK_UNQUALIFIED",
                          "error_class":type(e).__name__,
                          "public_source_contract_reason":str(e)[:150] if isinstance(e,ValueError) else "I/O failure",
                          "chemical_identifications":0}))
        raise SystemExit(2)
