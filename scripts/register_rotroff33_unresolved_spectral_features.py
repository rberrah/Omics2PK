#!/usr/bin/env python3
"""Register real original Rotroff GC-TOF source feature labels only, NO patient rows.

Original original publisher supplementary ZIP, SHA-confirmed workbooks and
all 33 original sample A predose profiles are qualified upstream; only
first-row chemical labels / numbered bins enter SQLite. DB stays ephemeral
in GitHub Actions, can be preserved locally for next-step author BinBase join.
"""
from __future__ import annotations
from urllib.request import Request,urlopen
import argparse
import json
import re

from scripts.audit_rotroff33_metformin_original_pk_omics import LIMIT,SUPPLEMENT
from omics2pk.rotroff33_clinical_e3 import measured_original_data,actual_original_tables
from omics2pk.spectral_knowledge_bank import (
    make_db,register_original_omic_feature,summarize_registered_features,
)
from scripts.audit_rotroff33_predose_feature_labels import UNRESOLVED

STUDY="ROTROFF_METFORMIN_2016_GC_TOF"
SOURCE="https://pmc.ncbi.nlm.nih.gov/articles/PMC4906013/"

def register_verified_rotroff_features(db,raw):
    # Direct original study is SHA-locked and exact subject+predose sample
    # crosswalk is checked, but NO source patient rows enter the registry.
    original_people,metadata=measured_original_data(raw)
    _,omics=actual_original_tables(raw)
    if len(original_people)!=33 or len(omics)!=2:
        raise ValueError("Original Rotroff source schema changed")
    first=omics[0][0]
    source_ids=list(range(1,len(original_people[0]["predose_omics"])+1))
    # Original 384 source labels are stored at columns 1..384;
    # never use a sample value, patient ID, or feature intensity here.
    if len(source_ids)!=384 or max(first)!=384:
        raise ValueError("Original 384 first-row feature schema changed")
    for col in source_ids:
        label=first.get(col,"").strip()
        if not label or len(label)>250:
            raise ValueError("Missing original source chemical/bin header")
        numeric=bool(UNRESOLVED.fullmatch(label))
        register_original_omic_feature(
            db,study=STUDY,feature_id=label,label=label,
            assay="GC_EI",specimen="human_plasma_predose",
            provenance_url=SOURCE,
            known_source_chemical_name=not numeric,
        )
    counts=summarize_registered_features(db,STUDY)
    if (counts["n_registered_source_features"]!=384 or
        counts["annotation_stages"].get("UNRESOLVED_BIN_ID_ONLY_NO_SPECTRUM")!=229 or
        counts["annotation_stages"].get("SOURCE_CHEMICAL_NAME_NOT_STANDARD_CONFIRMED")!=155):
        raise ValueError("Original source 155/229 identity counts changed")
    return {**counts,
            "source_original_omics_SHA256":metadata["original_untargeted_metabolomic_workbook_SHA256"],
            "source_original_clinical_SHA256":metadata["original_clinical_workbook_SHA256"],
            "n_original_predose_study_subjects":33,
            "original_spectra_or_retention_indices_extracted":False,
            "no_clinical_or_person_level_source_data_exported":True,
            "status":"ORIGINAL_ROTROFF_384_FEATURE_LABELS_REGISTERED_SPECTRAL_IDENTITIES_STILL_UNRESOLVED"}

def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument("--db",default=":memory:",help="Local SQLite; CI uses :memory:")
    args=p.parse_args(argv)
    with urlopen(Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK spectral original Rotroff feature-label-only bank"}),timeout=90) as resp:
        raw=resp.read(LIMIT+1)
    db=make_db(args.db)
    try:return register_verified_rotroff_features(db,raw)
    finally:db.close()

if __name__=="__main__":
    try:print(json.dumps(main(),sort_keys=True,indent=2))
    except (OSError,ValueError) as err:
        print(json.dumps({"status":"ABSTAIN_ORIGINAL_ROTROFF_SOURCE_NOT_QUALIFIED",
                          "error_type":type(err).__name__,"unknown_bins_identified":0}))
        raise SystemExit(2)
