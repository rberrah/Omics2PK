#!/usr/bin/env python3
"""Ephemeral official public Rotroff 2016 n33 source → n32 heldout PK comparison.

Only aggregate metrics and publisher XLSX SHA emitted.
"""
from __future__ import annotations
import json
from urllib.request import Request,urlopen
from scripts.audit_rotroff33_metformin_original_pk_omics import LIMIT,SUPPLEMENT
from omics2pk.rotroff33_clinical_e3 import E3SourceError
from omics2pk.rotroff33_pathway_e3 import run_source_real_experiment

def main():
    with urlopen(Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK-2026-11 scientific-validated-predose-feature-PK-ablation"}),timeout=90) as response:
        data=response.read(LIMIT+1)
    return run_source_real_experiment(data)
if __name__=="__main__":
    try:print(json.dumps(main(),sort_keys=True,indent=2))
    except (OSError,E3SourceError,ValueError) as e:
        print(json.dumps({"status":"ABSTAIN_ORIGINAL_PAIRED_PK_PATHWAY_SOURCE",
                          "error_class":type(e).__name__,
                          "quantitative_PK_gain_claimed":False},sort_keys=True))
        raise SystemExit(2)
