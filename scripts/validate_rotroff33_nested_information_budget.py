#!/usr/bin/env python3
"""Original Rotroff source 32-person XGB nested minimal omic budgets.

All patient-level source files ephemeral, aggregate only to stdout.
"""
from __future__ import annotations
import json
from urllib.request import Request,urlopen

from scripts.audit_rotroff33_metformin_original_pk_omics import LIMIT,SUPPLEMENT
from omics2pk.rotroff33_clinical_e3 import E3SourceError
from omics2pk.rotroff33_nested_information_budget import run_original_source

def main():
    with urlopen(Request(SUPPLEMENT,headers={
        "User-Agent":"Omics2PK 2026-11-04 real Rotroff nested omics 80/90% information gain-budget research"
    }),timeout=90) as reply:
        source=reply.read(LIMIT+1)
    return run_original_source(source)

if __name__=="__main__":
    try:
        print(json.dumps(main(),indent=2,sort_keys=True))
    except (ValueError,E3SourceError,OSError) as exc:
        print(json.dumps({"status":"ABSTAIN_ORIGINAL_NESTED_OMICS_BUDGET_SOURCE",
                          "error_class":type(exc).__name__,
                          "validated_clinical_improvement_claimed":False}))
        raise SystemExit(2)
