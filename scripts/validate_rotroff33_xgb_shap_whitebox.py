#!/usr/bin/env python3
"""Run verified public Rotroff 2016 n33 source: XGB Tree SHAP and fold-safe white-box.

Source workbook remains ephemeral. Output aggregate-only features/accuracy.
"""
from __future__ import annotations
import json
from urllib.request import Request,urlopen
from scripts.audit_rotroff33_metformin_original_pk_omics import LIMIT,SUPPLEMENT
from omics2pk.rotroff33_clinical_e3 import E3SourceError
from omics2pk.rotroff33_xgb_shap_whitebox import run_real_source

def main():
    with urlopen(Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK-2026-11 real clinical metformin PK XGBoost and TreeSHAP audit"}),timeout=90) as response:
        payload=response.read(LIMIT+1)
    return run_real_source(payload)

if __name__=="__main__":
    try:print(json.dumps(main(),indent=2,sort_keys=True))
    except (E3SourceError,ValueError,OSError) as error:
        print(json.dumps({"status":"ABSTAIN_ORIGINAL_METFORMIN_HUMAN_XGB_PK_SOURCE",
                          "exception_class":type(error).__name__,
                          "clinical_superiority_claimed":False}))
        raise SystemExit(2)
