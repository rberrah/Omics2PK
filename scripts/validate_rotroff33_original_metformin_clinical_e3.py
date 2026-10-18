"""Original publisher n33 Rotroff predose omics → actually measured n32 metformin PK.

Single ephemeral original supplementary acquisition, original workbook SHA
verification, exact true 33 biological subject-key matching, source A-only
measurement, and 32 outer donor heldout folds per AUClast/renal clearance.
Only aggregate source provenance and model diagnostics are printed.
"""
from __future__ import annotations
import json
from urllib.request import Request,urlopen

from scripts.audit_rotroff33_metformin_original_pk_omics import LIMIT,SUPPLEMENT
from omics2pk.rotroff33_clinical_e3 import (
    E3SourceError, measured_original_data,evaluate_e3_original
)


def main():
    with urlopen(Request(SUPPLEMENT,headers={
        "User-Agent":"Omics2PK-real-public-n33-predose-GCTOF-metformin-PK-E3-benchmark/1.0"
    }),timeout=65) as res:
        original=res.read(LIMIT+1)
    observations,provenance=measured_original_data(original)
    return evaluate_e3_original(observations,provenance)


if __name__=="__main__":
    try:
        print(json.dumps(main(),indent=2,sort_keys=True))
    except (OSError,E3SourceError,ValueError) as e:
        print(json.dumps({
            "status":"SOURCE_QUALIFICATION_OR_NUMERIC_E3_VALIDATION_ABSTAIN",
            "error_type":type(e).__name__,
            "genuine_original_clinical_e3_benchmark_reproduced":False,
        },indent=2,sort_keys=True))
        raise SystemExit(2)
