#!/usr/bin/env python3
"""Local-only authorized multi-CYP RNA/protein/activity benchmark.

No original donor rows are printed or committed. The input CSV remains local;
stdout contains aggregate source QC and held-donor-out results only.
"""
from __future__ import annotations

import argparse
import json

from omics2pk.multicyp_rna_protein_activity_source_gate import (
    MultiCYPSourceNotQualified,
    analyze_authorized_source_csv,
)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("source_csv")
    parser.add_argument("--report",default=None)
    args=parser.parse_args()
    try:
        result=analyze_authorized_source_csv(args.source_csv)
    except MultiCYPSourceNotQualified as exc:
        result={
            "status":"ABSTAIN_ORIGINAL_MULTI_CYP_SOURCE_NOT_QUALIFIED",
            "reason":str(exc),
            "original_individual_rows_or_predictions_exported":False,
        }
        print(json.dumps(result,indent=2,sort_keys=True))
        raise SystemExit(2)
    payload=json.dumps(result,indent=2,sort_keys=True)+"\n"
    print(payload,end="")
    if args.report:
        from pathlib import Path
        Path(args.report).write_text(payload,encoding="utf-8")


if __name__=="__main__":
    main()
