#!/usr/bin/env python3
"""One-shot local research ZIP exchange readiness audit. Prints only aggregate QC."""
import argparse
import json
from pathlib import Path

from omics2pk.multiomics_pk_readiness import inspect_pk_readiness


def main():
    parser=argparse.ArgumentParser(description="Check MultiOmics original exchange readiness for measured PK/NCA and future PopPK")
    parser.add_argument("exchange_dir",type=Path,help="Unzipped research dossier exchange/ folder (local only)")
    args=parser.parse_args()
    print(json.dumps(inspect_pk_readiness(args.exchange_dir),indent=2,sort_keys=True))


if __name__=="__main__":
    main()
