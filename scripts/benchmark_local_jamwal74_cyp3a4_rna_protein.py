"""Local authorized Jamwal74 liver donors original mRNA, protein and observed MDZ activity."""
from __future__ import annotations
import argparse
import json

from omics2pk.cyp3a4_rna_protein_mdz_source_gate import (
    SourceNotQualified,analyze_authorized_source_csv
)


if __name__=="__main__":
    args=argparse.ArgumentParser(
        description="Fail-closed same-original-liver CYP3A4 transcript/protein/midazolam activity ablation"
    )
    args.add_argument("verified_original_74_donor_csv",
                      help="Research-authorized original patient rows (never commit or upload)")
    p=args.parse_args()
    try:
        report=analyze_authorized_source_csv(p.verified_original_74_donor_csv)
    except (SourceNotQualified,OSError) as err:
        raise SystemExit("ABSTAIN: "+str(err))
    print(json.dumps(report,sort_keys=True,indent=2))
