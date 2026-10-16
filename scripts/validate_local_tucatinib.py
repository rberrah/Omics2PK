"""Local ONLY: requires provenance-checked donor-aligned Table-1 CSV.
No original published donor rows included in the public code repository.
"""
import argparse
import json
from omics2pk.tucatinib_local_validation import run_local_csv, TucatinibEvidenceError

def main():
    p=argparse.ArgumentParser(description="Original donor-specific tucatinib apparent HLM clearance benchmark (local-only)")
    p.add_argument("source_table_csv")
    a=p.parse_args()
    try:
        report=run_local_csv(a.source_table_csv)
    except (TucatinibEvidenceError,OSError,UnicodeError) as exc:
        p.exit(2,"INVALID_OR_UNAVAILABLE: "+str(exc)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=="__main__":main()
