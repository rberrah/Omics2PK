"""Local-only donor-table qualification and exploratory NA-SMX benchmark."""
import argparse
import json
from omics2pk.smx_nat_local_validation import NATSourceError, run_local_csv

def main():
    parser=argparse.ArgumentParser(
        description="Source-verified NAT1/2 vs observed NA-SMX formation (local source table only)")
    parser.add_argument("original_public_table_csv",
                        help="Locally curated original Lanphier et al. Table 1. Never commit these source donor rows.")
    args=parser.parse_args()
    try:
        result=run_local_csv(args.original_public_table_csv)
    except (OSError,ValueError,UnicodeError,NATSourceError) as exc:
        parser.exit(2,"ABSTAIN: "+str(exc)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
