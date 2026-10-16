"""Local-only genuine published matched-donor Table 1 benchmark."""
import argparse
import json
from omics2pk.matched_tki_metabolism import MatchedMetabolismError, from_local_csv

if __name__ == "__main__":
    p=argparse.ArgumentParser(description="Actual 21-donor tucatinib/imatinib microsomal metabolism")
    p.add_argument("tucatinib_table_csv")
    p.add_argument("imatinib_table_csv")
    args=p.parse_args()
    try:
        result=from_local_csv(args.tucatinib_table_csv,args.imatinib_table_csv)
    except (MatchedMetabolismError,OSError,ValueError,UnicodeError) as exc:
        p.exit(2,"ABSTAIN: "+str(exc)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))
