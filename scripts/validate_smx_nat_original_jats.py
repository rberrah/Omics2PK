"""Source-gated benchmark from legitimately acquired original local article JATS XML.

Original donor values remain strictly in memory, never logged or committed.
"""
import argparse
import json
from omics2pk.smx_nat_original_jats import run_original_jats
from omics2pk.smx_nat_local_validation import NATSourceError


if __name__=="__main__":
    p=argparse.ArgumentParser(description="Local-only source-verified PMC12861138 Table 1 NAT/SMX evaluation")
    p.add_argument("original_article_jats_xml",
                   help="Private original publisher/PMC JATS article XML, DOI 10.1016/j.dmd.2025.100195")
    a=p.parse_args()
    try:
        report=run_original_jats(a.original_article_jats_xml)
    except (OSError,ValueError,NATSourceError) as exc:
        p.exit(2,"SOURCE_ABSTAIN: "+str(exc)+"\n")
    print(json.dumps(report,sort_keys=True,indent=2))
