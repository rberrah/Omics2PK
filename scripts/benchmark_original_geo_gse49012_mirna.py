"""Ephemeral live acquisition and verification of GEO GSE49012 original samples.

Fetch the official DOI/PMID-anchored GEO family SOFT. It contains per-sample
miRNA microarray intensities and independently measured microsomal CYP3A
activity in each GEO sample's original characteristics. All individual rows
remain in memory. Only source SHA-256 and aggregate heldout results are logged.
"""
from __future__ import annotations

from hashlib import sha256
import json
import re
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from omics2pk.geo_cyp3a_mirna import (
    GEO49012Error, MAX_SOFT_BYTES, compare_geo_mirna_vs_cyp3a_activity,
    parse_original_geo_soft,
)

ORIGINAL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49012/soft/GSE49012_family.soft.gz"


def run_original():
    req=Request(ORIGINAL,headers={
        "User-Agent":"Omics2PK original GEO CYP3A donor-heldout research/1.0"})
    with urlopen(req,timeout=50) as response:
        raw=response.read(MAX_SOFT_BYTES+1)
    samples=parse_original_geo_soft(raw)
    # Probe-identity-only metadata is safe to report. It contains no donor or
    # measured expression/phenotype values and helps verify assay aliases.
    feature_ids=sorted({k for r in samples for k in r["intensity"]})
    relevant=[k for k in feature_ids if re.search(
        r"(?:miR|mir)[-_ ]?(?:27b|206|103|107|21|130a|142)(?:\\b|[-_*])",
        k, re.IGNORECASE)]
    result=compare_geo_mirna_vs_cyp3a_activity(samples)
    result["source_n_unique_molecular_feature_IDs"]=len(feature_ids)
    result["source_molecular_feature_ID_examples"]=feature_ids[:15]
    result["source_candidate_miRNA_feature_IDs_only"]=relevant[:100]
    result["original_GEO_family_SOFT_SHA256"]=sha256(raw).hexdigest()
    result["original_GEO_public_source_url"]=ORIGINAL
    result["raw_GEO_allele_or_activity_or_miRNA_rows_exported"]=False
    result["live_source_fulltext_and_sample_join_verified"]=True
    return result


if __name__=="__main__":
    try:
        print(json.dumps(run_original(),sort_keys=True,indent=2))
    except (HTTPError,URLError,OSError,TimeoutError,GEO49012Error) as exc:
        print(json.dumps({
            "status":"ABSTAIN_GSE49012_ORIGINAL_SOURCE_UNAVAILABLE_OR_UNQUALIFIED",
            "reason_category":type(exc).__name__,
            "actual_original_GEO_CYP3A_activity_benchmark_reproduced":False,
        },sort_keys=True))
        raise SystemExit(2)
