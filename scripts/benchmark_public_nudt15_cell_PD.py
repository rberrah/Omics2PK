#!/usr/bin/env python3
"""Actual original MaveDB NUDT15 protein stability vs thiopurine cell survival.

Source set immutable URN, currently pin original content SHA256. Refuse any
unexpected source mutation or incomplete variant mapping. No row export.
"""
import hashlib
import json
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

from omics2pk.public_nudt15_variant_pd_benchmark import (
    URN,EXPECTED_SOURCE_SHA256,parse_original_mavedb,
    benchmark_nudt15_original,NUDT15PublicSourceError,
)

URL="https://api.mavedb.org/api/v1/score-sets/"+URN+"/scores"


def replay():
    try:
        with urlopen(Request(URL,headers={
             "User-Agent":"Omics2PK-public-NUDT15-20261105/1.0",
             "Accept":"text/csv"}),timeout=40) as r:
            b=r.read(5_000_001)
    except (HTTPError,URLError,OSError,TimeoutError) as exc:
        raise NUDT15PublicSourceError("SOURCE_UNAVAILABLE_NO_SYNTHETIC_SUBSTITUTION") from exc
    if hashlib.sha256(b).hexdigest()!=EXPECTED_SOURCE_SHA256:
        raise NUDT15PublicSourceError("ORIGINAL_MAVEDB_SCORE_SET_CONTENT_SHA256_CHANGED")
    data,meta=parse_original_mavedb(b,strict_original=True)
    out=benchmark_nudt15_original(data)
    out.update(meta)
    out["original_public_source_SHA256"]=EXPECTED_SOURCE_SHA256
    out["source_fetched_live_no_variant_table_committed"]=True
    return out


if __name__=="__main__":
    try:
        print(json.dumps(replay(),indent=2,sort_keys=True))
    except NUDT15PublicSourceError as exc:
        print(json.dumps({
            "status":"ABSTAIN_PUBLIC_SOURCE_NOT_QUALIFIED",
            "reason_category":str(exc)[:220],
            "source_replay_validated":False,
        }))
        raise SystemExit(2)
