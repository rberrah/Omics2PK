#!/usr/bin/env python3
"""Replay an immutable public CYP2C9 genetic-variant activity+abundance table.

GitHub original blob is SOURCE-VERIFIED via git SHA-1; network failure must
ABSTAIN with exit code 2. We never persist public large source tables,
individual variants, per-variant losses, or clinical patient information.
"""
from __future__ import annotations
import base64
import hashlib
import json
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen

from omics2pk.public_cyp2c9_variant_function_benchmark import (
    SOURCE_BLOB_SHA,SOURCE_REPOSITORY,PublicVariantDMSUnavailable,
    parse_public_original,benchmark_public_original,describe_known_class_controls,
)

URL=f"https://api.github.com/repos/{SOURCE_REPOSITORY}/git/blobs/{SOURCE_BLOB_SHA}"
MAX_RESPONSE_BYTES=5_000_000


def qualify_original_public_blob(response):
    if not isinstance(response,(bytes,bytearray)) or len(response)>MAX_RESPONSE_BYTES:
        raise PublicVariantDMSUnavailable("Wrong public original GitHub response or oversized blob")
    try:
        parsed=json.loads(response)
    except (ValueError,UnicodeError) as exc:
        raise PublicVariantDMSUnavailable("Original upstream GitHub blob not JSON") from exc
    if (parsed.get("sha")!=SOURCE_BLOB_SHA
            or parsed.get("encoding")!="base64"
            or not isinstance(parsed.get("content"),str)):
        raise PublicVariantDMSUnavailable("Original upstream immutable blob source identity changed")
    try:
        raw=base64.b64decode(parsed["content"],validate=False)
    except (ValueError,TypeError) as exc:
        raise PublicVariantDMSUnavailable("Original upstream base64 cannot be decoded") from exc
    calculated=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if calculated!=SOURCE_BLOB_SHA or parsed.get("size")!=len(raw):
        raise PublicVariantDMSUnavailable("Original git blob content cryptographic provenance failed")
    return raw


def run_original():
    request=Request(URL,headers={"User-Agent":"Omics2PK-public-known-functional-variant-replay/1.0",
                                 "Accept":"application/vnd.github+json"})
    try:
        with urlopen(request,timeout=35) as response:
            raw_response=response.read(MAX_RESPONSE_BYTES+1)
    except (HTTPError,URLError,OSError,TimeoutError) as exc:
        raise PublicVariantDMSUnavailable(
            "OFFICIAL_GITHUB_UPSTREAM_SOURCE_UNAVAILABLE_NO_REPLAY") from exc
    source=qualify_original_public_blob(raw_response)
    rows=parse_public_original(source,strict_original=True)
    result=benchmark_public_original(rows)
    result.update(describe_known_class_controls(source,strict_original=True))
    result["source_fetched_via_immutable_GitHub_git_blob_API"]=True
    result["source_bytes_and_individual_scores_not_committed_or_exported"]=True
    return result


if __name__=="__main__":
    try:
        print(json.dumps(run_original(),sort_keys=True,indent=2))
    except PublicVariantDMSUnavailable as exc:
        print(json.dumps({
            "status":"SOURCE_UNAVAILABLE_OR_FAILED_QUALIFICATION_ABSTAIN",
            "reason_type":type(exc).__name__,
            "original_source_reproduced":False,
        }))
        raise SystemExit(2)
