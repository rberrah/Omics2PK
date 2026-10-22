#!/usr/bin/env python3
"""Audit official Zhang 2016 JPET supplemental-material routes.

Source: DOI 10.1124/jpet.116.233635, PMID 27189963, n=100 human liver
microsomes with CYP mRNA, protein content and probe activities.

This audit is metadata-only. It never parses or exports donor-level values.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.error
import urllib.request
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin,urlparse

DOI="10.1124/jpet.116.233635"
PMID="27189963"
PII="S0022356524193494"
LEGACY_SUPPLEMENT_URL=(
    "https://jpet.aspetjournals.org/lookup/suppl/doi:"+DOI+"/-/DC1"
)
CURRENT_ARTICLE_URL=(
    "https://jpet.aspetjournals.org/article/"
    "S0022-3565%2824%2919349-4/abstract"
)
MAX_BYTES=3_000_000
APPROVED_HOSTS={
    "jpet.aspetjournals.org",
    "www.sciencedirect.com",
    "sciencedirect.com",
    "linkinghub.elsevier.com",
    "ars.els-cdn.com",
}
SUPPLEMENT_EXTENSIONS={".pdf",".xls",".xlsx",".doc",".docx",".csv",".tsv",".zip"}


class _Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag.lower()=="a":
            href=dict(attrs).get("href")
            if isinstance(href,str) and href.strip():
                self.links.append(href.strip())


def _safe_url_shape(url):
    p=urlparse(url)
    return {
        "scheme":p.scheme,
        "host":p.hostname,
        "path":p.path[:500],
        "extension":Path(p.path).suffix.lower() or None,
    }


def _candidate_links(base_url,html):
    parser=_Links(); parser.feed(html)
    out=[]; seen=set()
    for href in parser.links:
        full=urljoin(base_url,href)
        p=urlparse(full)
        if p.scheme not in {"http","https"} or p.hostname not in APPROVED_HOSTS:
            continue
        low=(p.path+"?"+p.query).lower()
        ext=Path(p.path).suffix.lower()
        if not (
            ext in SUPPLEMENT_EXTENSIONS
            or "supp" in low or "mmc" in low or "additional" in low
            or "attachment" in low or "/dc1" in low
        ):
            continue
        shape=_safe_url_shape(full)
        key=(shape["host"],shape["path"],shape["extension"])
        if key not in seen:
            seen.add(key); out.append(shape)
    return out[:50]


def inspect_payload(final_url,content_type,payload):
    info={
        "final_url":_safe_url_shape(final_url),
        "content_type":content_type.split(";",1)[0].strip().lower() or None,
        "size_bytes":len(payload),
        "sha256":hashlib.sha256(payload).hexdigest(),
        "donor_level_content_parsed":False,
    }
    if payload.startswith(b"%PDF-"):
        info["payload_kind"]="pdf"
    elif payload.startswith(b"PK\x03\x04"):
        info["payload_kind"]="zip_or_office"
    elif payload.startswith(b"\xd0\xcf\x11\xe0"):
        info["payload_kind"]="legacy_office"
    else:
        sample=payload[:4096].lstrip().lower()
        if sample.startswith(b"<!doctype html") or sample.startswith(b"<html"):
            info["payload_kind"]="html"
            text=payload.decode("utf-8","replace")
            info["candidate_supplement_links"]=_candidate_links(final_url,text)
            plain=" ".join(text.lower().split())
            info["mentions_supplement"]=("supplement" in plain)
            info["mentions_doi"]=(DOI.lower() in plain)
        else:
            info["payload_kind"]="unknown"
    return info


def fetch_public(url,*,timeout=30,limit=MAX_BYTES):
    if url not in {LEGACY_SUPPLEMENT_URL,CURRENT_ARTICLE_URL}:
        raise ValueError("Only preregistered publisher routes may be fetched")
    req=urllib.request.Request(url,headers={
        "User-Agent":"Omics2PK-source-audit/0.1 (+https://github.com/rberrah/Omics2PK)",
        "Accept":"text/html,application/pdf,application/octet-stream,*/*",
    })
    try:
        with urllib.request.urlopen(req,timeout=timeout) as response:
            final=response.geturl()
            if urlparse(final).hostname not in APPROVED_HOSTS:
                raise ValueError("Publisher route redirected to unapproved host")
            declared=response.headers.get("Content-Length")
            if declared and int(declared)>limit:
                raise ValueError("Publisher response exceeds audit byte limit")
            payload=response.read(limit+1)
            if len(payload)>limit:
                raise ValueError("Publisher response exceeded audit byte limit")
            return {
                "state":"RETRIEVED",
                **inspect_payload(final,response.headers.get("Content-Type",""),payload),
            }
    except urllib.error.HTTPError as exc:
        return {
            "state":"HTTP_ERROR",
            "http_status":exc.code,
            "requested_url":_safe_url_shape(url),
            "final_url":_safe_url_shape(exc.geturl() or url),
            "donor_level_content_parsed":False,
        }
    except urllib.error.URLError as exc:
        return {
            "state":"NETWORK_ERROR",
            "error_type":type(exc.reason).__name__,
            "requested_url":_safe_url_shape(url),
            "donor_level_content_parsed":False,
        }


def audit(fetch=fetch_public):
    legacy=fetch(LEGACY_SUPPLEMENT_URL)
    article=fetch(CURRENT_ARTICLE_URL)
    candidates=[]
    for source in (legacy,article):
        candidates.extend(source.get("candidate_supplement_links",[]))
    unique=[]; seen=set()
    for item in candidates:
        key=(item.get("host"),item.get("path"),item.get("extension"))
        if key not in seen:
            seen.add(key); unique.append(item)
    direct=next(
        (x for x in (legacy,article)
         if x.get("state")=="RETRIEVED"
         and x.get("payload_kind") in {"pdf","zip_or_office","legacy_office"}),
        None,
    )
    status=(
        "OFFICIAL_SUPPLEMENT_FILE_ROUTE_RETRIEVED" if direct else
        "OFFICIAL_SUPPLEMENT_LINKS_DISCOVERED" if unique else
        "ABSTAIN_OFFICIAL_SUPPLEMENT_FILE_NOT_RETRIEVED"
    )
    return {
        "status":status,
        "source_doi":DOI,
        "source_pmid":PMID,
        "source_pii":PII,
        "published_n_human_livers":100,
        "article_reports_same_cohort_cyp_mrna_protein_activity":True,
        "article_reports_supplemental_material":True,
        "legacy_publisher_route":legacy,
        "current_publisher_route":article,
        "candidate_official_supplement_links":unique,
        "candidate_link_count":len(unique),
        "source_rows_acquired":False,
        "rna_protein_activity_crosswalk_qualified":False,
        "model_fitted":False,
        "next_gate":(
            "Retrieve a legitimate original publisher/author supplement, then "
            "qualify donor IDs, assay units and missingness before modeling."
        ),
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    report=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
