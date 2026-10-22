#!/usr/bin/env python3
"""Audit the open-access Zhang/POR companion supplement (PMC4986620).

The companion DMD paper explicitly states that mRNA, protein and activity of
the same 10 CYPs were measured in the same set of 100 HLMs as Zhang et al.
2016. This audit asks only whether its official Europe PMC supplement contains
source material capable of proving the donor-level crosswalk.

No original donor IDs or numerical assay values are emitted.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen

from scripts.audit_all_hepatic_rna_protein_original_public_sources import inspect_package

DOI="10.1124/dmd.116.069849"
PMCID="PMC4986620"
SOURCE=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
MAX_BYTES=20_000_000


def audit(fetch=None):
    if fetch is None:
        def fetch(url):
            req=Request(url,headers={
                "User-Agent":"Omics2PK-original-source-audit/1.0",
                "Accept":"application/octet-stream,*/*",
            })
            with urlopen(req,timeout=35) as response:
                blob=response.read(MAX_BYTES+1)
                ctype=response.headers.get("Content-Type","").split(";")[0]
            if len(blob)>MAX_BYTES:
                raise ValueError("Europe PMC supplement exceeds byte limit")
            return blob,ctype
    report={
        "source_doi":DOI,
        "source_pmcid":PMCID,
        "official_source":SOURCE,
        "same_100_HLM_panel_as_primary_Zhang_publication_level":True,
        "same_panel_statement_is_not_a_donor_crosswalk":True,
        "source_rows_acquired":False,
        "rna_protein_activity_crosswalk_qualified":False,
        "model_fitted":False,
    }
    try:
        blob,ctype=fetch(SOURCE)
        report["response_content_type"]=ctype
        report["response_bytes"]=len(blob)
        report["response_sha256"]=sha256(blob).hexdigest()
        structure=inspect_package(blob)
        report["source_structure"]=structure
        report["access"]=(
            "OFFICIAL_COMPANION_SUPPLEMENT_INSPECTED_STRUCTURE_ONLY"
            if structure.get("archive_usable")
            else "OFFICIAL_ENDPOINT_NOT_NUMERICAL_ARCHIVE"
        )
        docs=structure.get("documents") or []
        candidate=any(
            (d.get("n_rough_original_paired_table_candidates",0)>0
             or d.get("n_original_tables_with_50plus_numeric_rows",0)>0
             or d.get("max_worksheet_numeric_rows",0)>=50
             or d.get("n_column_header_roles",0)>=3)
            for d in docs
        )
        report["structural_candidate_for_donor_level_source"]=bool(candidate)
    except (HTTPError,URLError,TimeoutError,OSError,ValueError) as exc:
        report["access"]="SOURCE_UNAVAILABLE_ABSTAIN"
        report["error_type"]=type(exc).__name__
        report["structural_candidate_for_donor_level_source"]=False
    report["next_gate"]=(
        "If a structural candidate exists, inspect source schema and stable donor IDs "
        "privately before any model; otherwise retain Zhang ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED."
    )
    return report


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    out=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
