#!/usr/bin/env python3
"""Availability and format audit of primary public PK/proteomics research sources.

This is NOT a pharmacokinetic fit or patient-level data join. Downloaded
supplements remain in RAM; output is a small metadata/status report only.
"""
from __future__ import annotations

import argparse
import io
import json
import re
import sys
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import date
from pathlib import Path
from urllib.parse import urlparse

ARTICLE = "10.1021/acs.molpharmaceut.1c00053"
WEGLER_PMID = "33739838"
WEGLER_PMCID = "PMC8041379"
YIN_ZENODO_RECORD = 15723817
MAX_BYTES = 9_000_000
EXPECTED_PMC_SUPPLEMENTS = {
    "wegler_si_methods_pdf": ("mp1c00053_si_001.pdf", "pdf"),
    "wegler_s1_proteomics": ("mp1c00053_si_002.xlsx", "xlsx"),
    "wegler_s2_bioinformatics": ("mp1c00053_si_003.xlsx", "xlsx"),
}
# Restrict all requests to explicit source URLs, not user supplied redirects.
SOURCES = {
    "wegler_si_methods_pdf": (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8041379/bin/mp1c00053_si_001.pdf", "pdf"
    ),
    "wegler_s1_proteomics": (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8041379/bin/mp1c00053_si_002.xlsx", "xlsx"
    ),
    "wegler_s2_bioinformatics": (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8041379/bin/mp1c00053_si_003.xlsx", "xlsx"
    ),
    "wegler_cloud_metadata": (
        "https://pmc-oa-opendata.s3.amazonaws.com/metadata/PMC8041379.1.json", "pmc_cloud"
    ),
    "yin_battail_metadata": (
        "https://zenodo.org/api/records/15723817", "zenodo"
    ),
}


def bounded_http(url: str, limit: int = MAX_BYTES) -> bytes:
    """Bounded and explicit public download; an upstream error is reported, not hidden."""
    parsed = urlparse(url)
    approved_static = url in {value[0] for value in SOURCES.values()}
    approved_s3 = (
        parsed.scheme == "https"
        and parsed.hostname == "pmc-oa-opendata.s3.amazonaws.com"
        and parsed.path.startswith("/PMC8041379.1/")
        and parsed.path.rsplit("/", 1)[-1] in {
            filename for filename, _ in EXPECTED_PMC_SUPPLEMENTS.values()
        }
        and not parsed.fragment
    )
    if not (approved_static or approved_s3):
        raise ValueError("Unapproved public data URL")
    request = urllib.request.Request(url, headers={
        "User-Agent": "Omics2PK-research-source-inventory/0.1 (+https://github.com/rberrah/Omics2PK)",
        "Accept": "application/octet-stream,application/json,*/*",
    })
    with urllib.request.urlopen(request, timeout=30) as response:
        # Do not follow to unrelated domains or untrusted filenames.
        final = response.geturl()
        if urlparse(final).hostname not in {"pmc.ncbi.nlm.nih.gov", "zenodo.org", "pmc-oa-opendata.s3.amazonaws.com"}:
            raise ValueError("Upstream redirected to an unapproved host")
        declared = response.headers.get("Content-Length")
        if declared and int(declared) > limit:
            raise ValueError("Source exceeds allowed in-memory size")
        payload = response.read(limit + 1)
    if len(payload) > limit:
        raise ValueError("Source exceeded in-memory size")
    return payload


def workbook_structure(data: bytes) -> dict:
    """Read only XLSX workbook sheet names; never load or print proteomic values."""
    if not zipfile.is_zipfile(io.BytesIO(data)):
        raise ValueError("Upstream file is not an XLSX ZIP")
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        names = archive.namelist()
        if "xl/workbook.xml" not in names:
            raise ValueError("XLSX is missing xl/workbook.xml")
        if len(names) > 3000:
            raise ValueError("Too many ZIP member paths")
        xml_bytes = archive.read("xl/workbook.xml")
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    root = ET.fromstring(xml_bytes)
    sheets = [s.attrib.get("name", "") for s in root.findall(".//m:sheet", ns)]
    if not sheets or any(not s for s in sheets):
        raise ValueError("Missing XLSX sheet names")
    return {"format": "xlsx", "n_sheets": len(sheets),
            "sheet_names": sheets[:35], "has_further_sheets": len(sheets) > 35}


def zenodo_metadata(data: bytes) -> dict:
    obj = json.loads(data.decode("utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("Zenodo malformed record")
    real_id = int(obj.get("id", -1))
    concept_doi = str(obj.get("conceptdoi") or "")
    # This cited DOI is a concept DOI: an API GET can resolve to a concrete
    # version. Keep the exact version ID + DOI in the metadata report.
    if real_id != YIN_ZENODO_RECORD and concept_doi != (
        "10.5281/zenodo.%s" % YIN_ZENODO_RECORD
    ):
        raise ValueError("Zenodo record mismatch: actual ID %s, DOI %s, concept DOI %s"
                         % (real_id, obj.get("doi"), concept_doi))
    metadata = obj.get("metadata", {})
    if not isinstance(metadata, dict):
        raise ValueError("Missing record metadata")
    raw_files = obj.get("files") or []
    if isinstance(raw_files, dict):
        raw_files = raw_files.get("entries", [])
        if isinstance(raw_files, dict):
            raw_files = list(raw_files.values())
    if not isinstance(raw_files, list):
        raise ValueError("Unexpected Zenodo file listing")
    names = [str(f.get("key", f.get("filename", ""))) for f in raw_files if isinstance(f, dict)]
    names = [x for x in names if x]
    return {
        "record_id": real_id,
        "cited_concept_record": YIN_ZENODO_RECORD,
        "concept_doi": concept_doi or None,
        "doi": obj.get("doi") or metadata.get("doi"),
        "title": str(metadata.get("title", ""))[:180],
        "n_public_file_entries": len(raw_files),
        "public_filenames": names[:75],
        "listing_truncated": len(names) > 75,
    }


def pmc_cloud_metadata(data: bytes) -> dict:
    """PMC's 2026 AWS distribution; validates article identity and lists files only."""
    obj = json.loads(data.decode("utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("PMC metadata is not an object")
    if str(obj.get("pmcid", "")).upper() != WEGLER_PMCID:
        raise ValueError("Wrong PMC source ID in cloud release")
    if str(obj.get("doi", "")).lower() != ARTICLE.lower():
        raise ValueError("Wrong DOI in PMC cloud release")
    items = obj.get("media_urls") or []
    if isinstance(items, dict):
        items = list(items.values())
    if not isinstance(items, list):
        raise ValueError("Invalid PMC supplementary media listing")
    urls = []
    for item in items:
        if isinstance(item, str):
            urls.append(item)
        elif isinstance(item, dict):
            candidate = item.get("url") or item.get("href")
            if isinstance(candidate, str):
                urls.append(candidate)
    supplement_filenames = [
        file for file in ("mp1c00053_si_001.pdf", "mp1c00053_si_002.xlsx", "mp1c00053_si_003.xlsx")
        if any(file in u.split("?")[0] for u in urls)
    ]
    return {
        "pmcid": WEGLER_PMCID,
        "version": obj.get("version"),
        "license_code": obj.get("license_code"),
        "is_pmc_openaccess": obj.get("is_pmc_openaccess"),
        "n_media_entries": len(urls),
        "target_media_url_shapes": [{"scheme": urlparse(u).scheme,
                                      "host": urlparse(u).hostname,
                                      "path": urlparse(u).path}
                                     for u in urls if any(
                                       filename in u for filename,_ in EXPECTED_PMC_SUPPLEMENTS.values())][:8],
        "listed_target_supplements": supplement_filenames,
        "source_distribution": "official PMC 2026 AWS metadata, no downloaded donor data",
    }


def validate_payload(kind: str, data: bytes) -> dict:
    if kind == "pdf":
        if not data.startswith(b"%PDF-"):
            raise ValueError("Expected PDF signature; possible HTML block/error response")
        return {"format": "pdf", "file_size_bytes": len(data)}
    if kind == "xlsx":
        return {**workbook_structure(data), "file_size_bytes": len(data)}
    if kind == "zenodo":
        return zenodo_metadata(data)
    if kind == "pmc_cloud":
        return pmc_cloud_metadata(data)
    raise ValueError("Unsupported source kind")


def audit(fetch=bounded_http) -> dict:
    results = {}
    cloud_payload = None
    for name, (url, fmt) in SOURCES.items():
        try:
            data = fetch(url)
            metadata = validate_payload(fmt, data)
        except (ValueError, OSError, urllib.error.URLError, TimeoutError,
                json.JSONDecodeError, ET.ParseError, zipfile.BadZipFile) as exc:
            # No raw response data or source-linked participant information is logged.
            results[name] = {
                "state": "UNAVAILABLE_OR_INVALID",
                "source_url": url,
                "error_type": type(exc).__name__,
                "detail": str(exc)[:240],
            }
        else:
            results[name] = {"state": "FORMAT_VERIFIED", "source_url": url,
                             "metadata_only": metadata}
            if name == "wegler_cloud_metadata":
                cloud_payload = data
    # The historical /bin/ paths now return 404 for this article. Only use
    # supplement object links *actually listed* in the validated 2026 NLM
    # metadata; never guess an S3 object path or donor-level file content.
    if cloud_payload is not None:
        source_obj = json.loads(cloud_payload.decode("utf-8"))
        listed = source_obj.get("media_urls") or []
        if isinstance(listed, dict):
            listed = list(listed.values())
        candidate_urls = []
        for entry in listed:
            candidate = entry if isinstance(entry, str) else (
                (entry.get("url") or entry.get("href")) if isinstance(entry, dict) else None
            )
            if isinstance(candidate, str):
                candidate_urls.append(candidate)
        for name, (filename, fmt) in EXPECTED_PMC_SUPPLEMENTS.items():
            matches = [u for u in candidate_urls if urlparse(u).path.endswith("/" + filename)]
            if len(matches) != 1:
                continue
            url = matches[0]
            try:
                payload = fetch(url)
                info = validate_payload(fmt, payload)
            except (ValueError, OSError, urllib.error.URLError, TimeoutError,
                    ET.ParseError, zipfile.BadZipFile) as exc:
                results[name]["cloud_fallback"] = {
                    "state": "UNAVAILABLE_OR_INVALID", "error_type": type(exc).__name__,
                    "detail": str(exc)[:120]
                }
            else:
                results[name] = {
                    "state": "FORMAT_VERIFIED",
                    "source_url": url,
                    "source_distribution": "verified PMC 2026 AWS metadata media_urls",
                    "metadata_only": info,
                }
    accessible = sum(item["state"] == "FORMAT_VERIFIED" for item in results.values())
    return {
        "project": "Omics2PK",
        "qualification_date_utc": date.today().isoformat(),
        "research_status": "SOURCE_INVENTORY_ONLY_NO_PAIRED_CLINT_TARGET_CONFIRMED",
        "wegler_doi": ARTICLE,
        "wegler_pmcid": WEGLER_PMCID,
        "wegler_n_donors_reported": 15,
        "wegler_substrates": ["bupropion", "diclofenac", "omeprazole", "bufuralol", "midazolam"],
        "yin_zenodo_record": YIN_ZENODO_RECORD,
        "format_verified_count": accessible,
        "source_count": len(SOURCES),
        "sources": results,
        "blocking_conditions": [
            "Do not equate proteomic source availability with access to donor-matched CLint.",
            "No verified donor-level endpoint source field, matrix mapping, assay units or HLM-vs-HH endpoint pairing.",
            "No evaluation or clinical inference until lawful reproducible feature/outcome join is proven.",
            "Yin/Battail is response/survival rather than PK and its original SWEET procedure is transductive.",
            "No controlled participant-level data or supplementary molecular matrices may be committed or uploaded.",
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--strict-sources", action="store_true")
    args = parser.parse_args()
    report = audit()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["research_status"],
                      "accessible": report["format_verified_count"],
                      "total": report["source_count"],
                      "source_states": {k: (v["state"] + ": " + v.get("error_type", "") + " " + v.get("detail", "")[:200]) for k, v in report["sources"].items()},
                      "wegler_cloud_public_metadata": report["sources"].get("wegler_cloud_metadata", {}).get("metadata_only", {}),
                      "cloud_fallback_diagnostics": {k: v.get("cloud_fallback") for k, v in report["sources"].items() if v.get("cloud_fallback")},
                      "yin_public_metadata": report["sources"].get("yin_battail_metadata", {}).get("metadata_only", {})},
                     indent=2))
    if args.strict_sources and report["format_verified_count"] != len(SOURCES):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
