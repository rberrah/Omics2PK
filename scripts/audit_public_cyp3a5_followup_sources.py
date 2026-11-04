"""Source-only public primary-PMC table/supplement inventory for CYP3A5 work.

NO individual subject/source rows are printed or persisted. This inventory
establishes whether an original published donor-level crosswalk actually
exists; article metadata or an image of a figure is NOT a machine-readable
numeric cohort.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from html import unescape
import json
import re
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from xml.etree import ElementTree as ET

SOURCES = {
    "grangeon": ("PMC12735975", "10.3390/ph18121789"),
    "bissada": ("PMC6815940", "10.1124/dmd.119.088823"),
    "topletz": ("PMC3876810", "10.1124/dmd.113.051094"),
}
MAX_BYTES = 8_000_000


class PublicSourceAuditError(ValueError):
    pass


def _get(url: str) -> bytes:
    with urlopen(Request(url,headers={"User-Agent":"Omics2PK-public-source-audit/1.0"}),timeout=25) as response:
        original=response.read(MAX_BYTES+1)
    if len(original)>MAX_BYTES:
        raise PublicSourceAuditError("External original fulltext over size limit")
    return original


def _simplify(value: object) -> str:
    return " ".join(str(value).split())


def _xml_metadata(data: bytes, doi: str) -> dict:
    root = ET.fromstring(data)
    dois = {x.text.strip() for x in root.findall(".//article-id")
            if x.attrib.get("pub-id-type")=="doi" and x.text}
    if doi not in dois:
        raise PublicSourceAuditError("DOI mismatch")
    tables = []
    for wrap in root.findall(".//table-wrap"):
        l=wrap.find("./label")
        label=_simplify("".join(l.itertext())) if l is not None else "UNLABELED"
        trs=wrap.findall(".//table//tr")
        distribution = sorted(Counter(len(row.findall("./td")) for row in trs).items())
        tables.append({"label":label,"rows_total":len(trs),
                       "data_cell_column_count_distribution":distribution,
                       "n_rows_with_3plus_data_cells":sum(
                           len(row.findall("./td"))>=3 for row in trs)})
    supplements = []
    for item in root.findall(".//supplementary-material"):
        label = item.find("./label")
        supplements.append({"label": _simplify("".join(label.itertext())) if label is not None else "UNLABELED",
                            "extension": _simplify(item.attrib.get("mimetype","")) + "/" +
                                         _simplify(item.attrib.get("mime-subtype","")),
                            "links_reported":len([x for x in item.iter()
                                                  if any(k.endswith("href") for k in x.attrib)])})
    return {"tables": tables, "supplementary_materials": supplements,
            "n_original_donor_numeric_matrix_verified":0,
            "no_public_patient_values_emitted":True}


def _html_metadata(data: bytes, doi: str) -> dict:
    try:
        from bs4 import BeautifulSoup
    except ImportError as exc:
        raise PublicSourceAuditError("Need beautifulsoup4 for original official PMC fallback") from exc
    soup=BeautifulSoup(data,"html.parser")
    if doi not in soup.get_text(" ",strip=True):
        raise PublicSourceAuditError("DOI not located in official article HTML")
    tables=[]
    for table in soup.find_all("table"):
        caption=table.find("caption")
        label = _simplify(caption.get_text(" ",strip=True)) if caption else "UNLABELED"
        tr=table.find_all("tr")
        counts=Counter(len(row.find_all("td",recursive=False)) for row in tr)
        tables.append({"label":label[:100], "rows_total":len(tr),
                       "data_cell_column_count_distribution":sorted(counts.items()),
                       "n_rows_with_3plus_data_cells":sum(
                           len(row.find_all("td",recursive=False))>=3 for row in tr)})
    supplements=[]
    for el in soup.find_all("a",href=True):
        href=el.get("href","")
        if re.search(r"(supplement|\-s00\d\.|\-supp\.)",href,re.I):
            extension=href.split("?")[0].rsplit(".",1)[-1].lower()
            if extension not in {"pdf","doc","docx","xlsx","xls","zip","csv"}:
                extension="webpage_or_unknown"
            supplements.append({"label":_simplify(el.get_text(" ",strip=True))[:80],
                                "extension":extension})
    return {"tables":tables[:30],
            "supplementary_materials":supplements[:35],
            "n_original_donor_numeric_matrix_verified":0,
            "no_public_patient_values_emitted":True}


def audit_primary_source(label: str) -> dict:
    if label not in SOURCES:
        raise PublicSourceAuditError("Unknown allowlisted primary source")
    pmcid, doi=SOURCES[label]
    endpoints = [
        ("official_EuropePMC_JATS",f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",_xml_metadata),
        ("official_PMC_HTML",f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/",_html_metadata),
    ]
    errors=[]
    for mode,url,parser in endpoints:
        try:
            original=_get(url)
            meta=parser(original,doi)
            return {"status":"ORIGINAL_ARTICLE_STRUCTURE_VERIFIED_NOT_DONOR_DATA",
                    "article":label,"pmcid":pmcid,"doi":doi,
                    "acquisition":mode,"source_sha256":sha256(original).hexdigest(),**meta}
        except (HTTPError,URLError,TimeoutError,OSError,ET.ParseError,
                PublicSourceAuditError) as err:
            errors.append(type(err).__name__)
    return {"status":"ABSTAIN_OFFICIAL_ARTICLE_UNAVAILABLE",
            "article":label,"pmcid":pmcid,"doi":doi,
            "source_error_categories":errors,
            "n_original_donor_numeric_matrix_verified":0,
            "no_public_patient_values_emitted":True}


if __name__=="__main__":
    outputs=[audit_primary_source(name) for name in SOURCES]
    print(json.dumps(outputs,indent=2,sort_keys=True))
    if any(r["status"]!="ORIGINAL_ARTICLE_STRUCTURE_VERIFIED_NOT_DONOR_DATA"
           for r in outputs):
        raise SystemExit(2)
