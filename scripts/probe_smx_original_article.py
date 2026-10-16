"""Source reachability audit: only PMC/EPMC original XML, no donor values emitted.

A remote 403, non-article HTML, missing DOI, or unsupported licence means
ABSTAIN. Never substitute figure-derived pseudo-measurements or sample sums.
"""
from __future__ import annotations

from json import dumps
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from xml.etree import ElementTree as ET

SOURCES = {
    "europepmc": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12861138/fullTextXML",
    "pmc": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12861138/?report=xml",
    "pmc_official_oai": (
        "https://pmc.ncbi.nlm.nih.gov/api/oai/v1/mh/?verb=GetRecord"
        "&identifier=oai%3Apubmedcentral.nih.gov%3A12861138&metadataPrefix=pmc"
    ),
}
DOI="10.1016/j.dmd.2025.100195"
MAX_BYTES=3_000_000


def inspect_xml(data):
    if len(data)>MAX_BYTES or b"<!DOCTYPE" in data[:10000] or b"<!ENTITY" in data[:10000]:
        raise ValueError("Unrecognized/oversize XML article")
    root=ET.fromstring(data)
    if root.tag.endswith("OAI-PMH"):
        errors=[n.attrib.get("code","unknown") for n in root.iter()
                if n.tag.split("}")[-1]=="error"]
        if errors:
            raise ValueError("PMC_OAI_FULL_TEXT_PERMISSION_OR_SOURCE_ERROR:"+",".join(errors))
        articles=[n for n in root.iter() if n.tag.split("}")[-1]=="article"]
        if len(articles)!=1:
            raise ValueError("PMC OAI did not provide one licensed full-text article")
        root=articles[0]
    if root.tag.split("}")[-1]!="article":
        raise ValueError("Not an NLM JATS article")
    ids={x.attrib.get("pub-id-type"): (x.text or "").strip()
         for x in root.iter() if x.tag.split("}")[-1]=="article-id"}
    if ids.get("doi")!=DOI or ids.get("pmc")!="12861138":
        raise ValueError("Unrecognized DOI/PMC article")
    tables=[]
    for table in root.iter() if table.tag.split("}")[-1]=="table-wrap":
        trows=list(table.iter("tr"))
        text=" ".join(table.itertext()).lower()
        tables.append({"n_data_rows":len(trows)-1,
                       "has_original_donor_and_NA_SMX":
                           "donor" in text and ("smx" in text or "sulfamethoxazole" in text)})
    if not any(t["n_data_rows"]>=26 and t["has_original_donor_and_NA_SMX"] for t in tables):
        raise ValueError("No valid original 26-donor source table")
    return {"n_tables":len(tables),
            "matched_26_donor_table":True,
            "source_article":ids,
            "no_original_values_exported":True}


def audit():
    report={}
    for name,url in SOURCES.items():
        try:
            with urlopen(Request(url,headers={"User-Agent":
                 "Omics2PK-public-research/1.0 (source-metadata-only)",
                 "Accept":"application/xml,text/xml;q=0.9"}),timeout=25) as resp:
                final=resp.geturl()
                if not any(final.startswith(base.rsplit("/",1)[0]) for base in SOURCES.values()):
                    raise ValueError("Unapproved source redirect")
                data=resp.read(MAX_BYTES+1)
            parsed=inspect_xml(data)
            report[name]={"status":"ORIGINAL_ARTICLE_TABLE_STRUCTURE_VERIFIED",**parsed}
        except (ValueError,ET.ParseError,HTTPError,URLError,TimeoutError,OSError) as exc:
            report[name]={"status":"SOURCE_UNAVAILABLE_ABSTAIN",
                          "error_type":type(exc).__name__,
                          "error":str(exc)[:160]}
    return report


if __name__ == "__main__":
    print(dumps(audit(),sort_keys=True,indent=2))
