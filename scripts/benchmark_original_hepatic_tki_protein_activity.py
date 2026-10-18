"""Ephemeral source acquisition of ORIGINAL printed HLM Table 1 for two drugs.

Original public Europe PMC JATS XML is queried at DOI-fixed allowlisted PMCIDs;
21 unique source HLM donor rows plus one pooled HLM are read in memory, never
printed or checked into git. All original names/measurements remain transient.

Crucial: external primary studies are on the SAME 21 donor samples. A published
Table 1 is not a new assay per donor, and imatinib non-detects are NOT zeroes.
"""
from __future__ import annotations
from hashlib import sha256
import json
import re
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET

from omics2pk.hepatic_protein_activity_fusion import from_original_source_rows

SOURCES={
    "tucatinib":("PMC13095627","10.1016/j.dmd.2025.100061"),
    "imatinib":("PMC12990092","10.1016/j.dmd.2025.100196"),
}
MAX_XML=5_000_000

class OfficialTableError(ValueError):pass

def _get_table(xml_bytes,pmcid,doi):
    if len(xml_bytes)>MAX_XML or not xml_bytes.lstrip().startswith(b"<?xml") and not xml_bytes.lstrip().startswith(b"<"):
        raise OfficialTableError("Original Europe PMC XML source unavailable/invalid")
    root=ET.fromstring(xml_bytes)
    dois={x.text.strip() for x in root.findall(".//article-id") if x.attrib.get("pub-id-type")=="doi" and x.text}
    if doi not in dois:
        raise OfficialTableError("Original document DOI identity mismatch")
    candidates=[]
    for wrap in root.findall(".//table-wrap"):
        label=" ".join("".join(wrap.find("label").itertext()).split()) if wrap.find("label") is not None else ""
        if re.match(r"Table\s*1\b",label,re.I):
            candidates.append(wrap)
    if len(candidates)!=1:
        raise OfficialTableError("Original source has no unique verified Table 1")
    table=candidates[0].find(".//table")
    if table is None:
        raise OfficialTableError("Original source Table 1 lacks machine-readable table markup")
    body=table.find(".//tbody")
    if body is None:
        raise OfficialTableError("Original source Table 1 has no independent donor rows")
    data=[]
    for row in body.findall(".//tr"):
        elements=row.findall("./td")
        if not elements:continue
        cells=[" ".join("".join(c.itertext()).replace("\xa0"," ").split()) for c in elements]
        data.append(cells)
    if len(data)!=22:
        raise OfficialTableError("Original primary Table 1 must have 21 individuals + 1 original pooled HLM")
    expect_columns=13 if pmcid=="PMC13095627" else 14
    if any(len(r)!=expect_columns for r in data):
        raise OfficialTableError("Original Table 1 cell structure changed")
    original_rows=[]
    for cells in data:
        if not re.fullmatch(r"\d+|22",cells[0]):
            raise OfficialTableError("Original source donor-ordinal layout changed")
        kind="pool" if cells[1].strip().casefold()=="pooled" else "individual"
        if kind=="individual" and not re.fullmatch(r"H\d{4}",cells[1]):
            raise OfficialTableError("Original source human biological donor ID format changed")
        def number(index):
            s=cells[index].replace(",","").strip().replace("−","-")
            if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?",s):
                raise OfficialTableError("Unexpected original measured nonpositive/missing source protein/activity")
            return s
        common={
            "donor_id":cells[1], "row_kind":kind,"source_doi":doi,
            "cyp2c8_protein_pmol_mg":number(6),
            "cyp3a_protein_pmol_mg":number(9),
            "cyp2c8_activity_pmol_min_mg":number(11),
            "cyp3a_activity_pmol_min_mg":number(10),
        }
        if pmcid=="PMC13095627":
            common["tucatinib_cl_app_ul_min_mg"]=number(12)
        else:
            cl=cells[12].strip()
            common["imatinib_cl_app_ul_min_mg"]="" if cl in ("-","–","—") else number(12)
            common["ndesmethyl_imatinib_peak_area_ratio"]=(
                number(13) if kind=="individual" else "1"
            )
        original_rows.append(common)
    return original_rows


def source_published_fusion():
    rows={}
    hashes={}
    for name,(pmcid,doi) in SOURCES.items():
        uri=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
        with urlopen(Request(uri,headers={"User-Agent":"Omics2PK-original-TKI-source-research-meta-only/1.0"}),timeout=50) as response:
            raw=response.read(MAX_XML+1)
        if len(raw)>MAX_XML:
            raise OfficialTableError("Original Europe PMC fulltext exceeds bounded source cap")
        rows[name]=_get_table(raw,pmcid,doi)
        hashes[name]=sha256(raw).hexdigest()
    analysis=from_original_source_rows(rows["tucatinib"],rows["imatinib"])
    analysis["original_publisher_ XML_source_SHA256"] = hashes
    analysis["real_actual_original_article_rows_transient_only"] = True
    return analysis


if __name__=="__main__":
    try:
        print(json.dumps(source_published_fusion(),indent=2,sort_keys=True))
    except (OfficialTableError,OSError,ET.ParseError,ValueError) as err:
        print(json.dumps({"status":"OFFICIAL_ORIGINAL_TWO_DRUG_SOURCE_UNAVAILABLE_ABSTAIN",
                          "error_type":type(err).__name__,
                          "observed_source_metabolism_multi_view_scientifically_validated":False},sort_keys=True))
        raise SystemExit(2)
