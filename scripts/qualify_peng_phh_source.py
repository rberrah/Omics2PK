"""Original-source QUALIFICATION of publicly available Peng PHH transcriptomics.

Never exports row-level RNA/functional measurements, only aggregate source
structure. Native open article / PMC AWS and GEO archives are ephemeral.
Not yet RNA-to-activity *performance* validation.
"""
from __future__ import annotations

from collections import Counter
import gzip
import hashlib
from io import BytesIO
import json
import re
import tarfile
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET
from zipfile import ZipFile

PMCID = "PMC12532551"
GEO = "GSE289636"
PMC_META = f"https://pmc-oa-opendata.s3.amazonaws.com/metadata/{PMCID}.1.json"
GEO_SOFT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289636/soft/GSE289636_family.soft.gz"
GEO_TAR = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289636/suppl/GSE289636_RAW.tar"
ORIGINAL_SUPP = "13619_2025_258_MOESM1_ESM.docx"
DOCX_NS = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
EXPECTED_BATCHES = ("025", "211", "789", "330", "409", "416")


def load(url: str, cap: int) -> bytes:
    with urlopen(Request(url, headers={"User-Agent":"Omics2PK-public-PHH-audit/0.2"}), timeout=40) as res:
        body = res.read(cap + 1)
        if len(body) > cap:
            raise ValueError("Original source exceeds declared safety size cap")
        return body


def parse_aws_metadata(blob: bytes) -> tuple[str, str]:
    record = json.loads(blob)
    if record.get("pmcid") != PMCID or str(record.get("version")) not in ("1", "PMC12532551.1"):
        raise ValueError("Unexpected article/version in official OA metadata")
    media = [s for s in record.get("media_urls", []) if isinstance(s,str) and ORIGINAL_SUPP in s]
    if len(media) != 1 or not media[0].startswith("s3://pmc-oa-opendata/PMC12532551.1/"):
        raise ValueError("Missing/ambiguous official supplement object URL")
    return media[0].replace("s3://pmc-oa-opendata/", "https://pmc-oa-opendata.s3.amazonaws.com/",1), str(record.get("license_code"))


def word_table_counts(blob: bytes) -> dict:
    with ZipFile(BytesIO(blob)) as z:
        raw = z.read("word/document.xml")
    doc = ET.fromstring(raw)
    data = []
    for t in doc.iter(DOCX_NS+"tbl"):
        rows = list(t.findall(DOCX_NS+"tr"))
        if not rows:
            continue
        blocks=[]
        for row in rows:
            cells = [" ".join("".join(n.itertext()).split()) for n in row.findall(DOCX_NS+"tc")]
            blocks.append(cells)
        all_content = " ".join(x for cells in blocks for x in cells)
        # Some tables are images or layout, not numerical assays
        data.append({
            "rows":len(rows),
            "first_three_row_cell_counts":[len(x) for x in blocks[:3]],
            "contains_cyp_marker": bool(re.search(r"CYP(?:2B6|2C9|2C19|2D6|3A4|3A5)",all_content,re.I)),
            "contains_metabolic_terms": bool(re.search(r"(?:clearance|substrate|metabolite|MFR|SCR)",all_content,re.I)),
            "donor_code_tokens": sorted({
                b for b in EXPECTED_BATCHES if re.search(r"\b(?:PHH[-_ ]?)?"+re.escape(b)+r"\b",all_content,re.I)
            }),
            "numeric_cell_count":sum(bool(re.fullmatch(r"-?\d+(?:\.\d+)?",cell)) for cells in blocks for cell in cells),
        })
    return {"table_count":len(data),"tables":data}


def geo_sample_metadata(gzipped: bytes) -> dict:
    text = gzip.decompress(gzipped).decode("utf-8", "replace")
    names = re.findall(r"^!Sample_title = ([^\r\n]+)",text,re.M)
    gsm = re.findall(r"^\^SAMPLE = (GSM\d+)",text,re.M)
    if len(names) != 6 or len(gsm) != 6:
        raise ValueError("GEO must expose exactly six original donor samples")
    decoded = []
    for name in names:
        mat=re.search(r"(?<!\d)(025|211|789|330|409|416)(?!\d)",name)
        if mat is None: raise ValueError("Unable to recover original GEO donor batch from title")
        decoded.append(mat.group(1))
    if set(decoded) != set(EXPECTED_BATCHES):
        raise ValueError("Unexpected/unmatched six GEO donors")
    return {"geo_accession":GEO, "n_geo_samples":len(gsm),
            "n_unique_gsm":len(set(gsm)), "donor_batches_verified":len(set(decoded)),
            "sample_title_batch_codes":sorted(decoded),
            "has_raw_supplement":bool(re.search(r"^!Series_supplementary_file = ",text,re.M)),
            "geo_family_sha256":hashlib.sha256(gzipped).hexdigest()}


def geo_tar_inventory(blob: bytes) -> dict:
    files = []
    with tarfile.open(fileobj=BytesIO(blob)) as t:
        for m in t.getmembers():
            if not m.isfile():
                raise ValueError("Unexpected nonregular GEO TAR object")
            files.append({"name":m.name.rsplit("/",1)[-1], "unpacked_bytes":m.size})
    batch_found = {b:sum(b in f["name"] for f in files) for b in EXPECTED_BATCHES}
    return {"n_files":len(files),"files":files,"batch_filename_presence":batch_found,
            "raw_tar_sha256":hashlib.sha256(blob).hexdigest()}


def audit_source() -> dict:
    meta=load(PMC_META,200_000)
    link, license_code = parse_aws_metadata(meta)
    supplement=load(link,12_000_000)
    docx=word_table_counts(supplement)
    geo=geo_sample_metadata(load(GEO_SOFT,10_000_000))
    raw=geo_tar_inventory(load(GEO_TAR,40_000_000))
    return {"status":"PUBLIC_SOURCE_STRUCTURE_VERIFIED_NO_RNA_PK_FIT",
            "article_pmcid":PMCID,
            "geo":geo, "geo_raw":raw, "supplement":docx,
            "supplement_license_code":license_code,
            "supplement_sha256":hashlib.sha256(supplement).hexdigest(),
            "interpretation":"Only source format, six donor GEO pairing and table structure: original metabolite/CYP outcome values and units require independent verification.",
            "no_individual_expression_or_functional_measurements_exported":True}


def main():
    print(json.dumps(audit_source(),indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
