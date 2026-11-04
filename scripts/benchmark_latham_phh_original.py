"""Reproducible, transient source reader for public Latham PHH Table 4.

Never store individual original donor rows, send them to workflows as artifacts,
or write identifiers to stdout. Official DOI-pinned Europe PMC JATS first,
official PMC HTML fallback. Refuse a changed table layout or incomplete source.
"""
from __future__ import annotations

from hashlib import sha256
import json
import re
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

from omics2pk.latham_phh_source_gate import (
    LathamPHHError, SOURCE_DOI, SOURCE_PMCID, benchmark_latham_phh_table4,
)

MAX_BYTES = 6_000_000
JATS = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{SOURCE_PMCID}/fullTextXML"
PMC = f"https://pmc.ncbi.nlm.nih.gov/articles/{SOURCE_PMCID}/"


class PublishedTableUnavailable(LathamPHHError):
    pass


def _string(cell):
    return " ".join("".join(cell.itertext()).replace("\xa0", " ").split())


def _get_table4_xml(source: bytes):
    if len(source) > MAX_BYTES:
        raise PublishedTableUnavailable("Official XML fulltext over size limit")
    try:
        root = ET.fromstring(source)
    except ET.ParseError as exc:
        raise PublishedTableUnavailable("Official XML malformed") from exc
    ids = {i.text.strip() for i in root.findall(".//article-id")
           if i.attrib.get("pub-id-type") == "doi" and i.text}
    if SOURCE_DOI not in ids:
        raise PublishedTableUnavailable("Unexpected published source DOI")
    wraps = []
    for wrap in root.findall(".//table-wrap"):
        label = wrap.find("label")
        if label is not None and re.fullmatch(r"Table\s*4\.?",
                                               _string(label), re.I):
            wraps.append(wrap)
    if len(wraps) != 1:
        raise PublishedTableUnavailable("Original Table 4 is absent/duplicated")
    table = wraps[0].find(".//table")
    tbody = table.find(".//tbody") if table is not None else None
    if tbody is None:
        raise PublishedTableUnavailable("Original Table 4 does not contain tbody")
    cells = [[_string(x) for x in tr.findall("./td")]
             for tr in tbody.findall("./tr")]
    return _strict_table4_cells(cells)


def _get_table4_html(source: bytes):
    if len(source) > MAX_BYTES * 3:
        raise PublishedTableUnavailable("Official HTML exceeds source size limit")
    try:
        from bs4 import BeautifulSoup
    except ImportError as exc:
        raise PublishedTableUnavailable("beautifulsoup4 required for official HTML fallback") from exc
    soup = BeautifulSoup(source, "html.parser")
    plain = soup.get_text(" ", strip=True)
    if SOURCE_DOI not in plain:
        raise PublishedTableUnavailable("Official HTML has unexpected source DOI")
    candidates = []
    for table in soup.find_all("table"):
        header = " ".join(th.get_text(" ", strip=True)
                          for th in table.find_all("th")).lower()
        if not all(k in header for k in ("donor", "cyp3a5", "imatinib", "cyp2c8")):
            continue
        rows = []
        for tr in table.find_all("tr"):
            cells = [c.get_text(" ", strip=True) for c in tr.find_all("td", recursive=False)]
            if cells:
                rows.append(cells)
        if len(rows) == 14 and all(len(r) == 18 for r in rows):
            candidates.append(rows)
    if len(candidates) != 1:
        raise PublishedTableUnavailable("Official HTML lacks unique 14-donor Table 4")
    return _strict_table4_cells(candidates[0])


def _strict_table4_cells(rows):
    if len(rows) != 14 or any(len(c) != 18 for c in rows):
        raise PublishedTableUnavailable("Table 4 source row or column layout changed")
    out = []
    for cells in rows:
        if not re.fullmatch(r"[A-Z]{3}", cells[0]):
            raise PublishedTableUnavailable("Table 4 source donor identity layout changed")
        def value(idx, *, missing=False):
            text = cells[idx].replace(",", "").replace("\u2212", "-").strip()
            if missing and text in {"-", "—", "–", "", "NA", "N/A"}:
                return None
            if not re.fullmatch(r"\d+(?:\.\d+)?", text):
                raise PublishedTableUnavailable("Unexpected original table number/missingness")
            return float(text)
        out.append({
            "source_doi": SOURCE_DOI,
            "donor_id": cells[0],
            "cyp2c8_activity": value(5),
            "cyp3a_probe_activity": value(6),
            "cyp2c8_diplotype": cells[8].replace(" ", ""),
            "cyp3a5_diplotype": cells[10].replace(" ", ""),
            "published_cyp3a5_category": cells[11].strip().upper(),
            "cyp2c8_protein": value(12),
            "cyp3a4_protein": value(13, missing=True),
            "cyp3a5_protein": value(14, missing=True),
            "cyp3a_total_protein": value(15, missing=True),
            "imatinib_apparent_cl": value(16),
            "n_desmethyl_imatinib_formation": value(17),
        })
    # Contract validation is mandatory BEFORE any model computation.
    from omics2pk.latham_phh_source_gate import _qualified
    _qualified(out)
    return out


def _read_official(url: str, max_bytes: int) -> bytes:
    request = Request(url, headers={"User-Agent": "Omics2PK-research-original-PMC/1.0"})
    with urlopen(request, timeout=30) as response:
        return response.read(max_bytes + 1)


def reproduce_original():
    errors = []
    for method, url, parser, limit in (
        ("europepmc_official_JATS", JATS, _get_table4_xml, MAX_BYTES),
        ("pmc_official_HTML", PMC, _get_table4_html, MAX_BYTES * 3),
    ):
        try:
            raw = _read_official(url, limit)
            records = parser(raw)
            result = benchmark_latham_phh_table4(records)
            result["original_public_fulltext_SHA256"] = sha256(raw).hexdigest()
            result["source_acquisition_method"] = method
            result["source_rows_retained_or_released"] = False
            return result
        except (HTTPError, URLError, TimeoutError, OSError, PublishedTableUnavailable,
                LathamPHHError) as exc:
            errors.append(type(exc).__name__)
    raise PublishedTableUnavailable(
        "Official DOI-pinned public PHH source not reproducible; abstain: " + ",".join(errors))


if __name__ == "__main__":
    try:
        print(json.dumps(reproduce_original(), indent=2, sort_keys=True))
    except (PublishedTableUnavailable, LathamPHHError) as exc:
        print(json.dumps({
            "status": "OFFICIAL_PH H_SOURCE_UNAVAILABLE_ABSTAIN".replace(" ", ""),
            "reason_category": type(exc).__name__,
            "original_human_evidence_reproduced": False,
        }))
        raise SystemExit(2)
