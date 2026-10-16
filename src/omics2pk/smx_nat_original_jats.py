"""LOCAL-ONLY parser for original Lanphier et al. 2025 PMC/JATS Table 1.

Original article: DOI 10.1016/j.dmd.2025.100195, PMC12861138.
No original donor rows included in repository. Input must be a legitimately
acquired original XML file on the user's machine. Cannot digitize figures.
"""
from __future__ import annotations

import hashlib
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from omics2pk.smx_nat_local_validation import (
    NATSourceError, SOURCE_DOI, evaluate, qualify_original,
)
from omics2pk.tucatinib_local_validation import _spearman

PMC="12861138"
MAX_XML=5_000_000
SOURCE_ASSOCIATIONS={
    "NAT1_activity_vs_NA_SMX":0.576,
    "NAT2_activity_vs_NA_SMX":0.459,
    "NAT2_peptide_DNTDLVEFK_vs_activity":0.823,
    "NAT2_peptide_EFLNSHLLPK_vs_activity":0.734,
}
OUTPUT_FIELDS=(
    "donor_id","row_kind","source_doi","na_smx_pmol_min_million_cells",
    "na_paba_pmol_min_million_cells","na_smz_pmol_min_million_cells",
    "nat2_peptide_eflnshllpk_pmol_mg","nat2_peptide_dntdlvefk_pmol_mg",
    "nat2_genotype_phenotype","age_years","sex","bmi",
)


def _compact(element):
    text=" ".join("".join(element.itertext()).split())
    return text.replace("−","-").replace("–","-").replace("—","-").strip()


def rows_from_original_jats(data: bytes):
    if len(data)>MAX_XML or len(data)<200 or b"<!DOCTYPE" in data.upper() or b"<!ENTITY" in data.upper():
        raise NATSourceError("Invalid or unsafe original article XML input")
    try:
        root=ET.fromstring(data)
    except ET.ParseError as exc:
        raise NATSourceError("Not XML or malformed original source") from exc
    if root.tag != "article":
        raise NATSourceError("Source is not an NLM JATS article")
    ids={item.attrib.get("pub-id-type"): (item.text or "").strip()
         for item in root.iter("article-id")}
    if ids.get("doi")!=SOURCE_DOI or ids.get("pmc")!=PMC:
        raise NATSourceError("Original DOI/PMC source identity does not match")
    candidates=[]
    for container in root.iter("table-wrap"):
        for table in container.findall(".//table"):
            rows=[]
            for tr in table.iter("tr"):
                cells=[_compact(c) for c in tr if c.tag in ("th","td")]
                if cells:
                    rows.append(cells)
            if rows and len(rows[0])>=13:
                heads=[re.sub(r"[^a-z0-9]","",h.lower()) for h in rows[0]]
                if "donor" in heads[0] and "nat2" in "".join(heads) and "nasmx" in "".join(heads):
                    candidates.append(rows)
    if len(candidates)!=1:
        raise NATSourceError("Expected exactly one Table 1 donor-specific quantitative source")
    table=candidates[0]
    if len(table)!=27:
        raise NATSourceError("Expected exact original header and 26 donor rows")
    heads=[re.sub(r"[^a-z0-9]","",h.lower()) for h in table[0]]
    expected={0:"donor",1:"gender",3:"age",4:"bmi",
              7:"nat2genotypeinferredphenotype",8:"nasmx",9:"nasmz",
              10:"napaba",11:"nat2proteinconcentration",12:"nat2proteinconcentration"}
    for idx, token in expected.items():
        if token not in heads[idx]:
            raise NATSourceError(f"Original quantitative Table 1 column {idx} schema changed")
    if "eflnshllpk" not in heads[11] or "dntdlvefk" not in heads[12]:
        raise NATSourceError("Source NAT2 peptide header assignment changed")
    rows=[]
    for cells in table[1:]:
        if len(cells)!=len(table[0]):
            raise NATSourceError("Ragged original primary Table 1 row")
        values={
            "donor_id":cells[0], "row_kind":"individual","source_doi":SOURCE_DOI,
            "sex":cells[1],"age_years":cells[3],"bmi":cells[4],
            "nat2_genotype_phenotype":cells[7],
            "na_smx_pmol_min_million_cells":cells[8],
            "na_smz_pmol_min_million_cells":cells[9],
            "na_paba_pmol_min_million_cells":cells[10],
            "nat2_peptide_eflnshllpk_pmol_mg":cells[11],
            "nat2_peptide_dntdlvefk_pmol_mg":cells[12],
        }
        rows.append(values)
    qualify_original(rows)
    return rows


def published_source_validation(rows):
    """Independent original in-paper Spearman checks; not fitted endpoints."""
    people=qualify_original(rows)
    def rho(fieldx,fieldy):
        paired=[(p[fieldx],p[fieldy]) for p in people
                if p[fieldx] is not None and p[fieldy] is not None]
        if len(paired)<10:
            raise NATSourceError("Too few original matched published correlation observations")
        return _spearman([a for a,b in paired],[b for a,b in paired]),len(paired)
    found={
        "NAT1_activity_vs_NA_SMX":rho("paba","y"),
        "NAT2_activity_vs_NA_SMX":rho("smz","y"),
        "NAT2_peptide_DNTDLVEFK_vs_activity":rho("peptide2","smz"),
        "NAT2_peptide_EFLNSHLLPK_vs_activity":rho("peptide1","smz"),
    }
    if any(abs(found[k][0]-ref)>0.035 for k,ref in SOURCE_ASSOCIATIONS.items()):
        raise NATSourceError(
            "Original in-paper NAT Spearman source checks fail; cannot publish outcomes")
    return {k:{"source_rho":ref,"recomputed_rho":round(found[k][0],6),
               "n_paired_original_donors":found[k][1]}
            for k,ref in SOURCE_ASSOCIATIONS.items()}


def run_original_jats(filename: str):
    original=Path(filename).read_bytes()
    if len(original)>MAX_XML:
        raise NATSourceError("Source too large")
    rows=rows_from_original_jats(original)
    checks=published_source_validation(rows)
    result=evaluate(rows)
    return {
        **result,
        "source_document_sha256":hashlib.sha256(original).hexdigest(),
        "published_source_associations_verified":checks,
        "input":"Local-only author original PMC12861138 JATS Table 1",
        "source_rows_exported":False,
        "not_patient_pk_or_ivive":True,
    }
