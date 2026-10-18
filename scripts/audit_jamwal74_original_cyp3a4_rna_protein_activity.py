"""Audit ORIGINAL 74-donor Jamwal et al. CYP3A4 mRNA, protein and MDZ activity release.

Original DOI 10.1021/acs.molpharmaceut.8b00159; NIH manuscript
PMCID PMC7065849. Article measured the three variables in SAME human livers,
but a paper's summary/figures are NOT an original usable 74-row donor table.

Only official Europe PMC supplementaryFiles, bounded in memory. All DOCX
tables and embedded spreadsheet parts are inspected for machine-readable
per-donor outcome/CYP abundance/mRNA source columns. Never emit original
human IDs, raw RNA, enzyme activity, or protein measurements. Any
ambiguity => NO same-donor quantitative validation claim.

A positive structural gate still requires manual source semantics and
same-74-donor values; figures, group means and simulated PBPK are excluded.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
import re
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile,BadZipFile

DOI="10.1021/acs.molpharmaceut.8b00159"
PMCID="PMC7065849"
SOURCE=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
DIRECT_SOURCE=f"https://pmc.ncbi.nlm.nih.gov/articles/{PMCID}/bin/NIHMS1566110-supplement-Supplementary_material.docx"
MAX_PACKAGE_BYTES=12_000_000
MAX_ENTRY_BYTES=6_000_000
WN="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
SHEET=("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}")

# Direct original per-person analysis requires distinct numeric x, y, z,
# source donor IDs and NOT aggregates/predicted PBPK for imaginary groups.
ROLES={
    "cyp3a4_mrna":re.compile(r"CYP\s*3A4.{0,32}(?:mRNA|RNA|transcript)|(?:mRNA|RNA).{0,24}CYP\s*3A4",re.I),
    "cyp3a4_protein":re.compile(r"CYP\s*3A4.{0,35}(?:protein|abundance|content)|(?:protein|abundance).{0,24}CYP\s*3A4",re.I),
    "mdz_hydroxylation":re.compile(r"(?:1.{0,2}hydroxy\s*midazolam|midazolam.{0,45}(?:hydroxyl|formation|metabolic|activity|clearance|CLint)|CYP\s*3A4.{0,30}activity)",re.I),
    "donor_key":re.compile(r"(?:donor|patient|subject|sample|liver)\s*(?:ID|number|no\.?|code)",re.I),
    "population_summary":re.compile(r"(?:mean|median|95%\s*CI|geometric mean|NASH|NAFL|steatosis|sex|gender|diabetic|n\s*=\s*24|n\s*=\s*26)",re.I),
    "simulation":re.compile(r"(?:predicted|simulated|PBPK|physiologically.based)",re.I),
}
REQUIRED=("cyp3a4_mrna","cyp3a4_protein","mdz_hydroxylation","donor_key")


def _float(s):
    try:
        n=float(str(s).strip().replace(",",""))
        return n == n and abs(n)<float("inf")
    except (TypeError,ValueError,OverflowError):return False


def _row_audit(rows):
    head=[s.strip() for row in rows[:min(len(rows),8)] for s in row if len(s)<250]
    headstr=" | ".join(head)
    hits={k:bool(p.search(headstr)) for k,p in ROLES.items()}
    n_numerical_rows=sum(sum(_float(s) for s in row)>=3 for row in rows[1:])
    likely_74=73<=len(rows)-1<=76
    # DO NOT assert original donor rows just from 74 displayed group records;
    # all three role headers AND explicit donor identifier are required.
    candidate=likely_74 and all(hits[k] for k in REQUIRED) and n_numerical_rows>=65 and not hits["simulation"]
    return {
        "original_table_n_rows_including_header":len(rows),
        "max_original_cells_in_row":max((len(row) for row in rows),default=0),
        "n_source_rows_with_3plus_numeric_cells":n_numerical_rows,
        "top_8_rows_feature_roles_present":hits,
        "strict_candidate_numeric_74_donor_RNA_protein_MDZ_matrix":candidate,
        "source_donor_or_raw_values_output":False,
    }


def inspect_original_docx(payload):
    if len(payload)>MAX_ENTRY_BYTES or not payload.startswith(b"PK"):
        raise ValueError("Invalid bounded original DOCX")
    summaries=[]
    parts=Counter()
    embedded=0
    with ZipFile(BytesIO(payload)) as archive:
        for info in archive.infolist():
            if info.is_dir():continue
            ext="."+info.filename.rsplit(".",1)[-1].lower() if "." in info.filename else "unknown"
            parts[ext]+=1
            if info.filename.startswith("word/embeddings/"):
                embedded+=1
            if info.filename != "word/document.xml":continue
            if info.file_size>MAX_ENTRY_BYTES:
                raise ValueError("Unbounded original DOCX main XML")
            xml=ET.fromstring(archive.read(info))
            for table in xml.iter(WN+"tbl"):
                rows=[]
                for tr in table.findall(WN+"tr"):
                    # DOCX nested table cells must not be represented twice.
                    row=[]
                    for cell in tr.findall(WN+"tc"):
                        row.append(" ".join(
                            t.text or "" for t in cell.iter(WN+"t")
                        )[:500])
                    rows.append(row)
                summaries.append(_row_audit(rows))
    return {
        "original_docx_SHA256":sha256(payload).hexdigest(),
        "n_original_docx_tables":len(summaries),
        "n_embedded_source_files":embedded,
        "docx_internal_file_extensions":dict(sorted(parts.items())),
        "table_schema_aggregates_only":summaries,
        "n_strict_source_74person_numeric_RNA_protein_activity_tables":sum(
            x["strict_candidate_numeric_74_donor_RNA_protein_MDZ_matrix"] for x in summaries
        ),
    }


def _xls_schema(payload):
    """XLSX screening only; no person-specific values output."""
    if len(payload)>MAX_ENTRY_BYTES or not payload.startswith(b"PK"):
        raise ValueError("Invalid bounded original XLSX")
    workbook_count=0
    roles={k:False for k in ROLES}
    with ZipFile(BytesIO(payload)) as z:
        for part in z.infolist():
            if part.is_dir() or part.file_size>MAX_ENTRY_BYTES:
                continue
            if not (part.filename.startswith("xl/worksheets/sheet") or part.filename=="xl/sharedStrings.xml"):
                continue
            if not part.filename.endswith(".xml"):continue
            xml=ET.fromstring(z.read(part))
            vals=[(v.text or "") for v in xml.iter() if v.tag.endswith("}t") and v.text]
            txt=" ".join(vals[:2000])[:100000]
            for name,rx in ROLES.items():
                roles[name] |= bool(rx.search(txt))
            if part.filename.startswith("xl/worksheets/sheet"):workbook_count+=1
    return {
        "original_xlsx_SHA256":sha256(payload).hexdigest(),
        "n_worksheets":workbook_count,
        "roles_in_original_sheet_text":roles,
        "qualified_as_verified_original_numeric_74donor_table":False,
        "source_values_and_subjects_exported":False
    }


def audit_official_original_supplement(pkg):
    if len(pkg)>MAX_PACKAGE_BYTES or not pkg.startswith(b"PK"):
        raise ValueError("Official original source unavailable or oversized")
    inventory=Counter()
    files=[]
    # NIHMS supplements may be a standalone OOXML DOCX, itself a ZIP.
    with ZipFile(BytesIO(pkg)) as probe:
        if "word/document.xml" in probe.namelist():
            original=inspect_original_docx(pkg)
            original["original_source_format"]=".docx"
            return {
                "primary_article_DOI":DOI,
                "original_study_PMID":"29792708",
                "original_study_n_unique_human_livers_reported":74,
                "original_official_publisher_supplement_SHA256":sha256(pkg).hexdigest(),
                "original_supplement_member_format_counts":{".docx":1},
                "source_excel_and_word_structural_examination":[original],
                "n_candidate_original_numeric_74person_RNA_protein_activity_tables":
                    original["n_strict_source_74person_numeric_RNA_protein_activity_tables"],
                "verified_individual_74donor_RNA_protein_activity_values_available":False,
                "actual_original_74donor_RNA_vs_multimodal_clearance_comparison_justified":False,
                "simulation_not_treated_as_clinical_observed_PK":True,
                "no_original_person_identifiers_or_measured_values_exported":True,
            }
    with ZipFile(BytesIO(pkg)) as archive:
        for info in archive.infolist():
            if info.is_dir():continue
            if ".." in info.filename.split("/"):
                raise ValueError("Unsafe original source member name")
            ext="."+info.filename.rsplit(".",1)[-1].lower() if "." in info.filename else "unknown"
            inventory[ext]+=1
            if ext not in (".docx",".xlsx"):continue
            if info.file_size>MAX_ENTRY_BYTES:
                raise ValueError("Source data entry too large")
            raw=archive.read(info)
            record=(inspect_original_docx(raw) if ext==".docx" else _xls_schema(raw))
            record["original_source_format"]=ext
            files.append(record)
    proven=sum(
        x.get("n_strict_source_74person_numeric_RNA_protein_activity_tables",0)
        for x in files
    )
    return {
        "primary_article_DOI":DOI,
        "original_study_PMID":"29792708",
        "original_study_n_unique_human_livers_reported":74,
        "original_official_publisher_supplement_SHA256":sha256(pkg).hexdigest(),
        "original_supplement_member_format_counts":dict(sorted(inventory.items())),
        "source_excel_and_word_structural_examination":files,
        "n_candidate_original_numeric_74person_RNA_protein_activity_tables":proven,
        "verified_individual_74donor_RNA_protein_activity_values_available":
            False,
        "actual_original_74donor_RNA_vs_multimodal_clearance_comparison_justified":False,
        "simulation_not_treated_as_clinical_observed_PK":True,
        "no_original_person_identifiers_or_measured_values_exported":True,
    }


def main():
    source_errors=[]
    for origin,url in (("Europe_PMC_publisher_supplementary",SOURCE),
                       ("official_PMC_original_NIHMS_DOCX",DIRECT_SOURCE)):
        try:
            with urlopen(Request(url,headers={
                "User-Agent":"Omics2PK-original-CYP3A4-mRNA-protein-MDZ-PMC7065849-source-qualification/1.0"
            }),timeout=60) as response:
                data=response.read(MAX_PACKAGE_BYTES+1)
                mimetype=response.headers.get("Content-Type","")[:90]
            try:
                answer=audit_official_original_supplement(data)
                answer["actual_original_source_delivery"]=origin
                return answer
            except (ValueError,BadZipFile,ET.ParseError) as err:
                source_errors.append({"endpoint":origin,
                                      "error_type":type(err).__name__,
                                      "safe_error":str(err)[:120],
                                      "n_downloaded_bytes":len(data),
                                      "content_type":mimetype})
        except (OSError,URLError,HTTPError) as err:
            source_errors.append({"endpoint":origin,"error_type":type(err).__name__})
    return {"source_status":"OFFICIAL_PMC7065849_SUPPLEMENT_UNAVAILABLE_ABSTAIN",
            "official_source_attempt_diagnostics_only":source_errors,
            "actual_original_74donor_RNA_vs_multimodal_clearance_comparison_justified":False}


if __name__=="__main__":
    out=main()
    print(json.dumps(out,indent=2,sort_keys=True))
    if out.get("source_status")=="OFFICIAL_PMC7065849_SUPPLEMENT_UNAVAILABLE_ABSTAIN":
        # Original source inaccessible is a scientific abstention, NOT code failure.
        # Actual source evidence remains false in the emitted audit JSON.
        print("SOURCE_GATED_ABSTENTION: no original donor-level scientific fit executed")
