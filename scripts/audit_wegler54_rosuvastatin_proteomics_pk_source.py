"""Source-only audit: Wegler et al. 2021 n54 clinical rosuvastatin–proteomics.

Study DOI 10.1002/cpt.2056, PMCID PMC7984432. Official open-access
publisher supplementary S1 was described as donor-matched concentrations
of OATP1B1/1B3/2B1/NTCP and allometrics. This audit tests whether those
rows ALSO contain true measured individual AUC/Cmax outcomes.

Prints ONLY source hashes, aggregate dimensions and coarse binary role flags.
Never prints clinical rows, donor identifiers, column names, observed AUCs or
protein concentrations. Source unavailability or metadata-only data ABSTAIN.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
import re
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from xml.etree import ElementTree as ET
from zipfile import ZipFile, BadZipFile

DOI = "10.1002/cpt.2056"
PMCID = "PMC7984432"
PUBLIC_SUPP = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
MAX_PACKAGE_BYTES = 24_000_000
MAX_WORKBOOK_BYTES = 12_000_000
XMLNS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
ROLES = {
    "oatp1b1": re.compile(r"\b(?:oatp\s*1b1|slco\s*1b1)\b",re.I),
    "oatp1b3": re.compile(r"\b(?:oatp\s*1b3|slco\s*1b3)\b",re.I),
    "oatp2b1": re.compile(r"\b(?:oatp\s*2b1|slco\s*2b1)\b",re.I),
    "ntcp": re.compile(r"\b(?:ntcp|slc\s*10\s*a1)\b",re.I),
    "auc": re.compile(r"(?<![a-z])(?:auc|area\s+under\s+(?:the\s+)?curve)(?![a-z])",re.I),
    "cmax": re.compile(r"(?<![a-z])c[\s_\-]*max(?![a-z])|peak\s+(?:plasma\s+)?concentration",re.I),
    "donor": re.compile(r"\b(?:donor|subject|participant|volunteer|patient)(?:\s*(?:no|number|id|code))?\b",re.I),
    "bodyweight": re.compile(r"(?:body\s*weight|weight\s*\(?kg)",re.I),
}
RELEVANT = ("oatp1b1","oatp1b3","oatp2b1","ntcp","auc","cmax")


def _cell_text(cell, shared):
    t = cell.get("t","")
    v = cell.find(XMLNS+"v")
    if t == "s":
        try: return shared[int(v.text)] if v is not None else ""
        except (IndexError, ValueError, TypeError): return ""
    if t == "inlineStr": return "".join(cell.itertext())
    return v.text if v is not None and v.text else ""


def inspect_xlsx_schema(contents):
    if len(contents)>MAX_WORKBOOK_BYTES or not contents.startswith(b"PK"):
        raise ValueError("Invalid bounded original XLSX")
    summaries=[]
    with ZipFile(BytesIO(contents)) as zipfile:
        shared=[]
        if "xl/sharedStrings.xml" in zipfile.namelist():
            if zipfile.getinfo("xl/sharedStrings.xml").file_size > MAX_WORKBOOK_BYTES:
                raise ValueError("Unsafe large sharedStrings")
            st=ET.fromstring(zipfile.read("xl/sharedStrings.xml"))
            shared=["".join(x.itertext()) for x in st.iter(XMLNS+"si")]
        sheets=sorted(x for x in zipfile.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",x))
        for sheet in sheets:
            if zipfile.getinfo(sheet).file_size > MAX_WORKBOOK_BYTES:
                raise ValueError("Unsafe large sheet")
            xml=ET.fromstring(zipfile.read(sheet))
            row_records=[]
            role_hits=Counter()
            rowwise_numeric_cols=Counter()
            header_candidates=[]
            for row in xml.iter(XMLNS+"row"):
                cells=[]
                for cell in row.iter(XMLNS+"c"):
                    ref=cell.get("r","")
                    col=re.match(r"[A-Z]+",ref)
                    if col:
                        cells.append((col.group(0),_cell_text(cell,shared)))
                row_records.append(cells)
                roles={role for _,txt in cells for role,expr in ROLES.items() if len(txt)<400 and expr.search(txt)}
                role_hits.update(roles)
                if "auc" in roles and "cmax" in roles:
                    header_candidates.append(len(row_records)-1)
            for i in header_candidates:
                headers=row_records[i]
                numeric_by_col=Counter()
                auc_cols=[c for c,s in headers if ROLES["auc"].search(s)]
                cmax_cols=[c for c,s in headers if ROLES["cmax"].search(s)]
                prot_cols=[c for c,s in headers if any(ROLES[k].search(s) for k in RELEVANT[:4])]
                for row in row_records[i+1:]:
                    values=dict(row)
                    for col in set(auc_cols+cmax_cols+prot_cols):
                        try:
                            val=float(values.get(col,""))
                            if val==val and abs(val)<float("inf"):numeric_by_col[col]+=1
                        except (ValueError,OverflowError): pass
                if auc_cols and cmax_cols:
                    sharedPK=max((min(numeric_by_col[a],numeric_by_col[b]) for a in auc_cols for b in cmax_cols if a!=b),default=0)
                    rowwise_numeric_cols["n_numeric_auc_cmax_below_same_headers_max"]=max(rowwise_numeric_cols["n_numeric_auc_cmax_below_same_headers_max"],sharedPK)
                    if prot_cols:
                        rowwise_numeric_cols["n_numeric_protein_below_same_headers_max"]=max(rowwise_numeric_cols["n_numeric_protein_below_same_headers_max"],max(numeric_by_col[c] for c in prot_cols))
            summary={
                "n_excel_rows":len(row_records),
                "n_excel_cells":sum(len(r) for r in row_records),
                "role_keyword_present_anywhere":{role:role_hits[role]>0 for role in ROLES},
                "auc_cmax_joint_header_candidate_rows":len(header_candidates),
                "max_numeric_auc_cmax_subject_records_under_joint_headers":rowwise_numeric_cols["n_numeric_auc_cmax_below_same_headers_max"],
                "max_numeric_protein_records_under_same_headers":rowwise_numeric_cols["n_numeric_protein_below_same_headers_max"],
                "measured_54_person_auc_cmax_and_proteome_matrix_verified":(
                    rowwise_numeric_cols["n_numeric_auc_cmax_below_same_headers_max"]>=54
                    and rowwise_numeric_cols["n_numeric_protein_below_same_headers_max"]>=54
                    and role_hits["donor"]>0
                ),
                "source_cell_values_or_subject_names_exported":False,
            }
            summaries.append(summary)
    return {"original_workbook_SHA256":sha256(contents).hexdigest(),
            "original_sheet_count":len(summaries),
            "source_sheets_schema_only":summaries}


def inspect_supplement_package(contents):
    if len(contents)>MAX_PACKAGE_BYTES or not contents.startswith(b"PK"):
        raise ValueError("Original supplementary source unavailable/not ZIP")
    collected=[]
    with ZipFile(BytesIO(contents)) as container:
        for item in container.infolist():
            if item.is_dir() or not item.filename.lower().endswith(".xlsx"):
                continue
            if item.file_size>MAX_WORKBOOK_BYTES or ".." in item.filename.split("/"):
                raise ValueError("Unsafe supplementary xlsx")
            collected.append(inspect_xlsx_schema(container.read(item)))
    return {
        "study":"Wegler rosuvastatin 54 patients paired measured liver transporter protein and clinical PK",
        "DOI":DOI, "PMCID":PMCID,
        "publisher_original_supplement_SHA256":sha256(contents).hexdigest(),
        "public_supplement_xlsx_count":len(collected),
        "xlsx_original_source_schema_only":collected,
        "any_original_54_person_auc_cmax_proteomics_table_qualified":any(
            s["measured_54_person_auc_cmax_and_proteome_matrix_verified"]
            for w in collected for s in w["source_sheets_schema_only"]
        ),
        "neither_absence_of_PK_from_article_nor_manual_source_proof_claimed":True,
        "no_individual_subject_values_or_keys_printed":True,
    }


def run():
    try:
        with urlopen(Request(PUBLIC_SUPP,headers={"User-Agent":"Omics2PK-original-Wegler54-rosuvastatin-proteomics-PK-qualification/1.0"}),timeout=60) as res:
            contents=res.read(MAX_PACKAGE_BYTES+1)
        if len(contents)>MAX_PACKAGE_BYTES:
            raise ValueError("Oversized source")
        return inspect_supplement_package(contents)
    except (OSError, HTTPError, URLError, BadZipFile,ValueError,ET.ParseError) as err:
        return {"DOI":DOI,"PMCID":PMCID,"source_status":"PUBLIC_ORIGINAL_SUPPLEMENT_UNAVAILABLE_ABSTAIN",
                "error_type":type(err).__name__,
                "any_original_54_person_auc_cmax_proteomics_table_qualified":False}


if __name__=="__main__":
    print(json.dumps(run(),indent=2,sort_keys=True))
