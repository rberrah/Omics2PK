"""Inspect ORIGINAL public metabolomics papers' supplements for PK-linked patient matrices.

Does NOT publish authors' rows, features, measurements or patient IDs. MDPI
open research supplements are fetched into memory and only aggregate container
metadata emitted. The Phapale Wiley full source is metadata-catalogued as a
candidate but NOT downloaded by unauthorized bypass.

The published studies are fundamentally different: He 2022 renal recipients
focus on TAC efficacy/PD n109; Burghelea 2022 n42 is high/low TAC *groups*;
Phapale 2010 healthy n29 measured PRE-DOSE urine metabolites and 0–72h TAC
serial plasma PK (a high-priority NCA/E3 candidate).
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import io
import json
import re
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile,BadZipFile

CAP=10_000_000
PER_FILE_CAP=12_000_000
SOURCES={
  "Burghelea_2022_Biomedicines": {
    "DOI":"10.3390/biomedicines10051157",
    "article_n_human_outlier_TAC":42,
    "article_total_screened":135,
    "article_omics_assay":"serum lipidomics UHPLC-MS, 336 annotated metabolites",
    "article_pk":"one observed 24-hour TAC trough, HIGH>8 vs LOW<5 ng/mL",
    "original_question":"discrimination of selected TAC high/low clinical group, not continuous NCA or PopPK",
    "supplement_url":"https://pmc.ncbi.nlm.nih.gov/articles/instance/9138871/bin/biomedicines-10-01157-s001.zip",
  },
  "He_2022_Molecules": {
    "DOI":"10.3390/molecules27113517",
    "article_n_human_outlier_TAC":109,
    "article_omics_assay":"LC-MS untargeted metabolomics",
    "article_pk":"Tacrolimus dose and TDM may be in clinical source, verify per patient; PD efficacy was outcome",
    "original_question":"metabolomic predictors of tacrolimus clinical efficacy, not verified serial C(t)",
    "supplement_url":"https://pmc.ncbi.nlm.nih.gov/articles/instance/9182099/bin/molecules-27-03517-s001.zip",
  },
}

def _fetch(url):
    try:
        with urlopen(Request(url,headers={"User-Agent":"Omics2PK-public-open-supplement-scientific-container-QC/1.0"}),timeout=20) as r:
            content=r.read(CAP+1)
            ctype=r.headers.get("Content-Type","").split(";",1)[0][:60]
        if len(content)>CAP:return None,{"status":"oversized_source_supplement"}
        if content[:2]!=b"PK":return None,{"status":"http_response_not_source_zip","content_type":ctype}
        return content,{"status":"ORIGINAL_PUBLIC_SUPPLEMENT_RETRIEVED","source_sha256":sha256(content).hexdigest(),
                        "n_original_compressed_bytes":len(content)}
    except HTTPError as e:return None,{"status":"SOURCE_HTTP_ACCESS_DENIED_OR_UNAVAILABLE","http_status":e.code}
    except (URLError,TimeoutError,OSError) as e:
        return None,{"status":"SOURCE_ENDPOINT_UNAVAILABLE","failure_class":type(e).__name__}

def _xlsx_tables(blob):
    """Only sheet labels + matrix row/column dimensions. Never print cell contents."""
    try:
        with ZipFile(io.BytesIO(blob)) as z:
            parts=z.namelist()
            sheet_files=sorted(x for x in parts if re.fullmatch(r"xl/worksheets/sheet[0-9]+\.xml",x))
            sheets=[]
            ns="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
            for name in sheet_files[:35]:
                if z.getinfo(name).file_size>PER_FILE_CAP:raise ValueError("Source worksheet oversized")
                root=ET.fromstring(z.read(name))
                rows=root.findall(".//"+ns+"sheetData/"+ns+"row")
                maxcol=max((len(r.findall(ns+"c")) for r in rows),default=0)
                sheets.append({"n_original_source_rows":len(rows),"maximum_nonempty_cells_in_any_row":maxcol})
            return {"n_original_worksheets":len(sheet_files),"worksheet_shapes_no_source_values":sheets}
    except (BadZipFile,KeyError,ValueError,ET.ParseError) as e:
        return {"status":"XLSX_SCHEMA_INACCESSIBLE_OR_INVALID","exception_class":type(e).__name__}

def _inspect_archive(b):
    with ZipFile(io.BytesIO(b)) as z:
        members=z.infolist()
        if len(members)>300 or sum(m.file_size for m in members)>250_000_000:
            return {"status":"ORIGINAL_ZIP_UNSAFE_SIZE_ABSTAIN"}
        if any(n.filename.startswith("/") or ".." in n.filename.replace("\\","/").split("/") for n in members):
            return {"status":"ORIGINAL_ZIP_UNSAFE_PATH_ABSTAIN"}
        extensions=Counter(n.filename.rsplit(".",1)[-1].lower() for n in members if not n.is_dir())
        workbooks=[]
        for member in members:
            if member.filename.lower().endswith((".xlsx",".xlsm")):
                if member.file_size>PER_FILE_CAP:
                    workbooks.append({"status":"individual_workbook_too_large"})
                else:
                    workbooks.append(_xlsx_tables(z.read(member.filename)))
        return {
          "status":"ORIGINAL_RESEARCH_ARCHIVE_STRUCTURE_ONLY",
          "n_original_zip_members":len(members),
          "original_file_types_count":dict(extensions),
          "n_contained_spreadsheets":len(workbooks),
          "spreadsheet_matrix_shapes_without_original_values":workbooks[:8],
          "original_patient_identifiers_or_values_exported":False,
          "PK_sample_to_metabolomics_subject_ID_join_verified":False,
        }

def run():
    summary={}
    for study,meta in SOURCES.items():
        blob,fetch=_fetch(meta["supplement_url"])
        result={"article_DOI":meta["DOI"],"article_source_design":meta["original_question"],
                "article_PK_information":meta["article_pk"],"article_omics_assay":meta["article_omics_assay"],
                "original_supplement_archive_url":meta["supplement_url"],
                "supplement_retrieval":fetch,
                "individual_NCA_AUC_Cmax_Tmax_dataset_qualified":False}
        if blob:
            try:result["original_supplement_structure"]=_inspect_archive(blob)
            except (BadZipFile,ValueError) as e:
                result["original_supplement_structure"]={"status":"INVALID_ARCHIVE_ABSTAIN","error_type":type(e).__name__}
        summary[study]=result
    summary["Phapale_2010_high_priority"]= {
      "DOI":"10.1038/clpt.2009.296",
      "original_n_human_volunteers":29,
      "predose_urine_metabolite_ions":1256,
      "original_serial_plasma_TAC_horizon_h":72,
      "source_supplement_count_reported":13,
      "original_matching_biological_unit":"same volunteers, predose omic to postdose PK",
      "individual_omics_x_PK_source_matrix_retrieved":False,
      "numeric_AUC_Cmax_Tmax_on_identical_29_subjects_qualified":False,
      "priority":"HIGHEST_NATIVE_NCA_E3_SOURCE_CANDIDATE",
      "required":"original 13 Wiley supplementary .doc + deidentified urine ion intensity matrix 29×1256 + PK doses, hours and individual concentrations and identifier crosswalk",
    }
    return {"status":"PUBLIC_PRIMARY_SOURCE_SUPPLEMENT_TRIAGE_NO_FABRICATED_PAIRED_PK",
            "sources":summary,"original_patient_or_metabolite_values_exported":False}

if __name__=="__main__":print(json.dumps(run(),ensure_ascii=False,indent=2,sort_keys=True))
