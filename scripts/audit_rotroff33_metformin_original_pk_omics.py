"""Source-only original PMC4906013 Rotroff 2016 human metformin E3 qualification.

Original study n=33, 384 GC-TOF metabolites at baseline and post-dose;
metformin concentration PK AUC/Cmax measured in the manuscript.

Read ALL official Europe PMC supplementary Excel workbooks transiently.
An actual 33-donor E3 dataset requires same-person source measurements:
one original baseline metabolomic row with >=20 numerical features AND
measured original AUC+Cmax, joined through a verified donor ID.

Only original aggregate dimensions, DOI and content digests emitted; never
source subject identifiers, metabolomic intensities or clinical PK values.
If source has only processed aggregate associations or independent unjoined
sheets, ABSTAIN instead of fitting artificial patient data.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
import re
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile,BadZipFile

STUDY_DOI="10.3389/fphar.2016.00135"
PMCID="PMC4906013"
SUPPLEMENT=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
LIMIT=24_000_000
MEMBER_LIMIT=10_000_000
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
PK={
    "metformin_auc":re.compile(r"(?<![a-z])(?:auc|area\s+under\s+(?:the\s+)?curve)(?![a-z])",re.I),
    "metformin_cmax":re.compile(r"(?<![a-z])c[\s_\-]*max(?![a-z])|peak\s+(?:plasma\s+)?concentration",re.I),
    "metformin_concentration":re.compile(r"(?:metformin|met)\s*(?:plasma\s*)?(?:conc|concentration)",re.I),
    "source_subject_id":re.compile(r"^(?:subject|participant|donor|patient|volunteer|individual)(?:\s*(?:id|number|code|no))?$",re.I),
    "baseline":re.compile(r"\b(?:predose|pre-dose|baseline|timepoint\s*A|time\s*point\s*A|pre\s*metformin)\b",re.I),
}
MAX_ROWS=6000
MAX_COLS=1500


def _cell_str(c,shared):
    typ=c.attrib.get("t","")
    v=c.find(NS+"v")
    if typ=="s":
        try: return shared[int(v.text)] if v is not None else ""
        except (IndexError,ValueError,TypeError):return ""
    if typ=="inlineStr":return "".join(c.itertext())
    return v.text if v is not None and v.text else ""


def _col_index(ref):
    token=re.match(r"([A-Z]+)",ref)
    if not token:return None
    result=0
    for letter in token.group(1):result=result*26+(ord(letter)-64)
    return result-1


def _as_float(value):
    try:
        a=float(value)
        return a if a==a and abs(a)<float("inf") else None
    except (ValueError,TypeError,OverflowError):return None


def analyze_rows(rows):
    if len(rows)>MAX_ROWS:
        raise ValueError("Original XLSX contains excessive rows")
    head=rows[:min(len(rows),12)]
    role_occurrence={role:sum(
        1 for row in head for txt in row.values()
        if len(txt)<220 and pattern.search(txt)
    ) for role,pattern in PK.items()}
    any_pk=role_occurrence["metformin_auc"] and role_occurrence["metformin_cmax"]
    candidate_pk=0
    candidate_highdim=0
    candidate_baseline=0
    subject_header=0
    for row_i,h in enumerate(head):
        idcols=[i for i,v in h.items() if PK["source_subject_id"].fullmatch(v.strip())]
        auc=[i for i,v in h.items() if PK["metformin_auc"].search(v)]
        cmax=[i for i,v in h.items() if PK["metformin_cmax"].search(v)]
        if idcols:subject_header+=1
        if not idcols:continue
        numeric_count=Counter()
        paired=0
        ids=set()
        for row in rows[row_i+1:]:
            sid=row.get(idcols[0],"").strip()
            if not sid or len(sid)>60:continue
            ids.add(sid)
            numeric=sum(_as_float(v) is not None for c,v in row.items() if c not in idcols)
            if numeric>=20:paired+=1
            for c in set(auc+cmax):
                if _as_float(row.get(c)) is not None:numeric_count[c]+=1
        if len(ids)>=30 and paired>=30:
            candidate_highdim=max(candidate_highdim,paired)
            if role_occurrence["baseline"]:
                candidate_baseline=max(candidate_baseline,paired)
        if len(ids)>=30 and auc and cmax:
            best=max((min(numeric_count[a],numeric_count[b]) for a in auc for b in cmax if a!=b),default=0)
            candidate_pk=max(candidate_pk,best)
    return {
        "n_original_worksheet_rows":len(rows),
        "n_original_worksheet_cells":sum(len(r) for r in rows),
        "max_numeric_values_in_one_row":max((sum(_as_float(v) is not None for v in r.values()) for r in rows),default=0),
        "role_keyword_present_in_top_12_rows":{k:bool(v) for k,v in role_occurrence.items()},
        "n_early_rows_with_subject_header":subject_header,
        "max_subjects_with_numeric_auc_and_cmax":candidate_pk,
        "max_subjects_with_20_plus_numeric_features":candidate_highdim,
        "max_subjects_with_20_plus_features_and_baseline_keyword":candidate_baseline,
        "self_contained_33_person_baseline_omics_and_pk_table":
            candidate_pk>=33 and candidate_baseline>=33,
        "source_identifiers_headers_and_values_exported":False,
    }


def source_xlsx_schema(data):
    if len(data)>MEMBER_LIMIT or not data.startswith(b"PK"):
        raise ValueError("Invalid original bounded supplementary XLSX")
    worksheets=[]
    with ZipFile(BytesIO(data)) as z:
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            if z.getinfo("xl/sharedStrings.xml").file_size>MEMBER_LIMIT:
                raise ValueError("Unsafe large original strings")
            st=ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared=["".join(x.itertext()) for x in st.iter(NS+"si")]
        for name in sorted(x for x in z.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",x)):
            if z.getinfo(name).file_size>MEMBER_LIMIT:
                raise ValueError("Unsafe large original worksheet")
            root=ET.fromstring(z.read(name))
            rows=[]
            for node in root.iter(NS+"row"):
                values={}
                for cell in node.iter(NS+"c"):
                    i=_col_index(cell.attrib.get("r",""))
                    if i is None or i>=MAX_COLS:
                        continue
                    values[i]=_cell_str(cell,shared)
                rows.append(values)
            worksheets.append(analyze_rows(rows))
    return {"original_file_sha256":sha256(data).hexdigest(),
            "n_original_worksheets":len(worksheets),
            "worksheets_aggregate_only":worksheets}


def audit_package(package):
    if len(package)>LIMIT or not package.startswith(b"PK"):
        raise ValueError("Original publisher package not a bounded ZIP")
    workbooks=[]
    extensions=Counter()
    with ZipFile(BytesIO(package)) as z:
        for item in z.infolist():
            if item.is_dir():continue
            if ".." in item.filename.split("/"):raise ValueError("Archive path traversal")
            ext="."+item.filename.rsplit(".",1)[-1].lower() if "." in item.filename else "unknown"
            extensions[ext]+=1
            if ext==".xlsx":
                if item.file_size>MEMBER_LIMIT:raise ValueError("Original data source too large")
                workbooks.append(source_xlsx_schema(z.read(item)))
    return {
        "primary_study_DOI":STUDY_DOI,
        "PMCID":PMCID,
        "article_reports_original_individual_patients":33,
        "article_reports_384_original_metabolites":True,
        "article_reports_measured_metformin_AUC_Cmax":True,
        "original_publisher_archive_sha256":sha256(package).hexdigest(),
        "original_supplement_extension_counts":dict(sorted(extensions.items())),
        "n_original_workbooks":len(workbooks),
        "workbook_schemas_no_source_values":workbooks,
        "n_source_workbooks_containing_selfcontained_33_patient_E3_table":sum(
            any(sh["self_contained_33_person_baseline_omics_and_pk_table"] for sh in w["worksheets_aggregate_only"])
            for w in workbooks),
        "actual_baseline_omics_pk_patient_crosswalk_qualified":False,
        "caveat":"Even positive structural screen requires scientific manual original row provenance/time/units validation; no cross-file arbitrary join",
        "no_clinical_subjects_or_original_PGx_omics_exported":True,
    }


def main():
    try:
        with urlopen(Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK-Rotroff33-original-supplement-crosswalk-audit/1.0"}),timeout=55) as stream:
            raw=stream.read(LIMIT+1)
        return audit_package(raw)
    except (OSError,ValueError,BadZipFile,ET.ParseError) as err:
        return {"source_status":"ORIGINAL_PUBLISHER_SUPPLEMENT_NOT_ACCESSIBLE_ABSTAIN",
                "error_type":type(err).__name__,
                "actual_baseline_omics_pk_patient_crosswalk_qualified":False}


if __name__=="__main__":
    print(json.dumps(main(),indent=2,sort_keys=True))
