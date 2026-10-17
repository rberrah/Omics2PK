"""Safely test original 5 processed CSV links embedded in author's pinned notebook.

Silveira et al. 2021 DOI 10.3389/fphar.2021.752960, n40 human rosuvastatin
predose LC-MS + genuine targeted PK outcomes. The author README links 5 CSV
datasets in their 21 MB notebook; direct urllib returned Google Drive HTML.
This independent fallback uses gdown's public-file large-download handling.

Reads only original *code* cells, original outputs never inspected or copied.
Prints only counts, semantic column categories and hashes, never opaque cloud
IDs/URLs, clinical rows, intensities, or patient identifiers. Writes temporary
files only in ephemeral GitHub Actions runtime, no artifacts uploaded.
An author-linked CSV alone is NOT proof of a valid human omics/PK pairing.
"""
from __future__ import annotations
from collections import Counter
import csv
from hashlib import sha256
from io import StringIO
import json
import math
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from urllib.parse import urlsplit,parse_qs
from urllib.request import Request,urlopen

COMMIT="dbf9182b0f72f8b290f623f7cc6560c5a5f99d73"
NOTEBOOK=("https://raw.githubusercontent.com/GustavoHBDuarte/"
          f"Metabolomics_PK_proj/{COMMIT}/Metabolomics_Elastic_net.ipynb")
URL_RE=re.compile(r"https://drive\.google\.com/[^\s'\"<>\\)\]}]+",re.I)
MAX_NOTEBOOK_BYTES=28_000_000
MAX_CSV_BYTES=16_000_000
MAX_ROWS=10000
MAX_COLS=35000

def author_code_drive_ids(notebook_bytes):
    if len(notebook_bytes)>MAX_NOTEBOOK_BYTES:raise ValueError("Notebook too large")
    nb=json.loads(notebook_bytes)
    if not isinstance(nb,dict) or not isinstance(nb.get("cells"),list):
        raise ValueError("Invalid original author notebook")
    ids=[]
    for cell in nb["cells"]:
        if cell.get("cell_type")!="code":continue
        source=cell.get("source",[])
        if not isinstance(source,(list,str)):continue
        code="".join(source) if isinstance(source,list) else source
        for x in URL_RE.finditer(code):
            link=x.group(0).rstrip(".,;")
            parsed=urlsplit(link)
            if parsed.scheme!="https" or parsed.hostname!="drive.google.com":continue
            q=parse_qs(parsed.query)
            possible=q.get("id",[""])[0]
            if not possible:
                m=re.search(r"/file/d/([A-Za-z0-9_-]{10,200})",parsed.path)
                possible=m.group(1) if m else ""
            if re.fullmatch(r"[A-Za-z0-9_-]{10,200}",possible):
                ids.append(possible)
    return ids

def classify_csv(binary):
    if len(binary)>MAX_CSV_BYTES or b"\x00" in binary[:200]:
        raise ValueError("Not a bounded original CSV")
    prefix=binary[:512].lstrip().lower()
    if prefix.startswith((b"<html",b"<!doctype")):
        raise ValueError("HTML Google Drive access interstitial")
    content=binary.decode("utf-8-sig")
    initial=content[:8000]
    try:delimiter=csv.Sniffer().sniff(initial,delimiters=",;\t").delimiter
    except csv.Error:delimiter=","
    reader=csv.reader(StringIO(content),delimiter=delimiter)
    header=next(reader,None)
    if not header or len(header)>MAX_COLS:
        raise ValueError("Invalid source header")
    nrows=0
    numeric=0
    nonempty=0
    row_lengths=Counter()
    for row in reader:
        nrows+=1
        if nrows>MAX_ROWS:raise ValueError("Too many original source rows")
        row_lengths[len(row)]+=1
        for val in row:
            v=val.strip()
            if not v:continue
            nonempty+=1
            try:
                if math.isfinite(float(v)):numeric+=1
            except ValueError:pass
    norm=[" ".join(h.casefold().replace("_"," ").split()) for h in header]
    return {
      "n_original_data_rows":nrows,
      "n_original_columns":len(header),
      "n_numeric_cells":numeric,
      "n_nonempty_cells":nonempty,
      "n_consistent_width_rows":row_lengths.get(len(header),0),
      "has_auc_related_column":any(re.search(r"\bauc\b|\barea under\b",h) for h in norm),
      "has_cmax_related_column":any(re.search(r"\bc\s*max\b|\bmaximum concentration\b",h) for h in norm),
      "has_subject_label_column":any(re.search(r"\b(subject|volunteer|participant|patient)\s*(id|code|number|no)\b",h) for h in norm),
      "has_drug_concentration_timepoint_column":any("rosuvastatin" in h and ("conc" in h or "time" in h) for h in norm),
      "original_header_values_exported":False,
      "original_individual_rows_exported":False,
    }

def public_gdown_inspect(nb):
    import gdown  # Optional, only required in the dedicated online source CI
    ids=author_code_drive_ids(nb)
    report={
        "author_original_notebook":"GustavoHBDuarte/Metabolomics_PK_proj",
        "original_notebook_commit":COMMIT,
        "source_links_in_original_code":len(ids),
        "distinct_source_files":len(set(ids)),
        "original_notebook_sha256":sha256(nb).hexdigest(),
        "qualified_same_donor_metabolomics_PK_join":False,
        "output_contains_human_or_google_drive_id":False,
        "results":[],
    }
    with TemporaryDirectory(prefix="omics2pk-author-source-") as temp:
        for index,drive_id in enumerate(ids):
            rec={"source_position_one_based":index+1}
            target=Path(temp)/f"source_{index}.csv"
            try:
                # No logging of ID, URL, source data, exception message.
                out=gdown.download(id=drive_id,output=str(target),quiet=True)
                if not out or not target.is_file():
                    rec["status"]="PUBLIC_GDRIVE_DOWNLOAD_UNAVAILABLE"
                elif target.stat().st_size>MAX_CSV_BYTES:
                    rec["status"]="SOURCE_TOO_LARGE"
                else:
                    data=target.read_bytes()
                    rec["sha256"]=sha256(data).hexdigest()
                    rec["bytes"]=len(data)
                    try:
                        rec["source_dimensions"]=classify_csv(data)
                        rec["status"]="ORIGINAL_AUTHOR_CSV_RETRIEVED"
                    except (ValueError,UnicodeError,csv.Error):
                        rec["status"]="SOURCE_NOT_A_USABLE_CSV"
            except Exception:
                rec["status"]="PUBLIC_GDRIVE_DOWNLOAD_UNAVAILABLE"
            report["results"].append(rec)
    return report

def run():
    with urlopen(Request(NOTEBOOK,headers={"User-Agent":"Omics2PK-rosuvastatin-author-csv-source-audit/1.0"}),timeout=80) as stream:
        notebook=stream.read(MAX_NOTEBOOK_BYTES+1)
    if len(notebook)>MAX_NOTEBOOK_BYTES:raise ValueError("Original notebook too large")
    return public_gdown_inspect(notebook)

if __name__=="__main__":
    try: report=run()
    except Exception as e:
        report={"source_status":"ORIGINAL_AUTHOR_NOTEBOOK_UNAVAILABLE",
                "error_type":type(e).__name__,
                "qualified_same_donor_metabolomics_PK_join":False}
    print(json.dumps(report,indent=2,sort_keys=True))
