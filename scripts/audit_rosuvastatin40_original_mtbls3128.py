"""Original rosuvastatin n40 predose LC-MS + quantitative AUC/Cmax source.

Read official MetaboLights MTBLS3128 study metadata and PMC8633954 supplements.
Printed output: source hashes and aggregate file/header structure only, never
observed volunteer IDs, metabolomic intensities or plasma PK values.
Finding manuscript-linked original mzMLs and a reference PK table is evidence
of data availability, NOT an exact join or a reproduced supervised model.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
import re
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen
from zipfile import BadZipFile,ZipFile
from xml.etree import ElementTree as ET

STUDY="MTBLS3128"
PMCID="PMC8633954"
DOI="10.3389/fphar.2021.752960"
SOURCES={
  "metabolights_official_metadata":f"https://www.ebi.ac.uk/metabolights/ws/studies/{STUDY}",
  "metabolights_official_file_inventory":f"https://www.ebi.ac.uk/metabolights/ws/studies/{STUDY}/files",
  "europepmc_official_publisher_supplements":f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles",
}
MAX=20_000_000
SS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

def original(url):
    with urlopen(Request(url,headers={"User-Agent":"Omics2PK-academic-MTBLS3128-original-source-audit/1.0",
         "Accept":"application/json,application/zip,application/octet-stream,*/*;q=.5"}),
         timeout=40) as conn:
        data=conn.read(MAX+1)
        if len(data)>MAX:
            raise ValueError("Original small-source metadata exceeds explicit size limit")
    return data

def json_structure(data):
    obj=json.loads(data)
    if not isinstance(obj,(dict,list)):
        raise ValueError("Invalid original MetaboLights study response")
    labels=[]
    lengths=Counter()
    def walk(node,depth=0):
        if depth>24:return
        if isinstance(node,dict):
            for key,value in node.items():
                lk=key.lower()
                if lk in ("studyidentifier","accession","studyid") and isinstance(value,str):
                    labels.append(value)
                if lk in ("filename","name","path") and isinstance(value,str):
                    ext=value.rsplit(".",1)[-1].lower() if "." in value else "no-ext"
                    lengths[ext]+=1
                if lk in ("assays","samples","sources","metaboliteassignments","datafiles","files",
                         "studyassays","studyfiles") and isinstance(value,list):
                    lengths["key_"+lk+"_entries"]+=len(value)
                walk(value,depth+1)
        elif isinstance(node,list):
            for val in node[:6000]:walk(val,depth+1)
    walk(obj)
    if labels and STUDY not in labels:
        raise ValueError("MetaboLights metadata identity mismatch")
    return {
       "root_type":type(obj).__name__,
       "study_ID_checked":STUDY in labels,
       "root_keys_schema_only":sorted(obj.keys())[:35] if isinstance(obj,dict) else [],
       "file_extensions_or_metadata_collection_counts":dict(sorted(lengths.items())),
       "no_human_subject_rows_exported":True,
    }

def sheet_structure(xlsx):
    # Inspect only worksheet dimensions and presence of PK-specific COLUMN LABELS.
    with ZipFile(BytesIO(xlsx)) as wb:
        ns=wb.namelist()
        sheets=[]
        names=[]
        if "xl/workbook.xml" not in ns:
            raise ValueError("No Excel workbook metadata")
        names=[n.attrib.get("name","")[:80] for n in ET.fromstring(wb.read("xl/workbook.xml")).iter(SS+"sheet")]
        shared=[]
        if "xl/sharedStrings.xml" in ns:
            root=ET.fromstring(wb.read("xl/sharedStrings.xml"))
            shared=[" ".join(e.itertext())[:100] for e in root.iter(SS+"si")]
        for path in sorted(x for x in ns if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",x)):
            root=ET.fromstring(wb.read(path))
            rows=list(root.iter(SS+"row"))
            cell_count=sum(1 for _ in root.iter(SS+"c"))
            header_terms=[]
            # Read ONLY column labels of FIRST row; immediately map to
            # coarse YES/NO categories, do not export labels/subject IDs.
            if rows:
                for c in rows[0].iter(SS+"c"):
                    t=c.attrib.get("t")
                    v=c.find(SS+"v")
                    if t=="s" and v is not None:
                        try:s=shared[int(v.text or "-1")]
                        except (ValueError,IndexError):s=""
                    elif t=="inlineStr":s=" ".join(c.itertext())
                    else:s=""
                    header_terms.append(s)
            txt=" ".join(header_terms).lower()
            sheets.append({
               "n_rows_including_header":len(rows),
               "n_cells":cell_count,
               "header_has_auc":bool(re.search(r"\bauc\b",txt)),
               "header_has_cmax":bool(re.search(r"c[\s_]*max",txt)),
               "header_has_subject_key":bool(re.search(r"(?:subject|volunteer|participant|patient|sample[\s_-]*id)",txt)),
            })
    return {"sheet_names":names,"n_sheets":len(names),"sheet_structure":sheets}

def zip_structure(data,depth=0):
    if len(data)>MAX or depth>2 or not data.startswith(b"PK"):
        raise ValueError("Invalid original ZIP archive")
    count=Counter()
    files=[]
    with ZipFile(BytesIO(data)) as z:
        for f in z.infolist():
            if f.is_dir():continue
            if len(files)>=120 or f.file_size>MAX or ".." in f.filename.split("/"):
                raise ValueError("Unsafe excessive original archive")
            fname=f.filename.rsplit("/",1)[-1]
            kind=fname.rsplit(".",1)[-1].lower() if "." in fname else "unknown"
            count[kind]+=1
            item={"basename":fname[:140],"size":f.file_size}
            if kind=="zip":
                item["nested_archive"]=zip_structure(z.read(f),depth+1)
            elif kind=="xlsx":
                item["sha256"]=sha256(z.read(f)).hexdigest()
                item["metadata"]=sheet_structure(z.read(f))
            files.append(item)
    return {"file_type_counts":dict(count),"original_file_structure":files}

def audit():
    r={
      "study_id":STUDY,
      "article_doi":DOI,
      "pmcid":PMCID,
      "published_n_subjects":40,
      "published_drug":"rosuvastatin",
      "published_biomarkers":"predose untargeted LC-MS serum + independently calibrated plasma PK AUC and Cmax",
      "true_individual_PK_omics_crosswalk_verified":False,
      "external_holdout_or_nested_cv_reproduced":False,
      "attempts":{},
    }
    for kind,url in SOURCES.items():
        try:
            raw=original(url)
            info=zip_structure(raw) if kind.startswith("europepmc") else json_structure(raw)
            r["attempts"][kind]={"access":"ORIGINAL_SOURCE_RETRIEVED",
                "sha256":sha256(raw).hexdigest(),"source_bytes":len(raw),
                "structure":info}
        except (HTTPError,URLError,TimeoutError,OSError,ValueError,BadZipFile,ET.ParseError,RecursionError) as err:
            r["attempts"][kind]={"access":"SOURCE_UNAVAILABLE_ABSTAIN",
                                  "error_type":type(err).__name__}
    return r

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True,ensure_ascii=False))
