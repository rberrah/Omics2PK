"""Privacy-safe structural crosswalk for 2021 US EPA n=100 donor RNA/enzymatic activity.

The original xlsx is downloaded in memory from the exact EPA DOI and never
committed or uploaded. We report ONLY counts / column names, not donor IDs or
raw measurements. Patient/subject-level identifiers remain local.
"""
from __future__ import annotations

from collections import Counter
from io import BytesIO
import json
import re
from urllib.error import HTTPError,URLError
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from scripts.audit_epa_liver_rna_activity import _get, MS, RELPKG, RELS


def _column_index(ref):
    letters=re.match(r"^[A-Z]+",ref)
    if letters is None: raise ValueError("Malformed spreadsheet cell reference")
    value=0
    for c in letters.group(0):
        value=value*26+ord(c)-64
    return value-1


def parse_epa_workbook(data):
    result={}
    with ZipFile(BytesIO(data)) as z:
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        root_rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        paths={n.attrib["Id"]:n.attrib.get("Target","") for n in root_rel.findall(RELPKG+"Relationship")}
        strings=[]
        if "xl/sharedStrings.xml" in z.namelist():
            sr=ET.fromstring(z.read("xl/sharedStrings.xml"))
            strings=["".join(t.itertext()) for t in sr.findall(MS+"si")]
        for s in wb.iter(MS+"sheet"):
            target=paths.get(s.attrib.get(RELS+"id"),"").lstrip("/")
            if not target.startswith("xl/"):target="xl/"+target
            tree=ET.fromstring(z.read(target))
            rows=[]
            for r in tree.findall(".//"+MS+"sheetData/"+MS+"row"):
                fields={}
                for c in r.findall(MS+"c"):
                    index=_column_index(c.attrib.get("r",""))
                    if index>250:raise ValueError("Unexpectedly wide source")
                    v=c.find(MS+"v")
                    if v is None:
                        ins=c.find(MS+"is")
                        value="".join(ins.itertext()) if ins is not None else ""
                    elif c.attrib.get("t")=="s":
                        value=strings[int(v.text)]
                    elif c.attrib.get("t")=="b":
                        value=bool(int(v.text))
                    elif c.attrib.get("t") in ("str","inlineStr"):
                        value=v.text or ""
                    else:
                        try:value=float(v.text)
                        except (ValueError,TypeError):value=None
                    fields[index]=value
                if fields:
                    rows.append([fields.get(j) for j in range(max(fields)+1)])
            if not rows:raise ValueError("Empty source worksheet")
            header=[str(x).strip() if x is not None else "" for x in rows[0]]
            result[s.attrib["name"]]=(header,rows[1:])
    return result


def _col_indices(header):
    lowered=[re.sub(r"\s+"," ",x).casefold() for x in header]
    K=[i for i,x in enumerate(lowered) if "kumc sample id" in x]
    X=[i for i,x in enumerate(lowered) if "xenotech" in x and "id" in x]
    cyps={x:i for i,x in enumerate(header) if re.match(r"^CYP[0-9A-Z*]+$",x,re.I)}
    return K,X,cyps


def source_crosswalk_report(wb):
    sheets={}
    keysets={}
    for name,(header,rows) in wb.items():
        K,X,cyps=_col_indices(header)
        keys=[]
        n_numeric_cyp=0
        n_rows_nonempty_id=0
        for row in rows:
            key=[]
            for col in (K[:1]+X[:1]):
                v=row[col] if col<len(row) else None
                if v is not None and str(v).strip():key.append(str(v).strip())
            if key:
                n_rows_nonempty_id+=1
                keys.append(tuple(key))
            if any(col<len(row) and isinstance(row[col],(int,float)) for col in cyps.values()):
                n_numeric_cyp+=1
        keysets[name]=set(keys)
        sheets[name]={
            "header_columns":header,
            "n_data_rows":len(rows),
            "n_unique_sample_key_tuples":len(set(keys)),
            "n_rows_with_sample_key":n_rows_nonempty_id,
            "n_sample_key_duplicates":len(keys)-len(set(keys)),
            "n_rows_with_numeric_CYP_field":n_numeric_cyp,
            "has_kumc_id":bool(K),
            "has_xenotech_id":bool(X),
            "cyp_feature_columns":list(cyps),
        }
    overlap={}
    if "Figure 2" in keysets and "Table 3" in keysets:
        f,t=keysets["Figure 2"],keysets["Table 3"]
        overlap={
            "n_unique_matching_exact_composite_keys":len(f&t),
            "n_fig2_keys_not_table3":len(f-t),
            "n_table3_keys_not_fig2":len(t-f),
        }
    return {"source":"US EPA 10.23719/1519197", "sheets":sheets,
            "source_identity_join":overlap,
            "individual_ids_or_expression_exported":False,
            "qualitative_status":"STRUCTURE_AND_KEY_OVERLAP_NOT_FUNCTIONAL_ENDPOINT_VALIDATION"}


if __name__=="__main__":
    try:
        result=source_crosswalk_report(parse_epa_workbook(_get()))
        print(json.dumps(result,sort_keys=True,indent=2))
    except (HTTPError,URLError,TimeoutError,OSError,ValueError) as e:
        print(json.dumps({"status":"SOURCE_UNAVAILABLE_OR_MALFORMED",
                          "reason_type":type(e).__name__},sort_keys=True))
        raise
