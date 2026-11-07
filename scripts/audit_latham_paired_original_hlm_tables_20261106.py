"""Original Latham 2025 donor-matched table structural audit (two actual drug studies).

Original *primary* openly accessible full-text XML tables via Europe PMC:
- Imatinib PMCID PMC12990092, DOI 10.1016/j.dmd.2025.100196
- Tucatinib PMCID PMC13095627, DOI 10.1016/j.dmd.2025.100061

Audit only donor-ID set intersection, column schema, source row counts,
clearance measure presence, source page structural provenance. Do NOT print
per-donor genetic, protein, activity, CLint or demographic measurements.
No statistical model on unvalidated original pair tables.
"""
from __future__ import annotations
import json
from collections import Counter
import xml.etree.ElementTree as ET
from urllib.error import HTTPError,URLError
from urllib.request import urlopen,Request

ORIGINAL_PMC={
 "imatinib":{"pmcid":"PMC12990092","doi":"10.1016/j.dmd.2025.100196"},
 "tucatinib":{"pmcid":"PMC13095627","doi":"10.1016/j.dmd.2025.100061"},
}
SOURCE_LIMIT=4_000_000

def _get_original_jats(pmcid):
    # EBI Europe PMC serves original NLM/JATS XML, not an extracted language
    # model answer. Accept JATS only and verify PMCID in article-id.
    url=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    req=Request(url,headers={
       "Accept":"application/xml,text/xml",
       "User-Agent":"Omics2PK-original-academic-reproducibility/1.0"
    })
    with urlopen(req,timeout=65) as response:
        payload=response.read(SOURCE_LIMIT+1)
    if len(payload)>SOURCE_LIMIT:
        raise ValueError("Original article full text XML oversized")
    root=ET.fromstring(payload)
    if root.tag!="article":
        raise ValueError("Original fulltext source is not JATS article XML")
    identifiers={x.attrib.get("pub-id-type",""):((x.text or "").strip())
                 for x in root.iter("article-id")}
    if identifiers.get("pmc")!=pmcid.removeprefix("PMC"):
        raise ValueError("Original source JATS PMCID mismatch")
    return root

def _text(element):
    return " ".join(" ".join(element.itertext()).split())

def _tables(root):
    records=[]
    for tw in root.iter("table-wrap"):
        id=tw.attrib.get("id","")
        title=_text(tw.find("caption")) if tw.find("caption") is not None else ""
        table=tw.find("table")
        if table is None:
            # An inline alternate external graphic does not establish numeric
            # same-donor data; do not invent a cohort from figure text.
            continue
        trs=table.findall(".//tr")
        rows=[]
        for tr in trs:
            cells=[_text(cell) for cell in tr if cell.tag in ("td","th")]
            if cells:rows.append(cells)
        records.append({
            "original_source_table_id":id,
            "table_caption_excerpt":title[:250],
            "n_numeric_table_rows_including_header":len(rows),
            "max_original_source_columns":max([len(r) for r in rows],default=0),
            "_source_rows_private":rows,
        })
    return records

def _donor_table(tables,drug):
    candidates=[t for t in tables
        if "donor" in t["table_caption_excerpt"].lower()
        and t["max_original_source_columns"]>=10
        and t["n_numeric_table_rows_including_header"]>=21]
    if len(candidates)!=1:
        raise ValueError(f"Expected one original {drug} human liver donor Table1")
    table=candidates[0]
    rows=table["_source_rows_private"]
    data=[r for r in rows if len(r)>=10 and r[0].strip().isdigit()]
    if len(data)!=21:
        raise ValueError(f"Expected exactly 21 source individual {drug} liver donors")
    donors=[r[1].strip() for r in data]
    if len(set(donors))!=21:
        raise ValueError("Duplicated source donor ID in original JATS table")
    for r in data:
        if not r[1].strip() or not r[0].isdigit():
            raise ValueError("Unqualified original donor source row")
    return table,donors,data

def run():
    sources={}
    origin={}
    for name,p in ORIGINAL_PMC.items():
        try:
            root=_get_original_jats(p["pmcid"])
            tables=_tables(root)
            original_table,donors,data=_donor_table(tables,name)
            # 21 donor IDs retained only in ephemeral memory to compute source
            # crosswalk, never printed, saved or committed.
            origin[name]=set(donors)
            missing=0
            outcome=[r[-2 if name=="imatinib" else -1] for r in data]
            for v in outcome:
                if v.strip() in ("","-","–","NA","N/A"):missing+=1
            sources[name]={
                "status":"ORIGINAL_PRIMARY_JATS_DONOR_NUMERIC_TABLE_QUALIFIED",
                "PMCID":p["pmcid"],"DOI":p["doi"],
                "source_original_table_ID":original_table["original_source_table_id"],
                "source_caption":original_table["table_caption_excerpt"],
                "n_true_original_individual_liver_donors":21,
                "n_original_table_columns":original_table["max_original_source_columns"],
                "n_observed_last_or_penultimate_outcome_col_nonmissing":21-missing,
                "n_absent_or_ND_outcome_cells":missing,
                "individual_donor_values_or_source_IDs_not_exported":True,
                "no_genotype_claim_for_HLM_source":True,
            }
        except (HTTPError,URLError,ET.ParseError,OSError,ValueError) as err:
            sources[name]={
                "status":"ABSTAIN_ORIGINAL_NUMERIC_JATS_SOURCE_NOT_QUALIFIED",
                "PMCID":p["pmcid"],
                "exception_class":type(err).__name__,
                "no_source_row_or_IDs_printed":True,
            }
    if len(origin)==2:
        shared=origin["imatinib"]&origin["tucatinib"]
        overlap={
          "status":"ORIGINAL_INDIVIDUAL_DONOR_ID_CROSSWALK_CHECKED",
          "n_samematched_original_donor_IDs":len(shared),
          "n_original_imatinib_donors":len(origin["imatinib"]),
          "n_original_tucatinib_donors":len(origin["tucatinib"]),
          "these_are_not_two_independent_biological_cohorts":len(shared)>0,
          "donor_IDs_not_distributed":True,
        }
    else:overlap={
        "status":"ABSTAIN_PAIR_CROSSWALK_UNTIL_BOTH_PRIMARY_TABLES_ACQUIRED",
        "donor_IDs_not_distributed":True,
    }
    return {
      "status":"ORIGINAL_LATHAM_PAIRED_HLM_SOURCE_ACQUISITION_GATE",
      "same_donor_drug_exposures_are_NOT_independent_validation":True,
      "these_are_measured_IN_VITRO_metabolism_not_human_patient_PK":True,
      "original_public_source_tables":sources,
      "cross_study_same_donor":overlap,
    }

if __name__=="__main__":
    print(json.dumps(run(),sort_keys=True,indent=2))
