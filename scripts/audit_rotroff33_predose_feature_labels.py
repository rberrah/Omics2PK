#!/usr/bin/env python3
"""Audit original public Rotroff 2016 predose feature annotation coverage.

Only original feature HEADER metadata and aggregate counts are exported,
never a person's PK, molecular intensities or original sample labels.
Scientific names do not prove MSI 1 identification; paper says 155/384
identified. This script verifies whether the publisher's original workbook
actually includes feature→name mapping suitable for biological pathways.
"""
from __future__ import annotations
import json,re
from urllib.request import Request,urlopen
from collections import Counter
from scripts.audit_rotroff33_metformin_original_pk_omics import LIMIT,SUPPLEMENT,_as_float
from omics2pk.rotroff33_clinical_e3 import measured_original_data,actual_original_tables

# No patient labels emitted: first row only plus source feature annotation candidates.
UNRESOLVED=re.compile(r"^\*?\d{4,}$|^unknown|^x-\d|^unidentified",re.I)
TARGETS={
  "urea_cycle":("ornithine","citrulline","arginine","urea"),
  "purine_pools":("inosine","hypoxanthine","xanthine","guanine","adenine","uridine"),
  "glycolysis_redox":("lactate","pyruvate","glutathione","glutamic acid","2-hydroxybutyrate"),
}
def limited_label(x):
    x=str(x or "").strip()
    # Never output person identifiers or sample names even if mislabeled as feature.
    if not x or len(x)>70 or re.fullmatch(r"[A-Za-z]{4}\d{4}-[ABC]",x):
        return None
    return x
def audit_original_labels(zipdata):
    observations,m=measured_original_data(zipdata)
    clinical,omics=actual_original_tables(zipdata)
    mol=omics[0]
    if len(observations)!=33 or len(omics)!=2 or len(mol)<100:
        raise ValueError("Original source original assay/sample structure changed")
    numeric_cols=sorted(j for j in set().union(*(set(r) for r in mol[1:])) if j>0
                        and sum(_as_float(row.get(j)) is not None for row in mol[1:])>=25)
    # Above screening is intentionally broad across all source rows, then
    # compare to the genuine official qualified pre-dose 384-feature count.
    # Some clinical workbook feature labels can be in a separate worksheet.
    original_meta=[]
    for j in numeric_cols:
        original_meta.append((j,limited_label(mol[0].get(j,""))))
    n_unlabeled=sum(not name for j,name in original_meta)
    n_numeric_headers=sum(bool(name and UNRESOLVED.fullmatch(name)) for j,name in original_meta)
    n_named=sum(bool(name and not UNRESOLVED.fullmatch(name)) for j,name in original_meta)
    source_col_count=len(numeric_cols)
    # Source column labels are public chemical annotations not private values.
    # Only print presence/absence of predetermined chemical labels, no raw list.
    matched={}
    for key,terms in TARGETS.items():
        matched[key]={term:sum(bool(name and name.lower()==term.lower()) for j,name in original_meta)
                      for term in terms}
    second=omics[1]
    summary={
      "source":"Rotroff 2016 doi:10.3389/fphar.2016.00135 original supplemental MS workbook",
      "n_original_participants":m["n_original_unique_source_patients"],
      "n_original_predose_qualified_features":m["n_source_predose_molecular_features"],
      "n_candidate_feature_columns_among_all_samples":source_col_count,
      "n_candidate_source_feature_headers_missing":n_unlabeled,
      "n_candidate_feature_headers_numeric_or_unknown":n_numeric_headers,
      "n_candidate_feature_headers_non_numeric_text":n_named,
      "second_omics_sheet_rows":len(second),
      "second_omics_sheet_max_columns":max(map(lambda row:max(row,default=0)+1,second)),
      "second_omics_sheet_numeric_cells_per_row_max":max(sum(_as_float(x) is not None for x in r.values()) for r in second),
      "first_omics_sheet_source_header_non_person_token_examples":[name for _,name in original_meta if name and not UNRESOLVED.fullmatch(name)][:12],
      "literal_exact_name_matches_predeclared_pathways":matched,
      "paper_reports_structurally_identified_count":155,
      "source_chemical_id_confidence_proven_from_this_audit":False,
      "n_original_observation_rows_or_PK_targets_exported":0,
      "model_fitted":False,
    }
    return summary
def main():
    with urlopen(Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK original-source metadata-only metabolite quality audit"}),timeout=90) as response:
        data=response.read(LIMIT+1)
    return audit_original_labels(data)
if __name__=="__main__":
    try:print(json.dumps(main(),sort_keys=True,indent=2))
    except Exception as e:
        print(json.dumps({"status":"ABSTAIN_SOURCE_ANNOTATIONS_UNVERIFIED",
                          "error_class":type(e).__name__},sort_keys=True))
        raise SystemExit(2)
