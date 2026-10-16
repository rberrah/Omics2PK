"""Audit original public GEO longitudinal human acetaminophen RNA sample metadata.

Downloads ONLY official NCBI series-matrix metadata up to the beginning
of the numeric probe table. Not full omics expression or PK concentrations.
Never prints source participant IDs, sample accession strings, or values.
"""
from __future__ import annotations

import gzip
import json
import re
from collections import Counter
from urllib.request import Request, urlopen

SERIES={
    "GSE15158":{
        "study":"Fannin 2010 single 4g APAP / placebo, 9 volunteers",
        "publication_doi":"10.1002/hep.23330",
        "expected_rna_after_dosing_h":[6,18,24,48,72,96],
        "actual_patient_pk_source":"published_figure_2_concentration_curves_not_original_tabular_PK_verified"
    },
    "GSE70784":{
        "study":"Bushel/Fannin 2016 APAP 4g/day for 7 days vs placebo",
        "publication_doi":"10.1038/tpj.2016.8",
        "expected_rna_after_dosing_h":None,
        "actual_patient_pk_source":"no_confirmed_individual_original_drug_concentration_source"
    },
}
MAX_METADATA_LINES=600
MAX_METADATA_UNCOMPRESSED_CHARS=4_000_000


def audit_geo(series):
    if series not in SERIES:
        raise ValueError("Only explicitly researched NCBI GEO Series IDs are allowed")
    folder=series[:-3]+"nnn"
    url=f"https://ftp.ncbi.nlm.nih.gov/geo/series/{folder}/{series}/matrix/{series}_series_matrix.txt.gz"
    lines=[]
    n_bytes=0
    headers=Counter()
    with urlopen(Request(url,headers={"User-Agent":"Omics2PK-source-qualification-geo/1.0"}),timeout=60) as f:
        with gzip.GzipFile(fileobj=f) as inp:
            for idx,raw in enumerate(inp):
                if idx > MAX_METADATA_LINES:
                    raise ValueError("GEO series matrix without bounded metadata header")
                n_bytes+=len(raw)
                if n_bytes>MAX_METADATA_UNCOMPRESSED_CHARS:
                    raise ValueError("Oversized GEO sample metadata header")
                s=raw.decode("utf-8-sig",errors="replace").rstrip("\n\r")
                if s.startswith("!series_matrix_table_begin"):break
                if not s.startswith("!"):continue
                key=s.split("\t",1)[0]
                headers[key]+=1
                # Only Series metadata and SAMPLE COUNT are preserved, never the
                # array of raw Sample IDs or subject characteristics.
                if key in ("!Series_title","!Series_summary"):
                    lines.append(s[:200])
            else:
                raise ValueError("GEO series matrix did not expose a probe-table marker")
    output={
        "geo_accession":series,
        "geo_matrix_url":url,
        "study":SERIES[series]["study"],
        "publication_doi":SERIES[series]["publication_doi"],
        "metadata_rows_read_before_expression":idx+1,
        "metadata_first_section_n_bytes":n_bytes,
        "n_source_sample_columns_declared":None,
        "n_sample_characteristic_fields":sum(n for name,n in headers.items() if name.startswith("!Sample_characteristics")),
        "has_sample_geo_accession_array":"!Sample_geo_accession" in headers,
        "has_sample_title_array":"!Sample_title" in headers,
        "has_processed_GEO_expression_matrix_table":True,
        "original_blood_RNA_at_serial_times":True,
        "study_SOP_source_access":"public_GEO_not_PHI_export",
        "matched_dose_and_true_individual_concentrations":False,
        "matched_numeric_human_PK_NCA_endpoints":False,
        "requires_original_patient_source_PK_and_temporal_crosswalk":True,
        "source_claim_scope":"PUBLIC_LONGITUDINAL_RNA_METADATA_NOT_PK_PREDICTION",
        "all_individual_GEO_identifiers_and_expression_values_redacted":True,
    }
    return output


if __name__=="__main__":
    print(json.dumps({accession:audit_geo(accession) for accession in SERIES},
                    indent=2,sort_keys=True))
