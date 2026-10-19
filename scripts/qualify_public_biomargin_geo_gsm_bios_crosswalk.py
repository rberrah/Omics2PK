"""GEO Series Matrix annotation search for exact BIOMARGIN author BIOS keys.

GEO expression columns may be GSM, laboratory alias, or other IDs.
No sample-order matching and no unverified fuzzy alias matching.
Only aggregate counts are emitted. Research source metadata read ephemeral.
"""
from __future__ import annotations
import csv
from gzip import GzipFile
from io import BytesIO, TextIOWrapper
import json
import re

from scripts.qualify_public_biomargin_131_original_source import (
    PublicBIOMARGINAuditError, bounded_download, audit_author_notebook,
    audit_public_geo_counts_gzip, NOTEBOOK, COUNTS, MAX_NOTEBOOK_BYTES,
    MAX_COUNTS_COMPRESSED_BYTES, SAMPLE_PATTERN,
)

SERIES_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175718/matrix/"
    "GSE175718_series_matrix.txt.gz"
)
MAX_MATRIX_COMPRESSED = 8_000_000


def inspect_series_matrix(payload, author_people, count_column_ids, *, production=True):
    if not payload.startswith(b"\x1f\x8b") or len(payload)>MAX_MATRIX_COMPRESSED:
        raise PublicBIOMARGINAuditError("GEO Series Matrix is not bounded original gzip")
    records = {}
    try:
        with GzipFile(fileobj=BytesIO(payload)) as z:
            with TextIOWrapper(z,encoding="utf-8-sig",newline="") as stream:
                for line in stream:
                    if line.startswith("!series_matrix_table_begin"):
                        break
                    if not line.startswith("!Sample_"):
                        continue
                    cells=next(csv.reader([line.rstrip("\r\n")],delimiter="\t"))
                    key=cells[0].lstrip("!")
                    if key not in ("Sample_geo_accession","Sample_title",
                                   "Sample_characteristics_ch1","Sample_description",
                                   "Sample_source_name_ch1"):
                        continue
                    values=[x.strip('"') for x in cells[1:]]
                    if key in records:
                        records[key].append(values)
                    else:
                        records[key]=[values]
    except (OSError,UnicodeError,EOFError,csv.Error) as exc:
        raise PublicBIOMARGINAuditError("Original GEO Series Matrix metadata malformed") from exc
    gsm_lists=records.get("Sample_geo_accession",[])
    if len(gsm_lists)!=1:
        raise PublicBIOMARGINAuditError("Original GEO sample accession metadata missing")
    gsms=gsm_lists[0]
    if len(gsms)!=len(set(gsms)) or (production and len(gsms)!=384):
        raise PublicBIOMARGINAuditError("GEO accession original sample count invalid")
    n=len(gsms)
    for key,lines in records.items():
        if any(len(line)!=n for line in lines):
            raise PublicBIOMARGINAuditError("GEO metadata field columns not aligned to samples")
    lookup={}
    for pos,gsm in enumerate(gsms):
        ids=set()
        for key in ("Sample_title","Sample_characteristics_ch1","Sample_description",
                    "Sample_source_name_ch1"):
            for row in records.get(key,[]):
                ids.update(SAMPLE_PATTERN.findall(row[pos]))
        if len(ids)>1:
            raise PublicBIOMARGINAuditError("GEO analytical sample matches multiple original people")
        lookup[gsm]=next(iter(ids)) if ids else None

    titles=records.get("Sample_title",[[]])[0]
    direct_accessions=set(count_column_ids)&set(gsms)
    direct_titles=set(count_column_ids)&set(titles)
    if len(titles)!=len(set(titles)):
        # Not an error: nonunique GEO titles must never act as a unique join.
        direct_titles=set()
    exact_author_people_from_metadata=set(lookup.values())&set(author_people)
    paired=[gsm for gsm in direct_accessions
            if lookup[gsm] in author_people]
    paired_title=[gsm for gsm,title in zip(gsms,titles)
                  if title in direct_titles and lookup[gsm] in author_people]
    return {
        "status": "GEO_SERIES_MATRIX_METADATA_AGGREGATE_QC",
        "n_source_geo_accessions":n,
        "n_geo_metadata_samples_with_any_BIOS_style_id":sum(v is not None for v in lookup.values()),
        "n_distinct_author_131_BIOS_subjects_detected_in_GEO_metadata":len(exact_author_people_from_metadata),
        "n_GEO_count_columns_exact_GSM_match":len(direct_accessions),
        "n_GEO_count_columns_exact_unique_title_match":len(direct_titles),
        "n_exact_three_way_GSM_count_and_author_BIOS_matches":len(paired),
        "n_exact_three_way_title_count_and_author_BIOS_matches":len(paired_title),
        "identifiable_external_join_without_fuzzy_or_row_order":
            bool(paired or paired_title),
        "source_table_rows_or_GSM_BIOS_mapping_never_published":True,
    }


def peek_counts_columns(payload):
    try:
        with GzipFile(fileobj=BytesIO(payload)) as z:
            with TextIOWrapper(z,encoding="utf-8-sig",newline="") as f:
                header=next(csv.reader(f))
    except (OSError,UnicodeError,csv.Error,StopIteration) as exc:
        raise PublicBIOMARGINAuditError("Cannot inspect bounded GEO source column header") from exc
    if len(header)!=385 or len(header)!=len(set(header)):
        raise PublicBIOMARGINAuditError("Expected 384 original GEO count columns")
    return set(header[1:])


def main():
    notebook=bounded_download(NOTEBOOK,max_bytes=MAX_NOTEBOOK_BYTES)
    authors,author_qc=audit_author_notebook(notebook)
    counts=bounded_download(COUNTS,max_bytes=MAX_COUNTS_COMPRESSED_BYTES)
    counts_qc=audit_public_geo_counts_gzip(counts,authors)
    matrix=bounded_download(SERIES_MATRIX,max_bytes=MAX_MATRIX_COMPRESSED)
    metadata_qc=inspect_series_matrix(matrix,authors,peek_counts_columns(counts))
    result={
        "cohort_author_original":author_qc,
        "geo_original_counts":counts_qc,
        "geo_sample_metadata":metadata_qc,
        "PK_outcomes_and_genotypes_still_not_qualified":True,
        "candidate_join_does_not_substitute_for_authenticated_original_subject_crosswalk":True,
    }
    print(json.dumps(result,indent=2,sort_keys=True))
    return result


if __name__=="__main__":
    main()
