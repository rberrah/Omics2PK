"""Verify REAL GSE9588 427-human-liver array matrix, never export raw expression.

Only GEO's official accession-pinned processed series matrix. This is an
input assay gate, **not** transcriptome→measured CYP activity validation:
the independently measured enzyme activity original phenotype table and
exact GEO GSM↔biological donor crosswalk are NOT available here. Never
predict enzyme clearance from a different study or invent activity labels.

Read the compressed matrix as a STREAM (never load 427×feature rows into
memory). Original actual donor IDs/FPKM/probe values are not printed, saved
or checked into this repository. GEO platform probes are microarray arrays,
not RNA-seq counts, and raw expression processing preceded this analysis.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
from io import BytesIO
import gzip
import json
import re
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9588/matrix/GSE9588_series_matrix.txt.gz"
MAX_UNCOMPRESSED=370_000_000
MAX_PROBES=130_000
N_EXPECTED=427


def audit_gse_matrix_stream(gzip_stream,*,source_label="GEO source"):
    sample_accessions=[]
    in_table=False
    n_probes=0
    n_non_numeric=0
    n_wrong_width=0
    n_lines=0
    n_bytes=0
    platforms=set()
    expression_headers_seen=False
    for raw in gzip_stream:
        n_bytes+=len(raw)
        n_lines+=1
        if n_bytes>MAX_UNCOMPRESSED or n_lines>MAX_PROBES+3000:
            raise ValueError("Original GEO decompressed matrix exceeds bounded source limits")
        if not in_table and raw.startswith(b"!Series_platform_id"):
            vals=raw.decode("utf-8",errors="replace").split("\t")[1:]
            platforms|={v.strip().strip('"') for v in vals if re.fullmatch(
                r'"?GPL\d{2,8}"?',v.strip())}
        if raw.startswith(b"!series_matrix_table_begin"):
            if in_table or expression_headers_seen:
                raise ValueError("Unexpected repeated GEO data table")
            in_table=True
            continue
        if not in_table:
            continue
        if raw.startswith(b"!series_matrix_table_end"):
            break
        fields=raw.rstrip(b"\r\n").split(b"\t")
        if not expression_headers_seen:
            if fields[0].strip(b'"')!=b"ID_REF":
                raise ValueError("GEO source does not provide ID_REF gene probe matrix")
            sample_accessions=[x.decode("utf-8").strip('"') for x in fields[1:]]
            if len(sample_accessions)!=N_EXPECTED or len(set(sample_accessions))!=N_EXPECTED:
                raise ValueError("Source GEO not exactly 427 unique sample columns")
            if not all(re.fullmatch(r"GSM\d{5,9}",x) for x in sample_accessions):
                raise ValueError("Source series matrix sample headers are not GEO GSM")
            expression_headers_seen=True
            continue
        n_probes+=1
        if len(fields)!=N_EXPECTED+1:
            n_wrong_width+=1
            if n_wrong_width>10:raise ValueError("GEO original matrix has incompatible sample/probe widths")
            continue
        # Bounded independent numeric checks at beginning/mid/end, NOT full
        # retention/copy of original source signal values.
        for j in (1,2,3,107,215,319,425,426,427):
            s=fields[j].strip().strip(b'"')
            try:
                x=float(s)
                if x!=x or abs(x)==float("inf"):n_non_numeric+=1
            except ValueError:
                if s.lower() not in (b"na",b"null",b""):n_non_numeric+=1
    if not expression_headers_seen or n_probes<1000 or n_wrong_width>0:
        raise ValueError("GEO expression matrix not source qualified: insufficient probe rows")
    return {
        "source":"NCBI GEO GSE9588 processed original human liver expression matrix",
        "type":"microarray expression matrix; NOT RNA-seq counts",
        "original_unique_GSM_sample_columns":len(sample_accessions),
        "n_original_numeric_probe_rows":n_probes,
        "n_original_probe_rows_with_wrong_width":n_wrong_width,
        "n_bad_sample_numeric_checks":n_non_numeric,
        "source_GPL_platform_labels_only":sorted(platforms),
        "original_full_dataset_values_or_GSM_labels_exported":False,
        "RNA_only_gene_expression_matrix_structurally_qualified":n_non_numeric==0,
        "original_398_CYP_activity_same_person_crosswalk_verified":False,
        "original_proteomics_same_person_matrix_available":False,
        "clinical_PK_validation_executed":False,
        "no_supervised_CYP_phenotype_label_in_this_source":True,
    }


def main():
    try:
        with urlopen(Request(URL,headers={
            "User-Agent":"Omics2PK-original-GEO-GSE9588-data-schema-audit/1.0"
        }),timeout=90) as response:
            length=response.headers.get("Content-Length")
            if length and int(length)>150_000_000:
                raise ValueError("Original GEO zipped file too large for source runner")
            with gzip.GzipFile(fileobj=response) as stream:
                report=audit_gse_matrix_stream(stream)
        return {"status":"OFFICIAL_GEO_GSE9588_MATRIX_CHECK_COMPLETE",**report}
    except (OSError,URLError,HTTPError,ValueError,EOFError) as exc:
        return {"status":"GSE9588_ORIGINAL_MATRIX_NOT_RETRIEVED_ABSTAIN",
                "error_type":type(exc).__name__,
                "RNA_only_gene_expression_matrix_structurally_qualified":False,
                "original_398_CYP_activity_same_person_crosswalk_verified":False}


if __name__=="__main__":
    print(json.dumps(main(),sort_keys=True,indent=2))
