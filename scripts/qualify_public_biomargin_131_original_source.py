"""Public BIOMARGIN/ROCKET MOFA2 source-level pre-join audit.

Never print or commit individual BIOS IDs or expression values. Download
public source material ephemerally; discard bytes after aggregate QC.
This is NOT a real TAC pharmacokinetic cohort or MOFA2 numerical reproduction.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from gzip import GzipFile
from hashlib import sha1, sha256
from io import BytesIO, TextIOWrapper
import json
import re
from urllib.request import Request, urlopen

AUTHOR_COMMIT = "afad658913857955abc9dbac160acecda22e92c5"
AUTHOR_NOTEBOOK_BLOB = "616915d949797900e211d1595b07401d64c6c93c"
NOTEBOOK = (
    "https://raw.githubusercontent.com/InsermRightLab/MOFA2-Tinel-et-al/"
    + AUTHOR_COMMIT + "/Notebook_MOFA2_model_training.ipynb"
)
GEO_GSE = "GSE175718"
COUNTS = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175718/"
    "suppl/GSE175718_counts.csv.gz"
)
ORIGINAL_MOFA_VIEWS = (
    "Blood_RNAseq_Step1_filter_VSTnorm_top5k.csv",
    "blood_mRNA_Step1_ct_top5k_170621.csv",
    "Blood_miRNA_Step1_ct_norm_exp_log_all751_130721.csv",
    "biopsy_mRNA_Step1_ct_top5k_170621.csv",
    "Biopsy_miRNA_Step1_ct_norm_exp_log_all752_130721.csv",
    "Urine_mRNA_Step1_norm_exp_log_010721.csv",
)
SAMPLE_PATTERN = re.compile(r"BIOS\.[0-9]{2}\.[0-9]{4}")
FILE_PATTERN = re.compile(r'read\.csv2?\(file="([^"]+\.csv)"')
MAX_NOTEBOOK_BYTES = 2_000_000
MAX_COUNTS_COMPRESSED_BYTES = 15_000_000


class PublicBIOMARGINAuditError(ValueError):
    pass


def bounded_download(url, *, max_bytes):
    req = Request(url, headers={"User-Agent": "Omics2PK-aggregate-public-omics-source-audit/1.0"})
    try:
        with urlopen(req, timeout=75) as response:
            payload = response.read(max_bytes + 1)
    except OSError as exc:
        raise PublicBIOMARGINAuditError("Public source unavailable; no substitute data used") from exc
    if not payload or len(payload) > max_bytes:
        raise PublicBIOMARGINAuditError("Public source payload unexpectedly large or empty")
    return payload


def audit_author_notebook(payload):
    if not isinstance(payload, bytes) or not payload:
        raise PublicBIOMARGINAuditError("Exact original author notebook bytes required")
    blob = sha1(b"blob " + str(len(payload)).encode() + b"\0" + payload).hexdigest()
    if blob != AUTHOR_NOTEBOOK_BLOB:
        raise PublicBIOMARGINAuditError("Original author Git blob identity changed")
    try:
        obj = json.loads(payload.decode("utf-8-sig"))
        code = "\n".join("".join(c.get("source", [])) for c in obj.get("cells", [])
                         if c.get("cell_type") == "code")
    except (UnicodeError, TypeError, ValueError) as exc:
        raise PublicBIOMARGINAuditError("Invalid original notebook structure") from exc
    order = re.search(r'samples\.order\s*=\s*c\(([\s\S]+?)\)\s*(?:\n|$)', code)
    if not order:
        raise PublicBIOMARGINAuditError("Original biological person order absent")
    original = SAMPLE_PATTERN.findall(order.group(1))
    if len(original) != 131 or len(set(original)) != 131:
        raise PublicBIOMARGINAuditError("Author's exact 131 independent BIOS subject list changed")
    sites = Counter(s.split(".")[1] for s in original)
    if len(sites) != 4 or sum(sites.values()) != 131:
        raise PublicBIOMARGINAuditError("Expected four original source center strata")
    files = FILE_PATTERN.findall(code)
    if len(files) != 6 or tuple(files) != ORIGINAL_MOFA_VIEWS:
        raise PublicBIOMARGINAuditError("Original MOFA molecular views or source filenames changed")
    # Return sample ID *set* only internally for exact GEO overlap. Never log it.
    return set(original), {
        "source": "original_author_training_notebook",
        "source_commit": AUTHOR_COMMIT,
        "source_git_blob": blob,
        "n_author_unique_original_biological_subjects": len(original),
        "n_original_centers": len(sites),
        "n_source_molecular_views": len(files),
        "source_molecular_view_classes": [
            "blood_RNAseq", "blood_mRNA_microarray", "blood_miRNA",
            "biopsy_mRNA_microarray", "biopsy_miRNA", "urine_targeted_mRNA",
        ],
        "n_complete_six_view_subjects_published_not_reverified_in_this_audit": 35,
        "author_input_processed_matrices_present_in_this_notebook": False,
        "no_drug_PK_linked_or_evaluated": True,
    }


def audit_public_geo_counts_gzip(payload, mofa_subjects, *, production=True):
    """Inspect released counts without retaining/publicizing the person matrix.

    GEO IDs may be GSM accessions, whereas original MOFA2 prepared VST
    matrices use BIOS IDs: ZERO direct match means ABSTAIN, never row-align.
    """
    if not isinstance(payload, bytes) or not payload.startswith(b"\x1f\x8b"):
        raise PublicBIOMARGINAuditError("Expected original gzip GEO count CSV")
    if len(payload) > MAX_COUNTS_COMPRESSED_BYTES:
        raise PublicBIOMARGINAuditError("Unexpected large GEO matrix")
    gene_rows = 0
    n_assays = None
    id_matches = None
    zero_gene = 0
    try:
        with GzipFile(fileobj=BytesIO(payload)) as gz:
            with TextIOWrapper(gz, encoding="utf-8-sig", newline="") as fh:
                reader = csv.reader(fh)
                header = next(reader)
                if len(header) != len(set(header)) or not 3 <= len(header) <= 650:
                    raise PublicBIOMARGINAuditError("GEO count column names invalid")
                n_assays = len(header) - 1
                if production and not 100 <= n_assays <= 500:
                    raise PublicBIOMARGINAuditError("GEO count assay count incompatible with study")
                id_matches = len(set(header[1:]) & mofa_subjects)
                for row in reader:
                    if len(row) != len(header):
                        raise PublicBIOMARGINAuditError("Ragged GEO gene/assay count matrix")
                    gene_rows += 1
                    if gene_rows > 150_000:
                        raise PublicBIOMARGINAuditError("GEO count file unexpectedly large")
                    if gene_rows <= 50:
                        # Numeric columns only, true source raw-count validation.
                        try:
                            values = [float(v) for v in row[1:]]
                        except ValueError as exc:
                            raise PublicBIOMARGINAuditError("GEO RNA measurements not numeric") from exc
                        if any((not v.is_integer() or v < 0) for v in values):
                            raise PublicBIOMARGINAuditError("GEO raw counts appear nonintegral/negative")
                        zero_gene += int(all(v == 0 for v in values))
    except (EOFError, OSError, UnicodeError, csv.Error, StopIteration) as exc:
        raise PublicBIOMARGINAuditError("Original GEO compressed count source unreadable") from exc
    if production and not 5_000 <= gene_rows <= 100_000:
        raise PublicBIOMARGINAuditError("Original GEO features not consistent with transcriptomics")
    if gene_rows < 1:
        raise PublicBIOMARGINAuditError("GEO RNA matrix contains no features")
    return {
        "geo_accession": GEO_GSE,
        "source_file": "GSE175718_counts.csv.gz",
        "original_compressed_sha256": sha256(payload).hexdigest(),
        "number_gene_feature_rows": gene_rows,
        "number_assay_columns": n_assays,
        "n_direct_exact_BIOS_column_matches": id_matches,
        "join_status": ("DIRECT_ID_INTERSECTION_FOUND_NOT_VALIDATED" if id_matches
                        else "CROSSWALK_REQUIRED_GEO_GSM_VS_AUTHOR_BIOS_NO_POSITION_JOIN"),
        "RNA_raw_counts_are_not_VST_normalized_author_MOFA_input": True,
        "true_patient_pk_and_CYP_genotype_not_presented_by_this_count_audit": True,
        "no_source_rows_or_biological_IDs_exported": True,
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-geo-counts", action="store_true",
                    help="Download real public GEO count source; never save it")
    args = ap.parse_args(argv)
    source = bounded_download(NOTEBOOK, max_bytes=MAX_NOTEBOOK_BYTES)
    subjects, report = audit_author_notebook(source)
    result = {
        "status": "ORIGINAL_PUBLIC_AUTHOR_NOTEBOOK_131_BIOLOGICAL_SUBJECTS_SOURCE_CONFIRMED",
        "original_author_notebook": report,
        "external_original_geo_counts_retrieved": False,
        "numerical_six_view_MOFA2_reproduction_claimed": False,
        "real_tacrolimus_PK_or_PGx_join_confirmed": False,
    }
    if args.include_geo_counts:
        geofile = bounded_download(COUNTS, max_bytes=MAX_COUNTS_COMPRESSED_BYTES)
        result["original_geo_counts"] = audit_public_geo_counts_gzip(geofile, subjects)
        result["external_original_geo_counts_retrieved"] = True
    print(json.dumps(result, sort_keys=True, indent=2))
    return result


if __name__ == "__main__":
    main()
