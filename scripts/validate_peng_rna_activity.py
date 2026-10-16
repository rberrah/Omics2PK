#!/usr/bin/env python3
"""Pre-registered small-N public human hepatocyte RNA->CYP assay pilot.

Downloads original CC-BY supplementary DOCX + GEO counts, validates byte hashes
and six donor identities, computes donor pseudobulk CYP RNA in memory, compares
direct CYP expression with measured paired assays using donor-held-out folds.
No raw source matrices, CYP donor records, cell barcodes or patient data exported.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
from io import BytesIO, TextIOWrapper
import json
import re
import tarfile
import xml.etree.ElementTree as ET
from zipfile import ZipFile

from omics2pk.phh_cyp_evidence import CYPS, SIX_DONORS, EvidenceError, held_out_biomarker_test, parse_published_functional_table
from scripts.qualify_peng_phh_source import (
    DOCX_NS, GEO_TAR, GEO_SOFT, PMC_META,
    load, parse_aws_metadata, geo_sample_metadata,
)

EXPECTED_SUPPLEMENT_SHA = "deeb835fe82d7e0e883447268354b9abb8f3ca4950639be3b68b802d3c6183bf"
EXPECTED_GEO_TAR_SHA = "7f114365ba5e977e685f9cee108def64deb0a31c96ea13d95e19c0b42c2b01e8"
EXPECTED_GEO_SOFT_SHA = "443635f9a9c527ee32b1e0d3559f91b2583304c36e78cba5c972d9681a72012a"
SOURCE_CITATION = "Peng et al 2025 Cell Regeneration DOI 10.1186/s13619-025-00258-6; GEO GSE289636; PMC12532551.1 S4/S5"

def original_functional_matrices(docx: bytes) -> tuple[dict,dict]:
    if hashlib.sha256(docx).hexdigest()!=EXPECTED_SUPPLEMENT_SHA:
        raise EvidenceError("Frozen original PMC functional assay supplement SHA changed")
    with ZipFile(BytesIO(docx)) as z:
        root=ET.fromstring(z.read("word/document.xml"))
    matched={}
    for tbl in root.iter(DOCX_NS+"tbl"):
        rows=[]
        for tr in tbl.findall(DOCX_NS+"tr"):
            rows.append([" ".join("".join(tc.itertext()).split()) for tc in tr.findall(DOCX_NS+"tc")])
        if not rows or not rows[0]:
            continue
        head=" ".join(rows[0][:2])
        if "Enzyme" not in head or not ("Substrate" in head or "Metabolite" in head):
            continue
        kind = "substrate_clearance" if "Substrate" in head else "metabolite_formation"
        if kind in matched: raise EvidenceError("Duplicate supplement functional endpoint")
        matched[kind] = parse_published_functional_table(rows,kind=kind)
    if set(matched)!={"substrate_clearance","metabolite_formation"}:
        raise EvidenceError("Both distinct measured original CYP panel assays required")
    return matched["substrate_clearance"],matched["metabolite_formation"]


def expression_from_public_geo(gse_tar: bytes) -> tuple[dict,dict]:
    if hashlib.sha256(gse_tar).hexdigest()!=EXPECTED_GEO_TAR_SHA:
        raise EvidenceError("Frozen original public GEO raw archive SHA changed")
    cpm={}
    source_summary={}
    with tarfile.open(fileobj=BytesIO(gse_tar),mode="r:") as tar:
        members=[m for m in tar.getmembers() if m.isfile()]
        if len(members)!=6:
            raise EvidenceError("Expected precisely six original .csv.gz member files")
        for entry in members:
            suffix=re.search(r"_row_counts_(025|211|789|330|409|416)\.csv\.gz$",entry.name)
            if suffix is None: raise EvidenceError("Unknown source scRNA counts filename/batch")
            batch=suffix.group(1)
            if batch in cpm: raise EvidenceError("Biological donor duplicated")
            gz=tar.extractfile(entry)
            if gz is None: raise EvidenceError("GEO compressed count matrix not readable")
            seq=TextIOWrapper(gzip.GzipFile(fileobj=gz),encoding="utf-8-sig",newline="")
            reader=csv.reader(seq)
            try: head=next(reader)
            except StopIteration: raise EvidenceError("Empty GEO transcriptome")
            if head and head[0].strip().casefold()=="barcode" and "graphclust" in head and "nUMI" in head:
                raise EvidenceError("GEO_SUPPLEMENT_IS_SEURAT_CELL_METADATA_NOT_TRANSCRIPT_COUNTS")
            n_cells=len(head)-1
            if n_cells<10 or n_cells>150000: raise EvidenceError("Expected cell-by-gene count matrix, not donor-level summary")
            n_genes=0
            total=0.0
            gene_counts={}
            for row in reader:
                if len(row)!=len(head): raise EvidenceError("Ragged cell count matrix")
                gene=row[0].strip().upper()
                if not gene: raise EvidenceError("Unnamed expression row")
                if gene in gene_counts: raise EvidenceError("Duplicate CYP genes")
                val=0.0
                for x in row[1:]:
                    if not x:
                        raise EvidenceError("Missing UMI count; no imputation")
                    try: v=float(x)
                    except ValueError: raise EvidenceError("Non-numeric expression count") from None
                    if not (0.0<=v<1e10): raise EvidenceError("Invalid single-cell count")
                    val+=v
                total+=val
                n_genes+=1
                if gene in CYPS:
                    gene_counts[gene]=val
            if total<=0: raise EvidenceError("No observed RNA UMIs in library")
            if n_genes<1000: raise EvidenceError("Unexpectedly few transcriptome gene rows")
            cpm[batch]={g:gene_counts.get(g,0.0)*1e6/total for g in CYPS}
            source_summary[batch]={
                "n_cells":n_cells,"n_genes":n_genes,
                "cyp_genes_found_in_rows":sorted(gene_counts),
                "cyp_genes_nonzero":sum(gene_counts.get(g,0)>0 for g in CYPS),
                "library_total_umi":round(total),
            }
    if set(cpm)!=SIX_DONORS:
        raise EvidenceError("GEO and six original functional donors fail crosswalk")
    return cpm,source_summary


def validate_public_phh() -> dict:
    document_url,license_code = parse_aws_metadata(load(PMC_META,200_000))
    if license_code.upper()!="CC BY":
        raise EvidenceError("Source reuse licence must remain CC BY")
    docx=load(document_url,12_000_000)
    scr,mfr=original_functional_matrices(docx)
    meta=geo_sample_metadata(load(GEO_SOFT,10_000_000))
    if meta["geo_family_sha256"]!=EXPECTED_GEO_SOFT_SHA:
        raise EvidenceError("Original GEO donor metadata SHA changed")
    raw=load(GEO_TAR,40_000_000)
    gene_cpm,structure=expression_from_public_geo(raw)
    report=[]
    for endpoint_name,endpoint in (("substrate_clearance",scr),("metabolite_formation",mfr)):
        for gene in CYPS:
            measured=endpoint[gene]
            x={donor:gene_cpm[donor][gene] for donor in sorted(SIX_DONORS)}
            n_obs=sum(measured[d] is not None for d in SIX_DONORS)
            if n_obs<4:
                report.append({"endpoint":endpoint_name, "cyp_gene":gene,
                               "n_observed_donors":n_obs,
                               "status":"ABSTAIN_TOO_FEW_SOURCE_MEASUREMENTS"})
                continue
            result=held_out_biomarker_test(x,measured)
            report.append({"endpoint":endpoint_name, "cyp_gene":gene,"status":"EXPLORATORY_LODO",**result})
    n_improved=sum(r.get("rna_improved_over_mean_baseline_in_this_cohort",False) for r in report)
    return {
        "status":"REAL_SOURCE_SIX_DONOR_EXPLORATORY_RNA_TO_IN_VITRO_CYP_FUNCTION",
        "source":SOURCE_CITATION,
        "n_source_unique_biological_donors":6,
        "source_license_code":license_code,
        "file_sha256":{"supplement":EXPECTED_SUPPLEMENT_SHA,"geo_raw":EXPECTED_GEO_TAR_SHA,"geo_soft":EXPECTED_GEO_SOFT_SHA},
        "donor_batch_crosswalk_exact":True,
        "gene_count_QC_aggregate":{
            "min_genes_per_donor":min(s["n_genes"] for s in structure.values()),
            "max_genes_per_donor":max(s["n_genes"] for s in structure.values()),
            "min_cells_per_donor":min(s["n_cells"] for s in structure.values()),
            "max_cells_per_donor":max(s["n_cells"] for s in structure.values()),
            "all_donors_have_six_cyp_gene_rows":all(len(s["cyp_genes_found_in_rows"])==6 for s in structure.values()),
            "n_total_cells_not_independent":sum(s["n_cells"] for s in structure.values()),
        },
        "n_prespecified_cyp_assay_comparisons":12,
        "n_comparisons_rna_lower_mae_than_mean":n_improved,
        "n_comparisons_abstained":sum(r["status"].startswith("ABSTAIN") for r in report),
        "individual_enzyme_results":report,
        "all_assay_units_are_source_specific":True,
        "limits":[
            "Original PHH biological preparation (not timed patient dosing/concentration); RNA and assay sample times need source confirmation",
            "Only six donors, LO(D)O trains five; no stable predictive performance or external generalization inference",
            "Protein activity/isoform expression cannot substitute directly for molecular PK parameters",
            "No data-trained feature selection, no SWEET, no high-dimensional machine learning, no clinical PK",
            "Both assay endpoints and all CYP share the same six donors: do not interpret 12 tests as 72 independent donors",
        ],
        "individual_donor_or_single_cell_values_exported":False,
    }

def main():
    try:
        report=validate_public_phh()
    except EvidenceError as exc:
        if str(exc)!="GEO_SUPPLEMENT_IS_SEURAT_CELL_METADATA_NOT_TRANSCRIPT_COUNTS":
            raise
        report={
            "status":"ABSTAIN_GEO_HAS_SEURAT_METADATA_NOT_EXPRESSION_COUNTS",
            "n_original_donors":6,
            "paired_functional_assay_source_verified":True,
            "raw_fastq_or_gene_counts_required":True,
            "missing_prerequisite":"Reconstruct CYP expression from original raw RNA reads; GEO supplementary files contain cell clusters, UMAP, PCA, nUMI and QC fields rather than gene-count matrices",
            "no_cpm_or_rna_to_cyp_model_fitted":True,
            "no_clinical_pk_claim":True,
        }
    print(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True))

if __name__=="__main__":
    main()
