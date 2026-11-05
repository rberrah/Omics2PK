#!/usr/bin/env python3
"""NCBI GEO primary hepatocyte CYP3A5 RNA/genotype source pairing audit.

Zhong et al. 2020 (PLoS Genet 10.1371/journal.pgen.1008662):
GSE124074 RNA-seq and GSE147628 source genotype microarrays, component
of GSE124076. Source-qualified sample identifiers may overlap,
but presence of a genotype microarray != accessible genotype call at
rs776746, and gene-level TMM != exon-junction RNA nor CYP3A5 function.
No original per-donor values, names, or genomic IDs are exported.
"""
import csv
import gzip
import hashlib
import io
import json
import re
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

SERIES={
 "RNA_seq_GSE124074":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124074/matrix/GSE124074_series_matrix.txt.gz",
 "genotyping_GSE147628":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147628/matrix/GSE147628_series_matrix.txt.gz",
}
class GEOQualificationError(ValueError):
    pass


def parse_geo_series_matrix(raw,*,series,max_uncompressed=35_000_000):
    if not isinstance(raw,(bytes,bytearray)) or not 20<len(raw)<10_000_000:
        raise GEOQualificationError("Bounded authentic GEO compressed source required")
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(raw),mode="rb") as h:
            expanded=h.read(max_uncompressed+1)
    except (OSError,EOFError) as exc:
        raise GEOQualificationError("Original GEO series metadata compression invalid") from exc
    if len(expanded)>max_uncompressed:
        raise GEOQualificationError("Original GEO metadata exceeded source size guard")
    try:
        lines=expanded.decode("utf-8",errors="strict").splitlines()
    except UnicodeDecodeError as exc:
        raise GEOQualificationError("Original GEO series metadata encoding invalid") from exc
    samples={}
    for line in lines:
        if not line.startswith("!Sample_"):
            continue
        parts=next(csv.reader([line],delimiter="\t",quotechar='"'))
        if len(parts)<2 or parts[0] not in {"!Sample_title","!Sample_geo_accession"}:
            continue
        if parts[0] in samples:
            raise GEOQualificationError("Duplicate original GEO sample annotation")
        samples[parts[0]]=parts[1:]
    title=samples.get("!Sample_title")
    gsm=samples.get("!Sample_geo_accession")
    if title is None or gsm is None or len(title)!=len(gsm) or len(set(gsm))!=len(gsm):
        raise GEOQualificationError("Original GEO sample title/GSM pairing missing or duplicated")
    if any(not re.fullmatch(r"GSM[0-9]{6,10}",z) for z in gsm):
        raise GEOQualificationError("Unexpected GEO original accession")
    # Strip title whitespace only, preserve original provenance.
    return {
        "source_series":series,
        "n_GEO_sample_records":len(gsm),
        "source_SHA256":hashlib.sha256(raw).hexdigest(),
        "n_unique_sample_titles":len(set(title)),
        "original_titles":title,
    }


def audit_original():
    parsed={}
    for name,url in SERIES.items():
        try:
            with urlopen(Request(url,headers={
                "User-Agent":"Omics2PK-original-GEO-variant-omics-source-audit/1.0"}),timeout=35) as response:
                raw=response.read(10_000_001)
        except (HTTPError,URLError,OSError,TimeoutError) as exc:
            raise GEOQualificationError(name+": OFFICIAL_SOURCE_UNAVAILABLE") from exc
        parsed[name]=parse_geo_series_matrix(raw,series=name)
    r=parsed["RNA_seq_GSE124074"]
    d=parsed["genotyping_GSE147628"]
    # This is a *metadata-only* overlap of source titles. It is NOT
    # the original 60-donor expression/genotype quality-controlled join.
    overlapping=set(r["original_titles"])&set(d["original_titles"])
    result={
        "status":"ORIGINAL_GEO_RNA_GENOTYPE_SAMPLE_METADATA_AUDIT_NOT_FUNCTIONAL_PK",
        "original_RNA_series":"GSE124074",
        "original_genotype_series":"GSE147628",
        "parent_superseries":"GSE124076",
        "original_RNA_GEO_sample_count":r["n_GEO_sample_records"],
        "original_genotyping_GEO_sample_count":d["n_GEO_sample_records"],
        "n_exact_matching_source_titles_before_genotype_QC":len(overlapping),
        "original_source_file_SHA256_by_series":{
            k:v["source_SHA256"] for k,v in parsed.items()
        },
        "original_RNA_expression_values_or_donor_identifiers_exported":False,
        "published_exonQTL_rs776746_CYP3A5_not_eGene_not_independently_recalculated":True,
        "original_rs776746_genotype_calls_qualified":False,
        "original_exon_specific_splice_expression_qualified":False,
        "original_CYP3A5_specific_activity_assay_qualified":False,
        "no_original_clinical_tacrolimus_PK_or_PD":True,
        "additional_source_needed":"original_rs776746_genotype_call_matrix_and_exon_quantification_per_matched_sample",
    }
    if result["original_RNA_GEO_sample_count"]!=60:
        raise GEOQualificationError("Unexpected original RNA Geno series sample count")
    if result["original_genotyping_GEO_sample_count"]!=96:
        raise GEOQualificationError("Unexpected original genotyping series sample count")
    return result


if __name__=="__main__":
    try:
        print(json.dumps(audit_original(),sort_keys=True,indent=2))
    except GEOQualificationError as exc:
        print(json.dumps({"status":"ORIGINAL_GEO_SOURCE_NOT_QUALIFIED_ABSTAIN",
                          "reason_category":str(exc)[:200]}))
        raise SystemExit(2)
