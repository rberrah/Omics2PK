#!/usr/bin/env python3
"""GTEx V8 ORIGINAL public bulk tissue median TPM: source-qualified CYP3A4/5.

This is ONLY a tissue-median expression screen, not within-donor concordance,
PBMC, CYP enzymatic activity, clinical drug PK, or blood-derived vesicle RNA.
Original GTEx V8 official open access sample-wise GCT is substantially larger;
a same-person blood-liver correlation needs the sample-level data and ID join.
"""
from __future__ import annotations
import csv,gzip,hashlib,io,json,math,urllib.request

SOURCE="https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
TISSUES=("Liver","Whole_Blood","Small_Intestine_Terminal_Ileum")
GENES=("CYP3A4","CYP3A5","FMO3","ALB","PTPRC")

class GTExSourceError(ValueError):pass

def extract_original_gct(stream, source_file_sha256):
    lines=(line.decode("utf-8") for line in gzip.GzipFile(fileobj=stream))
    try:
        version=next(lines).strip()
        shape=next(lines).strip().split()
        header=next(csv.reader([next(lines)],delimiter="\t"))
    except (StopIteration, UnicodeError, OSError) as e:
        raise GTExSourceError("Missing original GCT source header") from e
    if not version.startswith("#1.2") or len(shape)!=2 or len(header)<53 or header[:2]!=["Name","Description"]:
        raise GTExSourceError("Not an original GTEx median expression GCT source")
    if any(tissue not in header for tissue in TISSUES):
        raise GTExSourceError("Required liver/blood/intestine original tissue labels absent")
    if "Whole_Blood"=="Cells_EBV-transformed_lymphocytes":
        raise GTExSourceError("Impossible internal ontology declaration")
    out={}
    nrows=0
    for line in lines:
        row=next(csv.reader([line],delimiter="\t"))
        if len(row)!=len(header):raise GTExSourceError("Ragged original expression source")
        nrows+=1
        gene=row[1].strip()
        if gene not in GENES:continue
        if gene in out:raise GTExSourceError("Duplicate source gene symbols; need Ensembl disambiguation")
        vals={}
        for tissue in TISSUES:
            raw=row[header.index(tissue)]
            try:val=float(raw)
            except (ValueError,TypeError) as e:raise GTExSourceError("Non-numeric tissue median") from e
            if not math.isfinite(val) or val<0:raise GTExSourceError("Nonfinite source median")
            vals[tissue]=round(val,7)
        out[gene]={"source_gencode_gene_id":row[0],
                   "median_TPM_by_original_GTEx_tissue":vals,
                   "blood_over_liver_median_TPM_ratio":round(vals["Whole_Blood"]/vals["Liver"],7) if vals["Liver"]>0 else None,
                   "blood_is_a_qualified_proxy_of_individual_hepatic_expression":False}
    if nrows<30000:raise GTExSourceError("Unqualified GTEx source gene count")
    if set(out)!=set(GENES):
        raise GTExSourceError("Missing gene(s) in original GTEx source: "+",".join(sorted(set(GENES)-set(out))))
    return {
      "status":"OFFICIAL_GTEX_V8_ORIGINAL_TISSUE_MEDIAN_EXPRESSION_VERIFIED_NOT_PAIRED_NOT_PK",
      "source_url":SOURCE,"source_sha256":source_file_sha256,
      "original_gct_version":version,"source_gene_rows_qualified":nrows,
      "source_tissues":list(TISSUES),
      "measured_specimen":"bulk post-mortem tissue RNA-seq; GTEx Whole_Blood is NOT PBMC",
      "gene_level_original_TPM":out,
      "same_donor_cross_tissue_correlation_measured":False,
      "donor_genotype_available_in_this_open_median_file":False,
      "CYP3A4_CYP3A5_splice_isoform_or_activity_measured":False,
      "observed_human_drug_AUC_or_CL_present":False,
      "no_source_expression_vectors_or_individual_donor_identifiers_exported":True,
    }

def main():
    req=urllib.request.Request(SOURCE,headers={"User-Agent":"Omics2PK-original-source GTEx V8 cross-tissue scientific pilot"})
    with urllib.request.urlopen(req,timeout=90) as response:
        body=response.read(45_000_001)
    if len(body)>=45_000_001 or len(body)<100000 or body[:2]!=b"\x1f\x8b":
        raise GTExSourceError("Expected bounded authentic compressed GTEx GCT")
    return extract_original_gct(io.BytesIO(body),hashlib.sha256(body).hexdigest())

if __name__=="__main__":
    try:print(json.dumps(main(),indent=2,sort_keys=True))
    except (OSError,GTExSourceError) as error:
        print(json.dumps({"status":"ABSTAIN_NO_AUTHENTIC_GTEX_SOURCE",
                          "error_type":type(error).__name__,
                          "http_code":getattr(error,"code",None),
                          "source_contract_reason":str(error)[:140] if isinstance(error,GTExSourceError) else "network_download_failed"},sort_keys=True))
        raise SystemExit(2)
