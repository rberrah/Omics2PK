"""Rerun eight related CYP functional assays from original EPA paired 100 livers."""
from __future__ import annotations
import hashlib
import json

from scripts.audit_epa_liver_rna_activity import _get, SOURCE_DOI
from scripts.inspect_epa_liver_crosswalk import parse_epa_workbook
from scripts.validate_epa_liver_rna_cyp import EXPECTED_SHA
from omics2pk.epa_liver_cyp_validation import EPALinkageError, read_paired_epa
from omics2pk.epa_eight_cyp_sensitivity import calculate_eight_gene_sensitivity


def original_real_epa():
    payload=_get()
    if hashlib.sha256(payload).hexdigest()!=EXPECTED_SHA:
        raise EPALinkageError("Original public EPA source checksum differs: abstain")
    original=parse_epa_workbook(payload)
    records,qc=read_paired_epa(original)
    result=calculate_eight_gene_sensitivity(records)
    result["source_doi"]=SOURCE_DOI
    result["original_document_sha256"]=EXPECTED_SHA
    result["original_paired_source_qc_all_eight"]=qc
    return result


if __name__=="__main__":
    print(json.dumps(original_real_epa(),sort_keys=True,indent=2,ensure_ascii=False))
