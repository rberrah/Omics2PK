"""Run original U.S. EPA 100-donor CYP RNA mismatched-label negative controls.

Official pinned EPA source only, loaded transiently in the Actions runner.
No raw rows, unique original subject IDs, concentrations or mRNA values are
printed, committed or uploaded. E2 enzyme activity, not patient PK.
"""
from __future__ import annotations
from hashlib import sha256
import json

from scripts.audit_epa_liver_rna_activity import _get, SOURCE_DOI, PAPER_DOI
from scripts.inspect_epa_liver_crosswalk import parse_epa_workbook
from scripts.validate_epa_liver_rna_cyp import EXPECTED_SHA
from omics2pk.epa_liver_cyp_validation import EPALinkageError, read_paired_epa
from omics2pk.epa_rna_negative_control import (
    PERMUTATIONS, SEED, evaluate_epa_permutation_negative_control,
)


def original_epa_source_negative_controls():
    payload=_get()
    if sha256(payload).hexdigest()!=EXPECTED_SHA:
        raise EPALinkageError("Official original EPA source SHA changed; ABSTAIN")
    matched,qc=read_paired_epa(parse_epa_workbook(payload))
    report=evaluate_epa_permutation_negative_control(
        matched,n_permutations=PERMUTATIONS,seed=SEED
    )
    report["original_verified_DOI"]=SOURCE_DOI
    report["source_primary_paper_DOI"]=PAPER_DOI
    report["original_document_sha256"]=EXPECTED_SHA
    report["source_donor_pair_qc_counts_only"]=qc
    return report


if __name__=="__main__":
    print(json.dumps(original_epa_source_negative_controls(),indent=2,sort_keys=True))
