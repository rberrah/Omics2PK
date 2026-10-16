"""Independent actual-source replication of EPA CYP held-donor-out models via sklearn.

Uses separately implemented scikit-learn StandardScaler/Ridge, not the
project's normal-equation solver. Original human liver records are ephemerally
loaded and NEVER output. This is a numerical/parity replication, NOT an
independent biological donor cohort.
"""
from __future__ import annotations
import hashlib
import json
import math

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from scripts.audit_epa_liver_rna_activity import _get
from scripts.inspect_epa_liver_crosswalk import parse_epa_workbook
from scripts.validate_epa_liver_rna_cyp import EXPECTED_SHA
from omics2pk.epa_liver_cyp_validation import (
    EXPECTED_GENES, EPALinkageError, read_paired_epa, validate_one_gene
)


def independent_ridge_loo(records):
    y=np.array([math.log1p(r[2]) for r in records],dtype=float)
    X=np.array([[r[3],r[4],math.log1p(r[1])] for r in records],dtype=float)
    ref_errors, rna_errors = [], []
    for held in range(len(records)):
        idx=np.arange(len(records))!=held
        true=records[held][2]
        for take, collector in (([0,1],ref_errors),([0,1,2],rna_errors)):
            scale=StandardScaler().fit(X[idx][:,take])
            train=scale.transform(X[idx][:,take])
            test=scale.transform(X[held:held+1][:,take])
            model=Ridge(alpha=1.0, fit_intercept=True, solver="cholesky")
            model.fit(train,y[idx])
            # expm1 of original model's log target; floor physical negative activity
            transformed=float(model.predict(test)[0])
            estimate=max(0.,math.expm1(min(transformed,30)))
            collector.append(abs(estimate-true))
    return {"age_sex":round(float(np.mean(ref_errors)),6),
            "age_sex_plus_rna":round(float(np.mean(rna_errors)),6)}


def main():
    original=_get()
    if hashlib.sha256(original).hexdigest()!=EXPECTED_SHA:
        raise EPALinkageError("Original EPA source SHA changed")
    data,qc=read_paired_epa(parse_epa_workbook(original))
    results={}
    for gene in EXPECTED_GENES:
        expected=validate_one_gene(data[gene])["mae_original_units"]
        independently_computed=independent_ridge_loo(data[gene])
        residuals={name:abs(independently_computed[name]-expected[name])
                   for name in independently_computed}
        if max(residuals.values())>0.00001:
            raise AssertionError(f"Independent sklearn replica differs for {gene}; no validation claim")
        results[gene]={
            "n_original_donors":len(data[gene]),
            "reference_mae_sklearn":independently_computed["age_sex"],
            "rna_adjusted_mae_sklearn":independently_computed["age_sex_plus_rna"],
            "max_disagreement_original_units":max(residuals.values()),
            "status":"INDEPENDENT_ESTIMATOR_PARITY_PASS",
        }
    print(json.dumps({
        "status":"ORIGINAL_SOURCE_EPA100_ALL_8_CYP_INDEPENDENT_ESTIMATOR_NUMERIC_PARITY",
        "scikit_learn_version":__import__("sklearn").__version__,
        "source_sha256":EXPECTED_SHA,
        "no_individual_source_rows_exported":True,
        "not_external_biological_replication":True,
        "n_validated_enzyme_pairs":len(results),
        "results":results
    },indent=2,sort_keys=True))


if __name__=="__main__":
    main()
