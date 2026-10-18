"""Exploratory source-annotated bile-acid MS1 proxy vs real midazolam CL/F.

This is neither a CYP3A enzyme activity measurement nor a confirmatory
biological PK biomarker. Nine GNPS library spectral matches may represent
adducts/isomers, not nine verified molecules. Predefined source-name panel,
never selected from PK outcomes. All real data supplied LOCAL ONLY.
"""
from __future__ import annotations
from hashlib import file_digest
from io import TextIOWrapper
from pathlib import Path
import argparse
import csv
import json
import re
from zipfile import ZipFile
import numpy as np

from scripts.validate_private_jarmusch13_multiomics_pipeline import (
    read_prepared_private_zip, ORIGINAL_GNPS_SHA256
)

PANEL_REGEX = re.compile(
    r"^(?:Spectral Match to )?(?:Cholic acid|Chenodeoxycholic acid|"
    r"Deoxycholic acid|Glycocholic acid|Taurocholic acid)\b",
    re.IGNORECASE,
)


class BileProxyError(ValueError):
    pass


def candidate_scans_from_spectral_matches(rows):
    """Extract putative bile acid spectral library scan IDs, not identities."""
    scans = set()
    for row in rows:
        if PANEL_REGEX.search(str(row.get("Compound_Name",""))):
            scan=str(row.get("#Scan#",""))
            if not scan.isdecimal():
                raise BileProxyError("Original MS2 spectral scan index missing")
            scans.add("GNPS_MS1:"+scan)
    return sorted(scans)


def inspect_original_annotations(source_path):
    path=Path(source_path)
    if not path.is_file():
        raise BileProxyError("Original GNPS source must be provided privately")
    with path.open("rb") as fd:
        if file_digest(fd,"sha256").hexdigest()!=ORIGINAL_GNPS_SHA256:
            raise BileProxyError("Unqualified GNPS source archive SHA")
    with ZipFile(path) as z:
        candidates=[m for m in z.namelist() if m.startswith("DB_result/") and m.endswith(".tsv")]
        if len(candidates)!=1:
            raise BileProxyError("Original source spectral library annotations absent")
        with z.open(candidates[0]) as fd:
            annotations=list(csv.DictReader(TextIOWrapper(fd,encoding="utf-8-sig"),
                                            delimiter="\t"))
    if len(annotations)!=346:
        raise BileProxyError("Original GNPS library-match source count changed")
    scans=candidate_scans_from_spectral_matches(annotations)
    if len(scans)!=9:
        raise BileProxyError("Source-derived exploratory bile acid panel changed")
    return scans


def compare_proxy(matrix, clinical, genotypes, pk, feature_ids, target_scans):
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import OneHotEncoder, StandardScaler

    names={f:i for i,f in enumerate(feature_ids)}
    if len(names)!=len(feature_ids) or not set(target_scans).issubset(names):
        raise BileProxyError("MS2 annotated ions do not match GNPS MS1 IDs")
    selected=[names[s] for s in target_scans]
    if not all(np.sum(matrix[:,i]>0)>=10 for i in selected):
        raise BileProxyError("Unexpected predose observation support for spectral panel")

    def predict(train,test,y_train):
        scaler=StandardScaler().fit(train)
        return float(np.exp(Ridge(alpha=1).fit(
            scaler.transform(train),np.log(y_train)
        ).predict(scaler.transform(test))[0]))

    outputs={}
    for endpoint, column in (("day1",0),("day9",1)):
        valid=np.isfinite(pk[:,column])
        X=matrix[valid]; C=clinical[valid]; G=genotypes[valid]; y=pk[valid,column]
        logx=np.log1p(X/np.maximum(X.sum(axis=1,keepdims=True),1)*1_000_000)
        errs={"clinical_pgx":[],"clinical_pgx_plus_putative_bile_ms1_score":[]}
        for held in range(len(y)):
            tr=np.arange(len(y))!=held;te=~tr
            encoder=OneHotEncoder(handle_unknown="ignore").fit(G[tr])
            btrain=np.column_stack((C[tr],encoder.transform(G[tr]).toarray()))
            btest=np.column_stack((C[te],encoder.transform(G[te]).toarray()))
            scaler=StandardScaler().fit(logx[tr][:,selected])
            # Fixed zero-outcome-information orientation/weights. Fitted
            # standardization uses ONLY outer train persons.
            score_tr=scaler.transform(logx[tr][:,selected]).mean(axis=1)[:,None]
            score_te=scaler.transform(logx[te][:,selected]).mean(axis=1)[:,None]
            a=predict(btrain,btest,y[tr])
            b=predict(np.column_stack((btrain,score_tr)),
                      np.column_stack((btest,score_te)),y[tr])
            errs["clinical_pgx"].append(abs(a-y[held]))
            errs["clinical_pgx_plus_putative_bile_ms1_score"].append(abs(b-y[held]))
        baseline=np.mean(errs["clinical_pgx"])
        omic=np.mean(errs["clinical_pgx_plus_putative_bile_ms1_score"])
        outputs[endpoint]={
            "n_independent_original_people":len(y),
            "clinical_pgx_MAE":round(float(baseline),6),
            "clinical_pgx_plus_bile_proxy_MAE":round(float(omic),6),
            "delta_MAE_incremental":round(float(omic-baseline),6),
            "independent_external_validation":False,
        }
    return outputs


def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument("--private-input-zip",required=True)
    ap.add_argument("--original-gnps-zip",required=True)
    args=ap.parse_args(argv)
    cohort,_=read_prepared_private_zip(args.private_input_zip)
    scans=inspect_original_annotations(args.original_gnps_zip)
    results=compare_proxy(
        cohort["matrix"],cohort["clinical"],cohort["genotypes"],cohort["pk"],
        cohort["feature_ids"],scans)
    print(json.dumps({
        "status":"SOURCE_MATCHED_MS2_BILE_ACID_CLASS_PROXY_EXPLORATORY_ONLY",
        "n_putatively_annotated_ms1_features":len(scans),
        "results":results,
        "spectral_identity_and_CYP3A_activity_unverified":True,
        "pk_measurement_units_unverified":True,
        "true_external_validation":False,
        "individual_source_data_never_exported":True,
    },indent=2,sort_keys=True))


if __name__=="__main__":
    main()
