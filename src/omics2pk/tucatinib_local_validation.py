"""Predeclared original Table-1 local-only human-liver CYP protein->clearance benchmark.

No downloaded publisher article data are included in this public code. Run only
against a legitimately acquired, donor-aligned source table. This is in-vitro
apparent microsomal clearance, NOT clinical plasma clearance.
"""
from __future__ import annotations

import csv
import hashlib
import math
import random
from statistics import mean
from pathlib import Path

SOURCE_DOI="10.1016/j.dmd.2025.100061"
COLS=("donor_id","row_kind","source_doi","tucatinib_cl_app_ul_min_mg",
      "cyp2c8_protein_pmol_mg","cyp3a_protein_pmol_mg",
      "cyp2c8_activity_pmol_min_mg","cyp3a_activity_pmol_min_mg")


class TucatinibEvidenceError(ValueError):
    """Unavailable or incompatible original source observations."""


def _positive(row,field):
    try: v=float(row[field])
    except (KeyError,TypeError,ValueError): raise TucatinibEvidenceError(f"Missing source value for {field}") from None
    if not math.isfinite(v) or v<=0:
        raise TucatinibEvidenceError(f"Measured values of {field} must be strictly positive")
    return v


def qualify(rows):
    if len(rows)!=22:
        raise TucatinibEvidenceError("Original Table 1 should have 21 individual donor HLM lots + 1 pooled lot")
    eligible=[]
    pooled=0
    donors=set()
    for row in rows:
        if not set(COLS).issubset(row):
            raise TucatinibEvidenceError("Missing source-aligned required column")
        if row["source_doi"].strip()!=SOURCE_DOI:
            raise TucatinibEvidenceError("Unexpected source DOI")
        lot=row["row_kind"].strip().lower()
        if lot=="pool":
            pooled+=1
            continue
        if lot!="individual":
            raise TucatinibEvidenceError("Pool/individual lot classification must be explicit")
        donor=row["donor_id"].strip()
        if not donor or donor in donors:
            raise TucatinibEvidenceError("Biological donor IDs must be unique and nonempty")
        donors.add(donor)
        y=_positive(row,"tucatinib_cl_app_ul_min_mg")
        prot=(_positive(row,"cyp2c8_protein_pmol_mg"),_positive(row,"cyp3a_protein_pmol_mg"))
        activity=(_positive(row,"cyp2c8_activity_pmol_min_mg"),_positive(row,"cyp3a_activity_pmol_min_mg"))
        eligible.append((donor,y,prot,activity))
    if len(eligible)!=21 or pooled!=1:
        raise TucatinibEvidenceError("Do not pseudo-replicate 150-donor pooled HLM as an individual")
    return sorted(eligible)



# Published source sanity checks: original 21 independent HLM donors,
# Spearman correlations, rounded to three decimals in the article.
# These are not a substitute for permissioned original source rows.
PUBLISHED_SPEARMAN = {
    "cyp2c8_protein": 0.904,
    "cyp3a_protein": 0.882,
    "cyp2c8_activity": 0.781,
    "cyp3a_activity": 0.907,
}


def _average_ranks(values):
    indices = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(indices):
        stop = start + 1
        while stop < len(indices) and values[indices[stop]] == values[indices[start]]:
            stop += 1
        average_rank = (start + 1 + stop) / 2.0
        for j in range(start, stop):
            ranks[indices[j]] = average_rank
        start = stop
    return ranks


def _spearman(x, y):
    a, b = _average_ranks(x), _average_ranks(y)
    xa, yb = mean(a), mean(b)
    cross = sum((t - xa) * (v - yb) for t, v in zip(a, b))
    norm = math.sqrt(sum((t - xa)**2 for t in a) * sum((v - yb)**2 for v in b))
    if norm == 0:
        raise TucatinibEvidenceError("Constant biomarker or clearance")
    return cross / norm


def verify_published_table_correlations(eligible, *, tolerance=0.012):
    """Fail closed on the four independent Spearman checks printed in Table 1 paper.

    These detect transcription/column-order errors in a *locally* curated table.
    The single-donor clearance mean from printed rounded numbers is deliberately
    NOT required to equal the article narrative mean (see results document).
    """
    if len(eligible) != 21:
        raise TucatinibEvidenceError("Published reference check requires 21 donors")
    y = [row[1] for row in eligible]
    vectors = {
        "cyp2c8_protein": [row[2][0] for row in eligible],
        "cyp3a_protein": [row[2][1] for row in eligible],
        "cyp2c8_activity": [row[3][0] for row in eligible],
        "cyp3a_activity": [row[3][1] for row in eligible],
    }
    found = {k: _spearman(v, y) for k, v in vectors.items()}
    if any(abs(found[k] - PUBLISHED_SPEARMAN[k]) > tolerance for k in found):
        raise TucatinibEvidenceError(
            "Original published CYP-clearance correlations mismatch: check source "
            "columns, donor pairing, transcription and exclusion of the 150-donor pool"
        )
    return {k: round(v, 6) for k, v in found.items()}

def _predict_ridge(train, x_test, alpha=1.0):
    """Fixed 2-variable ridge with training-only centering/scaling, no fitting in test fold."""
    m=len(train)
    x_mu=[mean(r[0][j] for r in train) for j in (0,1)]
    x_sd=[math.sqrt(mean((r[0][j]-x_mu[j])**2 for r in train)) for j in (0,1)]
    x_sd=[v if v>1e-12 else 1.0 for v in x_sd]
    y_mu=mean(r[1] for r in train)
    zz=[((r[0][0]-x_mu[0])/x_sd[0],(r[0][1]-x_mu[1])/x_sd[1]) for r in train]
    a=alpha+sum(z[0]*z[0] for z in zz)
    b=sum(z[0]*z[1] for z in zz)
    d=alpha+sum(z[1]*z[1] for z in zz)
    c=sum(z[0]*(v[1]-y_mu) for z,v in zip(zz,train))
    e=sum(z[1]*(v[1]-y_mu) for z,v in zip(zz,train))
    det=a*d-b*b
    if det<=0:raise TucatinibEvidenceError("Singular ridge design")
    beta1=(d*c-b*e)/det
    beta2=(-b*c+a*e)/det
    p=y_mu+beta1*(x_test[0]-x_mu[0])/x_sd[0]+beta2*(x_test[1]-x_mu[1])/x_sd[1]
    return max(0.0,p)


def score_donors(eligible):
    if len(eligible)!=21: raise TucatinibEvidenceError("Need 21 independent donor lots")
    ys=[r[1] for r in eligible]
    losses={"mean":[],"proteins":[],"cyp_activity":[]}
    for i in range(21):
        tr=[r for j,r in enumerate(eligible) if j!=i]
        base=mean(r[1] for r in tr)
        preds={
            "mean":base,
            "proteins":_predict_ridge([(r[2],r[1]) for r in tr],eligible[i][2]),
            "cyp_activity":_predict_ridge([(r[3],r[1]) for r in tr],eligible[i][3])
        }
        for k,p in preds.items():
            losses[k].append(abs(p-ys[i]))
    delta=[losses["proteins"][i]-losses["mean"][i] for i in range(21)]
    rng=random.Random(20261016)
    bs=sorted(mean(delta[rng.randrange(21)] for _ in range(21)) for _ in range(1500))
    report={
        "n_unique_human_donors":21,
        "pooled_hlm_excluded":True,
        "source":"2025 tucatinib original 21-donor HLM Table 1 (manually sourced; no source values redistributed)",
        "outcome":"apparent in-vitro HLM tucatinib metabolic clearance, uL/mg protein/min",
        "protocol":"21 donor-held-out folds, two-predictor ridge alpha1 fixed, train-only scaling",
        "baseline_mae_ul_min_mg":round(mean(losses["mean"]),6),
        "protein_mae_ul_min_mg":round(mean(losses["proteins"]),6),
        "enzyme_activity_mae_ul_min_mg":round(mean(losses["cyp_activity"]),6),
        "protein_minus_baseline_mae_delta":round(mean(delta),6),
        "paired_donor_bootstrap_95pct_delta_mae":[round(bs[37],6),round(bs[1462],6)],
        "external_validation":False,
        "experimental_clint_not_patient_pk":True,
    }
    return report


def run_local_csv(filename: str) -> dict:
    path=Path(filename)
    original=path.read_bytes()
    if len(original)>500_000:
        raise TucatinibEvidenceError("Unexpected size for 22 original donor rows")
    rows=list(csv.DictReader(original.decode("utf-8-sig").splitlines()))
    eligible=qualify(rows)
    report=score_donors(eligible)
    report["source_table_spearman_reproduced"] = verify_published_table_correlations(eligible)
    report["mean_table_clearance_original_rounded_values"] = round(
        mean(row[1] for row in eligible), 6
    )
    report["article_narrative_mean_reported"] = 9.24
    report["paper_table_vs_narrative_mean_agree"] = abs(
        report["mean_table_clearance_original_rounded_values"] - 9.24
    ) < 0.1
    report["curated_local_input_sha256"]=hashlib.sha256(original).hexdigest()
    report["no_source_patient_or_donor_rows_exported"]=True
    return report
