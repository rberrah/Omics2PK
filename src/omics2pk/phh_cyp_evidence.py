"""Exploratory six-donor PHH RNA to experimentally observed CYP activity.

All functions are pure and use no clinical data. Small n (6): descriptive
group-held-out comparisons ONLY, not pharmacokinetic clinical validation.
"""
from __future__ import annotations

from math import isfinite, log1p, sqrt
import re
from statistics import mean
from typing import Mapping, Sequence

CYPS = ("CYP1A2","CYP2B6","CYP2C9","CYP2C19","CYP2D6","CYP3A4")
SIX_DONORS = frozenset(("025","211","789","330","409","416"))

class EvidenceError(ValueError):
    """Measured-outcome / biological-unit linkage is not source-qualified."""


def extract_mean(value: str) -> float | None:
    """Only a clearly uncensored leading mean, such as '12.4 ± 1.3'."""
    s = str(value).strip().replace(",", ".").replace("−","-")
    if not s or s.startswith(("<",">","≤","≥","~")):
        return None
    v=re.match(r"^([+]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(?:\s*(?:±|\+/-|\(|$|[A-Za-z]))",s)
    if not v:
        return None
    x=float(v.group(1))
    return x if isfinite(x) and x >= 0 else None


def parse_published_functional_table(rows: Sequence[Sequence[str]], *, kind: str) -> dict:
    """Only tables whose original six source PHH codes and CYP rows match."""
    if kind not in ("substrate_clearance","metabolite_formation"):
        raise EvidenceError("Endpoint type must be prespecified")
    if len(rows) != 7 or any(len(r)!=8 for r in rows):
        raise EvidenceError("Original 6-CYP, 6-donor table structure must match 7x8")
    if rows[0][0].strip() != "Enzyme":
        raise EvidenceError("No explicit enzyme header")
    if kind == "substrate_clearance" and "Substrate" not in rows[0][1]:
        raise EvidenceError("Substrate missing from original endpoint heading")
    if kind == "metabolite_formation" and "Metabolite" not in rows[0][1]:
        raise EvidenceError("Metabolite missing from original endpoint heading")
    donor_cols=[]
    for head in rows[0][2:]:
        m=re.fullmatch(r"PHH(025|211|789|330|409|416)",head.strip(),re.I)
        if not m: raise EvidenceError("Unrecognized PHH source donor label")
        donor_cols.append(m.group(1))
    if len(set(donor_cols))!=6 or set(donor_cols)!=SIX_DONORS:
        raise EvidenceError("Six donor IDs are not in exact one-to-one correspondence")
    endpoints={}
    for row in rows[1:]:
        enzyme=row[0].strip().upper()
        if enzyme not in CYPS or enzyme in endpoints:
            raise EvidenceError("Unexpected or repeated CYP endpoint")
        endpoints[enzyme]={d:extract_mean(row[i+2]) for i,d in enumerate(donor_cols)}
    if set(endpoints)!=set(CYPS):
        raise EvidenceError("Original substrate panel incomplete")
    return endpoints


def pearson(x: Sequence[float], y: Sequence[float]) -> float | None:
    if len(x)<3 or len(x)!=len(y):
        return None
    xx,yy=mean(x),mean(y)
    a=sum((p-xx)*(q-yy) for p,q in zip(x,y))
    b=sqrt(sum((p-xx)**2 for p in x)*sum((q-yy)**2 for q in y))
    return round(a/b,6) if b>1e-12 else None


def held_out_biomarker_test(
    features: Mapping[str,float], outcomes: Mapping[str,float | None],
    *, fixed_ridge_alpha: float=1.0,
) -> dict:
    """Per donor LODO, feature standardization fit on training donors only."""
    if fixed_ridge_alpha!=1.0: raise EvidenceError("No post-hoc penalty tuning at n=6")
    donors=sorted(d for d,v in outcomes.items() if v is not None and d in features)
    if len(donors)<4 or len(donors)>6 or not set(donors)<=SIX_DONORS:
        raise EvidenceError("Require 4 to 6 source-matched donor outcomes")
    if any(not isfinite(features[d]) or features[d]<0 for d in donors):
        raise EvidenceError("RNA counts/CPM must be observed nonnegative values")
    xs=[log1p(features[d]) for d in donors]
    ys=[float(outcomes[d]) for d in donors]
    if any(not isfinite(y) or y<0 for y in ys):
        raise EvidenceError("CYP metabolic outcomes must be measured and nonnegative")
    base_errors,model_errors=[],[]
    for i in range(len(donors)):
        xx=[x for j,x in enumerate(xs) if j!=i]
        yy=[y for j,y in enumerate(ys) if j!=i]
        avg_x,avg_y=mean(xx),mean(yy)
        std_x=sqrt(mean([(z-avg_x)**2 for z in xx]))
        base=avg_y
        if std_x>1e-12:
            z=[(x-avg_x)/std_x for x in xx]
            beta=sum(a*(b-avg_y) for a,b in zip(z,yy))/(sum(a*a for a in z)+fixed_ridge_alpha)
            estimated=avg_y+beta*((xs[i]-avg_x)/std_x)
        else:
            estimated=avg_y
        # Physiological endpoint nonnegative; same rule predeclared in every fold.
        estimated=max(0.0,estimated)
        base_errors.append(abs(base-ys[i]))
        model_errors.append(abs(estimated-ys[i]))
    mae0=mean(base_errors)
    mae1=mean(model_errors)
    return {
        "n_independent_donors":len(donors),
        "n_donor_heldout_folds":len(donors),
        "rna_feature":"log1p(pseudobulk-CYP-CPM)",
        "method":"one-feature ridge, alpha fixed 1, train-only standardization",
        "baseline":"training-donor arithmetic mean of measured assay endpoint",
        "donor_expression_activity_pearson_r_exploratory":pearson(xs,ys),
        "baseline_mae_original_assay_units":round(mae0,6),
        "rna_mae_original_assay_units":round(mae1,6),
        "delta_mae_rna_minus_baseline":round(mae1-mae0,6),
        "rna_improved_over_mean_baseline_in_this_cohort":mae1<mae0,
        "interpretation":"n<=6, exploratory in-vitro donor-held-out; no population PK/external or clinical validation",
    }
