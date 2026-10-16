"""Reproducible human liver CYP qPCR mRNA -> measured microsomal enzyme activity.

Only original official EPA linked worksheets are used. Donor rows are loaded
temporarily and never exported. n=100 original biological donors, eight
biochemically defined CYP pairs; not an individual clinical PK dataset.
"""
from __future__ import annotations
from collections import Counter
from math import expm1, isfinite, log1p, sqrt
from random import Random
from statistics import mean
from typing import Mapping,Sequence

from scripts.inspect_epa_liver_crosswalk import _col_indices

EXPECTED_GENES=(
    "CYP1A2","CYP2B6","CYP2C9","CYP2C19",
    "CYP2D6","CYP2E1","CYP3A4","CYP4A11",
)
N_DONORS=100

class EPALinkageError(ValueError):
    """Identity or measurements cannot support a real donor-held-out test."""


def _lookup(header,label):
    normalized=[" ".join(str(x).lower().split()) for x in header]
    return next((i for i,x in enumerate(normalized) if label in x),None)


def _sample_key(row,k,x):
    if k>=len(row) or x>=len(row):
        return None
    a=row[k];b=row[x]
    if a is None or b is None or not str(a).strip() or not str(b).strip():
        return None
    return (str(a).strip(),str(b).strip())


def _numeric(x):
    if isinstance(x,bool):return None
    try:v=float(x)
    except (ValueError,TypeError):return None
    return v if isfinite(v) else None


def _panel(header,rows):
    k,x,cyps=_col_indices(header)
    if len(k)!=1 or len(x)!=1:
        raise EPALinkageError("Both original sample IDs must be present and unambiguous")
    age=_lookup(header,"age")
    sex=_lookup(header,"gender")
    if age is None or sex is None:
        raise EPALinkageError("Expected demographic controls absent from original source")
    lookup={}
    for row in rows:
        key=_sample_key(row,k[0],x[0])
        if key is None:continue
        if key in lookup:raise EPALinkageError("Duplicated original biological subject key")
        lookup[key]=row
    if len(lookup)!=N_DONORS:
        raise EPALinkageError("Expected exactly 100 biological donors in original source")
    return lookup,cyps,age,sex


def read_paired_epa(wb):
    if set(wb)!={"Figure 1","Figure 2","Table 3"}:
        raise EPALinkageError("Original three-sheet EPA source identity changed")
    activity,acols,aage,asex=_panel(*wb["Figure 2"])
    mrna,mcols,mage,msex=_panel(*wb["Table 3"])
    if set(activity)!=set(mrna):
        raise EPALinkageError("Actual donor identities for mRNA and activity do not match exactly")
    if not set(EXPECTED_GENES).issubset(acols) or not set(EXPECTED_GENES).issubset(mcols):
        raise EPALinkageError("Expected 8 RNA/activity pairs not present in source")
    measurements={}
    qc={}
    for gene in EXPECTED_GENES:
        records=[]
        miss={"missing_rna":0,"missing_activity":0,"nonpositive_rna":0,"negative_activity":0,"missing_age_or_sex":0}
        for key in sorted(activity):
            rA=activity[key];rM=mrna[key]
            rna=_numeric(rM[mcols[gene]]) if mcols[gene]<len(rM) else None
            act=_numeric(rA[acols[gene]]) if acols[gene]<len(rA) else None
            age=_numeric(rM[mage]) if mage<len(rM) else None
            sexval=str(rM[msex]).strip().lower() if msex<len(rM) else ""
            if age is None or not (0<=age<=120) or sexval not in ("male","female"):
                miss["missing_age_or_sex"]+=1;continue
            if rna is None:miss["missing_rna"]+=1;continue
            if rna<=0:miss["nonpositive_rna"]+=1;continue
            if act is None:miss["missing_activity"]+=1;continue
            if act<0:miss["negative_activity"]+=1;continue
            # Age/sex agree in the ORIGINAL two source sheets when available
            age2=_numeric(rA[aage]) if aage<len(rA) else None
            sex2=str(rA[asex]).strip().lower() if asex<len(rA) else ""
            if age2 is None or age2!=age or sex2!=sexval:
                raise EPALinkageError("Demographic contradiction between original two sheets")
            records.append((key,float(rna),float(act),float(age),1.0 if sexval=="male" else 0.0))
        measurements[gene]=records
        qc[gene]={"n_complete_independent_donors":len(records),**miss}
    return measurements,qc


def _solve(matrix, vector):
    """Gaussian pivoting; ridge SPD positive definite, no external ML dependency."""
    n=len(vector)
    a=[list(row)+[float(vector[i])] for i,row in enumerate(matrix)]
    for p in range(n):
        best=max(range(p,n),key=lambda r:abs(a[r][p]))
        if abs(a[best][p])<1e-10:raise EPALinkageError("Singular feature/penalty matrix")
        a[p],a[best]=a[best],a[p]
        scale=a[p][p]
        for j in range(p,n+1):a[p][j]/=scale
        for i in range(n):
            if i==p:continue
            fac=a[i][p]
            for j in range(p,n+1):a[i][j]-=fac*a[p][j]
    return [a[i][n] for i in range(n)]


def _ridge(train_x,train_y,test_x,alpha=1.0):
    if alpha!=1:raise EPALinkageError("No posthoc tuning allowed in small study")
    d=len(train_x[0])
    mu=[mean(row[j] for row in train_x) for j in range(d)]
    sd=[sqrt(mean((row[j]-mu[j])**2 for row in train_x)) for j in range(d)]
    sd=[s if s>1e-10 else 1 for s in sd]
    zz=[[(row[j]-mu[j])/sd[j] for j in range(d)] for row in train_x]
    avg=mean(train_y)
    mat=[[sum(row[i]*row[j] for row in zz)+(1.0 if i==j else 0.0)
          for j in range(d)] for i in range(d)]
    rhs=[sum(z[i]*(y-avg) for z,y in zip(zz,train_y)) for i in range(d)]
    beta=_solve(mat,rhs)
    zt=[(test_x[j]-mu[j])/sd[j] for j in range(d)]
    return avg+sum(b*z for b,z in zip(beta,zt))


def validate_one_gene(records, *, seed=20261016):
    if len(records)<60 or len(records)>100 or len(set(x[0] for x in records))!=len(records):
        raise EPALinkageError("Need >=60 unique source donors, no cell/assay pseudoreplication")
    errs={k:[] for k in ("reference_logmean","age_sex","rna_alone","age_sex_plus_rna")}
    raw_y=[x[2] for x in records]
    for hold in range(len(records)):
        train=[x for i,x in enumerate(records) if i!=hold]
        truth=records[hold][2]
        yt=[log1p(x[2]) for x in train]
        xref=(records[hold][3],records[hold][4])
        plans={
            "reference_logmean":mean(yt),
            "age_sex":_ridge([[x[3],x[4]] for x in train],yt,list(xref)),
            "rna_alone":_ridge([[log1p(x[1])] for x in train],yt,[log1p(records[hold][1])]),
            "age_sex_plus_rna":_ridge(
                [[x[3],x[4],log1p(x[1])] for x in train],
                yt,[*xref,log1p(records[hold][1])]),
        }
        for name,p in plans.items():
            pred=max(0.0,expm1(min(p,30)))  # no extreme extrapolation to infinity
            errs[name].append(abs(pred-truth))
    results={k:round(mean(v),6) for k,v in errs.items()}
    delta=[a-b for a,b in zip(errs["age_sex_plus_rna"],errs["age_sex"])]
    rng=Random(seed)
    boot=sorted(mean(delta[rng.randrange(len(delta))] for _ in delta)
                for _ in range(1500))
    return {
        "n_independent_original_donors":len(records),
        "n_donor_heldout_folds":len(records),
        "outcome":"enzyme activity, pmol/min/mg protein, original scale",
        "input":"original CYP transcript abundance relative to housekeeping qPCR",
        "models":"fixed ridge(alpha=1) log1p activity; log1p RNA; train-only standardization",
        "mae_original_units":results,
        "delta_mae_adjusted_plus_rna_minus_age_sex":round(mean(delta),6),
        "paired_donor_bootstrap_95pct_delta":[round(boot[37],6),round(boot[1462],6)],
        "RNA_reduces_adjusted_MAE_in_original_cohort":mean(delta)<0,
        "RNA_increment_clinical_PGx_or_PopPK":False,
    }
