"""HLC original matching RNA+germline SNP experimental ablations (not protein).

The original SNP feature metadata encodes dbSNP numeric IDs *without* 'rs'.
Only source-verified rs776746 (CYP3A5) and rs1799853 (CYP2C9)
are used. Genotype alleles are CC/CT/TT and dosage is a literal T count;
DO NOT call T the CYP3A5 expressor allele without source-strand QC.
The genotype source has 195 original person keys, but only 178 overlap
original liver expression source, and 171-172 have complete activity covariates.
"""
from __future__ import annotations
import csv
from zipfile import ZipFile

import numpy as np

from omics2pk.hlc_real_cyp_validation import (
    load_original_hlc,OriginalHLCQualificationError
)

TARGETS=(
    ("CYP3A4M_activity","776746",("CYP3A4","CYP3A5")),
    ("CYP3A4T_activity","776746",("CYP3A4","CYP3A5")),
    ("CYP2C9_activity","1799853",("CYP2C9",)),
)


def _load_cyp_genotypes(genotype_zip):
    import pandas as pd
    with ZipFile(genotype_zip) as z:
        feature=pd.read_csv(z.open("features.txt"),sep="\t",dtype=str)
        wanted=feature.loc[feature["dbSNP_rsID"].isin(("776746","1799853")),
                           ["dbSNP_rsID","feature_id"]]
        if len(wanted)!=2 or len(set(wanted["dbSNP_rsID"]))!=2:
            raise OriginalHLCQualificationError("Original named PGx rs loci unavailable")
        lookup=dict(zip(wanted["feature_id"],wanted["dbSNP_rsID"]))
        reader=csv.reader((s.decode("utf-8-sig",errors="replace")
                           for s in z.open("genotype.txt")),delimiter="\t")
        header=next(reader)
        subjects=header[1:]
        selected={}
        for line in reader:
            if line[0] in lookup:
                key=lookup[line[0]]
                if key in selected:
                    raise OriginalHLCQualificationError("Duplicate SNP from original source")
                selected[key]=line[1:]
    if len(subjects)!=195 or len(set(subjects))!=195 or set(selected)!=set(lookup.values()):
        raise OriginalHLCQualificationError("Original genotype source identity/missing SNP invalid")
    return subjects,selected


def _small_cv(x,y,seed):
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import KFold
    pred=np.empty(len(y))
    for train,test in KFold(5,shuffle=True,random_state=seed).split(y):
        model=make_pipeline(SimpleImputer(strategy="median"),StandardScaler(),
                            Ridge(alpha=1.0))
        model.fit(x[train],np.log1p(y[train]))
        pred[test]=np.maximum(0.,np.expm1(
            np.clip(model.predict(x[test]),-20.,20.)))
    return float(np.abs(pred-y).mean())


def benchmark_original_hlc_rna_plus_pgx(expression_zip,phenotype_zip,genotype_zip):
    r=load_original_hlc(expression_zip,phenotype_zip)
    gid,snps=_load_cyp_genotypes(genotype_zip)
    expr_ids={k:i for i,k in enumerate(r["sample_ids"])}
    common=sum(k in expr_ids for k in gid)
    if common!=178:
        raise OriginalHLCQualificationError("Observed genotype/expression biological key intersection changed")
    rows={}
    for target,rs,genes in TARGETS:
        yp=r["outcomes"][target]
        probe_indices=[int(k) for g in genes for k in r["gene_rows"][g]]
        models={"age_sex":[],"age_sex_plus_RNA":[],"age_sex_plus_SNP":[],
                "age_sex_plus_RNA_plus_SNP":[]}
        x=[]
        y=[]
        missing_genotypes=0
        for i,source_id in enumerate(gid):
            if source_id not in expr_ids:continue
            index=expr_ids[source_id]
            if not np.isfinite(r["clinical"][index]).all() or not np.isfinite(yp[index]):
                continue
            if yp[index]<0:continue
            alleles=snps[rs][i]
            dosage={"CC":0.,"CT":1.,"TT":2.,"NA":float("nan")}.get(alleles)
            if dosage is None:
                raise OriginalHLCQualificationError("Unexpected real source genotype allele encoding")
            missing_genotypes+=not np.isfinite(dosage)
            x.append([*r["clinical"][index],*r["rna"][index,probe_indices],dosage])
            y.append(yp[index])
        xx=np.asarray(x,dtype=float)
        yy=np.asarray(y,dtype=float)
        if len(yy) not in (171,172):
            raise OriginalHLCQualificationError("Original matching PGx subset size unexpectedly changed")
        nRNA=len(probe_indices)
        candidate={
            "age_sex":xx[:,:2],
            "age_sex_plus_RNA":xx[:,:2+nRNA],
            "age_sex_plus_SNP":xx[:,[0,1,2+nRNA]],
            "age_sex_plus_RNA_plus_SNP":xx,
        }
        for seed in (20261018,20261019,20261020):
            for name,a in candidate.items():
                models[name].append(_small_cv(a,yy,seed))
        rows[target]={
            "source_CYP_probe_activity":target,
            "dbSNP_rsID":"rs"+rs,
            "n_original_complete_source_people":len(yy),
            "n_SNP_missing_calls_fold_median_imputed":int(missing_genotypes),
            "n_original_mRNA_probe_features":nRNA,
            "model_MAE_mean_three_independent_fixed_CV_partitions":
                {name:round(float(np.mean(v)),6) for name,v in models.items()},
            "original_SNP_original_expression_and_CYP_assay_subject_rows_not_emitted":True,
        }
    return {
        "status":"HLC_AUTHENTIC_RNA_AND_GENOMICS_TWO_VARIANTS_E2_MEASURED_CYP_ACTIVITY",
        "n_genotyped_source_people":195,
        "n_exact_original_genotype_and_expression_people":178,
        "n_activity_endpoints":3,
        "source_SNP_allele_strand_function_not_confirmed":True,
        "source_proteomic_CYP_abundance_not_measured_here":True,
        "clinical_plasma_PK_outcome_not_measured_here":True,
        "all_normalization_and_genotype_imputation_train_fold_only":True,
        "source_donor_identifiers_genotypes_assays_or_per_person_predictions_exported":False,
        "results":rows,
    }
