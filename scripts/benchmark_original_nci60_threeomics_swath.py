#!/usr/bin/env python3
"""Original-source NCI-60 three-view doxorubicin GI50: research-only comparison.

Guo et al iScience 2019 DOI 10.1016/j.isci.2019.10.059, PMC6889472
original Table S1E (60 cell lines, 3171 proteins), cross-source matched
against the existing fixed mixOmics+IntLIM 47-line source. Molecular assays
come from DIFFERENT laboratories / culture batches. This is in-vitro PD,
not PK or paired same-biological-specimen multi-omics.
"""
from __future__ import annotations
import io,hashlib,json,urllib.request
from pathlib import Path
import numpy as np
from openpyxl import load_workbook
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder
from benchmark_public_nci60_twoomics import load_dataset,fit_view,predict,SEED

DOI="10.1016/j.isci.2019.10.059"
URL="https://ars.els-cdn.com/content/image/1-s2.0-S2589004219304407-mmc2.xlsx"
PIN="e9741865cfb12cced0a471f534966e2bcd2f9c6979f591f454daf9c0206779cb"
PREFIX_TO_CLASS={"BR":"BREAST","CNS":"CNS","CO":"COLON","LC":"NSCLC","LE":"LEUKEMIA","ME":"MELANOMA","OV":"OVARIAN","PR":"PROSTATE","RE":"RENAL"}
CAPS={"rna":300,"met":120,"prot":300}
ARMS={
"tumor_only":(),
"rna_only":("rna",),"metabolomics_only":("met",),"proteomics_only":("prot",),
"tumor_rna":("rna",),"tumor_metabolomics":("met",),"tumor_proteomics":("prot",),
"tumor_rna_metabolomics":("rna","met"),"tumor_rna_proteomics":("rna","prot"),
"tumor_proteomics_metabolomics":("prot","met"),
"tumor_threeomics":("rna","met","prot")
}
class SourceNotQualified(ValueError):pass
def canon(s):return "".join(ch for ch in str(s).upper() if ch.isalnum())
def distinct_map(keys,label):
    out={}
    for i,name in enumerate(keys):
        k=canon(name)
        if not k or k in out:raise SourceNotQualified("duplicate or empty "+label+" identity")
        out[k]=i
    return out
def load_official_proteome(return_gene_symbols=False):
    req=urllib.request.Request(URL,headers={"User-Agent":"Omics2PK published supplemental source audit"})
    with urllib.request.urlopen(req,timeout=90) as f:
        blob=f.read(45_000_001)
    sha=hashlib.sha256(blob).hexdigest()
    if sha!=PIN:raise SourceNotQualified("Original publisher supplement SHA256 mismatch: "+sha)
    book=load_workbook(io.BytesIO(blob),read_only=True,data_only=True)
    try:
        if "E" not in book.sheetnames:raise SourceNotQualified("Publisher Table S1E absent")
        ws=book["E"]
        if ws.max_row!=3172 or ws.max_column!=68:
            raise SourceNotQualified("Publisher protein Table S1E shape drift")
        rows=ws.iter_rows(values_only=True)
        header=next(rows)
        expected=("protein accession number","gene symbol","numPep","R2","best_mscore","numNA","transition_group_id","average protein intensity")
        if tuple(str(x) for x in header[:8])!=expected:
            raise SourceNotQualified("Original SWATH feature/source-unit schema drift")
        names,groups=[],[]
        for name in header[8:]:
            name=str(name)
            parts=name.split("_",1)
            if len(parts)!=2 or parts[0] not in PREFIX_TO_CLASS:raise SourceNotQualified("Unrecognized original tissue-prefixed cell-line label")
            names.append(canon(parts[1]));groups.append(PREFIX_TO_CLASS[parts[0]])
        if len(names)!=60 or len(set(names))!=60:raise SourceNotQualified("Original SWATH individual cell-line labels ambiguous")
        proteins=[]
        gene_symbols=[]
        values=[]
        accessions=set()
        for row in rows:
            accession=str(row[0]).strip()
            if not accession or accession in accessions:raise SourceNotQualified("Duplicate original SwissProt accession")
            accessions.add(accession)
            vals=[]
            for v in row[8:]:
                if v is None or str(v).strip() in ("","NA"):
                    vals.append(np.nan)
                else:
                    try:val=float(v)
                    except (TypeError,ValueError):raise SourceNotQualified("Non-numeric source protein intensity")
                    if val<0 or not np.isfinite(val):raise SourceNotQualified("Invalid original protein intensity")
                    vals.append(val)
            if len(vals)!=60:raise SourceNotQualified("Incorrect original SWATH specimen columns")
            proteins.append(accession);gene_symbols.append(str(row[1]).strip());values.append(vals)
        if len(proteins)!=3171:raise SourceNotQualified("Original protein count not 3171")
        mat=np.asarray(values,dtype=float).T
        if np.isnan(mat).mean()>0.25:raise SourceNotQualified("Unexpected >25% original protein missingness")
        # Table E is original averaged raw intensity, not already log10.
        mat=np.log10(1+mat)
        result=(tuple(names),tuple(groups),mat,sha)
        if return_gene_symbols:
            return result+(tuple(gene_symbols),)
        return result
    finally:book.close()
def run(ids,classes,y,rna,met,prot_names,prot_classes,prot):
    if len(ids)!=47 or len(set(ids))!=47:raise SourceNotQualified("Original fixed 47-line source gate failed")
    base=distinct_map(ids,"RNA/met/GI50")
    source=distinct_map(prot_names,"SWATH proteome")
    common=[k for k in base if k in source]
    if len(common)<35:raise SourceNotQualified("Too few exact, original same-lineage three-view cell-line IDs: "+str(len(common)))
    indices=np.array([base[k] for k in common])
    pindices=np.array([source[k] for k in common])
    c=np.asarray(classes)[indices]
    if any(c[i]!=prot_classes[pindices[i]] for i in range(len(c))):
        raise SourceNotQualified("Cross-source tissue lineage mismatch")
    arrays={"rna":rna[indices],"met":met[indices],"prot":prot[pindices]}
    y=y[indices];n=len(y)
    pred={name:np.full((3,n),np.nan) for name in ARMS}
    for repeat in range(3):
        kfold=KFold(n_splits=5,shuffle=True,random_state=SEED+repeat)
        for tr,te in kfold.split(np.arange(n)):
            folds={}
            for key,values in arrays.items():
                train,test,*_=fit_view(values[tr],values[te],CAPS[key])
                folds[key]=(train,test)
            encoder=OneHotEncoder(handle_unknown="ignore",sparse_output=False)
            base_train=encoder.fit_transform(c[tr,None]);base_test=encoder.transform(c[te,None])
            for name,views in ARMS.items():
                a=np.concatenate(([base_train] if name.startswith("tumor_") else [])+[folds[v][0] for v in views],axis=1)
                b=np.concatenate(([base_test] if name.startswith("tumor_") else [])+[folds[v][1] for v in views],axis=1)
                pred[name][repeat,te]=predict(a,b,y[tr])
    if any(not np.isfinite(p).all() for p in pred.values()):raise SourceNotQualified("Incomplete heldout predictions")
    loss={name:np.abs(p-y[None,:]).mean(axis=0) for name,p in pred.items()}
    scores={name:round(float(arr.mean()),6) for name,arr in loss.items()}
    boot=np.random.default_rng(SEED).integers(0,n,size=(1999,n))
    contrasts={}
    for against in ("tumor_only","tumor_rna","tumor_metabolomics","tumor_proteomics","tumor_rna_metabolomics","tumor_rna_proteomics","tumor_proteomics_metabolomics"):
        d=loss["tumor_threeomics"]-loss[against]
        lo,hi=np.quantile(d[boot].mean(axis=1),[.025,.975])
        contrasts["three_minus_"+against]={"delta_MAE":round(float(d.mean()),6),"descriptive_95pct_paired_loss_bootstrap":[round(float(lo),6),round(float(hi),6)]}
    return {
      "status":"EXPLORATORY_NCI60_ORIGINAL_CELL_LINE_MATCHED_CROSS_LAB_PD_NOT_CLINICAL_PK",
      "n_matched_unique_cell_lines":n,"n_RNA_features_original":rna.shape[1],
      "n_metabolites_original":met.shape[1],"n_proteins_original":prot.shape[1],
      "source_protein_DOI":DOI,"source_protein_original_table_SHA256":PIN,
      "same_original_cell_line_id":True,"same_contemporaneous_specimen":False,
      "response":"original processed mixOmics doxorubicin in-vitro GI50; NOT clearance, PK, or AUC",
      "repeated_folds":3,"folds":5,"independent_biological_observation_count":n,
      "protocol":"same n/folds for all arms; fixed Ridge alpha10; train-only feature selection cap RNA300/met120/protein300; train-only PCA3 per modality; original protein intensities log10(1+x)",
      "mae_same_source_GI50_units":scores,
      "paired_descriptive_contrasts":contrasts,
      "interpretation_limits":"Original RNA/metabolite/protein assays were performed in different studies/batches; no matched molecular specimen, no nested model selection, no external cell-line validation, no pharmacokinetic claim.",
      "original_line_profiles_or_GI50_exported":False}
def main():
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--data",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    opts=ap.parse_args()
    ids,c,y,_,_,r,m=load_dataset(opts.data)
    pn,pc,p,hashsum=load_official_proteome()
    report=run(ids,c,y,r,m,pn,pc,p)
    opts.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=="__main__":main()
