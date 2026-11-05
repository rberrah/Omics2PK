#!/usr/bin/env python3
"""Original-data NCI-60 pathway-coverage audit and train-frozen exploratory ablations.

This compares all Metabolon source features with ONLY SOURCE-NAMED candidates,
unknown/isobar/artefact features, and six predeclared biochemistry hypothesis
panels. Panels do not prove biochemical identity or clinically predictive PK.
Do not choose a 'best' post-hoc pathway or claim metabolome = pathway coverage.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder
from benchmark_public_nci60_twoomics import load_dataset,fit_view,predict,SEED
from benchmark_original_nci60_threeomics_swath import load_official_proteome,distinct_map,PIN
from omics2pk.nci60_biochemical_coverage import PATHWAY_PANELS,audit_labels,status_of_original_label

REPO_SOURCE="Mathelab/NCI60_GeneMetabolite_Data@4c140708c18c41910458414a626b1a51d3dd2a85"
CHEM_PANELS=tuple(PATHWAY_PANELS)
PRIMARY=("tumor_only","tumor_rna","tumor_protein","tumor_metabolome_all",
         "tumor_metabolome_named","tumor_metabolome_unresolved",
         "tumor_rna_metabolome_all","tumor_rna_metabolome_named",
         "tumor_rna_protein","tumor_rna_protein_metabolome_all",
         "tumor_rna_protein_metabolome_named")
def read_source_annotated_metabolites(path):
    with path.open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=["metabnames","id"]:
            raise ValueError("Original IntLIM metabolite annotation source columns changed")
        records=list(reader)
    if any(r["metabnames"]!=r["id"] for r in records):
        raise ValueError("Source original label and the only provided id disagree")
    return tuple(r["id"] for r in records)

def mode_names():
    return PRIMARY+tuple("tumor_path_"+p for p in CHEM_PANELS)+tuple("tumor_rna_path_"+p for p in CHEM_PANELS)

def compare_original(ids,classes,y,rna,met,met_ids,prot_names,prot_classes,prot,protein_genes,
                     meta_labels):
    original_set=set(meta_labels)
    if set(met_ids)!=original_set or len(met_ids)!=len(meta_labels):
        raise ValueError("NCI original metadata labels not identical to measured metabolite feature columns")
    meta=audit_labels(meta_labels)
    named_idx=[i for i,n in enumerate(met_ids) if status_of_original_label(n)=="named_candidate_identity_unverified"]
    unresolved_idx=[i for i,n in enumerate(met_ids) if status_of_original_label(n)!="named_candidate_identity_unverified"]
    if len(named_idx)!=111 or len(unresolved_idx)!=169:
        raise ValueError("Source named/unnamed source-level count mismatch")
    src=distinct_map(prot_names,"original SWATH proteome")
    base=distinct_map(ids,"original RNA/metabolomics/GI50")
    common=[k for k in base if k in src]
    if len(common)!=47:raise ValueError("Retain exact original 47 source-concordant cell lines for primary inference")
    a=np.asarray([base[k] for k in common],dtype=int)
    b=np.asarray([src[k] for k in common],dtype=int)
    cc=np.asarray(classes)[a]
    if any(cc[i]!=prot_classes[b[i]] for i in range(len(common))):
        raise ValueError("Source lineage mismatch")
    rna=rna[a];met=met[a];prot=prot[b];y=y[a]
    original_rna_names=set(meta.get("unused_rna_genes",[]))
    # Array ID annotations in exported CSV are original exact gene symbols.
    # Do not infer that high-throughput probes cover isoforms or genotype.
    # Coverage counts ignore case differences but not synonyms or aliases.
    coverage={}
    direct_drug_genes=("TOP2A","TOP2B","ABCB1","ABCC1","ABCG2","CBR1","CBR3")
    all_genes=set(protein_genes)
    # Actual original RNA feature strings passed by caller in separate tuple.
    for name,p in PATHWAY_PANELS.items():
        coverage[name]={"n_nominated_source_named_metabolites":len(p["metabolites"]),
                        "source_names_not_MSI_level_1_proven":True,
                        "missing_sentinel_analytes":list(p["not_measured_key_metabolites"]),
                        "mechanism_hypothesis_only":p["limitation"]}
    panel_idx={k:[met_ids.index(n) for n in p["metabolites"]] for k,p in PATHWAY_PANELS.items()}
    def views(train,test):
        rv=fit_view(rna[train],rna[test],300)[:2]
        pv=fit_view(prot[train],prot[test],300)[:2]
        allm=fit_view(met[train],met[test],120)[:2]
        named=fit_view(met[train][:,named_idx],met[test][:,named_idx],120)[:2]
        unresolved=fit_view(met[train][:,unresolved_idx],met[test][:,unresolved_idx],120)[:2]
        panels={k:fit_view(met[train][:,ii],met[test][:,ii],len(ii))[:2]
                for k,ii in panel_idx.items()}
        return {"rna":rv,"protein":pv,"metabolome_all":allm,
                "metabolome_named":named,"metabolome_unresolved":unresolved,"panels":panels}
    names=mode_names()
    preds={k:np.full((3,47),np.nan) for k in names}
    for repeat in range(3):
        for train,test in KFold(5,shuffle=True,random_state=SEED+repeat).split(np.arange(47)):
            o=views(train,test)
            encoder=OneHotEncoder(handle_unknown="ignore",sparse_output=False)
            basetrain=encoder.fit_transform(cc[train,None])
            basetest=encoder.transform(cc[test,None])
            def model(*keys,subpath=None):
                a=[basetrain]; b=[basetest]
                for k in keys:
                    a.append(o[k][0]); b.append(o[k][1])
                if subpath is not None:
                    a.append(o["panels"][subpath][0]);b.append(o["panels"][subpath][1])
                return predict(np.concatenate(a,axis=1),np.concatenate(b,axis=1),y[train])
            options={
             "tumor_only":model(),
             "tumor_rna":model("rna"),
             "tumor_protein":model("protein"),
             "tumor_metabolome_all":model("metabolome_all"),
             "tumor_metabolome_named":model("metabolome_named"),
             "tumor_metabolome_unresolved":model("metabolome_unresolved"),
             "tumor_rna_metabolome_all":model("rna","metabolome_all"),
             "tumor_rna_metabolome_named":model("rna","metabolome_named"),
             "tumor_rna_protein":model("rna","protein"),
             "tumor_rna_protein_metabolome_all":model("rna","protein","metabolome_all"),
             "tumor_rna_protein_metabolome_named":model("rna","protein","metabolome_named"),
            }
            for panel in CHEM_PANELS:
                options["tumor_path_"+panel]=model(subpath=panel)
                options["tumor_rna_path_"+panel]=model("rna",subpath=panel)
            for name,val in options.items():
                preds[name][repeat,test]=val
    if any(not np.isfinite(v).all() for v in preds.values()):
        raise ValueError("Nonfinite / missing original OOF predictions")
    loss={name:np.abs(mat-y[None,:]).mean(axis=0) for name,mat in preds.items()}
    mae={key:round(float(v.mean()),6) for key,v in loss.items()}
    # Replay prior source-qualified original benchmark to 6-decimal invariance
    prior={"tumor_only":.342488,"tumor_rna":.334666,"tumor_protein":.348451,
           "tumor_metabolome_all":.355840,"tumor_rna_metabolome_all":.351537,
           "tumor_rna_protein":.341053,"tumor_rna_protein_metabolome_all":.358124}
    for name,original in prior.items():
        if abs(mae[name]-original)>0.000005:
            raise ValueError("Original previously validated 3omics metrics failed exact replay: "+name)
    draw=np.random.default_rng(SEED).integers(0,47,size=(1999,47))
    pairs=[
      ("named_minus_all","tumor_metabolome_named","tumor_metabolome_all"),
      ("unresolved_minus_named","tumor_metabolome_unresolved","tumor_metabolome_named"),
      ("rna_named_minus_rna","tumor_rna_metabolome_named","tumor_rna"),
      ("rna_all_minus_rna","tumor_rna_metabolome_all","tumor_rna"),
      ("rna_protein_named_minus_rna_protein","tumor_rna_protein_metabolome_named","tumor_rna_protein"),
      ("rna_protein_all_minus_rna_protein","tumor_rna_protein_metabolome_all","tumor_rna_protein"),
      ("three_named_minus_three_all","tumor_rna_protein_metabolome_named","tumor_rna_protein_metabolome_all"),
    ]
    pairs += [("rna_"+name+"_minus_rna","tumor_rna_path_"+name,"tumor_rna") for name in CHEM_PANELS]
    delta={}
    for label,new,ref in pairs:
        d=loss[new]-loss[ref]
        lo,hi=np.quantile(d[draw].mean(axis=1),[.025,.975])
        delta[label]={"delta_MAE":round(float(d.mean()),6),
                      "descriptive_95pct_paired_heldout_loss_CI":[round(float(lo),6),round(float(hi),6)]}
    return {
      "status":"EXPLORATORY_ORIGINAL_NCI60_47_SOURCE_MATCHED_PATHWAY_QC_SUCCESS_NOT_CLINICAL_PK",
      "n_source_concordant_line_ids":47,"original_metabolite_provenance":REPO_SOURCE,
      "proteome_official_DOI":"10.1016/j.isci.2019.10.059",
      "proteome_original_source_sha256":PIN,
      "metabolite_annotation":meta,
      "source_matched_PCA_Ridge_model_MAE":mae,
      "descriptive_paired_loss_deltas":delta,
      "original_genotype_measurements_available":False,
      "same_biological_specimen_omics":False,
      "original_donor_or_cell_line_outcomes_exported":False,
      "interpretation":"Pathways manually nominated and not KEGG/Reactome-enzyme-validated; name does not prove MSI level 1; no pretreatment timing or direct doxorubicinol/target protein activity measured; fixed post-existing cross-study GI50 benchmark exploratory with multiple pathways; no universal omic ranking.",
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data",required=True,type=Path)
    ap.add_argument("--annotation",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    o=ap.parse_args()
    ids,classes,y,rna_ids,met_ids,rna,met=load_dataset(o.data)
    named=read_source_annotated_metabolites(o.annotation)
    prot_names,prot_classes,prot,sha,genes=load_official_proteome(return_gene_symbols=True)
    result=compare_original(ids,classes,y,rna,met,met_ids,
                            prot_names,prot_classes,prot,genes,named)
    result["source_biochemical_coverage"]={}
    for pathway,p in PATHWAY_PANELS.items():
        result["source_biochemical_coverage"][pathway]={
            "n_source_named_metabolites":len(p["metabolites"]),
            "n_original_rna_genes_source_named":sum(g in rna_ids for g in p["genes"]),
            "n_original_proteome_gene_symbols_source_named":sum(g in genes for g in p["genes"]),
            "requested_gene_symbols_present_in_both":[g for g in p["genes"] if g in rna_ids and g in genes],
            "missing_important_metabolites":list(p["not_measured_key_metabolites"]),
        }
    # On-target drug mechanisms versus metabolic proxies: evidence is the source
    # column lists, NOT proof of functional enzyme/transporter activity.
    direct=("TOP2A","TOP2B","ABCB1","ABCC1","ABCG2","CBR1","CBR3")
    result["direct_doxorubicin_gene_targets_and_metabolizing_transport_candidates"]={
        "predeclared_hypothesis_genes":direct,
        "RNA_measured_exact_gene_symbols":[g for g in direct if g in rna_ids],
        "protein_gene_symbols_present_in_SWATH":[g for g in direct if g in genes],
        "measured_DOX_or_doxorubicinol_as_metabolite_feature":any("doxorubicin" in n.lower() for n in named),
        "actual_TOP2_activity_or_ABCB1_transporter_function_measured":False,
    }
    o.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=="__main__":main()
