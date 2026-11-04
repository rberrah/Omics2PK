"""Layer-specific biological priors vs assay-matched negative controls conditional on somatic mutation.

Public n47 NCI-60 original cross-study RNA, SWATH proteins, named candidate
metabolites, CellMiner deleterious somatic mutations and doxorubicin GI50.
No function, molecule structure, germline diplotype, or human PK inferred.
"""
from __future__ import annotations

import numpy as np
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

from omics2pk.mechanism_matched_control import matched_source_control_panels
from omics2pk.nci60_biochemical_coverage import PATHWAY_PANELS

RNA_PRIOR=("TOP2A","TOP2B","ABCB1","ABCC1","ABCG2","CBR1","CBR3")
PROTEIN_PRIOR=("TOP2B","ABCB1","CBR1")
MET_PRIOR=PATHWAY_PANELS["glutathione_redox"]["metabolites"]
PRIOR={"rna":RNA_PRIOR,"prot":PROTEIN_PRIOR,"met":MET_PRIOR}
SEED=20261104


class MutationMechanismError(ValueError):
    pass


def _subset_prior_indices(all_feature_names, mode):
    names=tuple(all_feature_names)
    if len(set(names))!=len(names) and mode!="prot":
        raise MutationMechanismError("Duplicate source molecular feature columns")
    selected=[]
    for name in PRIOR[mode]:
        locations=[i for i,s in enumerate(names) if s==name]
        if len(locations)!=1:
            raise MutationMechanismError(
                "ABSTAIN: source "+mode+" missing or ambiguous predeclared molecular feature "+name)
        selected.append(locations[0])
    return tuple(selected)


def _predict(tr,te,y,lineage,mutation,assays,indices_by_layer,*,interaction=False):
    encoder=OneHotEncoder(sparse_output=False,handle_unknown="ignore")
    x=encoder.fit_transform(lineage[tr,None]);z=encoder.transform(lineage[te,None])
    msc=StandardScaler().fit(mutation[tr])
    x=np.column_stack((x,msc.transform(mutation[tr])))
    z=np.column_stack((z,msc.transform(mutation[te])))
    groups=mutation[:,0].astype(int)
    for mode,idx in indices_by_layer.items():
        a=assays[mode][tr][:,list(idx)]
        b=assays[mode][te][:,list(idx)]
        imp=SimpleImputer(strategy="median",keep_empty_features=True)
        xa=imp.fit_transform(a);xb=imp.transform(b)
        if xa.shape[1]!=len(idx) or not np.isfinite(xa).all() or not np.isfinite(xb).all():
            raise MutationMechanismError("Incomplete original source molecular prior")
        sc=StandardScaler().fit(xa)
        va=sc.transform(xa);vb=sc.transform(xb)
        x=np.column_stack((x,va));z=np.column_stack((z,vb))
        if interaction:
            prior_rate=groups[tr].mean()
            x=np.column_stack((x,va*(groups[tr,None]-prior_rate)))
            z=np.column_stack((z,vb*(groups[te,None]-prior_rate)))
    pred=np.asarray(Ridge(alpha=10.0).fit(x,y[tr]).predict(z),dtype=float)
    if not np.isfinite(pred).all():raise MutationMechanismError("Invalid true GI50 heldout prediction")
    return pred


def compare_source_priors_by_mutation(ids,lineage,y,mutation,mutation_labels,
                                      source_names,assays,*,repeats=3,controls=32):
    ids=tuple(ids);n=len(ids);y=np.asarray(y,dtype=float)
    line=np.asarray(lineage);mut=np.asarray(mutation,dtype=float)
    if (n<30 or len(set(ids))!=n or len(line)!=n or y.shape!=(n,) or
            not np.isfinite(y).all() or mut.shape!=(n,len(mutation_labels))
            or not mutation_labels or mutation_labels[0]!="mut_TP53"
            or not np.isfinite(mut).all() or not np.isin(mut[:,0],[0,1]).all()
            or tuple(assays)!=("rna","prot","met")
            or tuple(source_names)!=("rna","prot","met")
            or any(np.asarray(assays[k]).shape!=(n,len(source_names[k])) for k in PRIOR)):
        raise MutationMechanismError("Original mutation matched omic source required")
    if (repeats<1 or repeats>4 or controls<2 or controls>100 or
            not isinstance(repeats,int) or not isinstance(controls,int)):
        raise MutationMechanismError("Predeclared repeat and control counts required")
    xs={k:np.asarray(v,dtype=float) for k,v in assays.items()}
    prior={k:_subset_prior_indices(source_names[k],k) for k in PRIOR}
    models=("mutation_only","rna_prior","prot_prior","met_prior",
            "all_three_priors","rna_tp53_interaction","prot_tp53_interaction",
            "met_tp53_interaction","all_three_tp53_interaction")
    predictions={k:np.full((repeats,n),np.nan) for k in models}
    null={k:np.full((controls,repeats,n),np.nan) for k in PRIOR}
    for rep in range(repeats):
        for fold,(tr,te) in enumerate(KFold(5,shuffle=True,random_state=SEED+rep).split(np.arange(n))):
            predictions["mutation_only"][rep,te]=_predict(tr,te,y,line,mut,xs,{})
            for key in PRIOR:
                predictions[key+"_prior"][rep,te]=_predict(tr,te,y,line,mut,xs,{key:prior[key]})
                predictions[key+"_tp53_interaction"][rep,te]=_predict(
                    tr,te,y,line,mut,xs,{key:prior[key]},interaction=True)
                panels=matched_source_control_panels(
                    xs[key][tr],prior[key],n_controls=controls,
                    pool_size=min(35,max(4,xs[key].shape[1]-len(prior[key])-1)),
                    seed=SEED+100*rep+fold)
                for j,idx in enumerate(panels):
                    null[key][j,rep,te]=_predict(tr,te,y,line,mut,xs,{key:idx})
            predictions["all_three_priors"][rep,te]=_predict(tr,te,y,line,mut,xs,prior)
            predictions["all_three_tp53_interaction"][rep,te]=_predict(
                tr,te,y,line,mut,xs,prior,interaction=True)
    if any(not np.isfinite(a).all() for a in
           list(predictions.values())+list(null.values())):
        raise MutationMechanismError("Incomplete original cell-line fold coverage")
    abs_err={k:np.abs(v-y[None,:]) for k,v in predictions.items()}
    mae={k:round(float(err.mean()),6) for k,err in abs_err.items()}
    background=float(abs_err["mutation_only"].mean())
    controls_result={}
    for key,vec in null.items():
        ctrl=np.abs(vec-y[None,None,:]).mean(axis=(1,2))
        lower,median,upper=np.quantile(ctrl,[0.05,0.5,0.95])
        label=key+"_prior"
        controls_result[key]={
            "n_predeclared_measured_source_features":len(prior[key]),
            "n_assay_matched_negative_control_pipelines":controls,
            "prior_MAE":mae[label],
            "prior_minus_mutation_only_delta_MAE":round(mae[label]-background,6),
            "matched_source_random_MAE_q05_q50_q95":[round(float(v),6) for v in (lower,median,upper)],
            "fraction_of_random_matched_panels_with_better_MAE_than_prior":
                round(float(np.mean(ctrl<mae[label])),6),
            "candidate_biological_identity_grade":(
                "source_gene_symbol_not_functional_protein" if key=="rna" else
                "source_SWATH_peptide_quantification_not_functional_protein" if key=="prot"
                else "candidate_original_metabolite_names_not_MSI1_confirmed"),
            "not_an_external_independent_panel_reproduction":True,
        }
    group_results={}
    for g in (0,1):
        mask=mut[:,0]==g
        if mask.sum()<8:
            group_results[str(g)]={"status":"ABSTAIN_TOO_FEW_ORIGINAL_SOMATIC_MUTATION_GROUP_LINES"}
            continue
        ref=float(abs_err["mutation_only"][:,mask].mean())
        group_results[str(g)]={
            "n_distinct_lines":int(mask.sum()),
            "mutation_only_MAE":round(ref,6),
            "prior_delta_MAE_negative_better":{
                key:round(float(err[:,mask].mean()-ref),6)
                for key,err in abs_err.items() if key!="mutation_only"},
            "subgroup_posthoc_descriptive_only":True,
        }
    return {
        "status":"REAL_SOURCE_SOMATIC_MUTATION_CONDITIONAL_RNA_PROTEIN_METABOLITE_PRIOR_NULL_GI50",
        "n_original_biologically_distinct_cell_lines":n,
        "n_repeated_correlated_test_predictions":n*repeats,
        "genetic_effect_interactions_refer_only_to_somatic_TP53_flag":True,
        "somatic_gene_flag_not_complete_germline_PGx_or_function_assay":True,
        "cross_laboratory_omics_not_same_biological_preparation":True,
        "mutations_and_all_assay_preprocessing_train_only":True,
        "original_human_drug_pharmacokinetics_assayed":False,
        "source_measured_layer_specific_predeclared_prior_genes_or_metabolites":
            {key:list(value) for key,value in PRIOR.items()},
        "fixed_models_MAE_original_in_vitro_GI50":mae,
        "matched_random_control_comparison_by_layer":controls_result,
        "exploratory_TP53_stratum_MAE":group_results,
        "source_original_measurements_or_cell_line_predictions_exported":False,
        "confirmatory_mechanistic_or_genotype_effect_claim":False,
    }
