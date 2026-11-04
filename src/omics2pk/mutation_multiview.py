"""Mutation-aware NCI-60 three-omic view comparison, original-source research.

NCI60 mutations: CellMiner gene-level potentially deleterious *somatic* calls,
not complete germline CYP3A5 genotype, nor mechanistic knockout evidence.
Doxorubicin GI50 is cell-line pharmacodynamic response, not human PK.
Omics arise from distinct laboratory studies/preparations.
"""
from __future__ import annotations

from collections import Counter
import numpy as np
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from omics2pk.modalities_information_budget import choose_minimum_modalities

SEED=20261104
CAPS={"rna":300,"prot":300,"met":120}
MODES={
    "mutation_only":(),
    "mutation_rna":("rna",),
    "mutation_protein":("prot",),
    "mutation_metabolite":("met",),
    "mutation_rna_protein":("rna","prot"),
    "mutation_rna_metabolite":("rna","met"),
    "mutation_protein_metabolite":("prot","met"),
    "mutation_threeomics":("rna","prot","met"),
}
PAIRED_INTERACTIONS=("mutation_rna","mutation_protein","mutation_metabolite","mutation_threeomics")


class MutationMultiViewError(ValueError):
    pass


def _view(train, test, cap):
    """Original assay preprocessing is fitted exclusively on TRAIN lines."""
    train=np.asarray(train,dtype=float);test=np.asarray(test,dtype=float)
    if (train.ndim!=2 or test.ndim!=2 or len(train)<10 or
            train.shape[1]!=test.shape[1] or np.isinf(train).any() or np.isinf(test).any()):
        raise MutationMultiViewError("Unqualified measured molecular source view")
    coverage=np.isfinite(train).mean(axis=0)
    good=np.flatnonzero(coverage>=0.5)
    if len(good)<3:
        raise MutationMultiViewError("Insufficient measured train-visible omics")
    x=train[:,good];v=test[:,good]
    imputer=SimpleImputer(strategy="median",keep_empty_features=True)
    xx=imputer.fit_transform(x);vv=imputer.transform(v)
    variances=np.var(xx,axis=0)
    selected=np.argsort(-variances,kind="stable")
    selected=selected[variances[selected]>1e-12][:cap]
    if len(selected)<3:raise MutationMultiViewError("No source-assayed varying molecular panel")
    scaler=StandardScaler()
    a=scaler.fit_transform(xx[:,selected])
    b=scaler.transform(vv[:,selected])
    # Same PCA3 per original molecular modality; train-only.
    pca=PCA(n_components=3,svd_solver="full")
    return pca.fit_transform(a),pca.transform(b)


def _build_fold(y, lineage, mutation, assays, train, test, *, interactions=False, include_mutation=True):
    enc=OneHotEncoder(handle_unknown="ignore",sparse_output=False)
    cl_tr=enc.fit_transform(np.asarray(lineage)[train,None])
    cl_te=enc.transform(np.asarray(lineage)[test,None])
    # No unseen mutation calls: source 0/1 flags plus mutation burden.
    if include_mutation:
        s=StandardScaler().fit(mutation[train])
        base_tr=np.column_stack((cl_tr,s.transform(mutation[train])))
        base_te=np.column_stack((cl_te,s.transform(mutation[test])))
    else:
        base_tr=cl_tr
        base_te=cl_te
    groups=mutation[:,0].astype(int) # original predeclared TP53 column
    views={name:_view(values[train],values[test],CAPS[name])
           for name,values in assays.items()}
    out={}
    for name,modes in MODES.items():
        tr=[base_tr];te=[base_te]
        for mode in modes:
            tr.append(views[mode][0]);te.append(views[mode][1])
            if interactions and name in PAIRED_INTERACTIONS:
                # Add variant×omic PCs; mutation status itself is already in base.
                # No outcome/heldout-based subgroup view switching.
                tr.append(views[mode][0]*(groups[train,None]-groups[train].mean()))
                te.append(views[mode][1]*(groups[test,None]-groups[train].mean()))
        tr=np.column_stack(tr);te=np.column_stack(te)
        pred=Ridge(alpha=10.0).fit(tr,y[train]).predict(te)
        if not np.isfinite(pred).all():
            raise MutationMultiViewError("Nonfinite true response heldout prediction")
        out[name]=pred
    return out


def _errors(pred, truth):
    err=np.abs(pred-truth[None,:])
    return {"mae":round(float(err.mean()),6),
            "rmse":round(float(np.sqrt(np.mean((pred-truth[None,:])**2))),6)}


def evaluate_mutation_modality_budget(ids, lineage, response, mutation_labels,
                                      mutation, assays, *,
                                      repeats=3, inner_splits=3,
                                      n_controls=0):
    """Fixed and nested selected omic views, with and without TP53×omic PCs.

    The outer cell line never influences the choice of omic or fold transforms.
    Selection can legitimately return zero additional omic measurements.
    """
    ids=tuple(ids); lineage=np.asarray(lineage); y=np.asarray(response,dtype=float)
    mut=np.asarray(mutation,dtype=float); labels=tuple(mutation_labels)
    n=len(ids)
    if (n<30 or len(set(ids))!=n or len(lineage)!=n or y.shape!=(n,) or
            not np.isfinite(y).all() or not all(isinstance(g,str) and g for g in lineage) or
            mut.shape!=(n,len(labels)) or len(set(labels))!=len(labels) or
            not labels or labels[0]!="mut_TP53" or
            not np.isfinite(mut).all() or
            not np.isin(mut[:,0],[0.,1.]).all() or
            not np.all(np.isin(mut[:,:-1],[0.,1.])) or
            not np.all(mut[:,-1]>=0) or
            tuple(assays)!=("rna","prot","met") or
            any(np.asarray(v).ndim!=2 or len(v)!=n for v in assays.values())):
        raise MutationMultiViewError("Original same-source mutation and 3-view assays required")
    if not isinstance(repeats,int) or isinstance(repeats,bool) or not 1<=repeats<=4:
        raise MutationMultiViewError("Replicate count must be predeclared")
    if inner_splits!=3:raise MutationMultiViewError("Three train-only inner folds fixed")
    if n_controls!=0:raise MutationMultiViewError("No fabricated patient mutation control panels")
    if n-np.sum(mut[:,0])<8 or np.sum(mut[:,0])<8:
        raise MutationMultiViewError("ABSTAIN TP53 genotype strata too small for comparison")
    views={key:np.asarray(v,dtype=float) for key,v in assays.items()}
    fixed={name:np.full((repeats,n),np.nan) for name in MODES}
    fixed_without_mutations={name:np.full((repeats,n),np.nan) for name in MODES}
    interacts={name:np.full((repeats,n),np.nan) for name in PAIRED_INTERACTIONS}
    adaptive={q:np.full((repeats,n),np.nan) for q in (0.8,0.9)}
    choices={q:{0:Counter(),1:Counter()} for q in (0.8,0.9)}
    fallback_counts=Counter()
    for rep in range(repeats):
        outer=KFold(n_splits=5,shuffle=True,random_state=SEED+rep)
        for fold,(tr,te) in enumerate(outer.split(np.arange(n))):
            inner_pred={name:np.full(len(tr),np.nan) for name in MODES}
            inner=KFold(n_splits=inner_splits,shuffle=True,random_state=SEED+rep*100+fold)
            for it,iv in inner.split(tr):
                preds=_build_fold(y,lineage,mut,views,tr[it],tr[iv])
                for name in MODES:inner_pred[name][iv]=preds[name]
            if any(not np.isfinite(v).all() for v in inner_pred.values()):
                raise MutationMultiViewError("Incomplete mutation-aware inner CV")
            chosen={}
            for q in (0.8,0.9):
                eligible={}
                for g in (0,1):
                    mask=(mut[tr,0]==g)
                    if np.sum(mask)<8:
                        fallback_counts[(q,g)]+=1
                        mask=np.full(len(tr),True)
                    losses={arm:float(np.mean(np.abs(vec[mask]-y[tr][mask])))
                            for arm,vec in inner_pred.items()}
                    eligible[g]=choose_minimum_modalities(
                        losses,MODES,fraction=q,reference="mutation_only")["arm"]
                    choices[q][g][eligible[g]]+=1
                chosen[q]=eligible
            outer_preds=_build_fold(y,lineage,mut,views,tr,te)
            no_mutation_preds=_build_fold(y,lineage,mut,views,tr,te,include_mutation=False)
            int_preds=_build_fold(y,lineage,mut,views,tr,te,interactions=True)
            for arm in MODES:
                fixed[arm][rep,te]=outer_preds[arm]
                fixed_without_mutations[arm][rep,te]=no_mutation_preds[arm]
            for arm in PAIRED_INTERACTIONS:interacts[arm][rep,te]=int_preds[arm]
            for q in (0.8,0.9):
                for ix in te:
                    adaptive[q][rep,ix]=outer_preds[chosen[q][int(mut[ix,0])]][np.flatnonzero(te==ix)[0]]
    if any(not np.isfinite(v).all() for v in
           list(fixed.values())+list(fixed_without_mutations.values())+list(interacts.values())+list(adaptive.values())):
        raise MutationMultiViewError("Incomplete original line-disjoint heldout runs")
    obs={"fixed":{arm:_errors(v,y) for arm,v in fixed.items()},
         "fixed_without_original_somatic_mutation_features":{
             arm:_errors(v,y) for arm,v in fixed_without_mutations.items()},
         "observed_mutation_incremental_delta_MAE_negative_is_better":{
             arm:round(float(np.abs(fixed[arm]-y[None,:]).mean()-
                             np.abs(fixed_without_mutations[arm]-y[None,:]).mean()),6)
             for arm in MODES},
         "TP53_interaction":{arm:_errors(v,y) for arm,v in interacts.items()},
         "nested_TP53_guided_80_90":{}}
    for q,v in adaptive.items():
        obs["nested_TP53_guided_80_90"][str(int(q*100))]={
            **_errors(v,y),
            "selected_by_training_TP53_stratum":{
                str(g):dict(choices[q][g]) for g in (0,1)},
            "low_train_stratum_global_fallback_outer_folds":{
                str(g):fallback_counts[(q,g)] for g in (0,1)}}
    strata={}
    for group in (0,1):
        sel=mut[:,0]==group
        if sel.sum()<8:continue
        reference=np.abs(fixed["mutation_only"][:,sel]-y[sel][None,:]).mean()
        strata[str(group)]={
            "n_unique_original_lines":int(sel.sum()),
            "mutation_only_MAE":round(float(reference),6),
            "fixed_modalities_delta_MAE_negative_is_better":{
                arm:round(float(np.abs(v[:,sel]-y[sel][None,:]).mean()-reference),6)
                for arm,v in fixed.items() if arm!="mutation_only"},
            "subgroup_results_exploratory_not_external_replicates":True}
    # Additional *prespecified measured mutation* strata, never an
    # automatically optimized genotype-specific omic selection. A rare
    # carrier group is a source-qualified abstention rather than a fake
    # mutation interaction. These are repeated, descriptive errors only.
    by_mutation={}
    for j,label in enumerate(labels[:-1]):
        carriers=int(np.sum(mut[:,j]==1))
        noncarriers=int(np.sum(mut[:,j]==0))
        if carriers<8 or noncarriers<8:
            by_mutation[label]={
                "status":"ABSTAIN_TOO_FEW_MUTATION_CARRIERS_OR_NONCARRIERS",
                "n_original_mutated":carriers,"n_original_nonmutated":noncarriers}
            continue
        per_group={}
        for g in (0,1):
            mask=mut[:,j]==g
            baseline=float(np.abs(fixed["mutation_only"][:,mask]-y[mask][None,:]).mean())
            per_group[str(g)]={
                "n_original_distinct_lines":int(mask.sum()),
                "mutation_only_MAE":round(baseline,6),
                "fixed_modality_delta_MAE_negative_better":{
                    arm:round(float(np.abs(vec[:,mask]-y[mask][None,:]).mean()-baseline),6)
                    for arm,vec in fixed.items() if arm!="mutation_only"}}
        by_mutation[label]={"status":"DESCRIPTIVE_SUBGROUP_NOT_MUTATION_EFFECT_TEST",
                            "conditional_original_test_errors":per_group}
    return {
        "source_observed_mutation_gene_subgroups":by_mutation,
        "status":"SOURCE_MUTATION_AWARE_THREEOMIC_CELL_LINE_GI50_RESEARCH_NOT_PK",
        "n_unique_biological_cell_lines":n,
        "n_repeated_not_independent_test_predictions":n*repeats,
        "mutation_source":"rcellminerData CellMiner gene-level predicted-deleterious somatic binary calls, NOT complete germline CYP3A5 diplotypes",
        "n_measured_mutation_variables":len(labels),
        "TP53_group_predeclared":True,
        "TP53_allele_dosage_or_protein_function_not_known":True,
        "imputation_variance_PCA_and_selection_train_only":True,
        "genotype_subgroup_model_chosen_using_outer_response":False,
        "no_patient_or_original_cell_line_feature_rows_exported":True,
        "same_contemporaneous_multiomics_specimens":False,
        "clinical_human_PK_or_CYP3A5_genotype_outcome_demonstrated":False,
        "model_fixed_Ridge_alpha":10.0,
        "original_in_vitro_GI50_observed_holdout":obs,
        "TP53_somatic_variant_subgroup_descriptive":strata,
        "no_externally_replicated_modality_by_genotype_interaction":True,
    }
