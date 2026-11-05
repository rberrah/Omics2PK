"""Auditable mutation → RNA → protein → metabolite evidence, NOT a causal path.

Same cell-line IDs across separate lab studies are not same-date
same-sample matched multiomics; CellMiner presence >0 is not genotype dosage,
loss of function, or pharmacokinetically functional enzyme activity.
Named metabolites in source are not authenticated MSI-level-1 identities.
"""
from __future__ import annotations
import numpy as np
from omics2pk.nci60_biochemical_coverage import PATHWAY_PANELS
from omics2pk.predictive_diagnostics import _average_tie_ranks, _corr

DIRECT=("TOP2A","TOP2B","ABCB1","ABCC1","ABCG2","CBR1","CBR3")
METS=PATHWAY_PANELS["glutathione_redox"]["metabolites"]


class CascadeAuditError(ValueError):
    pass


def _rankspearman(x,y):
    # Biochemical assay comparisons cannot be computed from non-measured
    # proteins or single-gene mutation events with zero variance.
    return _corr(_average_tie_ranks(x),_average_tie_ranks(y))


def _r(x):
    return None if x is None else round(float(x),6)


def audit_original_measured_cascade(ids, y, mutations, mutation_labels,
                                    source_names, assays):
    """Observed feature-level source coverage + descriptive empirical coherence.

    Pre-specified 7 doxorubicin genes and the original 6 redox source names.
    No posthoc feature selection, borrowed proteins, abundance/activity
    equivalence, or interpretation of correlations as causal effects.
    """
    ids=tuple(ids)
    y=np.asarray(y,dtype=float)
    mut=np.asarray(mutations,dtype=float)
    labels=tuple(mutation_labels)
    n=len(ids)
    if (n<25 or len(set(ids))!=n or y.shape!=(n,) or
            not np.isfinite(y).all() or mut.shape!=(n,len(labels))
            or not np.isfinite(mut).all() or
            not labels or labels[0]!="mut_TP53" or
            tuple(source_names)!=("rna","prot","met") or
            tuple(assays)!=("rna","prot","met") or
            any(np.asarray(assays[k]).shape!=(n,len(source_names[k]))
                for k in ("rna","prot","met"))):
        raise CascadeAuditError("Original same cell-line source and TP53 mutation status required")
    tp53=mut[:,0]
    if not np.isin(tp53,[0,1]).all():
        raise CascadeAuditError("Original TP53 status is not an observed event presence")
    if np.sum(tp53==0)<8 or np.sum(tp53==1)<8:
        raise CascadeAuditError("Insufficient TP53 source-stratum representation")
    gene_results={}
    for gene in DIRECT:
        found={}
        for mode in ("rna","prot"):
            indices=[i for i,g in enumerate(source_names[mode]) if g==gene]
            if len(indices)>1:
                raise CascadeAuditError("Ambiguous original measured gene symbol: "+gene)
            found[mode]=indices[0] if indices else None
        m="mut_"+gene
        mutation_index=labels.index(m) if m in labels else None
        result={
            "RNA_measured_in_original_source":found["rna"] is not None,
            "SWATH_protein_measured_in_original_source":found["prot"] is not None,
            "mutation_event_gene_assayed_in_CellMiner":mutation_index is not None,
            "mutation_allele_dosage_or_function_known":False,
            "original_RNA_isoform_or_protein_active_site_verified":False,
            "pathway_drug_substrate_or_transport_function_measured":False,
        }
        if mutation_index is not None:
            v=mut[:,mutation_index]
            if not np.isin(v,[0,1]).all():
                raise CascadeAuditError("Unqualified original source mutation binary value")
            mutated=int(np.sum(v==1))
            result["n_original_line_mutation_event_present"]=mutated
            result["mutation_event_subgroup_comparison_status"]=(
                "DESCRIPTIVE_ONLY_NOT_GERMLINE_VARIANT_EFFECT" if 8<=mutated<=n-8
                else "ABSTAIN_INSUFFICIENT_MUTATION_CARRIERS_FOR_SUBGROUP")
        if found["rna"] is not None:
            rna=np.asarray(assays["rna"][:,found["rna"]],dtype=float)
            clean=np.isfinite(rna)
            result["RNA_n_source_observed"]=int(clean.sum())
            result["RNA_vs_GI50_spearman_cross_study_exploratory"]=_r(
                _rankspearman(rna[clean],y[clean]) if clean.sum()>=10 else None)
        if found["prot"] is not None:
            protein=np.asarray(assays["prot"][:,found["prot"]],dtype=float)
            clean=np.isfinite(protein)
            result["protein_n_source_observed"]=int(clean.sum())
            result["protein_vs_GI50_spearman_cross_study_exploratory"]=_r(
                _rankspearman(protein[clean],y[clean]) if clean.sum()>=10 else None)
        if found["rna"] is not None and found["prot"] is not None:
            rna=np.asarray(assays["rna"][:,found["rna"]],dtype=float)
            protein=np.asarray(assays["prot"][:,found["prot"]],dtype=float)
            clean=np.isfinite(rna)&np.isfinite(protein)
            result["n_original_RNA_and_SWATH_both_observed"]=int(clean.sum())
            if clean.sum()>=10:
                result["RNA_vs_SWATH_abundance_spearman_all_source_matched_lines"]=_r(
                    _rankspearman(rna[clean],protein[clean]))
                result["RNA_vs_SWATH_abundance_spearman_TP53_strata"]={
                    str(int(g)):_r(_rankspearman(rna[clean & (tp53==g)],
                                                protein[clean & (tp53==g)])
                                    if np.sum(clean & (tp53==g))>=8 else None)
                    for g in (0,1)}
                result["RNA_vs_SWATH_cross_lab_correlation_NOT_functional_protein"]=True
            else:
                result["RNA_vs_SWATH_status"]="ABSTAIN_FEWER_THAN_10_COOBSERVED_ORIGINAL_LINES"
        gene_results[gene]=result
    original_met=set(source_names["met"])
    present=sorted(set(METS)&original_met)
    if len(present)!=len(METS):
        raise CascadeAuditError("A predeclared original candidate redox source label missing")
    return {
        "status":"SOURCE_VERIFIED_OBSERVABLE_CASCADE_COVERAGE_NOT_CAUSAL_PATHWAY",
        "n_biologically_distinct_cell_lines":n,
        "TP53_mutation_signal_n_present":int(np.sum(tp53==1)),
        "biological_interpretation_nodes":gene_results,
        "redox_metabolite_source_named_candidates":list(METS),
        "n_redox_metabolites_source_named":len(METS),
        "GSSG_NADPH_redox_flux_or_drug_metabolite_assay_completed":False,
        "source_assays_same_time_same_biological_aliquot":False,
        "original_drug_transporter_flux_or_active_protein_assayed":False,
        "source_metabolites_MSI_level1_structure_verified":False,
        "direct_transcript_to_active_protein_or_metabolite_mechanism_proven":False,
        "reported_correlations_not_out_of_fold_predictive_attributions":True,
        "patient_or_cell_line_source_measurements_exported":False,
    }
