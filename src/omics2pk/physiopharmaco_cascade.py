"""Scientific layered PK mechanism graph: model signal != verified drug mechanism.

NetPK-style levels:
 original assay feature -> identity (only if independently confirmed) ->
 enzyme/transporter/pathway -> PK parameter -> observed PK endpoint.
STRING v12.5 associations (functional, physical, directed regulatory)
are separate evidence edges; never interpreted as drug transport or causation.
Unassigned GC-TOF bins remain original source nodes, not masses/PTMs.
"""
from __future__ import annotations
from math import isfinite
from urllib.parse import quote

STRING_DOC="https://string-db.org/help/api/"
STRING_VERSION="12.5"
MODEL_ATTRIBUTION="model_predictive_association_not_causal"
PK_SOURCES={
 "metformin":"https://pmc.ncbi.nlm.nih.gov/articles/PMC5709222/",
 "tacrolimus":"https://pmc.ncbi.nlm.nih.gov/articles/PMC4119065/",
}
KNOWN={
 "metformin":{
   "transporters":{
      "SLC22A1":("hepatic_uptake","liver","distribution_and_tissue_exposure"),
      "SLC22A2":("renal_uptake","kidney","renal_clearance"),
      "SLC47A1":("renal_efflux","kidney","renal_clearance"),
      "SLC47A2":("renal_efflux","kidney","renal_clearance"),
      "SLC29A4":("intestinal_uptake_candidate","intestine","bioavailability"),
   },
   "endpoints":("AUClast_0_24","CLr"),
   "transport_claims_not_reduced_to_Pgp":True,
 },
 "tacrolimus":{
   "transporters":{
      "ABCB1":("intestinal_efflux","intestine","bioavailability"),
      "CYP3A4":("metabolism","liver_and_intestine","clearance"),
      "CYP3A5":("metabolism","liver_and_intestine","clearance"),
   },
   "endpoints":("C0","AUC","CL_over_F"),
   "erythrocyte_partition_required":True,
 },
}
ORIGINS={"metabolomics","transcriptomics","proteomics","proteoform"}
CHEMICAL_GRADES={"source_name_only","putative_ms2","confirmed_standard"}

class PKCascadeError(ValueError):pass

def _feature(rec):
    if not isinstance(rec,dict):
        raise PKCascadeError("Original source evidence must be mapping")
    source=rec.get("source_id")
    omic=rec.get("omic")
    if not isinstance(source,str) or not 0<len(source)<=128 or omic not in ORIGINS:
        raise PKCascadeError("Original measured feature ID and layer required")
    if rec.get("is_original_measured") is not True:
        raise PKCascadeError("Cannot build measured feature from inferred protein/chemical")
    if omic=="metabolomics":
        grade=rec.get("chemical_grade","source_name_only")
        if grade not in CHEMICAL_GRADES:
            raise PKCascadeError("Original metabolite chemical identity grade required")
        chemical_id=rec.get("chemical_id")
        if chemical_id is not None and grade=="source_name_only":
            raise PKCascadeError("No chemical ID auto-assigned to original unknown GC-TOF bin")
        if grade=="confirmed_standard" and not chemical_id:
            raise PKCascadeError("Chemical confirmation requires source verified molecule ID")
    else:
        if rec.get("chemical_id"):
            raise PKCascadeError("A protein/RNA cannot be a small-molecule ChEBI object")
    return source,omic

def build_layered_pk_network(drug, measured_features, *,
                             molecular_pk_endpoint=None,string_interactions=()):
    """Pure offline graph with separate research/canonical/literature evidence.

    Input to this function is aggregate assay feature identities only; never
    require subject intensity matrix, clinical PK or patient metadata.
    """
    if drug not in KNOWN:
        raise PKCascadeError("No verified physiological PK reference for drug")
    if not isinstance(measured_features,(list,tuple)) or len(measured_features)>10000:
        raise PKCascadeError("List of measured source signals required")
    source_evidence=PK_SOURCES[drug]
    known=KNOWN[drug]
    nodes={}
    edges=[]
    def add_node(kind,key,**metadata):
        id=kind+":"+str(key)
        if id not in nodes:
            nodes[id]={"id":id,"layer":kind,"label":str(key),**metadata}
        return id
    def edge(src,dst,relation,level,source):
        if src not in nodes or dst not in nodes:
            raise PKCascadeError("Edges need verified source nodes")
        edges.append({"from":src,"to":dst,"relation":relation,"evidence_level":level,
                      "source_url":source,
                      "causal_or_functional_effect_demonstrated_in_this_cohort":False})
    drug_id=add_node("drug",drug,identity="original_drug_name")
    endpoints=known["endpoints"]
    if molecular_pk_endpoint is not None and molecular_pk_endpoint not in endpoints:
        raise PKCascadeError("PK endpoint must be a drug-specific original measured target")
    for endpoint in endpoints:
        end=add_node("pk_endpoint",endpoint,
                     observed_on_selected_source_patients=False)
        edge(drug_id,end,"drug_PK_endpoint_of_interest",
             "conceptual_endpoint_not_original_data",source_evidence)
    gene_nodes={}
    pk_nodes={}
    for gene,(role,tissue,pk_effect) in known["transporters"].items():
        g=add_node("gene_or_protein",gene,
                   biological_role_hypothesis=role,organ=tissue,
                   individual_function_assayed=False)
        pk=pk_nodes.setdefault(pk_effect,add_node("pk_process",pk_effect,organ=tissue))
        edge(g,pk,"literature_drug_specific_pathway",
             "literature_mechanistic_prior_not_patient_evidence",source_evidence)
        gene_nodes[gene]=g
    measured_gene_symbols=set()
    for record in measured_features:
        source,omic=_feature(record)
        ident=add_node("measured_"+omic,source,
                       assay_matrix=record.get("assay_matrix","source_not_specified"),
                       is_measured_feature=True)
        endpoint=molecular_pk_endpoint
        if record.get("has_model_attribution") is True:
            if endpoint is None:
                raise PKCascadeError("Model attribution requires actually observed PK endpoint label")
            edge(ident,"pk_endpoint:"+endpoint,MODEL_ATTRIBUTION,
                 "model_association_not_chemistry",record.get("model_run_url",""))
        if omic=="metabolomics":
            identity=record.get("chemical_id")
            if identity:
                chem=add_node("chemical_candidate",identity,
                              structure_source_grade=record["chemical_grade"])
                edge(ident,chem,"original_assay_to_structural_candidate",
                     record["chemical_grade"],record.get("chemical_source_url",""))
            # no edge to transporter without independently authenticated
            # structure plus source-functional substrate relation.
        else:
            symbol=record.get("gene_symbol")
            if symbol:
                if not isinstance(symbol,str) or not symbol.isupper() or len(symbol)>25:
                    raise PKCascadeError("Gene symbol must be source-qualified")
                measured_gene_symbols.add(symbol)
                protein=add_node("gene_or_protein",symbol,
                                 individual_function_assayed=False)
                edge(ident,protein,"original_assay_mapping_to_gene",
                     "measured_expression_or_protein_not_activity",
                     record.get("assay_provenance_url",""))
    if not isinstance(string_interactions,(list,tuple)):
        raise PKCascadeError("STRING graph edges must be a list")
    for e in string_interactions:
        if not isinstance(e,dict) or e.get("source") not in measured_gene_symbols or e.get("target") not in measured_gene_symbols:
            raise PKCascadeError("STRING can link only originally measured source genes")
        score=e.get("score")
        if isinstance(score,bool) or not isinstance(score,(int,float)) or not isfinite(score) or not 0<=score<=1:
            raise PKCascadeError("STRING confidence score must be [0,1]")
        if e.get("network_type") not in ("functional","physical","regulatory"):
            raise PKCascadeError("STRING functional association / physical / regulatory type missing")
        edge("gene_or_protein:"+e["source"],"gene_or_protein:"+e["target"],
             "STRING_"+e["network_type"]+"_association_confidence_"+str(round(score,4)),
             "STRING_database_prior_not_PK_causality",STRING_DOC)
    return {
        "status":"SOURCE_TYPED_PK_CASCADE_MECHANISTIC_HYPOTHESES_ONLY",
        "drug":drug,
        "source_physiology":source_evidence,
        "STRING_reference_version":STRING_VERSION,
        "nodes":list(nodes.values()),
        "edges":edges,
        "n_original_measured_source_features":len(measured_features),
        "n_automatically_identified_unknown_metabolites":0,
        "n_validated_omic_to_mechanistic_PK_links":0,
        "no_clinical_function_inferred_from_STRING_SHAP_or_feature_correlations":True,
        "PK_benefit_requires_actual_external_PK_validation":True,
    }
