"""Conservative source-specific annotation and hypothesis panels for NCI-60.

Original 280 metabolite labels: Mathelab/NCI60_GeneMetabolite_Data
@4c140708c18c41910458414a626b1a51d3dd2a85, fData.metab.csv.

Chemical names are *not* automatically high-confidence MS/MS identifications:
the source has no mapped InChIKey/HMDB/ChEBI, reference spectra or authentic
standards here. A named label is a candidate mapping, never level-1 proof.

Pathway panels below are PRESPECIFIED biochemical *hypotheses* only.
They are not a computationally certified KEGG/Reactome gene–reaction–metabolite
mapping; the panel membership is explicit and auditable. Do not conflate
doxorubicin treatment-derived analytes with endogenous pre-treatment assays.
"""
from __future__ import annotations
from collections import Counter

PATHWAY_PANELS={
 "glutathione_redox": {
  "metabolites":("cysteine","glutamic acid","glycine","gamma-glu-cys",
                 "glutathione, reduced","5-oxoproline"),
  "genes":("GCLC","GCLM","GSS","GSR","GPX1","GPX4","GSTP1","SLC7A11","ABCC1","NFE2L2"),
  "not_measured_key_metabolites":("glutathione, oxidized (GSSG)","NADPH","NADP+"),
  "limitation":"GSH abundance alone cannot infer GSH/GSSG redox state or flux",
 },
 "central_carbon_glycolysis_ppp": {
  "metabolites":("glucose-6-phosphate","3-phospho-d-glycerate",
                 "phosphoenolpyruvate","lactate","alpha-D-ribose 5-phosphate",
                 "sn-Glycerol 3-phosphate"),
  "genes":("HK1","HK2","G6PD","PGD","TALDO1","TKT","LDHA","PKM","PFKP"),
  "not_measured_key_metabolites":("6-phosphogluconate","NADPH","NADP+"),
  "limitation":"Levels without isotopic flux cannot identify PPP flux or NADPH regeneration",
 },
 "tca_energy": {
  "metabolites":("citric acid","alpha-keto-glutarate","succinate",
                 "fumaric acid","malic acid"),
  "genes":("CS","IDH1","IDH2","OGDH","SDHA","SDHB","FH","MDH1","MDH2"),
  "not_measured_key_metabolites":("acetyl-CoA","isocitrate"),
  "limitation":"Static pool sizes do not establish flux or mitochondrial respiration",
 },
 "nucleotide_pools": {
  "metabolites":("adenine","guanine","hypoxanthine","inosine",
                 "guanosine","guanosine 5'- monophosphate","uracil","uridine"),
  "genes":("TOP2A","TOP2B","RRM1","RRM2","TYMS","IMPDH1","IMPDH2","PNP"),
  "not_measured_key_metabolites":("dATP","dGTP","dCTP"),
  "limitation":"Doxorubicin DNA/topoisomerase injury is not directly measured by these pools",
 },
 "one_carbon_methylation": {
  "metabolites":("methionine","S-(5'-adenosyl)-L-methionine",
                 "S-adenosyl-l-homocysteine","serine","glycine",
                 "5'-s-methyl-5'-thioadenosine"),
  "genes":("MAT2A","AHCY","SHMT1","SHMT2","MTHFR","MTR","MTRR"),
  "not_measured_key_metabolites":("5-methyl-tetrahydrofolate","tetrahydrofolate"),
  "limitation":"SAM/SAH concentration ratio is not enzymatic methylation flux",
 },
 "polyamine_metabolism": {
  "metabolites":("ornithine","Putrescine","Spermidine",
                 "5'-s-methyl-5'-thioadenosine"),
  "genes":("ODC1","AMD1","SRM","SMS","SAT1"),
  "not_measured_key_metabolites":("spermine",),
  "limitation":"Polyamine turnover/transport cannot be inferred from a partial static pool",
 },
}

def status_of_original_label(name: str)->str:
    if not isinstance(name,str) or not name.strip():
        raise ValueError("Missing original source metabolite label")
    s=name.strip()
    if s.lower().startswith("isobar "):
        return "unresolved_isobar_group"
    if s.lower().startswith("x-") and s[2:3].isdigit():
        return "unidentified_X_feature"
    if "artifact" in s.lower() or "retired" in s.lower():
        return "possible_artifact_or_retired"
    return "named_candidate_identity_unverified"

def audit_labels(original_labels):
    values=list(original_labels)
    if len(values)!=280 or len(set(values))!=280:
        raise ValueError("Require exactly 280 unique original IntLIM/NCI60 source labels")
    classes=Counter(status_of_original_label(x) for x in values)
    named={x for x in values if status_of_original_label(x)=="named_candidate_identity_unverified"}
    if len(named)!=111:
        raise ValueError("Original source candidate-metabolite annotation count drift: "+str(len(named)))
    # Each proposed metabolite must be the literal original source label.
    panels={}
    for key,p in PATHWAY_PANELS.items():
        metabolites=p["metabolites"]
        if len(metabolites)<3 or len(metabolites)!=len(set(metabolites)):
            raise ValueError("Malformed explicit biochemical hypothesis panel "+key)
        missing=[x for x in metabolites if x not in named]
        if missing:
            raise ValueError("Panel "+key+" references unknown, ambiguous or absent labels: "+repr(missing))
        expected_absent=[x for x in p["not_measured_key_metabolites"] if x in values]
        if expected_absent:
            raise ValueError("Panel absent-analyte assumption drift: "+repr(expected_absent))
        panels[key]={"n_source_named_metabolites":len(metabolites),
                     "source_original_metabolite_labels":list(metabolites),
                     "genes_or_proteins_hypothesized":list(p["genes"]),
                     "important_unobserved_analytes":list(p["not_measured_key_metabolites"]),
                     "not_claimed_pathway_completeness":True}
    result={
     "source_total_features":len(values),
     "source_annotation_categories":dict(classes),
     "candidate_named_features_not_structurally_confirmed":len(named),
     "source_has_HMDB_ChEBI_InChIKey_or_MSMS_confidence_in_this_metadata":False,
     "source_unknown_isobar_artifact_count":len(values)-len(named),
     "source_original_feature_names_used_as_identity_keys_only":True,
     "hypothesis_panels_not_Reactome_or_KEGG_verified":panels,
    }
    return result
