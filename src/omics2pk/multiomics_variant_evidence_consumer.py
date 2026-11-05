"""Read-only optional MultiOmics variant evidence sidecar research consumer.

V1 MultiOmics downstream exchange is UNCHANGED. The sidecar is
a *qualitative annotation*; it cannot be silently promoted to
a numerical per-patient PK feature, missing genotype, or activity.
This adapter does not guess drug/organ from a gene name alone.
"""
from __future__ import annotations
from collections.abc import Mapping

SCHEMA="variant_omic_evidence_sidecar_v1"
READOUTS={
    "total_RNA":("transcriptomics","total_gene"),
    "canonical_splice_RNA":("transcriptomics","canonical_splice"),
    "aberrant_splice_RNA":("transcriptomics","aberrant_splice"),
    "protein_total":("proteomics","total_protein"),
    "protein_full_length":("proteomics","full_length"),
    "protein_variant_specific":("proteomics","variant_specific"),
    "endogenous_metabolite":("metabolomics","metabolite_source_peak_or_identity"),
    "measured_substrate_activity":("functional_activity","substrate_probe"),
}
HYPOTHESES={
    "HYPOTHESIS_POSTTRANSCRIPTIONAL_TRANSLATION_OR_PROTEIN_INSTABILITY",
    "HYPOTHESIS_CATALYTIC_OR_LOCALIZATION_OR_SUBSTRATE_BINDING_DEFECT",
    "HYPOTHESIS_ABERRANT_SPLICING_DESPITE_TOTAL_RNA",
    "HYPOTHESIS_CANONICAL_TO_ABERRANT_SPLICE_SHIFT",
    "HYPOTHESIS_PRESERVED_ABUNDANCE_LOST_MISSENSE_FUNCTION",
    "HYPOTHESIS_NMD_ESCAPE_OR_ALTERNATIVE_ISOFORM_POSSIBLE",
    "UNRESOLVED_SPLICE_CONSEQUENCE_REQUIRES_ISOFORM_ASSAY",
    "NO_UNIQUE_VARIANT_MECHANISM_ESTABLISHED_FROM_AVAILABLE_READOUTS",
}
ALLOWED_STATES={"reference_like","decreased","increased","absent","unknown"}
class OmicsVariantEvidenceHandoffError(ValueError):
    pass


def qualify_variant_evidence_sidecar(
        report,*, drug,endpoint,orig_PGx_pk_join_verified=False):
    """Interpret independent original assay hints, never invent PK gains.

    Input report from MultiOmics research-only
    reconcile_variant_omic_evidence(), not a V1 exchange manifest.
    """
    from .physiopharmaco_cascade import KNOWN
    if not isinstance(report,Mapping):
        raise OmicsVariantEvidenceHandoffError("Missing original sidecar research annotation")
    if report.get("schema_version")!=SCHEMA:
        raise OmicsVariantEvidenceHandoffError("Unknown MultiOmics optional sidecar schema version")
    if report.get("status")!="QUALIFIED_ORIGINAL_ASSAY_DESCRIPTIVE_HYPOTHESES_NOT_CAUSAL_PROOF":
        raise OmicsVariantEvidenceHandoffError("Abstained/unmatched source observations cannot assert mechanisms")
    for key in ("causal_variant_mechanism_proven_by_this_report",
                "protein_decay_observed_directly",
                "in_vivo_drug_PK_or_PD_improvement_estimated"):
        if report.get(key) is not False:
            raise OmicsVariantEvidenceHandoffError("Sidecar illegally asserts measured clinical or causal function")
    if (report.get("assays_share_source_biological_unit_tissue_and_time") is not True
            or report.get("can_export_as_optional_separate_annotation_not_exchange_V1_schema_change") is not True):
        raise OmicsVariantEvidenceHandoffError("Unqualified original assay context or destructive schema claim")
    if not isinstance(drug,str) or drug not in KNOWN or endpoint not in KNOWN[drug]["endpoints"]:
        raise OmicsVariantEvidenceHandoffError("Original drug-specific physiological endpoint not curated")
    gene=report.get("gene")
    if not isinstance(gene,str) or gene not in KNOWN[drug]["transporters"]:
        raise OmicsVariantEvidenceHandoffError("No source-qualified pharmacology gene for original drug")
    if report.get("genotype_coverage_state") not in {
            "original_variant_called","original_no_variant_detected_in_qualified_panel",
            "untyped","partial"}:
        raise OmicsVariantEvidenceHandoffError("Genetic original call status not qualified")
    if not isinstance(orig_PGx_pk_join_verified,bool):
        raise OmicsVariantEvidenceHandoffError("Original same-subject PK+PGx join must be separately declared")
    h=report.get("hypotheses_to_test")
    if not isinstance(h,list) or not h or any(x not in HYPOTHESES for x in h):
        raise OmicsVariantEvidenceHandoffError("Unknown sidecar causal mechanistic assertion")
    original=report.get("original_per_assay_signals_measured_not_inferred")
    if not isinstance(original,Mapping) or not original or any(
          k not in READOUTS or v not in ALLOWED_STATES for k,v in original.items()):
        raise OmicsVariantEvidenceHandoffError("Unsupported original assayed layer or misused state")
    assays=[{"readout":k,"modality":READOUTS[k][0],
             "measurement_specificity":READOUTS[k][1],
             "original_source_level_signal_not_numeric_pk_predictor":v}
            for k,v in sorted(original.items())]
    return {
        "status":"QUALIFIED_HYPOTHESIS_PRIOR_ONLY_PATIENT_PK_EFFECT_NOT_ESTABLISHED",
        "drug":drug,"endpoint":endpoint,
        "mechanism_gene":gene,
        "curated_drug_process":KNOWN[drug]["transporters"][gene],
        "source_variant_consequence":report.get("genomic_consequence"),
        "genotype_coverage":report["genotype_coverage_state"],
        "observed_omic_layers_not_fitted_as_PK_predictors":assays,
        "mechanistic_hypotheses_for_independent_test":sorted(set(h)),
        "original_PK_PGx_omic_linkage_separately_verified":orig_PGx_pk_join_verified,
        "no_original_numeric_omic_matrix_embedded":True,
        "eligible_as_patient_level_PK_predictor":False,
        "any_superiority_vs_PGx_only_proven":False,
        "need_original_per_subject_timed_omic_matrices":True,
        "need_nested_outcome_holdout_with_zero_omic_arm":True,
        "do_not_infer_variant_effect_or_protein_decay_from_SHAP":True,
        "original_MultiOmics_DOWNSTREAM_EXCHANGE_V1_stays_unchanged":True,
    }
