"""Prospective, assay-aware physiopharmacological omic acquisition triage.

No mutation-to-activity percentage, no learned PK gain, no inference of
functional protein from total RNA / peptide intensity / nonspecific metabolite.
Metadata-only planning, independent of heldout patient outcomes. For actual
clinical PK prediction retain original same-patient train/heldout comparators.

Primary literature: AMP/CAP CPIC consensus CYP3A5 mechanisms,
https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/ ;
NMD canonical EJC rule and exceptions,
https://pmc.ncbi.nlm.nih.gov/articles/PMC10045457/ ;
Drug/organ process reference: physiopharmaco_cascade.KNOWN.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from urllib.parse import urlparse

from .physiopharmaco_cascade import KNOWN, PK_SOURCES
from .functional_knowledge import CYP3A5_ALLELES
from .cyp3a_activity_phenotype import ASSAY_ROLES


class OmicTriageError(ValueError):
    """Reject contradictory or underspecified structural/source assumptions."""


METHODS = {
    "RNA_GENE_BULK": ("transcriptomics", "total_gene_abundance"),
    "RNA_3PRIME_TAG": ("transcriptomics", "total_gene_abundance"),
    "RNA_EXON_SHORT_READ": ("transcriptomics", "exon_abundance"),
    "RNA_SPLICE_JUNCTION": ("transcriptomics", "splice_junction"),
    "RNA_LONG_READ_ISOFORM": ("transcriptomics", "isoform_resolved"),
    "PROTEIN_BOTTOM_UP_SHARED_PEPTIDE": ("proteomics", "unresolved_protein_abundance"),
    "PROTEIN_TARGETED_PEPTIDE": ("proteomics", "targeted_peptide"),
    "PROTEIN_TOP_DOWN_INTACT": ("proteomics", "intact_proteoform"),
    "METABOLITE_MS1_UNTARGETED": ("metabolomics", "source_chemical_peak"),
    "METABOLITE_MSMS_TARGETED": ("metabolomics", "identified_molecule"),
    "METABOLITE_ISOTOPE_FLUX": ("metabolomics", "metabolic_flux"),
    "ACTIVITY_CYP3A_SHARED_PROBE": ("functional_activity", "shared_cyp3a"),
    "ACTIVITY_CYP3A5_SELECTIVE_PROBE": ("functional_activity", "cyp3a5_selective"),
    "ACTIVITY_TRANSPORT_SUBSTRATE": ("functional_activity", "measured_transport"),
}
MATRIX_TISSUES = {
    "liver", "intestine", "kidney", "whole_blood", "blood", "plasma",
    "urine", "hepatocytes", "microsomes", "cell_line", "unknown",
}
VARIANT_CONSEQUENCES = {
    "premature_stop", "frameshift", "splice_disrupting", "missense",
    "synonymous", "regulatory", "unknown",
}
GENOTYPE_STATUSES = {"complete", "partial", "missing"}
PEPTIDE_REGIONS = {"n_terminal", "c_terminal", "spans_variant", "unknown"}
IDENTITIES = {"unassigned", "source_name_only", "putative_ms2", "confirmed_standard"}
ROLES = ("direct_measured_function", "variant_or_isoform_resolution",
         "protein_or_isoform_abundance", "regulatory_expression_hypothesis",
    "total_rna_proxy", "pathway_proxy",
         "nonspecific_signal")
LEVELS = {
    "direct_measured_function": 1,
    "variant_or_isoform_resolution": 2,
    "protein_or_isoform_abundance": 3,
    "regulatory_expression_hypothesis": 3,
    "total_rna_proxy": 4,
    "pathway_proxy": 5,
    "nonspecific_signal": 6,
}


def _nonblank(value, name):
    if not isinstance(value, str) or not value.strip() or len(value) > 256:
        raise OmicTriageError(name + ": explicit short nonempty string required")
    return value.strip()


def _source_url(value, name):
    url = _nonblank(value, name)
    parse = urlparse(url)
    if parse.scheme != "https" or not parse.hostname:
        raise OmicTriageError(name + ": identifiable HTTPS source required")
    return url


def _classical_nmd_rule(variant):
    """Canonical EJC-dependent NMD *hypothesis* on one specified transcript.

    Stop coordinate is 0-based transcript position in spliced transcript;
    junction coordinate is last exon-exon junction in SAME transcript.
    Distance ~50–55 nt is not a hard deterministic biological threshold.
    """
    if variant["consequence"] not in {"premature_stop", "frameshift", "splice_disrupting"}:
        return "NOT_APPLICABLE_NO_PREMATURE_STOP_CLAIM"
    stop = variant.get("premature_stop_nt")
    last = variant.get("last_exon_junction_nt")
    if stop is None and last is None:
        return "UNKNOWN_NO_TRANSCRIPT_RESOLVED_STOP_AND_LAST_JUNCTION"
    if stop is None or last is None:
        raise OmicTriageError("NMD geometry must include BOTH premature_stop_nt and last_exon_junction_nt")
    if (not isinstance(stop, int) or isinstance(stop, bool) or stop < 0
            or not isinstance(last, int) or isinstance(last, bool) or last < 0):
        raise OmicTriageError("NMD coordinates must be nonnegative integers in spliced transcript")
    if not variant.get("transcript_id"):
        raise OmicTriageError("NMD geometry requires qualified original transcript accession")
    if stop >= last:
        return "CLASSICAL_EJC_NMD_NOT_PREDICTED_LAST_EXON_OR_AFTER_FINAL_JUNCTION"
    distance = last - stop
    if distance >= 55:
        return "CLASSICAL_EJC_NMD_PLAUSIBLE_NOT_CONFIRMED"
    if distance <= 50:
        return "CLASSICAL_EJC_NMD_NOT_PREDICTED_WITHIN_LAST_50NT"
    return "NMD_BOUNDARY_51_TO_54NT_INDETERMINATE"


def _variants(entries):
    if not isinstance(entries, (list, tuple)) or len(entries) > 1000:
        raise OmicTriageError("Explicit source-verified variant event list required")
    result = []
    for item in entries:
        if not isinstance(item, Mapping):
            raise OmicTriageError("Variant event must be a mapping")
        gene = _nonblank(item.get("gene"), "variant gene").upper()
        consequence = item.get("consequence")
        if consequence not in VARIANT_CONSEQUENCES:
            raise OmicTriageError("Variant consequence type must be qualified, not guessed from a score")
        call = item.get("call_status")
        if call not in GENOTYPE_STATUSES:
            raise OmicTriageError("Variant call status complete/partial/missing required")
        if call == "missing" and consequence != "unknown":
            raise OmicTriageError("Missing variant cannot carry an asserted functional consequence")
        if call != "missing":
            _source_url(item.get("source_url"), "variant source")
            _nonblank(item.get("source_variant_id"), "original variant identifier")
        if item.get("allele_state", "unresolved") not in {
                "heterozygous", "biallelic", "hemizygous", "unresolved"}:
            raise OmicTriageError("Variant haplotype and allele dosage not established by gene mutation flag")
        transcript = item.get("transcript_id")
        if transcript is not None and (not isinstance(transcript, str)
                or not (transcript.startswith("NM_") or transcript.startswith("ENST"))):
            raise OmicTriageError("Use qualified NM_ or ENST transcript accession")
        star = item.get("star_allele")
        if star is not None:
            if star == "*1":
                raise OmicTriageError("Reference *1 is a haplotype call, not a causal mutation event")
            if gene != "CYP3A5" or star not in CYP3A5_ALLELES:
                raise OmicTriageError("Only explicitly source-classified CYP3A5 *1/*3/*6/*7 allowed")
            expected = {
                "*1": ("unknown", "regulatory"),
                "*3": ("splice_disrupting",),
                "*6": ("splice_disrupting",),
                "*7": ("frameshift",),
            }[star]
            if consequence not in expected:
                raise OmicTriageError("Star-allele claim inconsistent with curated CYP3A5 variant consequence")
            if call != "complete":
                raise OmicTriageError("No star-allele mechanism inferred from partial or missing call")
        rec = dict(item)
        rec["gene"] = gene
        rec["nmd_rule"] = _classical_nmd_rule(rec)
        result.append(rec)
    return result


def _method_assay(entry):
    if not isinstance(entry, Mapping):
        raise OmicTriageError("Assay must be a metadata record")
    name = _nonblank(entry.get("assay_id"), "source assay ID")
    method = entry.get("method")
    if method not in METHODS:
        raise OmicTriageError("Unrecognized assay technology: explicit platform required")
    tissue = entry.get("tissue")
    if tissue not in MATRIX_TISSUES:
        raise OmicTriageError("Original tissue/compartment not qualified")
    gene = entry.get("gene")
    candidates = entry.get("candidate_genes")
    if METHODS[method][0] != "metabolomics":
        gene = _nonblank(gene, "assay gene").upper()
        if candidates is not None:
            raise OmicTriageError("Candidate-enzyme assignment is only for metabolite assays")
    else:
        if gene is not None:
            raise OmicTriageError("A metabolite cannot be assigned to a unique gene by name alone")
        if not isinstance(candidates, (list, tuple)) or any(
                not isinstance(g, str) or not g.strip() for g in candidates):
            raise OmicTriageError("Metabolite requires explicit candidate genes (may be empty)")
        gene = None
        candidates = [g.strip().upper() for g in candidates]
    if entry.get("status") not in {"proposed", "original_measured"}:
        raise OmicTriageError("Assay requires proposed or original_measured provenance")
    if entry["status"] == "original_measured":
        _source_url(entry.get("source_url"), "original assay source")
        _nonblank(entry.get("original_feature_id"), "measured source feature ID")
    if entry.get("collection_phase") not in {"predose", "postdose", "unknown"}:
        raise OmicTriageError("Assay original sample timing must be explicit")
    if method in {"RNA_SPLICE_JUNCTION", "RNA_LONG_READ_ISOFORM"}:
        _nonblank(entry.get("transcript_id"), "isoform-resolved RNA target")
    if method == "PROTEIN_TARGETED_PEPTIDE":
        if entry.get("peptide_region") not in PEPTIDE_REGIONS:
            raise OmicTriageError("Targeted protein assay needs peptide location relative to variant")
    if method == "PROTEIN_TOP_DOWN_INTACT":
        if entry.get("intact_identity_verified") not in (True, False):
            raise OmicTriageError("Intact-proteoform assertion requires explicit verification boolean")
        if entry.get("intact_identity_verified") is True:
            _source_url(entry.get("original_spectrum_url"), "Original proteoform spectrum and isoform evidence")
    if method.startswith("METABOLITE_"):
        if entry.get("chemical_identity_grade") not in IDENTITIES:
            raise OmicTriageError("Metabolite structure/identity confidence must be declared")
        if method == "METABOLITE_MS1_UNTARGETED" and entry["chemical_identity_grade"] == "confirmed_standard":
            raise OmicTriageError("Untargeted MS1 alone does not confirm chemical structure with reference standard")
    if method.startswith("ACTIVITY_"):
        if not isinstance(entry.get("assay_validated"), bool):
            raise OmicTriageError("Enzyme/transporter probe requires an explicit validation boolean")
        if method == "ACTIVITY_TRANSPORT_SUBSTRATE":
            _nonblank(entry.get("substrate_drug"), "transported substrate drug")
        if method == "ACTIVITY_CYP3A5_SELECTIVE_PROBE":
            original_ref=_source_url(entry.get("selectivity_evidence_url"), "selective enzyme assay method evidence")
            if (entry.get("probe_assay_code") != "T5_N_OXIDATION"
                    or original_ref != "https://doi.org/"+ASSAY_ROLES["T5_N_OXIDATION"]["reference"]):
                raise OmicTriageError("CYP3A5 selective-probe name/reference must exactly match curated T5 source method")
        if method == "ACTIVITY_CYP3A_SHARED_PROBE":
            code=entry.get("probe_assay_code")
            if code not in {"MIDAZOLAM_1_HYDROXYLATION", "TESTOSTERONE_6B_HYDROXYLATION"}:
                raise OmicTriageError("Shared CYP3A requires documented substrate-specific probe identity")
            if (_source_url(entry.get("probe_evidence_url"), "shared CYP3A source")
                    != "https://doi.org/"+ASSAY_ROLES[code]["reference"]):
                raise OmicTriageError("Shared CYP3A substrate probe reference not qualified")
    return {
        **entry, "assay_id":name, "method":method, "gene":gene,
        "candidate_genes":candidates,
        "modality":METHODS[method][0],
    }


def _is_relevant_tissue(drug, gene, tissue, endpoint):
    role, organ, _ = KNOWN[drug]["transporters"][gene]
    required = organ.split("_and_")
    aliases = {"hepatocytes":"liver", "microsomes":"liver"}
    actual = aliases.get(tissue, tissue)
    if endpoint in {"CLr"} and "renal" not in role and organ != "kidney":
        return "DRUG_PATHWAY_NOT_DIRECT_FOR_RENAL_CLEARANCE"
    if actual in required:
        return "DIRECT_DRUG_RELEVANT_TISSUE"
    if actual in {"plasma", "blood", "whole_blood", "urine"}:
        return "CIRCULATING_PROXY_TISSUE_NOT_DIRECT_ORGAN_ACTIVITY"
    return "ORGAN_OR_COMPARTMENT_MISMATCH"


def _rule_assay(drug, endpoint, a, variants):
    method = a["method"]
    roles = []
    evidence = []
    issues = []
    target_genes = set(a["candidate_genes"] or ()) if a["gene"] is None else {a["gene"]}
    relevant = [g for g in sorted(target_genes) if g in KNOWN[drug]["transporters"]]
    if not relevant:
        return {
            "status": "ABSTAIN_NOT_IN_DRUG_SPECIFIC_MECHANISM_REGISTRY",
            "source_gene_or_pathway_hypothesis":sorted(target_genes),
            "assay_id":a["assay_id"], "method":method, "modality":a["modality"],
            "mechanism_confirmed_for_subject":False,
            "eligible_for_outcome_ablation":False,
            "why": ["NO_VERIFIED_DRUG_SPECIFIC_PATHWAY_LINK_FOR_THIS_ASSAY"],
        }
    tissue = [_is_relevant_tissue(drug,g,a["tissue"],endpoint) for g in relevant]
    if any(t != "DIRECT_DRUG_RELEVANT_TISSUE" for t in tissue):
        issues.extend(sorted(set(t for t in tissue if t != "DIRECT_DRUG_RELEVANT_TISSUE")))
    # An uptake/transporter process in another organ cannot automatically
    # become a direct mechanistic covariate for a different PK endpoint.
    if all(t == "DRUG_PATHWAY_NOT_DIRECT_FOR_RENAL_CLEARANCE" for t in tissue):
        return {
            "status":"ABSTAIN_WRONG_DRUG_PROCESS_FOR_ENDPOINT",
            "assay_id":a["assay_id"],"method":method,"modality":a["modality"],
            "source_gene_or_pathway_hypothesis":relevant,
            "mechanism_confirmed_for_subject":False,
            "eligible_for_outcome_ablation":False,
            "why":["RENAL_CLEARANCE_NOT_INFERRED_FROM_HEPATIC_UPTAKE_PROXY"],
        }
    assoc = [v for v in variants if v["gene"] in relevant and v["call_status"] != "missing"]
    nmd = {v["nmd_rule"] for v in assoc}
    if assoc and any(v["nmd_rule"] == "CLASSICAL_EJC_NMD_PLAUSIBLE_NOT_CONFIRMED" for v in assoc):
        evidence.append("POTENTIAL_TRANSCRIPT_DECAY_NOT_CONFIRMED")
    if any(v["nmd_rule"].startswith("NMD_BOUNDARY") for v in assoc):
        issues.append("STOP_NEAR_NMD_BOUNDARY_NEEDS_FUNCTIONAL_ASSAY")
    if any(v["nmd_rule"].startswith("UNKNOWN_NO_TRANSCRIPT") for v in assoc):
        issues.append("PREMATURE_TERMINATION_NMD_GEOMETRY_UNRESOLVED")
    if method in {"RNA_GENE_BULK","RNA_3PRIME_TAG","RNA_EXON_SHORT_READ"}:
        if assoc and any(v["consequence"]=="regulatory" for v in assoc) and not any(
                v["consequence"] in {"splice_disrupting","premature_stop","frameshift"}
                for v in assoc):
            roles.append("regulatory_expression_hypothesis")
            evidence.append("RNA_CAN_REFLECT_REGULATORY_EXPRESSION_CHANGE_NOT_ACTIVITY")
        else:
            roles.append("total_rna_proxy")
        if assoc and any(v["consequence"]=="missense" for v in assoc):
            issues.append("MISSENSE_CAN_IMPAIR_PROTEIN_FUNCTION_WITH_UNCHANGED_RNA_ABUNDANCE")
        if assoc and any(v["consequence"] in {"splice_disrupting","premature_stop","frameshift"} for v in assoc):
            issues.append("TOTAL_RNA_CANNOT_DISTINGUISH_FUNCTIONAL_FROM_ABERRANT_TRANSCRIPT")
        if any(v.get("allele_state", "unresolved")=="heterozygous" for v in assoc):
            issues.append("WILDTYPE_ALLELE_MAY_STILL_PRODUCE_FUNCTIONAL_TRANSCRIPT")
        if method == "RNA_3PRIME_TAG":
            issues.append("3PRIME_TAG_MAY_PERSIST_DESPITE_UPSTREAM_STOP")
        if method == "RNA_EXON_SHORT_READ":
            issues.append("SHORT_READ_EXON_COUNTS_DO_NOT_RESOLVE_FULL_ISOFORM")
    elif method in {"RNA_SPLICE_JUNCTION", "RNA_LONG_READ_ISOFORM"}:
        roles.append("variant_or_isoform_resolution")
        if not assoc:
            issues.append("NO_QUALIFIED_VARIANT_TO_RESOLVE")
        if assoc and any(v.get("transcript_id") not in (None,a["transcript_id"]) for v in assoc):
            issues.append("ASSAY_TRANSCRIPT_MAY_NOT_MATCH_VARIANT_TRANSCRIPT")
        evidence.append("TRANSCRIPT_STRUCTURE_NOT_PROTEIN_ACTIVITY")
    elif method == "PROTEIN_BOTTOM_UP_SHARED_PEPTIDE":
        roles.append("nonspecific_signal")
        issues.append("SHARED_PEPTIDES_MAY_INCLUDE_TRUNCATED_OR_PARALOGOUS_PROTEIN")
    elif method == "PROTEIN_TARGETED_PEPTIDE":
        roles.append("protein_or_isoform_abundance")
        if a["peptide_region"] == "n_terminal":
            issues.append("N_TERMINAL_PEPTIDE_CAN_PERSIST_IN_TRUNCATED_INACTIVE_PROTEIN")
        elif a["peptide_region"] == "c_terminal":
            evidence.append("C_TERMINAL_PEPTIDE_MAY_HELP_DISTINGUISH_TRUNCATION_IF_COVERAGE_VERIFIED")
        elif a["peptide_region"] == "unknown":
            issues.append("PEPTIDE_REGION_UNKNOWN_CANNOT_VERIFY_FULL_LENGTH_PROTEIN")
        else:
            evidence.append("VARIANT_REGION_PEPTIDE_REQUIRES_ALLELE_SPECIFIC_SPECTRUM")
        if assoc and any(v["nmd_rule"] != "NOT_APPLICABLE_NO_PREMATURE_STOP_CLAIM" for v in assoc):
            issues.append("PEPTIDE_SIGNAL_NOT_DIRECT_ENZYME_ACTIVITY")
    elif method == "PROTEIN_TOP_DOWN_INTACT":
        roles.append("variant_or_isoform_resolution" if a["intact_identity_verified"] else "protein_or_isoform_abundance")
        if not a["intact_identity_verified"]:
            issues.append("INTACT_PROTEOFORM_ID_NOT_ORIGINALLY_VERIFIED")
        else:
            evidence.append("SOURCE_ATTESTED_INTACT_PROTEOFORM_ID_REQUIRES_ORIGINAL_SCAN_AUDIT")
        issues.append("INTACT_PROTEOFORM_ABUNDANCE_NOT_ENZYME_FUNCTION")
    elif method.startswith("METABOLITE_"):
        if len(relevant) > 1:
            issues.append("METABOLITE_IS_NOT_SPECIFIC_TO_ONE_CYP_OR_TRANSPORTER")
        if len(relevant) != len(target_genes):
            issues.append("METABOLITE_HAS_ADDITIONAL_UNQUALIFIED_CANDIDATE_GENES")
        grade=a["chemical_identity_grade"]
        if grade in {"source_name_only","unassigned"}:
            issues.append("SMALL_MOLECULE_IDENTITY_UNVERIFIED")
        elif grade == "putative_ms2":
            issues.append("PUTATIVE_MSMS_IDENTITY_NOT_STANDARD_CONFIRMED")
        if method == "METABOLITE_ISOTOPE_FLUX":
            roles.append("pathway_proxy")
            evidence.append("TRACER_FLUX_REQUIRES_SOURCE_ISOTOPE_AND_PATHWAY_SPECIFICITY")
        else:
            roles.append("pathway_proxy" if len(relevant)==1 and grade=="confirmed_standard"
                         else "nonspecific_signal")
            issues.append("METABOLITE_ABUNDANCE_NOT_SPECIFIC_ENZYME_OR_TRANSPORT_RATE")
    elif method == "ACTIVITY_CYP3A5_SELECTIVE_PROBE":
        if a["gene"]!="CYP3A5" or drug!="tacrolimus":
            issues.append("SELECTIVE_CYP3A5_READOUT_NOT_THIS_DRUG_PROCESS")
        if a["assay_validated"] and "SELECTIVE_CYP3A5_READOUT_NOT_THIS_DRUG_PROCESS" not in issues:
            roles.append("direct_measured_function")
        else:
            roles.append("nonspecific_signal")
            issues.append("SELECTIVE_ENZYME_ASSAY_NOT_VALIDATED_FOR_THIS_TARGET")
        evidence.append("SELECTIVE_CYP3A5_ACTIVITY_PROXY_NOT_ACTUAL_TACROLIMUS_CLEARANCE")
    elif method == "ACTIVITY_CYP3A_SHARED_PROBE":
        roles.append("pathway_proxy")
        issues.append("SHARED_CYP3A4_CYP3A5_ACTIVITY_NOT_ISOFORM_SPECIFIC")
        if not a["assay_validated"]:
            issues.append("PROBE_NOT_VALIDATED")
    elif method == "ACTIVITY_TRANSPORT_SUBSTRATE":
        match=a["substrate_drug"].casefold()==drug
        if a["tissue"]=="microsomes":
            issues.append("MICROSOME_FRACTION_LACKS_INTACT_TRANSPORTER_CELL_CONTEXT")
            match=False
        if a["assay_validated"] and match and len(relevant)==1:
            roles.append("direct_measured_function")
        else:
            roles.append("pathway_proxy")
            issues.append("TRANSPORT_ASSAY_NOT_SPECIFICALLY_VALIDATED_FOR_INDEX_DRUG_GENE")
    else:
        raise OmicTriageError("Unsupported modality rule path")
    role=min(roles,key=lambda r:LEVELS[r])
    tissue_quality=("DIRECT_PHYSIOLOGICAL_TARGET_ORGAN" if all(
                    t=="DIRECT_DRUG_RELEVANT_TISSUE" for t in tissue)
                    else "INDIRECT_CIRCULATING_SURROGATE_NOT_ORGAN_FUNCTION" if all(
                    t in {"DIRECT_DRUG_RELEVANT_TISSUE",
                          "CIRCULATING_PROXY_TISSUE_NOT_DIRECT_ORGAN_ACTIVITY"}
                    for t in tissue)
                    else "NONMATCHING_OR_PARTIALLY_UNSUPPORTED_ORGAN")
    tissue_priority={
        "DIRECT_PHYSIOLOGICAL_TARGET_ORGAN":1,
        "INDIRECT_CIRCULATING_SURROGATE_NOT_ORGAN_FUNCTION":2,
        "NONMATCHING_OR_PARTIALLY_UNSUPPORTED_ORGAN":3,
    }[tissue_quality]
    if tissue_priority>1:
        issues.append("METHODOLOGICALLY_VALID_ASSAY_IS_NOT_DIRECT_DRUG_ORGAN_FUNCTION")
    if (a["collection_phase"]=="postdose" or
            a["collection_phase"]=="unknown"):
        issues.append("NOT_QUALIFIED_PREDOSE_PREDICTOR")
    if a["status"]=="proposed":
        issues.append("NOT_YET_MEASURED_PROSPECTIVE_ASSAY")
    if any(t=="ORGAN_OR_COMPARTMENT_MISMATCH" for t in tissue):
        issues.append("NON_TARGET_ORGAN_READOUT")
    if a["status"]=="original_measured" and a["collection_phase"]=="predose":
        # This is only a *candidate*: original source subject joins,
        # biological specimen timing, outcome independence, and batch/feature
        # preprocessing MUST be independently checked by PK benchmark gates.
        eligible=True
    else:
        eligible=False
    return {
        "assay_id":a["assay_id"], "method":method, "modality":a["modality"],
        "source_gene_or_pathway_hypothesis":relevant,
        "drug_process_labels":[{
            "gene":g,
            "role":KNOWN[drug]["transporters"][g][0],
            "organ":KNOWN[drug]["transporters"][g][1],
            "pk_process":KNOWN[drug]["transporters"][g][2],
        } for g in relevant],
        "variant_NMD_statuses":sorted(nmd),
        "candidate_explanatory_role":role,
        "methodological_mechanistic_tissue_relevance":tissue_quality,
        "tissue_relevance_priority_not_expected_prediction_gain":tissue_priority,
        "heuristic_acquisition_priority_tier_not_predictive_rank":LEVELS[role],
        "status":("CANDIDATE_FOR_FOLD_SAFE_EVALUATION" if eligible else
                  "PROSPECTIVE_MEASUREMENT_NOT_YET_ELIGIBLE" if a["status"]=="proposed" else
                  "ABSTAIN_PREDOSE_LEAKAGE_OR_TIMING_UNRESOLVED"),
        "eligible_for_outcome_ablation":eligible,
        "evidence_observations":sorted(set(evidence)),
        "limitations":sorted(set(issues)),
        "methodology_is_not_evidence_of_patient_specific_PK_gain":True,
        "caller_asserted_source_identifiers_not_independently_authenticated_here":True,
        "mechanism_confirmed_for_subject":False,
    }


def plan_physio_omic_acquisition(*, drug, endpoint, variants=(), assays=()):
    """Pre-outcome ranked measurement design; not 'best omic' learned from data.

    The discrete ordering ranks *mechanistic interpretability of proposed
    assay type under supplied source-qualified variants*, not expected PK
    predictive improvement. Report ties and per-gene/organ coverage.
    """
    if not isinstance(drug, str) or drug not in KNOWN:
        raise OmicTriageError("Drug requires independently curated organ/process references")
    if endpoint not in KNOWN[drug]["endpoints"]:
        raise OmicTriageError("Drug-specific original PK endpoint must be explicit")
    variants=_variants(variants)
    if not isinstance(assays, (list, tuple)) or not assays or len(assays)>1000:
        raise OmicTriageError("At least one original or proposed assay descriptor is required")
    parsed=[_method_assay(a) for a in assays]
    names=[a["assay_id"] for a in parsed]
    if len(set(names))!=len(names):
        raise OmicTriageError("Original assay identifiers must be unique")
    vgenes={v["gene"] for v in variants}
    if vgenes-set(KNOWN[drug]["transporters"]):
        raise OmicTriageError("Cannot assign uncurated variant gene to drug pathway")
    outputs=[_rule_assay(drug,endpoint,a,variants) for a in parsed]
    outputs.sort(key=lambda a:(
        not a["status"].startswith(("CANDIDATE","PROSPECTIVE")),
        a.get("tissue_relevance_priority_not_expected_prediction_gain",99),
        a.get("heuristic_acquisition_priority_tier_not_predictive_rank",999),
        a["assay_id"]))
    grouped={}
    for item in outputs:
        if item["status"].startswith(("CANDIDATE","PROSPECTIVE")):
            for gene in item["source_gene_or_pathway_hypothesis"]:
                grouped.setdefault(gene,[]).append(item["assay_id"])
    available_modalities=sorted({
        a["modality"] for a in outputs
        if a.get("eligible_for_outcome_ablation") is True
        and a["modality"] in {"transcriptomics","proteomics","metabolomics"}
    })
    # A source-method declaration alone DOES NOT certify a same-donor
    # multiomic matrix. These are prospective arms only: individual
    # intersections must be source-verified and benchmarked independently.
    from itertools import combinations
    proposed_arms=[
        {"name":"clinical_plus_genotype_when_complete","modalities":[]}
    ] + [
        {"name":"+".join(layers), "modalities":list(layers)}
        for k in range(1,len(available_modalities)+1)
        for layers in combinations(available_modalities,k)
    ]
    required_clinical_covariates=(
        ["dose","dose_time","sampling_time","formulation","hematocrit",
         "interaction_exposure","original_CYP3A4_CYP3A5_genotype"]
        if drug=="tacrolimus"
        else ["dose","dose_time","sampling_time","original_renal_function",
              "measured_original_renal_transport_genotype_when_available"]
    )
    return {
        "prospective_fold_safe_evaluation_plan":{
            "candidate_arms_not_fitted_models":proposed_arms,
            "actual_same_subject_modalities_intersection_verified":False,
            "required_clinical_covariates":required_clinical_covariates,
            "baseline":"clinical_PopPK_plus_authentic_PGx_when_source_eligible",
            "requires_source_timed_dose_concentration_and_PGx":True,
            "selection":"TRAIN_ONLY_nested_CV_allow_zero_omics_and_control_correlated_SHAP",
            "test":"patient_disjoint_heldout_PK_MAE_RMSE_log_error_bias_and_coverage",
            "requires_predefined_prediction_endpoint_and_metric":True,
            "requires_transport_substrate_and_assay_specificity":True,
            "source_sampling_timing_and_outcome_leakage_independent_gate":True,
            "external_validation_required":True,
        },
        "status":"RESEARCH_PRE_OUTCOME_ASSAY_AWARE_PHYSIOPHARMACO_DESIGN_NOT_PK_VALIDATION",
        "drug":drug,
        "index_PK_endpoint":endpoint,
        "drug_specific_mechanism_source_url":PK_SOURCES[drug],
        "variant_source_consequence_interpretations":[{
            "gene":v["gene"],"source_variant_id":v.get("source_variant_id"),
            "star_allele":v.get("star_allele"),
            "call_status":v["call_status"],
            "consequence":v["consequence"],"transcript_id":v.get("transcript_id"),
            "allele_state":v.get("allele_state","unresolved"),
            "canonical_NMD_rule_not_function_proof":v["nmd_rule"],
        } for v in variants],
        "assay_candidates_by_interpretability_not_expected_model_performance":outputs,
        "genes_with_any_relevant_assay_by_source":grouped,
        "candidate_physio_assays_require_biological_validity_not_SHAP_only":True,
        "missing_measured_drug_process_genes":sorted(
            set(KNOWN[drug]["transporters"])-set(grouped)),
        "methodological_priorities_not_patient_PK_gain_estimates":True,
        "no_quantitative_activity_or_effect_sizes_invented":True,
        "no_PopPK_CYP_CL_or_dose_inferred":True,
        "no_false_CYP3A5_star_allele_from_gene_level_stop_call":True,
        "original_outcome_not_used_for_assay_prioritization":True,
        "next_validation_contract":(
            "First confirm original same-subject same-landmark and predose "
            "assay matrices/PGx/PK. Then compare clinical+PGx+PopPK vs "
            "each actually measured omic and their interactions in nested "
            "donor-heldout folds, including a zero-omics arm; test "
            "MAE/RMSE and clinical exposure metric predeclared, with "
            "external replication and biological negative controls."
        ),
    }
