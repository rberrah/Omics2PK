"""Research-only CYP3A5 functional-evidence triage, never an activity estimator.

Input observations are qualitative assertions supplied by authenticated assays.
No genotype is inferred from total RNA, metabolites, or missing variants. This
module does not access clinical data, estimate clearance, or recommend dosing.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


class FunctionalEvidenceError(ValueError):
    """Malformed molecular evidence; no silent biological inference."""


KNOWLEDGE_VERSION = "cyp3a5-functional-v0.1.0-20261103"
SOURCES = {
    "cpic": "https://cpicpgx.org/guidelines/guideline-for-tacrolimus-and-CYP3A5/",
    "pharmvar": "https://www.pharmvar.org/gene/cyp3a5",
    "gene_focus": "https://doi.org/10.1002/cpt.2563",
    "amp_cap": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/",
}
# CPIC/PharmVar categories, NOT numeric percent activity or patient-specific
# residual activity. *1 requires positive evidence of a qualified allele call;
# absence of a tested *3 locus does not establish *1 in an incomplete panel.
CYP3A5_ALLELES = {
    "*1": {"function": "normal", "variant": "reference functional haplotype",
           "mechanism": "functional CYP3A5 expression possible"},
    "*3": {"function": "no_function", "variant": "rs776746 / NM_000777.5:c.219-237A>G",
           "mechanism": "cryptic intron-3 acceptor, aberrant splicing and premature termination"},
    "*6": {"function": "no_function", "variant": "rs10264272 / NM_000777.5:c.624G>A",
           "mechanism": "synonymous exon-7 change disrupts splicing and causes premature termination"},
    "*7": {"function": "no_function", "variant": "rs41303343 / NM_000777.5:c.1035dup",
           "mechanism": "single-nucleotide insertion, frameshift and premature termination"},
}

_TARGETS = {
    "transcriptomics": {"total_gene", "canonical_splice_junction", "aberrant_splice_junction"},
    "proteomics": {"pan_cyp3a", "cyp3a5_full_length", "cyp3a5_unresolved_peptide"},
    "activity": {"cyp3a_shared", "cyp3a5_specific"},
    "metabolomics": {"metabolite_abundance", "tracer_flux"},
}
_SIGNALS = {"detected", "not_detected", "above_validated_reference", "below_validated_reference"}


def _parse_diplotype(diplotype: str | None, call_status: str) -> tuple[list[str], str | None]:
    if call_status not in {"complete", "partial", "missing"}:
        raise FunctionalEvidenceError("call_status must be complete, partial or missing")
    if call_status == "missing":
        if diplotype is not None:
            raise FunctionalEvidenceError("A missing genotype call cannot have a diplotype")
        return [], "GENOTYPE_MISSING"
    if not isinstance(diplotype, str) or not diplotype.strip():
        raise FunctionalEvidenceError("A supplied genotype requires an explicit diplotype")
    alleles = [s.strip().upper() for s in diplotype.split("/")]
    if len(alleles) != 2 or not all(a.startswith("*") and len(a) > 1 for a in alleles):
        raise FunctionalEvidenceError("Diplotype must be two explicit star alleles, e.g. *1/*3")
    if call_status == "partial":
        return alleles, "PARTIAL_ALLELE_COVERAGE"
    if any(a not in CYP3A5_ALLELES for a in alleles):
        return alleles, "UNRESOLVED_ALLELE_FUNCTION"
    return alleles, None


def assess_cyp3a5_evidence(
    *,
    diplotype: str | None = None,
    call_status: str = "missing",
    observations: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Interpret assay *limitations* and qualitative biological discordance.

    No arbitrary genotype activity scores, percentage losses, RNA-expression
    thresholds or weights are assigned. Each observation has ``omic``,
    ``target``, ``signal`` and optionally ``tissue``. Activity observations
    additionally require an explicit ``assay_validated`` boolean, while
    metabolomic observations require a non-empty ``candidate_enzymes`` list.
    Input contains no person identifiers; output contains only rule codes.
    """
    alleles, genotype_gate = _parse_diplotype(diplotype, call_status)
    if isinstance(observations, (str, bytes)) or not isinstance(observations, Sequence):
        raise FunctionalEvidenceError("observations must be a sequence of assay mappings")
    if len(observations) > 10000:
        raise FunctionalEvidenceError("Too many assay observations")

    warnings: list[str] = []
    def warn(code: str) -> None:
        if code not in warnings:
            warnings.append(code)

    phenotype = None
    if genotype_gate is not None:
        warn(genotype_gate)
    else:
        n_normal = sum(CYP3A5_ALLELES[a]["function"] == "normal" for a in alleles)
        phenotype = {0: "nonexpresser", 1: "intermediate_expresser", 2: "expresser"}[n_normal]

    discordance = False
    for item in observations:
        if not isinstance(item, Mapping):
            raise FunctionalEvidenceError("Each observation must be a mapping")
        omic = item.get("omic")
        target = item.get("target")
        signal = item.get("signal")
        if omic not in _TARGETS or target not in _TARGETS[omic]:
            raise FunctionalEvidenceError("Unsupported omic/target; do not infer assay specificity")
        if signal not in _SIGNALS:
            raise FunctionalEvidenceError("Unknown qualitative signal")
        if omic != "activity" and signal in {"above_validated_reference", "below_validated_reference"}:
            raise FunctionalEvidenceError("Reference-normalized activity labels require activity assays")
        if not isinstance(item.get("tissue"), str) or not item["tissue"].strip():
            warn("TISSUE_UNSPECIFIED")
        elif item["tissue"].strip().casefold() not in {"liver", "intestine", "hepatic", "intestinal"}:
            warn("NON_TARGET_TISSUE_EXTRAPOLATION")
        if omic == "transcriptomics":
            if target == "total_gene":
                warn("TOTAL_RNA_NOT_FUNCTIONAL_TRANSCRIPT")
            else:
                warn("SPLICE_JUNCTION_NOT_ENZYME_ACTIVITY")
        elif omic == "proteomics":
            if target != "cyp3a5_full_length":
                warn("PROTEIN_ASSAY_NOT_FULL_LENGTH_SPECIFIC")
            else:
                warn("FULL_LENGTH_PROTEIN_NOT_ACTIVITY")
        elif omic == "metabolomics":
            enzymes = item.get("candidate_enzymes")
            if isinstance(enzymes, (str, bytes)) or not isinstance(enzymes, (list, tuple)) or not enzymes:
                raise FunctionalEvidenceError("Metabolite requires candidate_enzymes (possibly 'unknown')")
            if any(not isinstance(x, str) or not x.strip() for x in enzymes):
                raise FunctionalEvidenceError("Invalid candidate enzyme ID")
            if len(set(enzymes)) != len(enzymes):
                raise FunctionalEvidenceError("Duplicate metabolite candidate enzyme")
            if len(enzymes) != 1 or enzymes[0].lower() == "unknown":
                warn("METABOLITE_PATHWAY_NOT_UNIQUE")
            if target == "metabolite_abundance":
                warn("METABOLITE_ABUNDANCE_NOT_FLUX")
            else:
                warn("TRACER_FLUX_NOT_AUTOMATICALLY_ENZYME_SPECIFIC")
        else:
            validated = item.get("assay_validated")
            if not isinstance(validated, bool):
                raise FunctionalEvidenceError("Activity assay requires explicit boolean assay_validated")
            if target != "cyp3a5_specific":
                warn("SHARED_CYP3A_ACTIVITY_NOT_CYP3A5_SPECIFIC")
            if not validated:
                warn("ACTIVITY_ASSAY_UNVALIDATED")
            # A *qualified* large specific activity signal should trigger
            # review rather than overwrite the CPIC classification. Mere
            # detection is compatible with low residual expression.
            if (validated and target == "cyp3a5_specific"
                    and signal == "above_validated_reference"
                    and phenotype == "nonexpresser"):
                discordance = True
                warn("GENOTYPE_FUNCTIONAL_ACTIVITY_DISCORDANCE_REVIEW")
    status = ("REVIEW_FUNCTIONAL_DISCORDANCE" if discordance else
              "ABSTAIN_GENOTYPE_UNQUALIFIED" if genotype_gate else
              "INTERPRETATION_ONLY" if observations else "GENOTYPE_ONLY")
    return {
        "status": status,
        "knowledge_version": KNOWLEDGE_VERSION,
        "gene": "CYP3A5",
        "diplotype_classification": phenotype,
        "allele_functions": (
            [CYP3A5_ALLELES[a]["function"] for a in alleles]
            if genotype_gate is None else None),
        "evidence_warnings": warnings,
        "n_assay_observations": len(observations),
        "percent_activity_loss": None,
        "estimated_enzyme_activity": None,
        "in_vivo_pk_gain_demonstrated": False,
        "clinical_dosing_use": False,
        "knowledge_sources": dict(SOURCES),
    }
