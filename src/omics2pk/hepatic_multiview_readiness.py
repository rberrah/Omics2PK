"""Fail-closed readiness gate for same-donor hepatic RNA/protein PK evidence.

This module does not fit a model. It answers one question first: are RNA,
protein and a drug-specific pharmacology outcome measured on the SAME original
human donors from one authenticated source?

It exists to prevent invalid cross-study fusion such as combining EPA RNA
donors with unrelated TKI/Wegler protein donors.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


class HepaticMultiviewReadinessError(ValueError):
    pass


@dataclass(frozen=True)
class HepaticMultiviewReadiness:
    status: str
    source_id: str
    outcome_kind: str
    n_rna: int
    n_protein: int
    n_outcome: int
    n_complete_case: int
    fit_allowed: bool
    missing_requirements: tuple[str, ...]


def _ids(values: Iterable[str] | None, name: str) -> set[str]:
    if values is None:
        return set()
    seq=[str(v).strip() for v in values]
    if any(not v for v in seq):
        raise HepaticMultiviewReadinessError(f"Blank biological donor ID in {name}")
    if len(set(seq)) != len(seq):
        raise HepaticMultiviewReadinessError(f"Duplicated biological donor ID in {name}")
    return set(seq)


def assess_same_donor_hepatic_multiview(
    *,
    source_id: str,
    rna_donors: Iterable[str] | None,
    protein_donors: Iterable[str] | None,
    outcome_donors: Iterable[str] | None,
    outcome_kind: str,
    source_rows_acquired: bool,
    human: bool = True,
) -> HepaticMultiviewReadiness:
    """Return eligibility for a true same-donor RNA+protein comparison.

    source_rows_acquired must mean authenticated individual rows, not article
    statements, scatterplot pixels, group summaries or row-order guesses.
    outcome_kind must be measured drug-specific metabolism or PK, not PBPK.
    """
    source_id=source_id.strip()
    if not source_id:
        raise HepaticMultiviewReadinessError("Source identifier is required")
    if not human:
        return HepaticMultiviewReadiness(
            "ABSTAIN_NONHUMAN_SOURCE",source_id,outcome_kind,0,0,0,0,False,
            ("human_same_donor_source",))
    allowed={"drug_specific_clint","drug_specific_metabolic_activity","clinical_pk"}
    if outcome_kind not in allowed:
        return HepaticMultiviewReadiness(
            "ABSTAIN_NO_DRUG_SPECIFIC_OUTCOME",source_id,outcome_kind,0,0,0,0,False,
            ("drug_specific_clint_or_activity_or_clinical_pk",))
    if not source_rows_acquired:
        return HepaticMultiviewReadiness(
            "ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED",source_id,outcome_kind,0,0,0,0,False,
            ("authenticated_individual_source_rows",))

    rna=_ids(rna_donors,"RNA")
    protein=_ids(protein_donors,"protein")
    outcome=_ids(outcome_donors,"outcome")
    missing=[]
    if not rna: missing.append("matched_RNA")
    if not protein: missing.append("matched_protein")
    if not outcome: missing.append("matched_drug_outcome")
    if missing:
        code=(
            "ABSTAIN_NO_MATCHED_RNA" if "matched_RNA" in missing else
            "ABSTAIN_NO_MATCHED_PROTEIN" if "matched_protein" in missing else
            "ABSTAIN_NO_MATCHED_DRUG_OUTCOME"
        )
        return HepaticMultiviewReadiness(
            code,source_id,outcome_kind,len(rna),len(protein),len(outcome),0,False,
            tuple(missing))

    common=rna & protein & outcome
    if not common:
        return HepaticMultiviewReadiness(
            "ABSTAIN_NO_SAME_DONOR_INTERSECTION",source_id,outcome_kind,
            len(rna),len(protein),len(outcome),0,False,
            ("same_biological_donors_across_all_views",))
    if common != rna or common != protein or common != outcome:
        return HepaticMultiviewReadiness(
            "ABSTAIN_INCOMPLETE_CASE_MISMATCH",source_id,outcome_kind,
            len(rna),len(protein),len(outcome),len(common),False,
            ("predeclared_complete_case_cohort",))

    return HepaticMultiviewReadiness(
        "READY_SAME_DONOR_RNA_PROTEIN_DRUG_OUTCOME",source_id,outcome_kind,
        len(rna),len(protein),len(outcome),len(common),True,())


def qualify_candidate_summary(
    source_id: str,
    *,
    has_rna: bool,
    has_protein: bool,
    has_drug_outcome: bool,
    individual_rows_acquired: bool,
) -> str:
    """Catalogue helper that never invents donor IDs."""
    if not source_id.strip():
        raise HepaticMultiviewReadinessError("Source identifier is required")
    if not individual_rows_acquired:
        return "ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED"
    if not has_rna:
        return "ABSTAIN_NO_MATCHED_RNA"
    if not has_protein:
        return "ABSTAIN_NO_MATCHED_PROTEIN"
    if not has_drug_outcome:
        return "ABSTAIN_NO_MATCHED_DRUG_OUTCOME"
    return "REQUIRES_EXPLICIT_DONOR_ID_INTERSECTION"
