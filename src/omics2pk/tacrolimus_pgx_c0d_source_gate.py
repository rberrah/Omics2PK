"""Fail-closed source gate for a tacrolimus clinical+PGx C0/D comparator.

Primary reference: Yang et al. 2022, Clin Transl Sci 15:2640-2651,
DOI 10.1111/cts.13388, GEO GSE196772.

The public GEO release contains genotype data, while the article reports that
patient-level clinical variables, tacrolimus dose and C0 were collected from
medical records. This gate therefore accepts a future legitimately obtained
same-subject clinical+genotype table, but it must never reconstruct individual
PK targets from published cohort summaries or SNP association tables.

Scope: stable tacrolimus trough/dose phenotype, not dense concentration-time
PopPK and not an Omics2PK superiority result.
"""
from __future__ import annotations

from collections import Counter
from math import isfinite, log
from typing import Mapping, Sequence


class TacrolimusPGxBaselineError(ValueError):
    pass


PUBLISHED_REFERENCE = {
    "doi": "10.1111/cts.13388",
    "geo": "GSE196772",
    "n_reported": 251,
    "outcome": "log(C0/(daily_dose/weight))",
    "published_in_sample_R2_percent": {
        "age": 11.03,
        "age_sex": 15.76,
        "age_sex_hct": 17.89,
        "age_sex_hct_rs776746": 35.12,
        "age_sex_hct_rs776746_rs75125371": 37.34,
    },
    "heldout_prediction_claim": False,
}


def _text(value: object, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise TacrolimusPGxBaselineError(f"{field} must be a non-empty string")
    return value.strip()


def _number(value: object, field: str, *, positive: bool = False) -> float:
    if isinstance(value, bool):
        raise TacrolimusPGxBaselineError(f"{field} must be numeric")
    try:
        out=float(value)
    except (TypeError,ValueError):
        raise TacrolimusPGxBaselineError(f"{field} must be numeric") from None
    if not isfinite(out):
        raise TacrolimusPGxBaselineError(f"{field} must be finite")
    if positive and out<=0:
        raise TacrolimusPGxBaselineError(f"{field} must be positive")
    return out


def _dosage(value: object, field: str) -> int:
    if isinstance(value,bool):
        raise TacrolimusPGxBaselineError(f"{field} must be allele dosage 0/1/2")
    try:
        numeric=float(value)
    except (TypeError,ValueError):
        raise TacrolimusPGxBaselineError(f"{field} must be allele dosage 0/1/2") from None
    if numeric not in (0.0,1.0,2.0):
        raise TacrolimusPGxBaselineError(f"{field} must be allele dosage 0/1/2")
    return int(numeric)


def qualify_tacrolimus_clinical_pgx_c0d(
    rows: Sequence[Mapping[str, object]],
    *,
    require_rs75125371: bool = False,
    require_stable_maintenance: bool = True,
) -> dict:
    """Validate one same-subject clinical+PGx+C0/D row per recipient.

    Required fields:
      subject_id, age_years, sex, weight_kg, hematocrit_fraction,
      c0_ng_ml, daily_dose_mg, stable_maintenance,
      rs776746_effect_allele, rs776746_effect_allele_dosage

    Optional second PGx layer:
      rs75125371_effect_allele, rs75125371_effect_allele_dosage

    Sex is source-normalized to one of {female,male}. Hematocrit is a fraction
    in [0.1,0.7], so values reported as percent/per-mille must be converted
    explicitly before this gate rather than guessed.

    The function derives the paper-style natural-log dose/weight normalized
    trough target and returns aggregate provenance only; individual values are
    not exported.
    """
    if not rows:
        raise TacrolimusPGxBaselineError("No patient rows")

    subjects=set()
    rs776746_alleles=set()
    rs75125371_alleles=set()
    sex_counts=Counter()
    genotype_counts=Counter()
    candidate_counts=Counter()
    candidate_rows=0
    outcomes=[]

    for row in rows:
        sid=_text(row.get("subject_id"),"subject_id")
        if sid in subjects:
            raise TacrolimusPGxBaselineError("Duplicate biological subject_id")
        subjects.add(sid)

        age=_number(row.get("age_years"),"age_years",positive=True)
        if not 18<=age<=100:
            raise TacrolimusPGxBaselineError("age_years outside adult range")

        sex=_text(row.get("sex"),"sex").casefold()
        if sex not in {"female","male"}:
            raise TacrolimusPGxBaselineError("sex must be female or male")
        sex_counts[sex]+=1

        weight=_number(row.get("weight_kg"),"weight_kg",positive=True)
        if not 25<=weight<=250:
            raise TacrolimusPGxBaselineError("weight_kg outside prespecified human range")

        hct=_number(row.get("hematocrit_fraction"),"hematocrit_fraction")
        if not 0.10<=hct<=0.70:
            raise TacrolimusPGxBaselineError(
                "hematocrit_fraction must be explicit fraction in 0.10..0.70")

        c0=_number(row.get("c0_ng_ml"),"c0_ng_ml",positive=True)
        daily=_number(row.get("daily_dose_mg"),"daily_dose_mg",positive=True)

        stable=row.get("stable_maintenance")
        if not isinstance(stable,bool):
            raise TacrolimusPGxBaselineError("stable_maintenance must be boolean")
        if require_stable_maintenance and not stable:
            raise TacrolimusPGxBaselineError(
                "All modeled rows must satisfy stable maintenance definition")

        allele=_text(row.get("rs776746_effect_allele"),"rs776746_effect_allele").upper()
        if allele not in {"A","C","G","T"}:
            raise TacrolimusPGxBaselineError("rs776746 effect allele must be A/C/G/T")
        rs776746_alleles.add(allele)
        dosage=_dosage(
            row.get("rs776746_effect_allele_dosage"),
            "rs776746_effect_allele_dosage")
        genotype_counts[dosage]+=1

        second_allele=row.get("rs75125371_effect_allele")
        second_dosage=row.get("rs75125371_effect_allele_dosage")
        if require_rs75125371 or second_allele is not None or second_dosage is not None:
            allele2=_text(second_allele,"rs75125371_effect_allele").upper()
            if allele2 not in {"A","C","G","T"}:
                raise TacrolimusPGxBaselineError(
                    "rs75125371 effect allele must be A/C/G/T")
            rs75125371_alleles.add(allele2)
            candidate_counts[_dosage(
                second_dosage,"rs75125371_effect_allele_dosage")]+=1
            candidate_rows+=1

        normalized=c0/(daily/weight)
        if not isfinite(normalized) or normalized<=0:
            raise TacrolimusPGxBaselineError("Derived C0/(dose/weight) invalid")
        outcomes.append(log(normalized))

    if len(rs776746_alleles)!=1:
        raise TacrolimusPGxBaselineError(
            "rs776746 effect-allele orientation must be constant across cohort")
    if rs75125371_alleles and len(rs75125371_alleles)!=1:
        raise TacrolimusPGxBaselineError(
            "rs75125371 effect-allele orientation must be constant across cohort")
    if rs75125371_alleles and candidate_rows!=len(rows):
        raise TacrolimusPGxBaselineError(
            "rs75125371 must be complete for every subject in the comparison cohort")
    if len(genotype_counts)<2:
        raise TacrolimusPGxBaselineError(
            "rs776746 cohort has no genotype variation for comparison")
    if require_rs75125371 and len(candidate_counts)<2:
        raise TacrolimusPGxBaselineError(
            "rs75125371 cohort has no genotype variation for comparison")

    return {
        "status":"TACROLIMUS_CLINICAL_PGX_C0D_SOURCE_GATE_ONLY",
        "n_unique_subjects":len(subjects),
        "n_female":sex_counts["female"],
        "n_male":sex_counts["male"],
        "rs776746_effect_allele":next(iter(rs776746_alleles)),
        "rs776746_dosage_counts":{
            str(k):genotype_counts[k] for k in sorted(genotype_counts)
        },
        "rs75125371_included":bool(rs75125371_alleles),
        "rs75125371_effect_allele":(
            next(iter(rs75125371_alleles)) if rs75125371_alleles else None
        ),
        "rs75125371_dosage_counts":{
            str(k):candidate_counts[k] for k in sorted(candidate_counts)
        },
        "outcome_definition":"ln(C0_ng_ml / (daily_dose_mg / weight_kg))",
        "stable_maintenance_required":require_stable_maintenance,
        "hematocrit_unit":"fraction",
        "published_reference":PUBLISHED_REFERENCE,
        "individual_subject_ids_or_values_exported":False,
        "not_dense_pk_or_popPK":True,
        "not_omics_validation":True,
        "published_R2_is_in_sample_explanatory_not_heldout_prediction":True,
    }
