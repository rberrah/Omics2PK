"""Fail-closed public Table-4 PHH imatinib pilot, not a clinical PK model.

14 separately identified hepatocyte donors; 13 have full proteomics.
The 14 in-vitro imatinib clearances are outcomes; *drug-incubation* formation
of N-desmethyl-imatinib is NEVER a baseline predictor. CYP3A midazolam probe
activity is shared CYP3A4/5, not CYP3A5-specific activity.
"""
from __future__ import annotations

import re
from math import isfinite, log
from statistics import mean
from typing import Mapping, Sequence

from omics2pk.functional_activity_benchmark import _ridge_predict
from omics2pk.functional_knowledge import assess_cyp3a5_evidence


class LathamPHHError(ValueError):
    pass


SOURCE_DOI = "10.1016/j.dmd.2025.100196"
SOURCE_PMCID = "PMC12990092"
PROTEIN_FIELDS = ("cyp2c8_protein", "cyp3a4_protein", "cyp3a5_protein", "cyp3a_total_protein")
NUMERIC_FIELDS = (
    "cyp2c8_activity", "cyp3a_probe_activity", "imatinib_apparent_cl",
    "n_desmethyl_imatinib_formation",
)


def _num(value: object, field: str) -> float:
    if isinstance(value, bool):
        raise LathamPHHError(field + ": invalid boolean")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        raise LathamPHHError(field + ": invalid source number") from None
    if not isfinite(number) or number <= 0:
        raise LathamPHHError(field + ": expected positive finite measured value")
    return number


def _qualified(rows: Sequence[Mapping[str, object]]) -> list[dict]:
    if len(rows) != 14:
        raise LathamPHHError("Expected exactly 14 PHH biological donors in Table 4")
    seen = set()
    out = []
    for raw in rows:
        if not isinstance(raw, Mapping):
            raise LathamPHHError("Invalid biological donor row")
        donor = raw.get("donor_id")
        if not isinstance(donor, str) or not re.fullmatch(r"[A-Z]{3}", donor) or donor in seen:
            raise LathamPHHError("Invalid/duplicate original PHH donor ID")
        seen.add(donor)
        if raw.get("source_doi") != SOURCE_DOI:
            raise LathamPHHError("Original article DOI mismatch")
        cyp2c8 = raw.get("cyp2c8_diplotype")
        if cyp2c8 not in {"*1/*1", "*1/*2", "*1/*3", "*1/*4"}:
            raise LathamPHHError("Unexpected CYP2C8 original diplotype")
        try:
            phenotype = assess_cyp3a5_evidence(
                diplotype=raw.get("cyp3a5_diplotype"), call_status="complete")
        except ValueError as exc:
            raise LathamPHHError("Invalid original CYP3A5 diplotype") from exc
        category = phenotype["diplotype_classification"]
        if category is None:
            raise LathamPHHError("Unresolved CYP3A5 source allele")
        class_name = {"nonexpresser": "PM", "intermediate_expresser": "IM",
                      "expresser": "NM"}[category]
        if raw.get("published_cyp3a5_category") != class_name:
            raise LathamPHHError("PGx classification disagrees with original source")
        x = {f: _num(raw.get(f), f) for f in NUMERIC_FIELDS}
        protein = []
        for f in PROTEIN_FIELDS:
            value = raw.get(f)
            protein.append(None if value is None else _num(value, f))
        if protein[0] is None:
            raise LathamPHHError("Original CYP2C8 protein cannot be absent")
        if any(v is None for v in protein[1:]):
            if protein[1:] != [None, None, None]:
                raise LathamPHHError("Incomplete correlated missing CYP3A proteins")
        elif abs(protein[1] + protein[2] - protein[3]) > 0.16:
            raise LathamPHHError("CYP3A4 + CYP3A5 protein does not reproduce total CYP3A")
        out.append({
            "donor_id": donor,
            "pgx": [float(cyp2c8 != "*1/*1"),
                    float({"PM": 0, "IM": 1, "NM": 2}[class_name])],
            "phenotype": class_name, "cyp2c8_diplotype": cyp2c8,
            "cyp3a5_diplotype": raw["cyp3a5_diplotype"],
            "protein": protein[:3] if protein[1] is not None else None,
            "activity": [x["cyp2c8_activity"], x["cyp3a_probe_activity"]],
            "outcome": x["imatinib_apparent_cl"],
            "formation": x["n_desmethyl_imatinib_formation"],
        })
    counts = {category: sum(r["phenotype"] == category for r in out)
              for category in ("NM", "IM", "PM")}
    if counts != {"NM": 1, "IM": 3, "PM": 10}:
        raise LathamPHHError("Original CYP3A5 1/3/10 genotype group distribution changed")
    if sum(r["protein"] is None for r in out) != 1:
        raise LathamPHHError("Expected 1 originally missing CYP3A protein panel")
    if sum(r["cyp2c8_diplotype"] == "*1/*1" for r in out) != 8:
        raise LathamPHHError("Original CYP2C8 genotype balance changed")
    if not any("*6" in r["cyp3a5_diplotype"] for r in out) or not any(
            "*7" in r["cyp3a5_diplotype"] for r in out):
        raise LathamPHHError("Original CYP3A5 *6/*7 molecular source control missing")
    return sorted(out, key=lambda x: x["donor_id"])


def benchmark_latham_phh_table4(rows: Sequence[Mapping[str, object]]) -> dict:
    """Leave-one-original-hepatocyte-donor-out Ridge on 13 complete cases.

    All four models use exactly the same 13 original donors, train-only
    scaling, fixed alpha=1, and log apparent in-vitro clearance. Not E3.
    Values and donor identifiers are never returned.
    """
    qualified = _qualified(rows)
    complete = [r for r in qualified if r["protein"] is not None]
    if len(complete) != 13:
        raise LathamPHHError("All four arms require the same 13 protein-complete donors")
    arms = ("pgx_only", "pgx_plus_proteins", "pgx_plus_probe_activities",
            "pgx_plus_proteins_and_probe_activities")
    losses = {name: [] for name in arms}
    for held in complete:
        train = [r for r in complete if r["donor_id"] != held["donor_id"]]
        if len(train) != 12:
            raise LathamPHHError("Invalid held-donor training partition")
        for name in arms:
            def vector(r):
                v = list(r["pgx"])
                if name in {"pgx_plus_proteins", "pgx_plus_proteins_and_probe_activities"}:
                    v.extend(r["protein"])
                if name in {"pgx_plus_probe_activities", "pgx_plus_proteins_and_probe_activities"}:
                    v.extend(r["activity"])
                return v
            guess = _ridge_predict([vector(r) for r in train],
                                   [log(r["outcome"]) for r in train], [vector(held)])[0]
            losses[name].append(abs(guess - held["outcome"]))
    scores = {name: round(mean(values), 6) for name, values in losses.items()}
    # Group-specific contrasts are DESCRIPTIVE fold-disjoint measurements,
    # not a discovered mutation-specific predictor: n=13 and PM dominates.
    # Original donor IDs stay private; never export per-person errors.
    strata = {
        "CYP3A5_nonexpresser_PM": [i for i,r in enumerate(complete)
                                  if r["phenotype"]=="PM"],
        "CYP3A5_IM_or_NM": [i for i,r in enumerate(complete)
                            if r["phenotype"]!="PM"],
        "CYP2C8_reference_star1_star1": [i for i,r in enumerate(complete)
                                        if r["cyp2c8_diplotype"]=="*1/*1"],
        "CYP2C8_variant_diplotype": [i for i,r in enumerate(complete)
                                    if r["cyp2c8_diplotype"]!="*1/*1"],
    }
    stratified={}
    for label,index in strata.items():
        if len(index)<4:
            stratified[label]={"n_original_donors":len(index),
                               "status":"ABSTAIN_SMALL_ORIGINAL_GENOTYPE_STRATUM"}
            continue
        stratified[label]={
            "n_original_donors":len(index),
            "status":"DESCRIPTIVE_LODO_ERROR_NOT_GENOTYPE_SPECIFIC_VALIDATION",
            "heldout_MAE_original_printed_scale":{
                arm:round(mean(losses[arm][i] for i in index),6)
                for arm in arms
            },
            "subgroup_no_model_or_feature_selected_using_holdout":True,
        }
    # Positive targeted peptide abundance and *shared* CYP3A probe activity in
    # CPIC nonexpressers must not be relabelled CYP3A5-specific catalysis.
    # Both assays can be technically correct without implying functional *3.
    pm = [r for r in qualified if r["phenotype"] == "PM"]
    pm_with_protein = [r for r in pm if r["protein"] is not None]
    n_pm_peptide = sum(r["protein"][2] > 0 for r in pm_with_protein)
    n_pm_shared_probe = sum(r["activity"][1] > 0 for r in pm)
    return {
        "status": "ORIGINAL_PUBLISHED_TABLE4_EXPLORATORY_IN_VITRO_PILOT_NOT_E3",
        "source_doi": SOURCE_DOI, "source_pmcid": SOURCE_PMCID,
        "source_rows": 14, "n_same_protein_complete_donors_all_arms": 13,
        "excluded_due_to_source_missing_protein": 1,
        "cyp3a5_phenotype_group_counts_original_14": {"NM": 1, "IM": 3, "PM": 10},
        "cyp3a5_PM_phenotype_protein_complete": len(pm_with_protein),
        "cyp3a5_PM_with_positive_targeted_peptide_quantitation": n_pm_peptide,
        "cyp3a5_PM_with_positive_shared_CYP3A_probe_activity": n_pm_shared_probe,
        "positive_peptide_or_shared_probe_is_not_CYP3A5_specific_function": True,
        "included_cyp3a5_star_6_and_star_7": True,
        "endpoint": "PHH_imatinib_apparent_in_vitro_clearance_not_human_plasma_PK",
        "endpoint_table4_unit": "microlitre_per_million_cells_per_minute",
        "source_prose_unit_discrepancy_not_resolved": True,
        "heldout_MAE_by_predeclared_arm_source_table_scale": scores,
        "genotype_stratified_MAE_exploratory":stratified,
        "allele_star6_or_star7_separate_effect_estimable":False,
        "no_RNA_or_endogenous_metabolomics_in_original_PHH_table4":True,
        "MAE_delta_vs_PGx_only": {k: round(v - scores["pgx_only"], 6)
                                 for k, v in scores.items() if k != "pgx_only"},
        "shared_CYP3A4_5_probe_not_CYP3A5_specific": True,
        "n_desmethyl_formation_is_outcome_not_predose_predictor": True,
        "outcome_not_untargeted_endogenous_metabolomics": True,
        "heldout_uncertainty_calibrated": False,
        "in_vivo_PopPK_or_tacrolimus_prediction_demonstrated": False,
        "individual_source_rows_or_donor_ids_exported": False,
    }
