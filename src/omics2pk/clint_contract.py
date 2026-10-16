"""Strict public-research donor/CLint/proteome input contract; never fits a PK model.

This module accepts only curated, documented local records from a legally
accessible dataset. Technical assay rows are not independent patients, and
HLM vs HH CLint units are NOT automatically interchangeable.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isfinite
import re
from typing import Mapping, Sequence

DRUG_TO_CYP = {
    "bupropion": "CYP2B6",
    "diclofenac": "CYP2C9",
    "omeprazole": "CYP2C19",
    "bufuralol": "CYP2D6",
    "midazolam": "CYP3A4",
}
SYSTEMS = {"HH", "HLM"}
_SHA = re.compile(r"^[0-9a-f]{64}$")


class CLintContractError(ValueError):
    """Input isn't sufficiently qualified for donor-appropriate analysis."""


@dataclass(frozen=True)
class Measurement:
    donor_id: str
    drug: str
    system: str
    clint: float
    clint_unit: str
    clint_basis: str
    measurement_status: str
    source_table: str
    source_sha256: str


@dataclass(frozen=True)
class Protein:
    donor_id: str
    system: str
    protein_id: str
    abundance: float
    abundance_unit: str
    source_table: str
    source_sha256: str


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise CLintContractError(f"{label}: must be nonempty text")
    return value.strip()


def _hash(value: str) -> str:
    if not _SHA.fullmatch(_text(value, "source_sha256")):
        raise CLintContractError("Source SHA-256 must be explicit and 64 lowercase hex digits")
    return value


def _number(value: object, label: str) -> float:
    try:
        n = float(value)
    except (ValueError, TypeError):
        raise CLintContractError(f"{label}: not numeric") from None
    if not isfinite(n) or n <= 0:
        raise CLintContractError(f"{label}: require finite strictly positive measured value")
    return n


def read_measurement(row: Mapping[str, object]) -> Measurement:
    m = Measurement(
        donor_id=_text(row.get("donor_id"), "donor_id"),
        drug=_text(row.get("drug"), "drug").lower(),
        system=_text(row.get("system"), "system").upper(),
        clint=_number(row.get("clint"), "clint"),
        clint_unit=_text(row.get("clint_unit"), "clint_unit"),
        clint_basis=_text(row.get("clint_basis"), "clint_basis").upper(),
        measurement_status=_text(row.get("measurement_status"), "measurement_status").upper(),
        source_table=_text(row.get("source_table"), "source_table"),
        source_sha256=_hash(row.get("source_sha256")),
    )
    if m.drug not in DRUG_TO_CYP:
        raise CLintContractError("Compound is not one of the five predeclared Wegler probes")
    if m.system not in SYSTEMS:
        raise CLintContractError("CLint system must explicitly be HH or HLM")
    if m.clint_basis not in {"UNBOUND", "TOTAL"}:
        raise CLintContractError("CLint basis must be UNBOUND or TOTAL, never inferred")
    if m.measurement_status != "MEASURED":
        raise CLintContractError("Inferred or imputed CLint values are not directly observed outcomes")
    return m


def read_protein(row: Mapping[str, object]) -> Protein:
    p = Protein(
        donor_id=_text(row.get("donor_id"), "donor_id"),
        system=_text(row.get("system"), "system").upper(),
        protein_id=_text(row.get("protein_id"), "protein_id"),
        abundance=_number(row.get("abundance"), "abundance"),
        abundance_unit=_text(row.get("abundance_unit"), "abundance_unit"),
        source_table=_text(row.get("source_table"), "source_table"),
        source_sha256=_hash(row.get("source_sha256")),
    )
    if p.system not in SYSTEMS:
        raise CLintContractError("Protein system must be HH or HLM (not liver homogenate)")
    return p


def qualify_pairing(
    measured_rows: Sequence[Mapping[str, object]],
    protein_rows: Sequence[Mapping[str, object]],
    *,
    target_system: str = "HH",
    feature_system: str = "HH",
) -> dict:
    """Return aggregate QC only; no donor-linked values or patient data are persisted."""
    if target_system not in SYSTEMS or feature_system not in SYSTEMS:
        raise CLintContractError("Systems must be explicit HH/HLM")
    measurements = [read_measurement(r) for r in measured_rows]
    proteins = [read_protein(r) for r in protein_rows]
    if not measurements or not proteins:
        raise CLintContractError("Both directly measured CLint and proteomic features are required")
    keys = [(m.donor_id, m.drug, m.system) for m in measurements]
    if len(keys) != len(set(keys)):
        raise CLintContractError("Duplicate donor/drug/system CLint rows: aggregate technical replicates with provenance")
    proteome_keys = [(p.donor_id, p.system, p.protein_id) for p in proteins]
    if len(proteome_keys) != len(set(proteome_keys)):
        raise CLintContractError("Duplicate donor/system/protein rows: technical replicates require documented aggregation")
    # Different assay systems and drugs can have fundamentally distinct CLint units.
    # Each drug/system must nevertheless be self-consistent; no silent conversions.
    units = {}
    for m in measurements:
        group = (m.drug, m.system)
        designation = (m.clint_unit, m.clint_basis)
        if group in units and units[group] != designation:
            raise CLintContractError("Inconsistent CLint unit/basis within same drug-system")
        units[group] = designation
    protein_units = {}
    for p in proteins:
        group = (p.system, p.protein_id)
        if group in protein_units and protein_units[group] != p.abundance_unit:
            raise CLintContractError("Inconsistent abundance units within protein/system")
        protein_units[group] = p.abundance_unit

    endpoint_donors = {m.donor_id for m in measurements if m.system == target_system}
    source_donors = {p.donor_id for p in proteins if p.system == feature_system}
    if not endpoint_donors:
        raise CLintContractError("Target assay system has no directly measured CLint")
    if not source_donors:
        raise CLintContractError("Requested proteomic compartment has no measured features")
    if endpoint_donors != source_donors:
        raise CLintContractError("Donor sets do not match exactly between CLint target and chosen proteome")
    if len(endpoint_donors) > 15:
        raise CLintContractError("Unexpected more than 15 original Wegler donors")
    per_drug = Counter(m.drug for m in measurements if m.system == target_system)
    # We don't assert complete 15x5; missing donor/drug endpoints are real QC.
    return {
        "status": "CURATED_DATA_STRUCTURE_ONLY_NO_PK_MODEL_FITTED",
        "n_independent_donors": len(endpoint_donors),
        "n_measured_endpoint_rows": sum(m.system == target_system for m in measurements),
        "n_protein_feature_rows": sum(p.system == feature_system for p in proteins),
        "n_target_drugs": len(per_drug),
        "endpoint_counts_per_drug": dict(sorted(per_drug.items())),
        "target_system": target_system,
        "feature_system": feature_system,
        "strictly_positive_measured_clint": True,
        "no_system_unit_conversion_automated": True,
        "donor_ids_excluded_from_report": True,
        "missingness_and_censoring_require_external_source_review": True,
    }


def leave_one_donor_out(measured_rows: Sequence[Mapping[str, object]]) -> tuple:
    """Pure split indices; excludes all assay systems and five probe drugs per donor."""
    rows = [read_measurement(r) for r in measured_rows]
    if len(rows) != len({(r.donor_id, r.drug, r.system) for r in rows}):
        raise CLintContractError("Duplicate measured outcome rows")
    donors = sorted({r.donor_id for r in rows})
    if len(donors) < 3:
        raise CLintContractError("At least 3 donors needed to define exploratory donor-held-out folds")
    return tuple(
        (tuple(i for i, row in enumerate(rows) if row.donor_id != donor),
         tuple(i for i, row in enumerate(rows) if row.donor_id == donor))
        for donor in donors
    )
