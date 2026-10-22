"""Private same-donor gate for original HLC Synapse expression and phenotype.

Preferred HLC route:
  syn88644 curatedExpressionLiver -> explicit source donor key
  syn4629  curatedPhenotype       -> same source donor key + CYP activity

The gate never joins by row order and never prints donor identifiers or
individual measurements. It returns an in-memory private join object for an
authorized local benchmark plus a separate aggregate-only summary.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
import re
from typing import Sequence

EXPRESSION_ENTITY = "syn88644"
PHENOTYPE_ENTITY = "syn4629"
CURRENT_EXPRESSION_ENTITY_SAMPLES = 465
HASH = re.compile(r"^[0-9a-f]{64}$")
ALLOWED_STATUS = {"observed", "blq", "missing", "not_measured"}

PHENOTYPE_FIELDS = (
    "donor_id",
    "phenotype_entity_id",
    "assay_id",
    "activity_value",
    "activity_unit",
    "activity_status",
    "site_id",
    "age",
    "sex",
)


class HLCSynapseJoinError(ValueError):
    pass


@dataclass(frozen=True)
class PrivateHLCActivityJoin:
    """Authorized in-memory join. Do not serialize to the public repository."""

    expression_indices: tuple[int, ...]
    donor_ids: tuple[str, ...]
    site_ids: tuple[str, ...]
    ages: tuple[float, ...]
    sexes: tuple[float, ...]
    activities: tuple[float, ...]
    assay_id: str
    activity_unit: str
    expression_sha256: str
    phenotype_sha256: str


def _text(value, label):
    text = str(value).strip()
    if not text or len(text) > 200:
        raise HLCSynapseJoinError(f"Invalid {label}")
    return text


def _numeric(value, label):
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise HLCSynapseJoinError(f"Invalid numeric {label}") from exc
    if not isfinite(number):
        raise HLCSynapseJoinError(f"Nonfinite {label}")
    return number


def join_original_hlc_synapse(
    expression_donor_ids: Sequence[str],
    phenotype_rows: Sequence[dict],
    *,
    assay_id: str,
    expression_sha256: str,
    phenotype_sha256: str,
) -> PrivateHLCActivityJoin:
    """Create an explicit donor-key join for one measured CYP endpoint."""
    assay_id = _text(assay_id, "assay_id")
    if not HASH.fullmatch(str(expression_sha256)) or not HASH.fullmatch(str(phenotype_sha256)):
        raise HLCSynapseJoinError("Original source SHA-256 values are required")
    expression_ids = tuple(_text(x, "expression donor ID") for x in expression_donor_ids)
    if (
        len(expression_ids) != CURRENT_EXPRESSION_ENTITY_SAMPLES
        or len(set(expression_ids)) != len(expression_ids)
    ):
        raise HLCSynapseJoinError(
            f"Require exactly {CURRENT_EXPRESSION_ENTITY_SAMPLES} unique donors "
            f"from current {EXPRESSION_ENTITY} expression entity"
        )
    if not phenotype_rows:
        raise HLCSynapseJoinError("Original phenotype rows are missing")

    expression_pos = {donor: i for i, donor in enumerate(expression_ids)}
    selected = {}
    units = set()
    for row in phenotype_rows:
        if set(row) != set(PHENOTYPE_FIELDS):
            raise HLCSynapseJoinError("Require exact HLC phenotype extraction schema")
        if _text(row["phenotype_entity_id"], "phenotype entity") != PHENOTYPE_ENTITY:
            raise HLCSynapseJoinError("Phenotype rows must originate from syn4629")
        if _text(row["assay_id"], "assay ID") != assay_id:
            continue
        donor = _text(row["donor_id"], "phenotype donor ID")
        if donor in selected:
            raise HLCSynapseJoinError("Duplicate donor for selected HLC CYP assay")
        status = _text(row["activity_status"], "activity status").lower()
        if status not in ALLOWED_STATUS:
            raise HLCSynapseJoinError("Invalid HLC activity measurement status")
        unit = _text(row["activity_unit"], "activity unit")
        units.add(unit)
        raw = str(row["activity_value"]).strip()
        if status != "observed":
            if raw:
                raise HLCSynapseJoinError(
                    "Non-observed CYP activity must not carry a numeric replacement"
                )
            continue
        activity = _numeric(raw, "activity")
        if activity < 0:
            raise HLCSynapseJoinError("Observed CYP activity cannot be negative")
        age = _numeric(row["age"], "age")
        if not 0 <= age <= 120:
            raise HLCSynapseJoinError("Age outside prespecified human range")
        sex_label = _text(row["sex"], "sex").lower()
        if sex_label not in {"female", "male"}:
            raise HLCSynapseJoinError("Sex must be female/male in original source extraction")
        site = _text(row["site_id"], "site")
        selected[donor] = {
            "activity": activity,
            "age": age,
            "sex": 1.0 if sex_label == "male" else 0.0,
            "site": site,
        }

    if len(units) != 1:
        raise HLCSynapseJoinError("Selected CYP activity must use one physical unit")
    matched = [donor for donor in expression_ids if donor in selected]
    if len(matched) < 60:
        raise HLCSynapseJoinError("Fewer than 60 same-donor expression/activity observations")
    sites = sorted({selected[d]["site"] for d in matched})
    if len(sites) < 3:
        raise HLCSynapseJoinError("Need at least three original HLC sites for primary validation")
    if any(sum(selected[d]["site"] == site for d in matched) < 5 for site in sites):
        raise HLCSynapseJoinError("Each original HLC site needs at least five matched donors")

    return PrivateHLCActivityJoin(
        expression_indices=tuple(expression_pos[d] for d in matched),
        donor_ids=tuple(matched),
        site_ids=tuple(selected[d]["site"] for d in matched),
        ages=tuple(float(selected[d]["age"]) for d in matched),
        sexes=tuple(float(selected[d]["sex"]) for d in matched),
        activities=tuple(float(selected[d]["activity"]) for d in matched),
        assay_id=assay_id,
        activity_unit=next(iter(units)),
        expression_sha256=str(expression_sha256),
        phenotype_sha256=str(phenotype_sha256),
    )


def summarize_private_hlc_join(join: PrivateHLCActivityJoin) -> dict:
    """Return aggregate source readiness only; never individual IDs/measurements."""
    site_counts = {
        site: join.site_ids.count(site) for site in sorted(set(join.site_ids))
    }
    return {
        "status": "QUALIFIED_DIRECT_HLC_SYNAPSE_SAME_DONOR_JOIN",
        "expression_entity_id": EXPRESSION_ENTITY,
        "phenotype_entity_id": PHENOTYPE_ENTITY,
        "expression_entity_expected_unique_donors": CURRENT_EXPRESSION_ENTITY_SAMPLES,
        "selected_assay_id": join.assay_id,
        "activity_unit": join.activity_unit,
        "n_same_donor_complete_observations": len(join.donor_ids),
        "n_original_sites": len(site_counts),
        "site_complete_case_counts": site_counts,
        "expression_source_sha256": join.expression_sha256,
        "phenotype_source_sha256": join.phenotype_sha256,
        "join_key": "explicit original HLC donor_id; never row order",
        "ready_for_site_disjoint_E2_benchmark": True,
        "GSE9588_crosswalk_required": False,
        "individual_donor_ids_or_activity_values_exported": False,
        "not_clinical_PK": True,
    }
