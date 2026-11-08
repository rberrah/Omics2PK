"""Validator for PUBLIC_DATA_EVIDENCE_REGISTRY.

The registry records which datasets may enter a confirmatory Omics2PK
analysis. This module enforces the registry's own rules so that an entry
cannot quietly claim provenance that was never checked, nor claim an
objective its endpoint class cannot support.

No network access. No patient-level values are read or returned: the
registry holds provenance metadata and aggregate counts only.
"""
from __future__ import annotations

import json
from pathlib import Path

REGISTRY_PATH = Path(__file__).resolve().parents[2] / "registry" / "PUBLIC_DATA_EVIDENCE_REGISTRY.json"

STATUSES = ("ACQUIRED", "QUALIFIED", "PARTIAL", "ABSTAIN")

ACCESS_TIERS = (
    "open_public",
    "public_registered",
    "author_request",
    "governed_application",
    "unavailable",
)

#: Tiers from which data may enter a confirmatory analysis in this phase.
#: ``public_registered`` was admitted by the ruling of 2026-11-07, recorded in
#: the registry ``amendments`` block, but only under the access-requirement
#: rule below: the tier name alone never settles admissibility.
CONFIRMATORY_TIERS = ("open_public", "public_registered")

#: How a repository gates a download, which is what actually decides
#: admissibility under the public-only rule. ``self_sign`` is a click-through
#: acceptance of data-use terms: it is granted by the act of accepting, with no
#: human deciding, so it is not an authorization request. ``managed_act`` is a
#: committee decision on a submitted application, which is.
ACCESS_REQUIREMENT_TYPES = ("none", "self_sign", "certified_user", "managed_act")

#: Requirement types a requester can satisfy unilaterally, with no third party
#: reviewing, deciding or being contacted.
SELF_SERVICE_REQUIREMENT_TYPES = ("none", "self_sign", "certified_user")

#: Tiers that must declare how the download is gated, because for them the
#: answer is not implied by the tier name.
TIERS_REQUIRING_ACCESS_REQUIREMENT = ("public_registered",)

ENDPOINT_CLASSES = (
    "pk_clinical",
    "functional_in_vitro",
    "pd_cell_line",
    "biomarker_only",
)

OBJECTIVES = ("H1", "H2", "H3", "H4", "H5")

#: H1 is the incremental-PK-prediction hypothesis. An endpoint that is not a
#: human pharmacokinetic observation cannot satisfy it, however good the model.
OBJECTIVE_REQUIRED_ENDPOINTS = {"H1": ("pk_clinical",)}

SHA256_LEN = 64


class RegistryError(ValueError):
    """Raised when the registry violates one of its own rules."""


def load_registry(path=None):
    """Read the registry JSON and return it as a dict."""
    path = Path(path) if path is not None else REGISTRY_PATH
    if not path.is_file():
        raise RegistryError(f"registry file not found: {path}")
    with path.open(encoding="utf-8") as handle:
        registry = json.load(handle)
    if not isinstance(registry, dict):
        raise RegistryError("registry root must be a JSON object")
    return registry


def _require(condition, message):
    if not condition:
        raise RegistryError(message)


def _check_hex_sha(value, where):
    _require(
        isinstance(value, str) and len(value) == SHA256_LEN,
        f"{where}: sha256 must be {SHA256_LEN} hex characters",
    )
    _require(
        all(char in "0123456789abcdef" for char in value),
        f"{where}: sha256 must be lowercase hexadecimal",
    )


def _check_artifacts(entry, key):
    artifacts = entry.get("local_artifacts") or []
    _require(isinstance(artifacts, list), f"{key}: local_artifacts must be a list")
    for artifact in artifacts:
        _require(isinstance(artifact, dict), f"{key}: each artifact must be an object")
        for field in ("filename", "sha256"):
            _require(field in artifact, f"{key}: artifact missing '{field}'")
        _check_hex_sha(artifact["sha256"], f"{key}/{artifact['filename']}")
    return artifacts


def validate_entry(entry):
    """Validate one dataset entry and return its key."""
    _require(isinstance(entry, dict), "each dataset entry must be an object")
    key = entry.get("key")
    _require(isinstance(key, str) and key, "every dataset entry needs a non-empty 'key'")

    status = entry.get("status")
    _require(status in STATUSES, f"{key}: status must be one of {STATUSES}, got {status!r}")

    tier = entry.get("access_tier")
    _require(tier in ACCESS_TIERS, f"{key}: access_tier must be one of {ACCESS_TIERS}, got {tier!r}")

    verified = entry.get("verified_this_session")
    _require(isinstance(verified, bool), f"{key}: verified_this_session must be a boolean")

    # A verified claim must say how it was verified; an unverified one must
    # name where it came from, so no provenance is asserted without a trail.
    if verified:
        _require(
            isinstance(entry.get("verification_method"), str) and entry["verification_method"].strip(),
            f"{key}: verified_this_session is true so verification_method is required",
        )
    else:
        _require(
            entry.get("inherited_from") or entry.get("note"),
            f"{key}: unverified entries must record inherited_from or a note",
        )

    # Abstaining is a scientific result and needs a stated reason.
    if status == "ABSTAIN":
        _require(
            isinstance(entry.get("abstain_reason"), str) and entry["abstain_reason"].strip(),
            f"{key}: ABSTAIN requires abstain_reason",
        )
        _require(
            not entry.get("objectives_supported"),
            f"{key}: an ABSTAIN entry must not claim supported objectives",
        )

    # Checksums are validated first, then the scientific claims, then the
    # bookkeeping requirement that bytes in hand carry a recorded checksum.
    artifacts = _check_artifacts(entry, key)

    objectives = entry.get("objectives_supported") or []
    _require(isinstance(objectives, list), f"{key}: objectives_supported must be a list")
    for objective in objectives:
        _require(objective in OBJECTIVES, f"{key}: unknown objective {objective!r}")

    endpoint = entry.get("endpoint")
    if endpoint is not None:
        _require(isinstance(endpoint, dict), f"{key}: endpoint must be an object")
        endpoint_class = endpoint.get("endpoint_class")
        _require(
            endpoint_class in ENDPOINT_CLASSES,
            f"{key}: endpoint_class must be one of {ENDPOINT_CLASSES}, got {endpoint_class!r}",
        )
        for objective, allowed in OBJECTIVE_REQUIRED_ENDPOINTS.items():
            if objective in objectives:
                _require(
                    endpoint_class in allowed,
                    f"{key}: claims {objective} but its endpoint_class {endpoint_class!r} "
                    f"is not in {allowed}; a non-PK endpoint cannot satisfy {objective}",
                )
    else:
        _require(
            not objectives,
            f"{key}: objectives_supported requires a declared endpoint",
        )

    # An entry may only be used confirmatorily from an admissible tier.
    if objectives and tier not in CONFIRMATORY_TIERS:
        _require(
            isinstance(entry.get("access_note"), str) and entry["access_note"].strip(),
            f"{key}: tier {tier!r} is not confirmatory-admissible, so an access_note "
            "explaining the ruling required is mandatory",
        )

    # How the download is gated decides admissibility, not the tier name. A
    # registered-access source must state the gate and prove it was inspected,
    # so that a committee-reviewed source can never be filed as self-service.
    requirement = entry.get("access_requirement")
    if tier in TIERS_REQUIRING_ACCESS_REQUIREMENT:
        _require(
            isinstance(requirement, dict),
            f"{key}: tier {tier!r} must declare an access_requirement object "
            "stating how the download is gated",
        )
        req_type = requirement.get("type")
        _require(
            req_type in ACCESS_REQUIREMENT_TYPES,
            f"{key}: access_requirement.type must be one of "
            f"{ACCESS_REQUIREMENT_TYPES}, got {req_type!r}",
        )
        _require(
            req_type in SELF_SERVICE_REQUIREMENT_TYPES,
            f"{key}: access_requirement.type {req_type!r} is a third-party decision, "
            f"so the tier must be 'governed_application', not {tier!r}",
        )
        _require(
            isinstance(requirement.get("verified_on"), str)
            and requirement["verified_on"].strip(),
            f"{key}: access_requirement must record verified_on, the date the gate "
            "was actually inspected",
        )
        _require(
            isinstance(requirement.get("evidence"), str) and requirement["evidence"].strip(),
            f"{key}: access_requirement must record evidence, how the gate type was "
            "established, since the admissibility ruling rests on it",
        )
    elif isinstance(requirement, dict):
        req_type = requirement.get("type")
        _require(
            req_type in ACCESS_REQUIREMENT_TYPES,
            f"{key}: access_requirement.type must be one of "
            f"{ACCESS_REQUIREMENT_TYPES}, got {req_type!r}",
        )
        if tier == "open_public":
            _require(
                req_type == "none",
                f"{key}: tier 'open_public' means no gate at all, so "
                f"access_requirement.type {req_type!r} contradicts it",
            )

    # A functional analogue of H1 must be labelled as such, never as H1.
    if entry.get("h1_functional_analogue"):
        endpoint_class = (entry.get("endpoint") or {}).get("endpoint_class")
        _require(
            endpoint_class in ("functional_in_vitro", "pd_cell_line"),
            f"{key}: h1_functional_analogue is only meaningful for a non-PK endpoint",
        )
        _require(
            "H1" not in objectives,
            f"{key}: an entry cannot claim both H1 and h1_functional_analogue",
        )

    # Claiming bytes in hand requires a checksum.
    if status in ("ACQUIRED", "QUALIFIED") and tier in CONFIRMATORY_TIERS:
        _require(
            artifacts or entry.get("inherited_from") or not entry.get("verified_this_session"),
            f"{key}: {status} on a confirmatory-admissible source verified here "
            "must record local_artifacts",
        )

    return key


KNOWLEDGE_SOURCE_ROLES = ("functional_network", "pathway_annotation",
                          "chemical_gene_interaction")


def validate_knowledge_source(entry):
    """Validate one prior-knowledge resource and return its key.

    A knowledge source is not a cohort: it has no participant, no outcome and
    no endpoint, so it can neither support nor refute an objective by itself.
    What must be true of it is narrower and is checked here. It has to be
    retrievable without a third-party decision, carry a verified checksum, and
    state its redistribution position explicitly -- with an unverified position
    treated as a refusal, so silence can never be read as permission.
    """
    for field in ("key", "label", "role", "url", "access_requirement",
                  "licence", "redistribution_permitted",
                  "committed_to_this_repository", "local_artifacts",
                  "used_for"):
        _require(field in entry, f"knowledge source missing field '{field}'")
    key = entry["key"]
    _require(isinstance(key, str) and key, "knowledge source needs a key")
    _require(entry["role"] in KNOWLEDGE_SOURCE_ROLES,
             f"{key}: role must be one of {KNOWLEDGE_SOURCE_ROLES}")

    requirement = entry["access_requirement"]
    _require(isinstance(requirement, dict), f"{key}: access_requirement must be an object")
    for field in ("type", "evidence", "verified_on"):
        _require(field in requirement, f"{key}: access_requirement missing '{field}'")
    _require(requirement["type"] in ACCESS_REQUIREMENT_TYPES,
             f"{key}: access requirement type must be one of {ACCESS_REQUIREMENT_TYPES}")
    _require(
        requirement["type"] in SELF_SERVICE_REQUIREMENT_TYPES,
        f"{key}: a knowledge source behind a third-party decision is not admissible",
    )

    _require(isinstance(entry["licence"], str) and entry["licence"],
             f"{key}: the licence must be stated, 'not verified' included")
    permitted = entry["redistribution_permitted"]
    _require(permitted in (True, False, None),
             f"{key}: redistribution_permitted must be true, false or null")
    _require(entry["committed_to_this_repository"] is False,
             f"{key}: source artefacts are not committed to this repository")
    if permitted is not True:
        _require(
            "redistribution_note" in entry,
            f"{key}: redistribution is unproven, so a note must say how that is handled",
        )
    _check_artifacts(entry, key)
    _require(isinstance(entry["used_for"], list) and entry["used_for"],
             f"{key}: 'used_for' must say what the resource is used for")
    _require(
        not entry.get("objectives_supported"),
        f"{key}: a knowledge source carries no endpoint and cannot support an objective",
    )
    return key


def validate_registry(path=None):
    """Validate the whole registry; return an aggregate summary.

    Raises :class:`RegistryError` on the first violation found.
    """
    registry = load_registry(path)
    for field in ("registry_name", "version", "frozen_at", "access_policy", "status_definitions"):
        _require(field in registry, f"registry missing top-level field '{field}'")
    _require(
        registry["registry_name"] == "PUBLIC_DATA_EVIDENCE_REGISTRY",
        "registry_name must be PUBLIC_DATA_EVIDENCE_REGISTRY",
    )

    datasets = registry.get("datasets")
    _require(isinstance(datasets, list) and datasets, "registry needs a non-empty 'datasets' list")

    keys = []
    for entry in datasets:
        keys.append(validate_entry(entry))
    duplicates = sorted({key for key in keys if keys.count(key) > 1})
    _require(not duplicates, f"duplicate dataset keys: {duplicates}")

    knowledge = registry.get("knowledge_sources") or []
    _require(isinstance(knowledge, list), "'knowledge_sources' must be a list")
    knowledge_keys = [validate_knowledge_source(entry) for entry in knowledge]
    knowledge_duplicates = sorted(
        {key for key in knowledge_keys if knowledge_keys.count(key) > 1})
    _require(not knowledge_duplicates,
             f"duplicate knowledge source keys: {knowledge_duplicates}")
    overlap = sorted(set(knowledge_keys) & set(keys))
    _require(not overlap, f"keys used for both a dataset and a knowledge source: {overlap}")
    knowledge_keys = sorted(knowledge_keys)

    by_status = {status: sorted(e["key"] for e in datasets if e["status"] == status) for status in STATUSES}
    by_tier = {tier: sorted(e["key"] for e in datasets if e["access_tier"] == tier) for tier in ACCESS_TIERS}

    confirmatory = sorted(
        e["key"]
        for e in datasets
        if e["access_tier"] in CONFIRMATORY_TIERS
        and e["status"] in ("ACQUIRED", "QUALIFIED")
        and e.get("objectives_supported")
    )
    pk_clinical = sorted(
        e["key"]
        for e in datasets
        if (e.get("endpoint") or {}).get("endpoint_class") == "pk_clinical"
        and e["status"] in ("ACQUIRED", "QUALIFIED", "PARTIAL")
    )
    h1_eligible = sorted(
        e["key"]
        for e in datasets
        if "H1" in (e.get("objectives_supported") or [])
        and e["access_tier"] in CONFIRMATORY_TIERS
        and e["status"] == "QUALIFIED"
    )
    # Admissible sources that are nonetheless gated. Reported separately so a
    # reader can see at a glance which conclusions depend on a source the
    # reviewer must accept terms to retrieve, rather than fetch by URL.
    gated_admissible = sorted(
        e["key"]
        for e in datasets
        if e["key"] in confirmatory
        and (e.get("access_requirement") or {}).get("type", "none") != "none"
    )

    return {
        "registry_name": registry["registry_name"],
        "version": registry["version"],
        "frozen_at": registry["frozen_at"],
        "n_datasets": len(datasets),
        "by_status": by_status,
        "by_access_tier": by_tier,
        "confirmatory_admissible": confirmatory,
        "gated_but_admissible": gated_admissible,
        "clinical_pk_endpoint_datasets": pk_clinical,
        "h1_eligible_datasets": h1_eligible,
        "verified_this_session": sorted(e["key"] for e in datasets if e.get("verified_this_session")),
        "unverified_inherited": sorted(e["key"] for e in datasets if not e.get("verified_this_session")),
        "knowledge_sources": knowledge_keys,
        "knowledge_sources_redistributable": sorted(
            e["key"] for e in knowledge
            if e["redistribution_permitted"] is True),
    }


if __name__ == "__main__":  # pragma: no cover - manual inspection helper
    print(json.dumps(validate_registry(), indent=2, ensure_ascii=False, sort_keys=True))
