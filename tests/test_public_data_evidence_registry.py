"""Contract tests for PUBLIC_DATA_EVIDENCE_REGISTRY.

These tests assert the governance rules of the registry, not model
performance. They must fail if someone adds a dataset that claims an
objective its endpoint cannot support, asserts provenance without a
verification trail, or marks an authorization-gated source as usable.
"""
from __future__ import annotations

import copy

import pytest

from omics2pk.public_data_evidence_registry import (
    CONFIRMATORY_TIERS,
    SELF_SERVICE_REQUIREMENT_TYPES,
    RegistryError,
    load_registry,
    validate_entry,
    validate_registry,
)


@pytest.fixture(scope="module")
def registry():
    return load_registry()


@pytest.fixture(scope="module")
def summary():
    return validate_registry()


def test_registry_validates(summary):
    assert summary["registry_name"] == "PUBLIC_DATA_EVIDENCE_REGISTRY"
    assert summary["n_datasets"] >= 10


def test_ariana_is_the_only_h1_eligible_public_dataset(summary):
    """The core Phase 1 finding, enforced as a test.

    If a second public cohort with a clinical PK endpoint is ever acquired,
    this test must be updated deliberately rather than drift silently.
    """
    assert summary["h1_eligible_datasets"] == ["ariana_metoprolol_47"]


def test_ariana_provenance_is_pinned(registry):
    entry = next(d for d in registry["datasets"] if d["key"] == "ariana_metoprolol_47")
    assert entry["access_tier"] == "open_public"
    assert entry["status"] == "QUALIFIED"
    assert entry["verified_this_session"] is True
    artifact = entry["local_artifacts"][0]
    assert artifact["filename"] == "Solanidine_metoprolol_ratios_PGx.xlsx"
    assert artifact["sha256"] == (
        "7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37"
    )
    assert entry["endpoint"]["observations_per_subject"] == 1
    assert entry["endpoint"]["endpoint_class"] == "pk_clinical"
    # The trough must never be relabelled as an exposure metric.
    assert set(entry["endpoint"]["must_not_be_called"]) >= {"AUC", "clearance"}


def test_single_observation_blocks_mapbe_and_nlme(registry):
    entry = next(d for d in registry["datasets"] if d["key"] == "ariana_metoprolol_47")
    caveats = " ".join(entry["caveats"]).lower()
    assert "circular" in caveats
    assert "not identifiable" in caveats
    # The metabolic ratio shares a denominator with the target.
    assert "validation target" in caveats


def test_hlc_is_functional_not_pk_and_cannot_claim_h1(registry):
    entry = next(d for d in registry["datasets"] if d["key"] == "hlc_syn4499")
    assert entry["endpoint"]["endpoint_class"] == "functional_in_vitro"
    assert "H1" not in entry["objectives_supported"]
    assert entry["access_tier"] == "public_registered"
    assert entry["access_note"]
    # Admitted on 2026-11-07, but the endpoint ruling is untouched by that:
    # becoming retrievable never turns an in vitro endpoint into clinical PK.
    assert entry["h1_functional_analogue"] is True


def test_hlc_gate_is_self_service_and_evidenced(registry):
    """The ruling rests on the gate being self-signed, so the entry must prove it."""
    entry = next(d for d in registry["datasets"] if d["key"] == "hlc_syn4499")
    requirement = entry["access_requirement"]
    assert requirement["type"] == "self_sign"
    assert requirement["type"] in SELF_SERVICE_REQUIREMENT_TYPES
    assert requirement["committee_review"] is False
    assert requirement["author_contact_required"] is False
    assert requirement["is_certified_user_required"] is False
    assert requirement["verified_on"] == "2026-11-07"
    # The evidence must name the inspected object, not merely assert the type.
    assert "SelfSignAccessRequirement" in requirement["evidence"]
    assert requirement["synapse_access_requirement_id"] == 9604913


def test_registered_tier_is_admissible_but_reported_as_gated(summary):
    """Admissible and ungated are different claims and must stay distinguishable."""
    assert "public_registered" in CONFIRMATORY_TIERS
    assert "hlc_syn4499" in summary["confirmatory_admissible"]
    assert summary["gated_but_admissible"] == ["hlc_syn4499"]
    # ARIANA is retrievable by URL, so it must never appear as gated.
    assert "ariana_metoprolol_47" not in summary["gated_but_admissible"]


def test_amendment_records_that_the_ruling_preceded_the_runs(registry):
    """A threshold or admissibility change made after seeing a result is invalid."""
    amendment = next(a for a in registry["amendments"] if a["version"] == "1.1.0")
    assert amendment["made_before_any_run_on_the_affected_data"] is True
    assert amendment["affected_datasets"] == ["hlc_syn4499"]
    assert "H1" in amendment["does_not_change"]


def test_a_committee_gate_cannot_be_filed_as_registered_access(registry):
    """The rule that keeps the ruling narrow: managed review stays excluded."""
    entry = copy.deepcopy(
        next(d for d in registry["datasets"] if d["key"] == "hlc_syn4499")
    )
    entry["access_requirement"]["type"] = "managed_act"
    with pytest.raises(RegistryError, match="third-party decision"):
        validate_entry(entry)


def test_registered_access_must_declare_its_gate(registry):
    entry = copy.deepcopy(
        next(d for d in registry["datasets"] if d["key"] == "hlc_syn4499")
    )
    del entry["access_requirement"]
    with pytest.raises(RegistryError, match="must declare an access_requirement"):
        validate_entry(entry)


def test_gate_type_must_be_evidenced_not_merely_asserted(registry):
    entry = copy.deepcopy(
        next(d for d in registry["datasets"] if d["key"] == "hlc_syn4499")
    )
    entry["access_requirement"]["evidence"] = "   "
    with pytest.raises(RegistryError, match="must record evidence"):
        validate_entry(entry)


def test_open_public_cannot_also_declare_a_gate(registry):
    entry = copy.deepcopy(
        next(d for d in registry["datasets"] if d["key"] == "ariana_metoprolol_47")
    )
    entry["access_requirement"] = {"type": "self_sign"}
    with pytest.raises(RegistryError, match="means no gate at all"):
        validate_entry(entry)


def test_hlc_retrieval_links_cover_every_local_archive(registry):
    """Provenance is only reproducible if each local file names its source entity."""
    entry = next(d for d in registry["datasets"] if d["key"] == "hlc_syn4499")
    files = entry["retrieval"]["files"]
    assert {f["local_filename"] for f in files} == {
        a["filename"] for a in entry["local_artifacts"]
    }
    for record in files:
        assert record["synapse_id"].startswith("syn")
        assert record["synapse_id"] in record["url"]


def test_abstained_sources_claim_no_objectives(registry):
    abstained = [d for d in registry["datasets"] if d["status"] == "ABSTAIN"]
    assert abstained, "the registry should record abstentions explicitly"
    for entry in abstained:
        assert not entry.get("objectives_supported")
        assert entry["abstain_reason"].strip()


def test_solanidine_replication_cohorts_are_all_abstained(registry):
    """No public cohort can replicate the ARIANA solanidine->PK finding."""
    expected = {"muller_cocktail_20", "wollmann_839", "kiiski_314"}
    by_key = {d["key"]: d for d in registry["datasets"]}
    for key in expected:
        assert by_key[key]["status"] == "ABSTAIN"
    finding = registry["phase1_findings"]["no_public_replication_cohort_for_the_solanidine_axis"]
    assert expected <= set(finding["checked_and_abstained"])


def test_endpoint_classes_are_never_pooled(registry):
    """Datasets must be separable by endpoint class for reporting."""
    classes = {
        d["key"]: (d.get("endpoint") or {}).get("endpoint_class")
        for d in registry["datasets"]
        if d.get("endpoint")
    }
    assert classes["ariana_metoprolol_47"] == "pk_clinical"
    assert classes["nci60"] == "pd_cell_line"
    assert classes["wegler_15"] == "functional_in_vitro"


def test_rejects_non_pk_endpoint_claiming_h1():
    bad = {
        "key": "bad_entry",
        "status": "QUALIFIED",
        "access_tier": "open_public",
        "verified_this_session": True,
        "verification_method": "fabricated for the test",
        "objectives_supported": ["H1"],
        "endpoint": {"name": "cell viability", "endpoint_class": "pd_cell_line"},
    }
    with pytest.raises(RegistryError, match="cannot satisfy H1"):
        validate_entry(bad)


def test_rejects_verified_claim_without_method():
    bad = {
        "key": "bad_entry",
        "status": "ACQUIRED",
        "access_tier": "open_public",
        "verified_this_session": True,
    }
    with pytest.raises(RegistryError, match="verification_method is required"):
        validate_entry(bad)


def test_rejects_abstain_without_reason():
    bad = {
        "key": "bad_entry",
        "status": "ABSTAIN",
        "access_tier": "author_request",
        "verified_this_session": False,
        "note": "placeholder",
    }
    with pytest.raises(RegistryError, match="requires abstain_reason"):
        validate_entry(bad)


def test_rejects_malformed_checksum():
    bad = {
        "key": "bad_entry",
        "status": "ACQUIRED",
        "access_tier": "open_public",
        "verified_this_session": True,
        "verification_method": "fabricated for the test",
        "local_artifacts": [{"filename": "x.csv", "sha256": "not-a-hash"}],
    }
    with pytest.raises(RegistryError, match="sha256"):
        validate_entry(bad)


def test_rejects_duplicate_keys(registry, tmp_path):
    import json

    doubled = copy.deepcopy(registry)
    doubled["datasets"].append(copy.deepcopy(doubled["datasets"][0]))
    path = tmp_path / "dup.json"
    path.write_text(json.dumps(doubled), encoding="utf-8")
    with pytest.raises(RegistryError, match="duplicate dataset keys"):
        validate_registry(path)


def test_every_entry_declares_a_verification_trail(registry):
    for entry in registry["datasets"]:
        if entry.get("verified_this_session"):
            assert entry["verification_method"].strip()
        else:
            assert entry.get("inherited_from") or entry.get("note")
