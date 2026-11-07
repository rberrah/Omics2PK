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
    # A tier that is not confirmatory-admissible must explain the ruling needed.
    assert "public_registered" not in CONFIRMATORY_TIERS
    assert entry["access_note"]


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
