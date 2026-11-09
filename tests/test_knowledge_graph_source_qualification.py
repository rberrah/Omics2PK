"""Contract tests for knowledge-source registration and qualification.

A prior-knowledge resource is not a cohort, and the rules that must hold for it
are different: it has to be retrievable without anyone's permission, it has to
carry a verified checksum, and its redistribution position has to be stated
rather than assumed. The point of these tests is that an unverified position is
treated as a refusal, so silence from a provider can never be read as consent.
"""
from __future__ import annotations

import copy
import gzip
import json
import unittest
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory

try:
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"pytest is required for this module: {exc}") from exc

from omics2pk.public_data_evidence_registry import (
    ACCESS_REQUIREMENT_TYPES,
    KNOWLEDGE_SOURCE_ROLES,
    RegistryError,
    load_registry,
    validate_knowledge_source,
    validate_registry,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
QUALIFICATION = REPO_ROOT / "docs" / "KNOWLEDGE_GRAPH_SOURCE_QUALIFICATION.json"

GOOD = {
    "key": "example_network",
    "label": "An example network",
    "role": "functional_network",
    "url": "https://example.org/network.txt.gz",
    "access_requirement": {
        "type": "none",
        "evidence": "Retrieved anonymously by HTTPS.",
        "verified_on": "2026-11-08",
    },
    "licence": "CC BY 4.0",
    "licence_verified": "licence page read on 2026-11-08",
    "redistribution_permitted": True,
    "committed_to_this_repository": False,
    "local_artifacts": [{"filename": "network.txt.gz", "bytes": 10, "sha256": "ab" * 32}],
    "used_for": ["edges of the mechanistic network"],
    "redistribution_note": "not committed regardless",
}


class TestKnowledgeSourceValidation(unittest.TestCase):
    def entry(self, **overrides):
        record = copy.deepcopy(GOOD)
        record.update(overrides)
        return record

    def test_a_well_formed_source_validates(self):
        self.assertEqual(validate_knowledge_source(self.entry()), "example_network")

    def test_every_required_field_is_required(self):
        for field in ("key", "label", "role", "url", "access_requirement",
                      "licence", "redistribution_permitted",
                      "committed_to_this_repository", "local_artifacts",
                      "used_for"):
            record = self.entry()
            del record[field]
            with self.assertRaisesRegex(RegistryError, field):
                validate_knowledge_source(record)

    def test_the_role_must_be_a_known_one(self):
        with self.assertRaisesRegex(RegistryError, "role must be"):
            validate_knowledge_source(self.entry(role="whatever"))
        for role in KNOWLEDGE_SOURCE_ROLES:
            validate_knowledge_source(self.entry(role=role))

    def test_a_source_behind_a_third_party_decision_is_refused(self):
        gated = self.entry(access_requirement={
            "type": "managed_act",
            "evidence": "a committee reviews each request",
            "verified_on": "2026-11-08",
        })
        with self.assertRaisesRegex(RegistryError, "third-party decision"):
            validate_knowledge_source(gated)

    def test_self_service_requirements_remain_admissible(self):
        for kind in ("none", "self_sign", "certified_user"):
            self.assertIn(kind, ACCESS_REQUIREMENT_TYPES)
            validate_knowledge_source(self.entry(access_requirement={
                "type": kind, "evidence": "self-service", "verified_on": "2026-11-08"}))

    def test_an_unknown_requirement_type_is_refused(self):
        with self.assertRaisesRegex(RegistryError, "requirement type"):
            validate_knowledge_source(self.entry(access_requirement={
                "type": "ask_nicely", "evidence": "x", "verified_on": "2026-11-08"}))

    def test_the_access_requirement_must_carry_evidence_and_a_date(self):
        for field in ("type", "evidence", "verified_on"):
            requirement = dict(GOOD["access_requirement"])
            del requirement[field]
            with self.assertRaisesRegex(RegistryError, field):
                validate_knowledge_source(self.entry(access_requirement=requirement))

    def test_an_unproven_redistribution_right_needs_an_explicit_note(self):
        record = self.entry(redistribution_permitted=None)
        del record["redistribution_note"]
        with self.assertRaisesRegex(RegistryError, "redistribution is unproven"):
            validate_knowledge_source(record)
        record["redistribution_note"] = "treated as refused; not committed"
        self.assertEqual(validate_knowledge_source(record), "example_network")

    def test_a_refused_redistribution_right_needs_the_same_note(self):
        record = self.entry(redistribution_permitted=False)
        del record["redistribution_note"]
        with self.assertRaisesRegex(RegistryError, "redistribution is unproven"):
            validate_knowledge_source(record)

    def test_redistribution_must_be_a_three_valued_decision(self):
        with self.assertRaisesRegex(RegistryError, "true, false or null"):
            validate_knowledge_source(self.entry(redistribution_permitted="maybe"))

    def test_the_licence_must_be_stated_even_when_unknown(self):
        with self.assertRaisesRegex(RegistryError, "licence must be stated"):
            validate_knowledge_source(self.entry(licence=""))
        validate_knowledge_source(self.entry(licence="not verified",
                                             redistribution_permitted=None))

    def test_a_committed_artefact_is_refused(self):
        with self.assertRaisesRegex(RegistryError, "not committed"):
            validate_knowledge_source(self.entry(committed_to_this_repository=True))

    def test_a_checksum_must_be_a_full_digest(self):
        with self.assertRaisesRegex(RegistryError, "sha256"):
            validate_knowledge_source(self.entry(
                local_artifacts=[{"filename": "n.gz", "sha256": "abc"}]))

    def test_a_knowledge_source_may_not_claim_to_support_an_objective(self):
        with self.assertRaisesRegex(RegistryError, "cannot support an objective"):
            validate_knowledge_source(self.entry(objectives_supported=["H2"]))

    def test_the_purpose_must_be_recorded(self):
        with self.assertRaisesRegex(RegistryError, "used_for"):
            validate_knowledge_source(self.entry(used_for=[]))


class TestTheCommittedRegistry(unittest.TestCase):
    def test_the_registry_validates_as_committed(self):
        summary = validate_registry()
        self.assertIn("knowledge_sources", summary)

    def test_the_four_sources_the_selector_needs_are_registered(self):
        # Checked by key rather than by count: the registry also holds the
        # annotation references that decide what a measured feature is, and a
        # count would have to be edited every time one is added.
        summary = validate_registry()
        self.assertGreaterEqual(len(summary["knowledge_sources"]), 4)
        entries = {entry["key"]: entry
                   for entry in load_registry()["knowledge_sources"]}
        for key, role in (("string_v12_human_links", "functional_network"),
                          ("string_v12_human_protein_info", "functional_network"),
                          ("reactome_pathways_gmt", "pathway_annotation"),
                          ("ctd_chemical_gene_interactions",
                           "chemical_gene_interaction")):
            self.assertIn(key, entries)
            self.assertEqual(entries[key]["role"], role, key)
        roles = {entry["role"] for entry in entries.values()}
        self.assertEqual(roles, set(KNOWLEDGE_SOURCE_ROLES))

    def test_the_annotation_references_behind_the_pools_are_registered(self):
        # They were used before being registered. The pool file carried their
        # checksums while the registry, which is the source of truth for public
        # data provenance, did not list what decides every feature's identity.
        entries = {entry["key"]: entry
                   for entry in load_registry()["knowledge_sources"]}
        for key, role in (("geo_gpl4372_platform_record", "probe_annotation"),
                          ("ncbi_human_gene_reference", "gene_reference"),
                          ("ncbi_gene2accession_human_subset", "gene_reference")):
            self.assertIn(key, entries)
            self.assertEqual(entries[key]["role"], role, key)
            self.assertTrue(entries[key]["used_for"], key)

    def test_the_registered_checksums_are_the_ones_the_pools_were_built_from(self):
        pools = json.loads(
            (REPO_ROOT / "registry" / "QUALIFIED_PROBE_POOLS.json")
            .read_text(encoding="utf-8"))
        registered = {artifact["sha256"]
                      for entry in load_registry()["knowledge_sources"]
                      for artifact in entry["local_artifacts"]}
        for field in ("platform_record", "human_gene_reference",
                      "accession_to_gene_map"):
            self.assertIn(pools["provenance"][field]["sha256"], registered, field)

    def test_no_registered_knowledge_source_is_committed(self):
        for entry in load_registry()["knowledge_sources"]:
            self.assertFalse(entry["committed_to_this_repository"], entry["key"])
            for artifact in entry["local_artifacts"]:
                self.assertFalse((REPO_ROOT / artifact["filename"]).exists(),
                                 msg=f"{artifact['filename']} must not be in the tree")

    def test_the_source_whose_terms_are_unverified_says_so(self):
        unverified = [entry for entry in load_registry()["knowledge_sources"]
                      if entry["redistribution_permitted"] is not True]
        self.assertTrue(unverified, "expected at least one unverified source")
        for entry in unverified:
            self.assertIn("not verified", entry["licence"].lower() + " "
                          + entry["licence_verified"].lower())
            self.assertIn("not permission", entry["redistribution_note"].lower())

    def test_the_amendment_records_the_finding_before_any_endpoint_was_read(self):
        # Found by version, not as "the last one": later amendments are expected
        # and must not be able to hide this finding by appending to the list.
        amendments = {entry["version"]: entry
                      for entry in load_registry()["amendments"]}
        self.assertIn("1.9.0", amendments)
        recorded = amendments["1.9.0"]
        self.assertTrue(recorded["made_before_any_run_on_the_affected_data"])
        joined = " ".join(recorded["corrections"]).lower()
        self.assertIn("cannot be given a calibrated significance threshold", joined)
        self.assertIn("ranking", joined)

    def test_the_registry_version_moved_with_the_amendment(self):
        registry = load_registry()
        self.assertEqual(registry["version"], registry["amendments"][-1]["version"])


class TestQualificationReport(unittest.TestCase):
    def report(self):
        return json.loads(QUALIFICATION.read_text(encoding="utf-8"))

    def test_the_report_exists_and_covers_four_sources(self):
        self.assertTrue(QUALIFICATION.is_file())
        self.assertEqual(len(self.report()["sources"]), 4)

    def test_the_report_states_that_no_endpoint_was_read(self):
        self.assertTrue(self.report()["no_endpoint_was_read"])

    def test_every_source_carries_a_verified_checksum_and_a_structure_check(self):
        for key, entry in self.report()["sources"].items():
            self.assertEqual(len(entry["bytes_verified"]["sha256"]), 64, key)
            self.assertGreater(entry["bytes_verified"]["bytes"], 0, key)
            self.assertTrue(entry["structure_verified"], key)
            self.assertFalse(entry["committed_to_this_repository"], key)

    def test_the_report_and_the_registry_agree_on_every_checksum(self):
        # The report covers the selector's graph sources; the registry also
        # holds the annotation references used elsewhere. Every checksum the
        # report states must be registered, and the report must state one for
        # each of its own sources -- a subset, not an equality, so adding a
        # source to the registry cannot break this and dropping one from the
        # registry still does.
        report = self.report()
        reported = {entry["bytes_verified"]["sha256"]
                    for entry in report["sources"].values()}
        registered = {artifact["sha256"]
                      for entry in load_registry()["knowledge_sources"]
                      for artifact in entry["local_artifacts"]}
        self.assertTrue(reported)
        self.assertLessEqual(reported, registered)
        self.assertEqual(len(reported), len(report["sources"]))

    def test_the_frozen_selector_parameters_are_recorded_with_the_sources(self):
        frozen = self.report()["frozen_selector_parameters"]
        self.assertEqual(frozen["string_min_combined_score"], 700)
        self.assertEqual(frozen["string_version"], "v12.0")


class TestQualificationScriptRefusesBadInput(unittest.TestCase):
    """The script must fail loudly on a changed release, not adapt to it."""

    def script(self):
        import importlib.util
        path = REPO_ROOT / "scripts" / "qualify_public_knowledge_graph_sources.py"
        spec = importlib.util.spec_from_file_location("qualify_kg", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_a_missing_file_is_reported_not_skipped(self):
        module = self.script()
        with TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(module.KnowledgeSourceQualificationError,
                                        "absent"):
                module.verify_bytes(Path(tmp))

    def test_a_file_with_the_wrong_bytes_is_refused(self):
        module = self.script()
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            for record in module.EXPECTED.values():
                (root / record["filename"]).write_bytes(b"not the real file")
            with self.assertRaisesRegex(module.KnowledgeSourceQualificationError,
                                        "not the file the registry records"):
                module.verify_bytes(root)

    def test_a_changed_network_header_is_refused(self):
        module = self.script()
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "info.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write("#a_different_layout\n")
            with self.assertRaisesRegex(module.KnowledgeSourceQualificationError,
                                        "header changed"):
                module.verify_string_info(path)

    def test_a_truncated_pathway_file_is_refused(self):
        module = self.script()
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "pathways.gmt.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("ReactomePathways.gmt", "Only one\tR-1\tAAA1\n")
            with self.assertRaisesRegex(module.KnowledgeSourceQualificationError,
                                        "only 1 pathways"):
                module.verify_reactome(path)

    def test_a_missing_archive_member_is_refused(self):
        module = self.script()
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "pathways.gmt.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("something_else.txt", "x")
            with self.assertRaisesRegex(module.KnowledgeSourceQualificationError,
                                        "member absent"):
                module.verify_reactome(path)

    def test_changed_interaction_columns_are_refused(self):
        module = self.script()
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "ixns.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write("# ChemicalName\tGeneSymbol\n")
                handle.write("x\ty\n")
            with self.assertRaisesRegex(module.KnowledgeSourceQualificationError,
                                        "columns changed|only "):
                module.verify_chemical_gene(path)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
