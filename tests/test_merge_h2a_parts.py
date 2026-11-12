"""Contract tests for the H2a batch merge.

Running ten endpoints in four parallel batches is a scheduling decision and must
leave no trace in the result. These tests check that the merge refuses anything
that would make it leave one, and that the void wrong-mechanism control is
declared rather than left to be spotted in a table of near-identical numbers.
"""
from __future__ import annotations

import copy
import importlib.util
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

try:
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"pytest is required for this module: {exc}") from exc

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "merge_h2a_parts.py"


def load_merger():
    spec = importlib.util.spec_from_file_location("merge_h2a", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MERGE = load_merger()


def endpoint_block(name, target="GENEA", wrong="GENEB"):
    return {
        "endpoint": name,
        "fold_design": "5_fold_x2_repeats",
        "repeats": 2,
        "n_donors": 397,
        "wrong_seed_label": "NEXT_activity",
        "first_fold_selection": {
            "netpk_target_seeded": {"1": [target]},
            "netpk_wrong_mechanism": {"1": [wrong]},
        },
    }


def part(endpoints, **overrides):
    body = {
        "protocol_version": "1.4.3",
        "endpoint_class": "functional_in_vitro",
        "frozen": {"seed": 20261108},
        "network": {"nodes": 16201},
        "universe": {"candidate_genes": 31862},
        "run_on": "2026-11-08",
        "endpoints": {name: endpoint_block(name) for name in endpoints},
    }
    body.update(overrides)
    return body


def write(root, name, body):
    path = Path(root) / name
    path.write_text(json.dumps(body), encoding="utf-8")
    return str(path)


ALL = MERGE.EXPECTED_ENDPOINTS


class TestTheMergeRefusesInconsistency(unittest.TestCase):
    def test_a_complete_set_merges(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", part(ALL[:5]))
            second = write(tmp, "b.json", part(ALL[5:]))
            merged = MERGE.merge([first, second])
        self.assertEqual(merged["endpoints_merged"], 10)
        self.assertEqual(merged["endpoints_missing"], [])
        self.assertEqual(len(merged["merged_from"]), 2)

    def test_an_incomplete_set_is_refused_by_default(self):
        with TemporaryDirectory() as tmp:
            only = write(tmp, "a.json", part(ALL[:5]))
            with self.assertRaisesRegex(MERGE.MergeError, "incomplete"):
                MERGE.merge([only])

    def test_an_incomplete_set_can_be_merged_deliberately_and_says_so(self):
        with TemporaryDirectory() as tmp:
            only = write(tmp, "a.json", part(ALL[:5]))
            merged = MERGE.merge([only], require_all=False)
        self.assertEqual(merged["endpoints_merged"], 5)
        self.assertEqual(len(merged["endpoints_missing"]), 5)

    def test_a_duplicated_endpoint_is_refused(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", part(ALL[:6]))
            second = write(tmp, "b.json", part(ALL[5:]))
            with self.assertRaisesRegex(MERGE.MergeError, "more than one batch"):
                MERGE.merge([first, second])

    def test_a_different_protocol_version_is_refused(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", part(ALL[:5]))
            second = write(tmp, "b.json", part(ALL[5:], protocol_version="1.4.2"))
            with self.assertRaisesRegex(MERGE.MergeError, "protocol_version"):
                MERGE.merge([first, second])

    def test_different_frozen_parameters_are_refused(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", part(ALL[:5]))
            second = write(tmp, "b.json", part(ALL[5:], frozen={"seed": 1}))
            with self.assertRaisesRegex(MERGE.MergeError, "frozen"):
                MERGE.merge([first, second])

    def test_a_different_network_or_universe_is_refused(self):
        for field, value in (("network", {"nodes": 999}),
                             ("universe", {"candidate_genes": 1})):
            with TemporaryDirectory() as tmp:
                first = write(tmp, "a.json", part(ALL[:5]))
                second = write(tmp, "b.json", part(ALL[5:], **{field: value}))
                with self.assertRaisesRegex(MERGE.MergeError, field):
                    MERGE.merge([first, second])

    def test_a_different_fold_design_is_refused(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", part(ALL[:5]))
            body = part(ALL[5:])
            for block in body["endpoints"].values():
                block["fold_design"] = "leave_one_centre_out"
            second = write(tmp, "b.json", body)
            with self.assertRaisesRegex(MERGE.MergeError, "fold designs"):
                MERGE.merge([first, second])

    def test_a_different_repeat_count_is_refused(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", part(ALL[:5]))
            body = part(ALL[5:])
            for block in body["endpoints"].values():
                block["repeats"] = 1
            second = write(tmp, "b.json", body)
            with self.assertRaisesRegex(MERGE.MergeError, "repeat counts"):
                MERGE.merge([first, second])

    def test_an_unexpected_endpoint_is_refused(self):
        with TemporaryDirectory() as tmp:
            body = part(ALL)
            body["endpoints"]["CYP9Z9_activity"] = endpoint_block("CYP9Z9_activity")
            only = write(tmp, "a.json", body)
            with self.assertRaisesRegex(MERGE.MergeError, "unexpected endpoints"):
                MERGE.merge([only])

    def test_no_batch_files_is_an_error(self):
        with self.assertRaisesRegex(MERGE.MergeError, "no batch files"):
            MERGE.merge([])

    def test_the_merge_records_how_it_was_assembled(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", part(ALL[:5]))
            second = write(tmp, "b.json", part(ALL[5:]))
            merged = MERGE.merge([first, second])
        self.assertIn("scheduling decision", MERGE.__doc__)
        self.assertIn("no endpoint appears twice", merged["merge_note"])
        self.assertEqual({entry["file"] for entry in merged["merged_from"]},
                         {"a.json", "b.json"})


class TestTheVoidControlIsDeclared(unittest.TestCase):
    """The frozen rotation can seed the wrong-mechanism arm with the right gene.

    Where a cohort measures one enzyme with two substrates, the next endpoint's
    label resolves to the same gene and the control compares an arm against
    itself. The numbers then look like "a wrong mechanism costs nothing", which
    is the opposite of what the control can support.
    """

    def test_a_different_gene_is_a_valid_control(self):
        audit = MERGE.audit_wrong_mechanism_control(
            {"A_activity": endpoint_block("A_activity", "GENEA", "GENEB")})
        self.assertEqual(audit["valid_count"], 1)
        self.assertEqual(audit["void_count"], 0)

    def test_the_same_gene_is_a_void_control_and_is_named(self):
        audit = MERGE.audit_wrong_mechanism_control(
            {"A_activity": endpoint_block("A_activity", "GENEA", "GENEA")})
        self.assertEqual(audit["void_count"], 1)
        self.assertEqual(audit["void_same_gene"][0]["endpoint"], "A_activity")
        self.assertEqual(audit["void_same_gene"][0]["gene"], "GENEA")

    def test_a_missing_selection_is_undetermined_not_assumed_valid(self):
        block = endpoint_block("A_activity")
        del block["first_fold_selection"]["netpk_wrong_mechanism"]
        audit = MERGE.audit_wrong_mechanism_control({"A_activity": block})
        self.assertEqual(audit["undetermined"], ["A_activity"])
        self.assertEqual(audit["valid_count"], 0)
        self.assertEqual(audit["void_count"], 0)

    def test_the_audit_refuses_to_let_a_void_control_be_read_as_a_result(self):
        audit = MERGE.audit_wrong_mechanism_control(
            {"A_activity": endpoint_block("A_activity", "GENEA", "GENEA")})
        text = audit["how_to_read"]
        self.assertIn("says nothing", text)
        self.assertIn("must not be read as evidence", text)
        self.assertIn("not re-run", text)

    def test_the_merge_carries_the_audit(self):
        with TemporaryDirectory() as tmp:
            body = part(ALL)
            body["endpoints"][ALL[0]] = endpoint_block(ALL[0], "GENEA", "GENEA")
            only = write(tmp, "a.json", body)
            merged = MERGE.merge([only])
        self.assertEqual(merged["wrong_mechanism_control_audit"]["void_count"], 1)
        self.assertEqual(merged["wrong_mechanism_control_audit"]["valid_count"], 9)


def h2d_part(endpoints, **overrides):
    """A batch shaped like an H2d report: no protocol version, no network."""
    body = {
        "analysis": "H2d: hybrid mechanistic selection",
        "endpoint_class": "functional_in_vitro",
        "generating_endpoint": "CYP2D6_activity",
        "validation_set_independence": {"dependence_index_not_a_replication_count": 3.55},
        "frozen": {"hybrid": {"oversample": 3}},
        "universe": {"candidate_genes": 22237,
                     "qualified_pool": {"pool_name": "A_qualified_human"}},
        "run_on": "2026-11-09",
        "endpoints": {name: {"endpoint": name,
                             "fold_design": "5_fold_x2_repeats",
                             "repeats": 2} for name in endpoints},
    }
    body.update(overrides)
    return body


class TestTheGateDoesNotDependOnWhichAnalysisRan(unittest.TestCase):
    """The named-field list could not see a field a report had just gained.

    The qualified pool lives inside `universe`, so a merge across two different
    pools would have passed a list that checked `universe` by name only if the
    rest of `universe` happened to agree. Every top-level field except the
    endpoint results and the run date is now compared.
    """

    def test_batches_of_another_analysis_merge_under_their_own_expected_set(self):
        nine = [e for e in ALL if e != "CYP2D6_activity"]
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", h2d_part(nine[:5]))
            second = write(tmp, "b.json", h2d_part(nine[5:]))
            merged = MERGE.merge([first, second], expected=tuple(nine))
        self.assertEqual(merged["endpoints_merged"], 9)
        self.assertEqual(merged["endpoints_missing"], [])

    def test_the_generating_endpoint_is_refused_when_it_is_not_expected(self):
        nine = [e for e in ALL if e != "CYP2D6_activity"]
        with TemporaryDirectory() as tmp:
            only = write(tmp, "a.json", h2d_part(nine + ["CYP2D6_activity"]))
            with self.assertRaisesRegex(MERGE.MergeError, "unexpected endpoints"):
                MERGE.merge([only], expected=tuple(nine))

    def test_two_different_qualified_pools_are_refused(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", h2d_part(["CYP2A6_activity"]))
            other = h2d_part(["CYP2B6_activity"])
            other["universe"]["qualified_pool"]["pool_name"] = "B_qualified_human_on_network"
            second = write(tmp, "b.json", other)
            with self.assertRaisesRegex(MERGE.MergeError, "universe"):
                MERGE.merge([first, second],
                            expected=("CYP2A6_activity", "CYP2B6_activity"))

    def test_a_field_absent_from_the_reference_is_still_compared(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", h2d_part(["CYP2A6_activity"]))
            other = h2d_part(["CYP2B6_activity"])
            other["an_extra_declaration"] = "only in this batch"
            second = write(tmp, "b.json", other)
            with self.assertRaisesRegex(MERGE.MergeError, "an_extra_declaration"):
                MERGE.merge([first, second],
                            expected=("CYP2A6_activity", "CYP2B6_activity"))

    def test_a_batch_without_the_frozen_block_is_refused(self):
        with TemporaryDirectory() as tmp:
            body = h2d_part(["CYP2A6_activity"])
            del body["frozen"]
            only = write(tmp, "a.json", body)
            with self.assertRaisesRegex(MERGE.MergeError, "declares no"):
                MERGE.merge([only], expected=("CYP2A6_activity",))

    def test_a_run_crossing_midnight_is_not_an_inconsistency(self):
        with TemporaryDirectory() as tmp:
            first = write(tmp, "a.json", h2d_part(["CYP2A6_activity"]))
            second = write(tmp, "b.json",
                           h2d_part(["CYP2B6_activity"], run_on="2026-11-10"))
            merged = MERGE.merge([first, second],
                                 expected=("CYP2A6_activity", "CYP2B6_activity"))
        self.assertEqual({entry["run_on"] for entry in merged["merged_from"]},
                         {"2026-11-09", "2026-11-10"})

    def test_an_unauditable_rotation_is_declared_absent_not_passed(self):
        with TemporaryDirectory() as tmp:
            only = write(tmp, "a.json", h2d_part(["CYP2A6_activity"]))
            merged = MERGE.merge([only], expected=("CYP2A6_activity",))
        audit = merged["wrong_mechanism_control_audit"]
        self.assertIn("not_applicable", audit)
        self.assertNotIn("valid_count", audit)
        self.assertIn("absent rather than reported as passed",
                      audit["not_applicable"])


class TestAgainstTheRealBatches(unittest.TestCase):
    """Runs only on the batch files that exist, so it cannot block a fresh clone."""

    def parts(self):
        return sorted((REPO_ROOT / "reports" / "data").glob("h2a_part_*.json"))

    def test_any_committed_batch_declares_the_frozen_specification(self):
        found = self.parts()
        if not found:
            self.skipTest("no batch files in the working tree")
        for path in found:
            body = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(body["endpoint_class"], "functional_in_vitro")
            self.assertEqual(body["frozen"]["seed"], 20261108)
            self.assertEqual(body["frozen"]["outer_repeats"], 2)
            self.assertGreaterEqual(body["frozen"]["random_replications"], 100)
            for result in body["endpoints"].values():
                self.assertEqual(result["repeats"], 2)
                self.assertEqual(result["folds"], 10)

    def test_committed_batches_merge_without_contradiction(self):
        found = self.parts()
        if not found:
            self.skipTest("no batch files in the working tree")
        merged = MERGE.merge([str(p) for p in found], require_all=False)
        self.assertGreater(merged["endpoints_merged"], 0)
        self.assertIn("wrong_mechanism_control_audit", merged)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
