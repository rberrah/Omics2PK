"""Synthetic-only tests for the direct HLC Synapse same-donor source gate."""
import random
import unittest

from omics2pk.hlc_synapse_same_donor_gate import (
    CURRENT_EXPRESSION_ENTITY_SAMPLES,
    HLCSynapseJoinError,
    join_original_hlc_synapse,
    summarize_private_hlc_join,
)


def fixture(n_matched=90):
    expression = tuple(f"HLC_EXPR_{i:03}" for i in range(CURRENT_EXPRESSION_ENTITY_SAMPLES))
    rows = []
    for i in range(n_matched):
        rows.append({
            "donor_id": expression[i],
            "phenotype_entity_id": "syn4629",
            "assay_id": "CYP3A4_midazolam",
            "activity_value": str(0.5 + i / 20),
            "activity_unit": "nmol/min/mg protein",
            "activity_status": "observed",
            "site_id": ("VANDERBILT" if i % 3 == 0 else ("PITTSBURGH" if i % 3 == 1 else "MERCK")),
            "age": str(20 + (i % 60)),
            "sex": "male" if i % 2 else "female",
        })
    return expression, rows


class HLCSynapseSameDonorGateTests(unittest.TestCase):
    def test_explicit_same_donor_join_is_ready_and_summary_is_aggregate_only(self):
        expression, rows = fixture()
        join = join_original_hlc_synapse(
            expression, rows,
            assay_id="CYP3A4_midazolam",
            expression_sha256="a" * 64,
            phenotype_sha256="b" * 64,
        )
        self.assertEqual(len(join.donor_ids), 90)
        self.assertEqual(len(join.expression_indices), 90)
        self.assertEqual(join.activity_unit, "nmol/min/mg protein")
        out = summarize_private_hlc_join(join)
        self.assertEqual(out["n_same_donor_complete_observations"], 90)
        self.assertEqual(out["n_original_sites"], 3)
        self.assertTrue(out["ready_for_site_disjoint_E2_benchmark"])
        self.assertFalse(out["GSE9588_crosswalk_required"])
        self.assertFalse(out["individual_donor_ids_or_activity_values_exported"])
        self.assertNotIn("HLC_EXPR_", str(out))
        self.assertNotIn("0.5", str(out))

    def test_row_order_is_irrelevant_because_join_uses_explicit_donor_key(self):
        expression, rows = fixture()
        shuffled = list(rows)
        random.Random(20261022).shuffle(shuffled)
        a = join_original_hlc_synapse(
            expression, rows, assay_id="CYP3A4_midazolam",
            expression_sha256="a" * 64, phenotype_sha256="b" * 64)
        b = join_original_hlc_synapse(
            expression, shuffled, assay_id="CYP3A4_midazolam",
            expression_sha256="a" * 64, phenotype_sha256="b" * 64)
        self.assertEqual(a.expression_indices, b.expression_indices)
        self.assertEqual(a.donor_ids, b.donor_ids)
        self.assertEqual(a.activities, b.activities)

    def test_nonobserved_activity_is_not_coerced_to_zero(self):
        expression, rows = fixture()
        rows[0]["activity_status"] = "blq"
        rows[0]["activity_value"] = ""
        join = join_original_hlc_synapse(
            expression, rows, assay_id="CYP3A4_midazolam",
            expression_sha256="a" * 64, phenotype_sha256="b" * 64)
        self.assertEqual(len(join.donor_ids), 89)
        rows[1]["activity_status"] = "missing"
        rows[1]["activity_value"] = "0"
        with self.assertRaisesRegex(HLCSynapseJoinError, "must not carry"):
            join_original_hlc_synapse(
                expression, rows, assay_id="CYP3A4_midazolam",
                expression_sha256="a" * 64, phenotype_sha256="b" * 64)

    def test_wrong_source_entity_duplicate_or_bad_sha_fails_closed(self):
        expression, rows = fixture()
        bad = [dict(x) for x in rows]
        bad[0]["phenotype_entity_id"] = "syn999"
        with self.assertRaisesRegex(HLCSynapseJoinError, "syn4629"):
            join_original_hlc_synapse(
                expression, bad, assay_id="CYP3A4_midazolam",
                expression_sha256="a" * 64, phenotype_sha256="b" * 64)
        dup = rows + [dict(rows[0])]
        with self.assertRaisesRegex(HLCSynapseJoinError, "Duplicate"):
            join_original_hlc_synapse(
                expression, dup, assay_id="CYP3A4_midazolam",
                expression_sha256="a" * 64, phenotype_sha256="b" * 64)
        with self.assertRaisesRegex(HLCSynapseJoinError, "SHA-256"):
            join_original_hlc_synapse(
                expression, rows, assay_id="CYP3A4_midazolam",
                expression_sha256="not-a-hash", phenotype_sha256="b" * 64)

    def test_requires_current_expression_entity_cardinality_and_three_sites(self):
        expression, rows = fixture()
        with self.assertRaisesRegex(HLCSynapseJoinError, "exactly 465"):
            join_original_hlc_synapse(
                expression[:-1], rows, assay_id="CYP3A4_midazolam",
                expression_sha256="a" * 64, phenotype_sha256="b" * 64)
        two_site = [dict(x) for x in rows]
        for i, row in enumerate(two_site):
            row["site_id"] = "A" if i % 2 else "B"
        with self.assertRaisesRegex(HLCSynapseJoinError, "at least three"):
            join_original_hlc_synapse(
                expression, two_site, assay_id="CYP3A4_midazolam",
                expression_sha256="a" * 64, phenotype_sha256="b" * 64)

    def test_too_few_same_donor_observations_abstains(self):
        expression, rows = fixture(n_matched=59)
        with self.assertRaisesRegex(HLCSynapseJoinError, "Fewer than 60"):
            join_original_hlc_synapse(
                expression, rows, assay_id="CYP3A4_midazolam",
                expression_sha256="a" * 64, phenotype_sha256="b" * 64)


if __name__ == "__main__":
    unittest.main()
