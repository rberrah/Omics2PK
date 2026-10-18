"""Synthetic-only original subject naming characters/timepoint separation."""
from __future__ import annotations
import unittest
from scripts.audit_rotroff33_original_sample_naming import (
    mask, structure_from_rows
)

class SourceCodebookPrivacyTests(unittest.TestCase):
    def test_no_original_identifier_embedded_in_patterns(self):
        s="SyntheticPerson23-A"
        self.assertNotIn(s,mask(s))
        self.assertEqual(mask("AB123-C"),"XX###-X")

    def test_triplet_prefixes_are_source_diagnostic_not_validated_join(self):
        rows=[{0:"Sample ID"}]
        for i in range(33):
            for part in "ABC": rows.append({0:f"PERSON{i:02d}-{part}"})
        output=structure_from_rows(rows)
        r=output["suffix_candidates_first_four_columns"][0]
        self.assertEqual(r["n_source_distinct_provisional_ABC_parents"],33)
        self.assertEqual(r["n_parents_with_all_three_A_B_C"],33)
        self.assertTrue(output["never_pair_automatically"])
        self.assertNotIn("PERSON",str(output))

if __name__=="__main__":unittest.main()
