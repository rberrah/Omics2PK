"""GEO 427-array original data streaming gate, fabricated matrices ONLY."""
from __future__ import annotations
from io import BytesIO
import unittest

from scripts.qualify_gse9588_real_human_liver_expression import (
    audit_gse_matrix_stream,N_EXPECTED
)


def fake_gse(n=427,nrows=1001):
    handle=BytesIO()
    handle.write(b"!Series_platform_id\t\"GPL4371\"\n")
    handle.write(b"!series_matrix_table_begin\n")
    handle.write(b"ID_REF\t"+b"\t".join(
        f"GSM{100000+i}".encode() for i in range(n))+b"\n")
    line=b"FICTITIOUS_PROBE\t"+b"\t".join(b"2.3" for _ in range(n))+b"\n"
    for _ in range(nrows):handle.write(line)
    handle.write(b"!series_matrix_table_end\n")
    handle.seek(0)
    return handle


class FabricatedGEOInputSourceTests(unittest.TestCase):
    def test_427_original_source_style_arrays_but_no_clinical_outcomes(self):
        result=audit_gse_matrix_stream(fake_gse())
        self.assertEqual(result["original_unique_GSM_sample_columns"],427)
        self.assertEqual(result["n_original_numeric_probe_rows"],1001)
        self.assertTrue(result["RNA_only_gene_expression_matrix_structurally_qualified"])
        self.assertFalse(result["original_398_CYP_activity_same_person_crosswalk_verified"])
        self.assertFalse(result["clinical_PK_validation_executed"])
        self.assertNotIn("GSM1000",str(result))

    def test_426_original_array_labels_denied(self):
        with self.assertRaisesRegex(ValueError,"427"):
            audit_gse_matrix_stream(fake_gse(n=426))

    def test_non_numeric_probe_rows_not_source_validated(self):
        h=fake_gse(nrows=999)
        with self.assertRaisesRegex(ValueError,"insufficient"):
            audit_gse_matrix_stream(h)

    def test_repeated_original_donor_GSM_denied(self):
        b=fake_gse().getvalue()
        changed=b.replace(b"GSM100001",b"GSM100000")
        with self.assertRaisesRegex(ValueError,"unique"):
            audit_gse_matrix_stream(BytesIO(changed))


if __name__=="__main__":unittest.main()
