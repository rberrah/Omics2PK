"""Local Synapse human liver source safety tests; strictly fictional files."""
from __future__ import annotations
import contextlib
import io
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from zipfile import ZipFile

from omics2pk.hlc_authorized_local_source import (
    LocalHLCSourceError,inspect_original_local_source,audit_local_hlc_sources
)


def make_mock_phenotype(path):
    path.write_text("donor_id,CYP3A4_activity,CYP3A4_protein\n"
                    "FICTIONAL_HUMAN_1,3.3,6.0\n"
                    "FICTIONAL_HUMAN_2,1.2,4.1\n",encoding="utf8")


class OriginalHLCPrivacyTests(unittest.TestCase):
    def test_tabular_original_source_scan_never_outputs_ids_or_values(self):
        with TemporaryDirectory() as t:
            path=Path(t)/"phenotype.csv"
            make_mock_phenotype(path)
            report=inspect_original_local_source(path,role="phenotype")
            self.assertEqual(report["original_Synapse_id"],"syn4629")
            self.assertEqual(report["structure"]["n_rows"],2)
            self.assertEqual(report["structure"]["n_columns"],3)
            self.assertTrue(report["structure"]["source_header_semantics_present"]["cyp"])
            self.assertFalse(report["same_donor_GEO_crosswalk_qualified"])
            self.assertFalse(report["CYP_model_fitted"])
            self.assertNotIn("FICTIONAL_HUMAN",str(report))
            self.assertNotIn("3.3",str(report))

    def test_synapse_phenotype_name_not_proof_assay_is_present(self):
        with TemporaryDirectory() as t:
            path=Path(t)/"curatedPhenotype"
            path.write_bytes(b"\x89HDF\r\nbinaryunknown")
            report=inspect_original_local_source(path,role="phenotype")
            self.assertEqual(report["format"],"binary")
            self.assertTrue(report["structure"]["requires_authorized_source_format_conversion"])
            self.assertFalse(report["original_numeric_CYP_activity_qualified"])

    def test_actual_archive_members_never_logged(self):
        with TemporaryDirectory() as t:
            path=Path(t)/"private.zip"
            with ZipFile(path,"w") as z:
                z.writestr("PERSON_RACYSAMPLE1234.csv",
                           "id,CYP3A4\nFICTIONAL_1,4\n")
            report=inspect_original_local_source(path,role="expression")
            self.assertEqual(report["structure"]["archive_members"],1)
            self.assertNotIn("PERSON_RACYSAMPLE",str(report))

    def test_no_source_uid_substitution_across_unmatched_cohorts(self):
        with TemporaryDirectory() as t:
            path=Path(t)/"curated.csv"
            make_mock_phenotype(path)
            report=audit_local_hlc_sources({"phenotype":path})
            self.assertFalse(report["GSE9588_GSM_to_original_Synapse_donor_ID_matched"])
            self.assertFalse(report["original_CYP_activity_per_donor_and_unit_verified"])
            self.assertFalse(report["training_or_model_PK_validation_executed"])

    def test_unexpected_source_role_rejected(self):
        with TemporaryDirectory() as t:
            file=Path(t)/"fake.csv"
            make_mock_phenotype(file)
            with self.assertRaises(LocalHLCSourceError):
                inspect_original_local_source(file,role="unrelated_study")

    def test_private_download_cli_uses_only_requested_synapse_ids(self):
        from scripts.download_authorized_synapse_hlc import main
        with TemporaryDirectory() as t:
            private=Path(t)/"private"/"hlc"
            private.mkdir(parents=True)
            phenotype=private/"source_data.csv"
            make_mock_phenotype(phenotype)
            class FakeSynapse:
                def get(self,synid,downloadLocation):
                    if synid=="syn4629":
                        return SimpleNamespace(path=str(phenotype))
                    raise PermissionError("FICTIONAL_SECRET_TOKEN_DO_NOT_LOG")
            fake=SimpleNamespace(login=lambda silent:FakeSynapse())
            output=io.StringIO()
            with patch.dict(sys.modules,{"synapseclient":fake}):
                with contextlib.redirect_stdout(output):
                    code=main(["--out",str(private)])
            self.assertEqual(code,2)
            self.assertIn("syn4629",output.getvalue())
            self.assertIn("syn88644",output.getvalue())
            self.assertNotIn("FICTIONAL_SECRET_TOKEN",output.getvalue())
            self.assertNotIn("FICTIONAL_HUMAN",output.getvalue())
            self.assertNotIn(str(private),output.getvalue())


if __name__=="__main__":
    unittest.main()
