"""Synthetic-only tests for donor/assay/CLint qualification, no public donor rows."""
import unittest

from omics2pk.clint_contract import (
    CLintContractError, leave_one_donor_out, qualify_pairing
)

HASH = "a" * 64
DRUGS = ("bupropion", "diclofenac", "omeprazole", "bufuralol", "midazolam")


def fixture(n=15):
    clint, protein = [], []
    for donor in range(n):
        donor_id = "SYNTH_DONOR_%02d" % donor
        for drug in DRUGS:
            for system in ("HH", "HLM"):
                clint.append(dict(
                    donor_id=donor_id, drug=drug, system=system,
                    clint=0.1 + donor * 0.02, clint_unit=(
                        "uL/min/million_cells" if system == "HH" else "uL/min/mg_protein"
                    ), clint_basis="UNBOUND", measurement_status="MEASURED",
                    source_table="synthetic_fixture", source_sha256=HASH
                ))
        for system in ("HH", "HLM"):
            for gene in ("CYP2B6", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4"):
                protein.append(dict(
                    donor_id=donor_id, system=system, protein_id=gene,
                    abundance=0.3 + donor * 0.01, abundance_unit="fmol/microgram_protein",
                    source_table="synthetic_fixture", source_sha256=HASH
                ))
    return clint, protein


class PairedMeasuredCLintContractTests(unittest.TestCase):
    def test_donor_held_out_all_drugs_and_both_systems(self):
        clint, prot = fixture()
        result = qualify_pairing(clint, prot)
        self.assertEqual(result["n_independent_donors"], 15)
        self.assertEqual(result["n_measured_endpoint_rows"], 75)
        self.assertEqual(result["n_target_drugs"], 5)
        self.assertNotIn("SYNTH_DONOR", repr(result))
        cv = leave_one_donor_out(clint)
        self.assertEqual(len(cv), 15)
        for train, test in cv:
            self.assertEqual(len(test), 10)
            self.assertEqual(len(train), 140)
            train_ids = {clint[i]["donor_id"] for i in train}
            test_ids = {clint[i]["donor_id"] for i in test}
            self.assertTrue(train_ids.isdisjoint(test_ids))

    def test_missing_protein_donor_must_not_be_silently_dropped(self):
        clint, prot = fixture()
        prot = [r for r in prot if r["donor_id"] != "SYNTH_DONOR_00"]
        with self.assertRaisesRegex(CLintContractError, "Donor sets do not match"):
            qualify_pairing(clint, prot)

    def test_no_mixed_units_within_same_system_and_drug(self):
        clint, prot = fixture()
        clint[-1]["clint_unit"] = "mL/min/kg"
        with self.assertRaisesRegex(CLintContractError, "Inconsistent CLint unit"):
            qualify_pairing(clint, prot)

    def test_no_pseudoreplicated_protein_rows(self):
        clint, prot = fixture()
        prot.append(dict(prot[0]))
        with self.assertRaisesRegex(CLintContractError, "technical replicates"):
            qualify_pairing(clint, prot)

    def test_no_pseudoreplicated_clearance_assays(self):
        clint, prot = fixture()
        clint.append(dict(clint[0]))
        with self.assertRaisesRegex(CLintContractError, "technical replicates"):
            qualify_pairing(clint, prot)

    def test_reject_nonmeasured_imputed_or_nonpositive_clearance(self):
        for key, value, expected in (
            ("measurement_status", "IMPUTED", "Inferred"),
            ("clint", 0, "strictly positive"),
            ("clint", float("nan"), "strictly positive"),
            ("clint_basis", "unknown", "must be"),
        ):
            with self.subTest(key=key, value=value):
                clint, prot = fixture()
                clint[0][key] = value
                with self.assertRaisesRegex(CLintContractError, expected):
                    qualify_pairing(clint, prot)

    def test_separate_feature_system_is_allowed_only_if_same_donors(self):
        clint, prot = fixture()
        result = qualify_pairing(clint, prot, target_system="HH", feature_system="HLM")
        self.assertEqual(result["feature_system"], "HLM")
        self.assertEqual(result["target_system"], "HH")

    def test_allows_explicit_missing_drug_measurements_without_inventing_values(self):
        clint, prot = fixture()
        clint = [r for r in clint if not (r["donor_id"] == "SYNTH_DONOR_01"
                                         and r["drug"] == "midazolam" and r["system"] == "HH")]
        result = qualify_pairing(clint, prot)
        self.assertEqual(result["endpoint_counts_per_drug"]["midazolam"], 14)
        self.assertEqual(result["n_independent_donors"], 15)


if __name__ == "__main__":
    unittest.main()
