"""Offline synthetic-only pharmacogene PD tests; no public participant/variant data."""
import csv
import io
import unittest

from omics2pk.public_nudt15_variant_pd_benchmark import (
    URN,_variant,NUDT15PublicSourceError,
    parse_original_mavedb,benchmark_nudt15_original,
)


def synthetic_rows():
    output=io.StringIO()
    cols=["accession","hgvs_nt","hgvs_splice","hgvs_pro","score",
          "abundance_score","sensitivity_score"]
    w=csv.DictWriter(output,fieldnames=cols)
    w.writeheader()
    for position in range(1,41):
        for alt in ("Arg","Gly"):
            high=position%2==0
            w.writerow({"accession":f"{URN}#{position}_{alt}",
                        "hgvs_nt":"","hgvs_splice":"",
                        "hgvs_pro":f"p.(Ala{position}{alt})",
                        "score":.8 if high else .2,
                        "abundance_score":.9 if high else .1,
                        "sensitivity_score":.85 if high else .15})
    return output.getvalue().encode()


class NUDT15OriginalSourceTests(unittest.TestCase):
    def test_known_protein_hgvs_3letter_vs_short_and_unsupported_rejected(self):
        self.assertEqual(_variant("p.(Arg139Cys)"),("R",139,"C"))
        self.assertEqual(_variant("p.Arg139Cys"),("R",139,"C"))
        self.assertEqual(_variant("p.R139C"),("R",139,"C"))
        self.assertEqual(_variant("p.(Gln6Glu)"),("Q",6,"E"))
        for invalid in ("p.Arg139Ter","p.Arg139Arg","p.Arg139fs","p.?",
                        "p.(Arg169Cys)","c.42G>A"):
            self.assertIsNone(_variant(invalid))

    def test_synthetic_PDxprotein_dual_assay_position_disjoint(self):
        measured,meta=parse_original_mavedb(synthetic_rows(),strict_original=False)
        self.assertEqual(meta["n_paired_missense_variant_ids"],80)
        self.assertEqual(meta["n_original_nonmissense_or_unparsed_variants"],0)
        r=benchmark_nudt15_original(measured)
        self.assertEqual(r["n_independent_aa_positions"],40)
        self.assertEqual(sum(r["n_independent_aa_positions_per_heldout_fold"]),40)
        self.assertEqual(r["status"],
                         "PUBLIC_NUDT15_VARIANT_MATCHED_CELL_LINE_PD_NOT_HUMAN_PATIENT_TOXICITY")
        self.assertTrue(r["functional_measure_not_same_as_protein_stability"])
        self.assertTrue(r["no_pretreatment_mRNA_or_protein_decay_rate_assayed"])
        self.assertTrue(r["CV_grouped_by_protein_position_not_random_variant"])
        self.assertTrue(r["no_model_uses_PD_score_as_omic_predictor"])
        self.assertNotIn("Ala1Arg",str(r))
        self.assertLess(r["held_position_out_methods"]["measured_protein_stability"]
                        ["heldout_MAE_thiopurine_tolerance_score"],
                        r["held_position_out_methods"]["training_mean"]
                        ["heldout_MAE_thiopurine_tolerance_score"])

    def test_original_row_count_abstains_instead_of_faking_source(self):
        with self.assertRaisesRegex(NUDT15PublicSourceError,"SOURCE_NOT_QUALIFIED|COHORT_NOT_QUALIFIED"):
            parse_original_mavedb(synthetic_rows(),strict_original=True)

    def test_duplicate_variant_source_or_damaged_csv_rejected(self):
        data=synthetic_rows()
        fragment=data.strip().splitlines()[-1]
        with self.assertRaisesRegex(NUDT15PublicSourceError,"duplicated"):
            parse_original_mavedb(data+fragment+b"\n",strict_original=False)
        with self.assertRaisesRegex(NUDT15PublicSourceError,"columns"):
            parse_original_mavedb(data.replace(b"abundance_score",b"bad_column_name"),
                                  strict_original=False)
        with self.assertRaises(NUDT15PublicSourceError):
            parse_original_mavedb(data.replace(b"0.85",b"not_a_number"),
                                  strict_original=False)

    def test_source_stability_and_function_are_distinct_phenotypes(self):
        data=synthetic_rows().replace(b",0.9,0.85",b",0.9,0.15",20)
        a,m=parse_original_mavedb(data,strict_original=False)
        r=benchmark_nudt15_original(a)
        self.assertEqual(r["n_stable_but_thiopurine_sensitive_research_thresholds_0_8_0_3"],20)
        self.assertTrue(r["not_direct_NUDT15_enzymatic_substrate_metabolite_measure"])


if __name__=="__main__":
    unittest.main()
