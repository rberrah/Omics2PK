"""Synthetic-only tests. No third-party activity/protein rows checked in."""
import base64
import csv
import io
import json
import unittest

from omics2pk.public_cyp2c9_variant_function_benchmark import (
    PublicVariantDMSUnavailable,parse_public_original,
    benchmark_public_original,SOURCE_BLOB_SHA,
)
from scripts.benchmark_public_cyp2c9_variant_activity_abundance import (
    qualify_original_public_blob,
)


def fixture():
    fields=("variant","class","position","start","end","activity_class",
            "activity_score","abundance_class","abundance_score")
    stream=io.StringIO()
    writer=csv.DictWriter(stream,fieldnames=fields)
    writer.writeheader()
    for pos in range(1,41):
        for b in ("C","G"):
            high=pos%2==0
            writer.writerow({
                "variant":f"A{pos}{b}","class":"missense","position":pos,
                "start":"A","end":b,
                "activity_class":"wt-like" if high else "decreased",
                "activity_score":str((.80 if high else .25)+.01*(pos%4)),
                "abundance_class":"wt-like" if high else "decreased",
                "abundance_score":str((.90 if high else .35)+.01*(pos%5)),
            })
    return stream.getvalue().encode()


class PublicKnownVariantTests(unittest.TestCase):
    def test_synthetic_residue_disjoint_CV_discriminates_protein_function(self):
        rows=parse_public_original(fixture(),strict_original=False)
        self.assertEqual(len(rows),80)
        observed=benchmark_public_original(rows)
        self.assertEqual(observed["n_individual_cyp2c9_missense_variants_with_both_scores"],80)
        self.assertEqual(observed["n_independent_protein_positions_grouped_for_holdout"],40)
        self.assertEqual(sum(observed["heldout_unique_positions_per_fold"]),40)
        self.assertEqual(observed["n_original_public_missense_decreased_activity_but_wt_like_abundance"],0)
        self.assertTrue(observed["no_specific_drug_metabolism_or_substrate_activity_measured"])
        self.assertFalse(observed["original_source_variant_and_activity_rows_released"])
        self.assertEqual(len(observed["heldout_metrics"]),4)
        self.assertLess(observed["heldout_metrics"]["protein_abundance_only"]
                        ["held_position_out_MAE_activity_score"],
                        observed["heldout_metrics"]["fold_training_mean"]
                        ["held_position_out_MAE_activity_score"])
        self.assertGreater(observed["heldout_metrics"]["protein_abundance_only"]
                           ["held_position_out_binary_activity_class_AUROC_not_PK_AUC"],.8)
        self.assertNotIn("A1C",str(observed))
        self.assertEqual(benchmark_public_original(rows),observed)

    def test_wrong_original_row_count_cannot_claim_real_replay(self):
        with self.assertRaisesRegex(PublicVariantDMSUnavailable,"source qualified"):
            parse_public_original(fixture(),strict_original=True)

    def test_duplicate_and_genomic_annotation_disagreement_fails_closed(self):
        body=fixture()
        last=body.strip().splitlines()[-1]
        duplicate=body+last+b"\n"
        with self.assertRaisesRegex(PublicVariantDMSUnavailable,"Duplicate"):
            parse_public_original(duplicate,strict_original=False)
        bad=body.replace(b"A1C",b"Q1C")
        with self.assertRaisesRegex(PublicVariantDMSUnavailable,"identity"):
            parse_public_original(bad,strict_original=False)
        nonsense=body.replace(b"missense",b"unknown_")
        with self.assertRaisesRegex(PublicVariantDMSUnavailable,"annotation"):
            parse_public_original(nonsense,strict_original=False)

    def test_source_blob_cryptographic_gate_not_just_server_label(self):
        raw=fixture()
        response=json.dumps({
            "sha":SOURCE_BLOB_SHA,"encoding":"base64",
            "size":len(raw),"content":base64.b64encode(raw).decode(),
        }).encode()
        with self.assertRaisesRegex(PublicVariantDMSUnavailable,"cryptographic"):
            qualify_original_public_blob(response)
        wrong=json.dumps({"sha":"wrong","content":"foo",
                          "encoding":"base64","size":0}).encode()
        with self.assertRaisesRegex(PublicVariantDMSUnavailable,"identity"):
            qualify_original_public_blob(wrong)

    def test_missing_independent_positions_abstains(self):
        rows=parse_public_original(fixture(),strict_original=False)
        with self.assertRaisesRegex(PublicVariantDMSUnavailable,"positions"):
            benchmark_public_original(rows[:20])


if __name__=="__main__":
    unittest.main()
