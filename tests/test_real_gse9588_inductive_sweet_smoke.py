"""All data fabricated, no original GEO samples or clinical donor labels."""
import importlib.util
import math
import unittest
from io import BytesIO

from scripts.smoke_real_gse9588_inductive_sweet import (
    TRAIN,PANEL,extract_train_variable_array_probes,real_gse_frozen_network_smoke
)

N=427

def toy(nrows=1100):
    payload=BytesIO()
    payload.write(b"!series_matrix_table_begin\n")
    payload.write(b"ID_REF\t"+b"\t".join(
        f"GSM{700000+i}".encode() for i in range(N))+b"\n")
    for j in range(nrows):
        nums=[str(round(3.0+(j%31)*0.1+0.9*math.sin((j+1)*(i+2)/47),6)).encode()
              for i in range(N)]
        payload.write(b"FAKE_PROBE\t"+b"\t".join(nums)+b"\n")
    payload.write(b"!series_matrix_table_end\n")
    payload.seek(0)
    return payload

class SourceFrozenReferenceTests(unittest.TestCase):
    def test_toy_427_sample_extracted_from_source_all_available(self):
        names,probes,n=extract_train_variable_array_probes(toy())
        self.assertEqual(len(names),N)
        self.assertEqual(len(probes),PANEL)
        self.assertEqual(n,1100)
        self.assertNotIn("GSM7000",repr([PANEL,n]))

    @unittest.skipUnless(importlib.util.find_spec("numpy"),"Optional numpy not installed")
    def test_toy_427_original_style_frozen_reference_does_not_claim_PK(self):
        ids,probes,n=extract_train_variable_array_probes(toy())
        result=real_gse_frozen_network_smoke(ids,probes,n)
        self.assertEqual(result["fixed_reference_training_GSM_samples"],TRAIN)
        self.assertEqual(result["unseen_GSM_samples_projected"],352)
        self.assertFalse(result["predictive_CYP_activity_test_executed"])
        self.assertFalse(result["equivalent_to_official_SWEET"])
        self.assertFalse(result["gene_symbol_annotation_provenance"])
        self.assertFalse(result["same_original_human_liver_donors_proven_independent"])
        self.assertNotIn("GSM700",str(result))

if __name__=="__main__":unittest.main()
