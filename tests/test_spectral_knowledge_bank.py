"""Fully synthetic spectra + real source identifiers, no fabricated biomedical identification."""
import sqlite3
import unittest

from omics2pk.spectral_biochemical_seed import seed_bibliography
from omics2pk.spectral_knowledge_bank import (
    REFERENCE_CATALOGUE, SpectralKnowledgeError, add_mechanism_link,
    add_reference, make_db, search_spectral_candidates, triage_unassigned_feature_id,
    register_original_omic_feature,summarize_registered_features
)

PEAKS=[[50.0,4.0],[57.0,9.0],[68.0,3.0],[73.0,10.0],
       [85.0,7.0],[103.0,1.0],[115.0,5.0],[128.0,3.0]]
def gc(source="massbank"):
    return {
      "library_source":source, "source_record_accession":"SYNTHETIC-NOT-REAL",
      "source_url":"https://example.org/TEST-NOT-AN-AUTHENTIC-RECORD",
      "instrument_method":"synthetic capillary GC-MS test",
      "reuse_permission":"explicitly_permitted",
      "acquisition":"GC_EI","ionization":"EI","derivatization":"TMS",
      "retention_index":1500.0,"retention_index_system":"FAME",
      "chemical_id":"CHEBI:16797","peaks":PEAKS
    }
def query():
    p=gc()
    return {k:p[k] for k in ("acquisition","ionization","derivatization",
                              "retention_index","retention_index_system","peaks")}

class SpectrumBankTests(unittest.TestCase):
    def setUp(self):
        self.db=make_db(":memory:")
    def tearDown(self):
        self.db.close()
    def find(self,q,**opts):
        return search_spectral_candidates(
            self.db,q,mass_tolerance_da=.03,retention_tolerance=40,**opts
        )
    def test_unresolved_original_study_features_can_be_indexed_without_mass(self):
        study="synthetic_study"
        status=register_original_omic_feature(self.db,study=study,feature_id="629905",
           label="629905",assay="GC_EI",specimen="human_plasma_predose",
           provenance_url="https://example.org/synthetic-source",
           known_source_chemical_name=False)
        self.assertEqual(status,"UNRESOLVED_BIN_ID_ONLY_NO_SPECTRUM")
        register_original_omic_feature(self.db,study=study,feature_id="citrulline",
           label="citrulline",assay="GC_EI",specimen="human_plasma_predose",
           provenance_url="https://example.org/synthetic-source",
           known_source_chemical_name=True)
        out=summarize_registered_features(self.db,study)
        self.assertEqual(out["n_registered_source_features"],2)
        self.assertEqual(out["annotation_stages"]["UNRESOLVED_BIN_ID_ONLY_NO_SPECTRUM"],1)
        self.assertEqual(out["n_transporter_or_enzyme_function_inferred"],0)
        with self.assertRaises(SpectralKnowledgeError):
            register_original_omic_feature(self.db,study=study,feature_id="629905",
               label="629905",assay="GC_EI",specimen="human_plasma_predose",
               provenance_url="https://example.org/synthetic-source",
               known_source_chemical_name=False)

    def test_real_chemical_id_is_not_spectral_identity_without_evidence(self):
        add_reference(self.db,gc())
        r=self.find(query())
        self.assertEqual(r["status"],"CANDIDATES_NOT_STRUCTURAL_CONFIRMATION")
        self.assertEqual(len(r["candidates"]),1)
        self.assertEqual(r["candidates"][0]["spectral_cosine"],1.0)
        self.assertEqual(r["candidates"][0]["n_matched_ions"],8)
        self.assertEqual(r["candidates"][0]["identification_grade"],"PUTATIVE_SPECTRAL_MATCH_NOT_CONFIRMED")
        self.assertFalse(r["organ_function_inferred"])
        self.assertFalse(r["chemical_identity_confirmed"])
        self.assertFalse(r["in_vivo_PK_validated"])
        self.assertEqual(r["candidates"][0]["mechanistic_associations_if_chemical_ID_is_confirmed"],[])
    def test_tissue_transport_relation_requires_explicit_citation_and_chemical_id(self):
        add_reference(self.db,gc())
        link={"chemical_id":"CHEBI:16797","gene":"SLC22A2",
              "relation":"reported_biomarker_association","organ":"human_kidney",
              "evidence_url":"https://pubmed.ncbi.nlm.nih.gov/40445153/","source_verified":True}
        with self.assertRaises(SpectralKnowledgeError):
            add_mechanism_link(self.db,{**link,"source_verified":False})
        with self.assertRaises(SpectralKnowledgeError):
            add_mechanism_link(self.db,{**link,"chemical_id":"629905"})
        add_mechanism_link(self.db,link)
        result=self.find(query())
        association=result["candidates"][0]["mechanistic_associations_if_chemical_ID_is_confirmed"]
        self.assertEqual(association[0]["gene"],"SLC22A2")
        self.assertEqual(association[0]["relation"],"reported_biomarker_association")
        self.assertFalse(result["candidates"][0]["organ_function_or_PK_change_demonstrated"])
    def test_ei_requires_ri_compatible_derivatization_and_method(self):
        add_reference(self.db,gc())
        for changes in ({"retention_index":1750.0},
                        {"retention_index_system":"Kovats"},
                        {"derivatization":"MOX-TMS"}):
            q={**query(),**changes}
            self.assertEqual(self.find(q)["candidates"],[])
    def test_gc_and_lc_are_not_interchangeable(self):
        add_reference(self.db,gc())
        q={**query(),"acquisition":"LC_MSMS","ionization":"ESI",
           "polarity":"positive","precursor_mz":150.0}
        self.assertEqual(self.find(q,precursor_tolerance_da=.01)["candidates"],[])
        with self.assertRaises(SpectralKnowledgeError):
            add_reference(self.db,{**gc(source="gnps")})
    def test_missing_spectrum_is_not_mass_or_transport(self):
        r=triage_unassigned_feature_id("*629905")
        self.assertIsNone(r["nominal_mass_da"])
        self.assertIsNone(r["candidate_ptm"])
        self.assertIsNone(r["candidate_transporter"])
        self.assertIsNone(r["candidate_enzyme"])
        self.assertEqual(self.find({"source_feature_id":"629905"})["status"],
                         "ABSTAIN_RAW_SPECTRAL_PEAKS_MISSING")
    def test_source_permission_identity_and_peak_verification(self):
        for changes in ({"reuse_permission":"unknown"}, {"chemical_id":"metformin"},
                        {"peaks":[[100,1]]},{"retention_index":None},
                        {"source_url":"http://example.org/record"}):
            with self.subTest(changes=changes):
                with self.assertRaises(SpectralKnowledgeError):
                    add_reference(self.db,{**gc(),**changes})
        add_reference(self.db,gc())
        with self.assertRaises(SpectralKnowledgeError):
            add_reference(self.db,gc())
    def test_curated_real_references_are_not_fabricated_spectra(self):
        outcome=seed_bibliography(self.db)
        self.assertEqual(outcome["n_bibliographic_links_seeded"],9)
        self.assertEqual(outcome["spectral_references_imported"],0)
        self.assertEqual(self.db.execute("SELECT COUNT(*) FROM mechanism_link").fetchone()[0],9)
        self.assertEqual(self.db.execute("SELECT COUNT(*) FROM spectral_reference").fetchone()[0],0)
        self.assertEqual(self.find(query())["candidates"],[])
        # Verifiable ChEBI links, but no automatic assignment of 229 Rotroff IDs.
        self.assertFalse(outcome["identifies_Rotroff_unknown_ids"])
        add_reference(self.db,gc())
        matched=self.find(query())
        links=matched["candidates"][0]["mechanistic_associations_if_chemical_ID_is_confirmed"]
        self.assertEqual({x["gene"] for x in links},
                         {"SLC22A2","SLC47A1","SLC47A2"})
        self.assertTrue(all(x["relation"]=="reported_biomarker_association" for x in links))

    def test_absent_mass_resolution_or_ri_tolerance_fails(self):
        add_reference(self.db,gc())
        with self.assertRaises(SpectralKnowledgeError):
            self.find(query(),minimum_cosine=1.2)
        with self.assertRaises(SpectralKnowledgeError):
            search_spectral_candidates(self.db,query(),mass_tolerance_da=5,retention_tolerance=40)
        with self.assertRaises(SpectralKnowledgeError):
            search_spectral_candidates(self.db,query(),mass_tolerance_da=.03)
    def test_reference_catalogue_distinguishes_protein_spectra(self):
        self.assertIn("PROTEOMICS_MSMS",REFERENCE_CATALOGUE["pride"]["techniques"])
        self.assertNotIn("GC_EI",REFERENCE_CATALOGUE["pride"]["techniques"])
        self.assertFalse(REFERENCE_CATALOGUE["rhea"]["spectral"])
        self.assertNotIn("LC_MSMS",REFERENCE_CATALOGUE["binbase"]["techniques"])

if __name__=="__main__":
    unittest.main()
