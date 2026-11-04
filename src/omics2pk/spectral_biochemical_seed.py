"""Small pre-reviewed clinical transporter/enzymatic reference SEEDS.

NO synthetic spectrum and NO prediction is introduced. Links are retrievable
ONLY after a chemical ID has independently passed the spectral candidate stage
(and still remain CONDITIONAL until chemical structure is confirmed).
Verified source means bibliographic link has been checked, NOT that a patient
sample was confirmed. Every molecule here has an external ChEBI reference.
"""
from __future__ import annotations
from .spectral_knowledge_bank import add_mechanism_link

# PubChem verified mapping: 1-methylnicotinamide CHEBI:16797;
# 1-methyladenosine CHEBI:16020; tacrolimus CHEBI:61049.
# Do not include medications among putative unknown Rotroff peaks.
_CURATED = (
    # Human renal OCT2/MATE endogenous biomarker axes, from clinical DDI study:
    # Nishii et al. 2025 doi:10.1016/j.xphs.2025.103776.
    ("CHEBI:16797", "SLC22A2", "reported_biomarker_association",
     "human_kidney", "https://pubmed.ncbi.nlm.nih.gov/40445153/"),
    ("CHEBI:16797", "SLC47A1", "reported_biomarker_association",
     "human_kidney", "https://pubmed.ncbi.nlm.nih.gov/40445153/"),
    ("CHEBI:16797", "SLC47A2", "reported_biomarker_association",
     "human_kidney", "https://pubmed.ncbi.nlm.nih.gov/40445153/"),
    ("CHEBI:16020", "SLC22A2", "reported_biomarker_association",
     "human_kidney", "https://pubmed.ncbi.nlm.nih.gov/40445153/"),
    ("CHEBI:16020", "SLC47A1", "reported_biomarker_association",
     "human_kidney", "https://pubmed.ncbi.nlm.nih.gov/40445153/"),
    ("CHEBI:16020", "SLC47A2", "reported_biomarker_association",
     "human_kidney", "https://pubmed.ncbi.nlm.nih.gov/40445153/"),
    # Tacrolimus: literature summaries of ABCB1 and CYP3A4/5 human PK.
    # Do not claim a fully specific CYP3A5 formation reaction or patient PK.
    ("CHEBI:61049", "ABCB1", "literature_reported_transporter_substrate",
     "human_intestine", "https://pmc.ncbi.nlm.nih.gov/articles/PMC4119065/"),
    ("CHEBI:61049", "CYP3A4", "literature_reported_enzyme_substrate",
     "human_liver_and_intestine", "https://pmc.ncbi.nlm.nih.gov/articles/PMC4119065/"),
    ("CHEBI:61049", "CYP3A5", "literature_reported_enzyme_substrate",
     "human_liver_and_intestine", "https://pmc.ncbi.nlm.nih.gov/articles/PMC4119065/"),
)


def seed_bibliography(db):
    """Nine independent literature assertions, zero spectral references."""
    for chemical_id, gene, relation, organ, reference in _CURATED:
        add_mechanism_link(db, {
            "chemical_id":chemical_id, "gene":gene, "relation":relation,
            "organ":organ, "evidence_url":reference, "source_verified":True,
        })
    return {"n_bibliographic_links_seeded":len(_CURATED),
            "spectral_references_imported":0,
            "patient_PK_function_proven":False,
            "tissue_specific_activity_estimated":False,
            "identifies_Rotroff_unknown_ids":False}
