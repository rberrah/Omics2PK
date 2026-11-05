"""Separate original protein / proteoform / posttranslational evidence registry.

This is NOT an EI-spectrum annotation algorithm. PRIDE original PXD refs,
protein IDs and original fragment/peptide evidence must be independently
authenticated. A modified albumin abundance is NOT by itself a drug-binding
assay or a patient PK clearance measurement.
"""
from __future__ import annotations
import re
import sqlite3
from .spectral_knowledge_bank import SpectralKnowledgeError

UNIPROT=re.compile(r"^[A-Z0-9]{6,10}$")
PXD=re.compile(r"^PXD[0-9]{6,}$")
VALID_ASSAYS={"BOTTOM_UP_PEPTIDE","TOP_DOWN_INTACT_PROTEOFORM","TARGETED_PROTEOMICS"}

def register_protein_evidence(db, *, source_accession, protein_uniprot,
                              gene, assay_type, specimen, modification,
                              modification_site, source_url, evidence_level,
                              original_spectrum_accession=None):
    """Metadata-only index. Never infer transporter kinetics or protein activity."""
    if not isinstance(source_accession,str) or not PXD.fullmatch(source_accession):
        raise SpectralKnowledgeError("Authentic PXD ProteomeXchange source ID required")
    if not isinstance(protein_uniprot,str) or not UNIPROT.fullmatch(protein_uniprot):
        raise SpectralKnowledgeError("Original UniProt protein accession required")
    if not isinstance(gene,str) or not re.fullmatch(r"[A-Z0-9-]{2,25}",gene):
        raise SpectralKnowledgeError("Original protein target gene required")
    if assay_type not in VALID_ASSAYS:
        raise SpectralKnowledgeError("Unknown proteomics method; no cross-matching with EI-MS")
    for key,value in (("specimen",specimen),("modification",modification),
                      ("site",modification_site),("source_url",source_url)):
        if not isinstance(value,str) or not value.strip() or len(value)>256:
            raise SpectralKnowledgeError("Explicit original proteoform metadata required: "+key)
    if not source_url.startswith("https://"):
        raise SpectralKnowledgeError("Proteomics source must be HTTPS original source")
    if evidence_level not in ("source_reported_unvalidated","spectral_assigned",
                              "standard_confirmed"):
        raise SpectralKnowledgeError("Protein/PTM evidence grade must be explicit")
    if evidence_level in ("spectral_assigned","standard_confirmed"):
        if not isinstance(original_spectrum_accession,str) or not original_spectrum_accession.strip():
            raise SpectralKnowledgeError("MS-confirmed protein/PTM requires the actual original scan/PSM accession")
    if modification_site!="unlocalized" and not re.fullmatch(r"(?:[A-Z][1-9][0-9]{0,5}|multiple_sites)",modification_site):
        raise SpectralKnowledgeError("Site must be a qualified residue-index or explicitly unlocalized")
    if evidence_level=="standard_confirmed" and modification_site=="unlocalized":
        raise SpectralKnowledgeError("Unlocalized modification cannot be marketed as site-confirmed")
    try:
        with db:
            db.execute("INSERT INTO protein_proteoform_reference VALUES(?,?,?,?,?,?,?,?,?)",
                       (source_accession,protein_uniprot,gene,assay_type,specimen,
                        modification,modification_site,source_url,evidence_level))
    except sqlite3.IntegrityError:
        raise SpectralKnowledgeError("Duplicate original proteoform reference") from None
    return {
      "status":"PROTEOFORM_OR_PROTEIN_SOURCE_METADATA_REGISTERED",
      "source_assay":assay_type,
      "proteoform_site_confirmed":evidence_level=="standard_confirmed",
      "protein_function_or_transport_rate_measured":False,
      "albumin_drug_binding_or_in_vivo_PK_change_established":False,
      "spectral_data_copied":False,
    }
