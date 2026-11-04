"""Offline spectral-evidence registry for mechanistic Omics2PK hypotheses.

Explicitly SEPARATES:
  (1) a real chemical spectrum from metabolomics GC-EI or LC-MS/MS,
  (2) tentative spectral matches / candidate chemical IDs,
  (3) literature-backed chemical -> transporter/enzyme relations,
  (4) downstream biological tissue/PK validation, NEVER inferred here.

References are user-supplied; this module does NOT download/scrape libraries,
grant redistribution rights, identify patient's omics, or claim an MS level-1 ID.
No patient ID, intensities by patient or clinical PK is accepted as a record.
"""
from __future__ import annotations

from math import isfinite, sqrt
from pathlib import Path
import json
import re
import sqlite3

REFERENCE_CATALOGUE = {
    "binbase": {
        "role": "GC-EI-MS plus derivatization-specific retention-index library",
        "url": "https://bindiscover.metabolomics.us",
        "spectral": True,
        "techniques": ["GC_EI"],
        "access": "Check source terms, per-spectrum access and retention-index method",
    },
    "massbank": {
        "role": "Experimental small-molecule spectra, library records and metadata",
        "url": "https://massbank.eu/",
        "spectral": True,
        "techniques": ["GC_EI", "LC_MSMS"],
        "access": "Redistribute only records with rights explicitly verified",
    },
    "gnps": {
        "role": "LC-MS/MS similarity-based metabolite candidate references",
        "url": "https://gnps.ucsd.edu/",
        "spectral": True,
        "techniques": ["LC_MSMS"],
        "access": "Library-specific permissions and provenance required",
    },
    "pride": {
        "role": "Proteomics peptides/protein/PTM raw evidence; NOT small-molecule EI spectra",
        "url": "https://www.ebi.ac.uk/pride/",
        "spectral": True,
        "techniques": ["PROTEOMICS_MSMS"],
        "access": "Use PXD original study accessions and protein/peptide identifiers",
    },
    "rhea": {
        "role": "Expert-curated chemical/transport reactions, ChEBI/UniProt links",
        "url": "https://www.rhea-db.org/",
        "spectral": False,
        "techniques": [],
        "access": "Reaction evidence must be sourced individually",
    },
    "uniprot": {
        "role": "Protein, enzyme, transporter and PTM annotations",
        "url": "https://www.uniprot.org/",
        "spectral": False,
        "techniques": [],
        "access": "Target identity and evidence-level provenance required",
    },
}
ACQUISITIONS = {"GC_EI", "LC_MSMS"}
RELATIONS = {"experimentally_transport_substrate", "experimentally_enzyme_substrate",
             "literature_reported_transporter_substrate", "literature_reported_enzyme_substrate",
             "reported_biomarker_association", "curated_reaction_participant"}
CHEMICAL_ID = re.compile(r"^(?:CHEBI:\d+|HMDB\d{5,}|INCHIKEY:[A-Z]{14}-[A-Z]{10}-[A-Z])$")
GENE = re.compile(r"^[A-Z][A-Z0-9-]{1,24}$")

class SpectralKnowledgeError(ValueError):
    pass


def _positive_number(value, name, *, zero=False):
    if isinstance(value, bool):
        raise SpectralKnowledgeError(name + " must be numeric")
    try:
        f = float(value)
    except (TypeError, ValueError, OverflowError):
        raise SpectralKnowledgeError(name + " must be numeric") from None
    if not isfinite(f) or (f < 0 if zero else f <= 0):
        raise SpectralKnowledgeError(name + " nonfinite or outside valid range")
    return f


def _peaklist(raw):
    if not isinstance(raw, list) or not 6 <= len(raw) <= 1024:
        raise SpectralKnowledgeError("Require 6..1024 actual spectrum ion pairs; source ID alone is insufficient")
    seen=set()
    peaks=[]
    for point in raw:
        if not isinstance(point, (list, tuple)) or len(point)!=2:
            raise SpectralKnowledgeError("Each spectral peak is an actual (m/z,intensity) pair")
        mz=_positive_number(point[0],"m/z")
        intensity=_positive_number(point[1],"peak intensity")
        if mz in seen:
            raise SpectralKnowledgeError("Duplicate m/z peaks require spectrum deconvolution upstream")
        seen.add(mz)
        peaks.append((mz,intensity))
    return sorted(peaks)


def validate_spectrum(record):
    if not isinstance(record, dict):
        raise SpectralKnowledgeError("Spectral record must be a mapping")
    acq=record.get("acquisition")
    if acq not in ACQUISITIONS:
        raise SpectralKnowledgeError("GC_EI and LC_MSMS only; do not compare EI with MS/MS or proteins")
    source=record.get("library_source")
    if source not in REFERENCE_CATALOGUE or acq not in REFERENCE_CATALOGUE[source]["techniques"]:
        raise SpectralKnowledgeError("Library reference does not document this spectral method")
    for key in ("source_record_accession","source_url","instrument_method"):
        if not isinstance(record.get(key),str) or not record[key].strip():
            raise SpectralKnowledgeError("Missing original reference provenance: "+key)
    if not record["source_url"].startswith("https://"):
        raise SpectralKnowledgeError("Reference URL must be an explicit HTTPS source")
    if record.get("reuse_permission")!="explicitly_permitted":
        raise SpectralKnowledgeError("Reference spectrum cannot be copied without approved reuse status")
    if not isinstance(record.get("chemical_id"),str) or not CHEMICAL_ID.fullmatch(record["chemical_id"]):
        raise SpectralKnowledgeError("Stable exact CHEBI/HMDB/InChIKey chemical identity is required for reference")
    if acq=="GC_EI":
        if record.get("ionization")!="EI":
            raise SpectralKnowledgeError("GC-EI requires electron ionization")
        if not isinstance(record.get("derivatization"),str) or not record["derivatization"].strip():
            raise SpectralKnowledgeError("GC derivatization must be recorded, including 'none'")
        _positive_number(record.get("retention_index"),"experimental GC retention index")
        if not isinstance(record.get("retention_index_system"),str) or not record["retention_index_system"].strip():
            raise SpectralKnowledgeError("GC retention index standard/method is required")
    else:
        if record.get("ionization")!="ESI":
            raise SpectralKnowledgeError("LC-MS/MS prototype only covers ESI")
        _positive_number(record.get("precursor_mz"),"true precursor m/z")
        if record.get("polarity") not in ("positive","negative"):
            raise SpectralKnowledgeError("LC-MS/MS polarity is required")
    _peaklist(record.get("peaks"))
    return record


def make_db(path):
    """Local SQLite only: no preloaded private/publisher mass spectra."""
    db=sqlite3.connect(str(path))
    db.executescript("""
        CREATE TABLE IF NOT EXISTS spectral_reference (
          source TEXT NOT NULL, accession TEXT NOT NULL, acquisition TEXT NOT NULL,
          payload TEXT NOT NULL, chemical_id TEXT NOT NULL,
          PRIMARY KEY(source,accession));
        CREATE INDEX IF NOT EXISTS spectrum_chemical ON spectral_reference(chemical_id);
        CREATE TABLE IF NOT EXISTS protein_proteoform_reference (
          source_accession TEXT NOT NULL, protein_uniprot TEXT NOT NULL,
          gene TEXT NOT NULL, assay_type TEXT NOT NULL,
          specimen TEXT NOT NULL, modification TEXT NOT NULL,
          site TEXT NOT NULL, source_url TEXT NOT NULL,
          source_evidence_level TEXT NOT NULL,
          PRIMARY KEY(source_accession,protein_uniprot,modification,site));
        CREATE TABLE IF NOT EXISTS omic_feature (
          study TEXT NOT NULL, source_feature_id TEXT NOT NULL,
          assay TEXT NOT NULL, specimen TEXT NOT NULL,
          source_label TEXT NOT NULL, annotation_stage TEXT NOT NULL,
          provenance_url TEXT NOT NULL,
          PRIMARY KEY(study,source_feature_id));
        CREATE INDEX IF NOT EXISTS omic_feature_stage ON omic_feature(annotation_stage);
        CREATE TABLE IF NOT EXISTS mechanism_link (
          chemical_id TEXT NOT NULL, gene TEXT NOT NULL, relation TEXT NOT NULL,
          organ TEXT NOT NULL, evidence_url TEXT NOT NULL, source_verified INTEGER NOT NULL,
          PRIMARY KEY(chemical_id,gene,relation,organ,evidence_url));
    """)
    return db


def add_reference(db,record):
    r=validate_spectrum(record)
    # Stored only when user declared reuse rights; not official source verified.
    try:
        with db:
            db.execute("INSERT INTO spectral_reference VALUES(?,?,?,?,?)",
                       (r["library_source"],r["source_record_accession"],r["acquisition"],
                        json.dumps(r,sort_keys=True),r["chemical_id"]))
    except sqlite3.IntegrityError:
        raise SpectralKnowledgeError("Same source accession already registered; no silent overwrite") from None


def add_mechanism_link(db,link):
    if not isinstance(link,dict):
        raise SpectralKnowledgeError("Mechanism must be a mapping")
    cid=link.get("chemical_id")
    if not isinstance(cid,str) or not CHEMICAL_ID.fullmatch(cid):
        raise SpectralKnowledgeError("Chemical structure ID required before a transporter hypothesis")
    gene=link.get("gene")
    if not isinstance(gene,str) or not GENE.fullmatch(gene):
        raise SpectralKnowledgeError("Human transporter/enzyme gene symbol required")
    if link.get("relation") not in RELATIONS:
        raise SpectralKnowledgeError("Undefined evidence relation")
    organ=link.get("organ")
    if not isinstance(organ,str) or not organ.strip():
        raise SpectralKnowledgeError("Tissue/organ context cannot be inferred")
    url=link.get("evidence_url")
    if not isinstance(url,str) or not url.startswith("https://"):
        raise SpectralKnowledgeError("Original mechanism reference required")
    if link.get("source_verified") is not True:
        raise SpectralKnowledgeError("Literature link remains unqualified without source verification")
    with db:
        db.execute("INSERT OR IGNORE INTO mechanism_link VALUES(?,?,?,?,?,?)",
                   (cid,gene,link["relation"],organ,url,1))


def _cosine(peaks_a,peaks_b,tolerance_da):
    # One-to-one peak pairing, normalized by ALL ions in each spectrum.
    # Greedy by intensity product; no reuse of an ion. For an experimental
    # implementation this is a transparent conservative approximate matcher,
    # not GNPS/MassBank identification/recalibration.
    a=_peaklist(peaks_a);b=_peaklist(peaks_b)
    pairs=[]
    for i,(ma,ia) in enumerate(a):
        for j,(mb,ib) in enumerate(b):
            if abs(ma-mb)<=tolerance_da:
                pairs.append((ia*ib,i,j))
    pairs.sort(reverse=True)
    useda=set();usedb=set();dot=0.0
    for product,i,j in pairs:
        if i not in useda and j not in usedb:
            dot+=product;useda.add(i);usedb.add(j)
    norm=sqrt(sum(x*x for _,x in a)*sum(x*x for _,x in b))
    return dot/norm,len(useda)


def search_spectral_candidates(db,query,*,mass_tolerance_da,retention_tolerance=None,
                               precursor_tolerance_da=None,min_matched_peaks=6,
                               minimum_cosine=0.7):
    """Candidate generation, NEVER an MSI level-1 verified identification."""
    if not isinstance(query,dict):
        raise SpectralKnowledgeError("Query spectrum metadata required")
    if not isinstance(query.get("peaks"),list):
        return {"status":"ABSTAIN_RAW_SPECTRAL_PEAKS_MISSING",
                "candidates":[],"organ_function_inferred":False}
    acq=query.get("acquisition")
    if acq not in ACQUISITIONS:
        raise SpectralKnowledgeError("Unsupported source ionization / assay family")
    q=_peaklist(query["peaks"])
    mz_tol=_positive_number(mass_tolerance_da,"mass tolerance")
    if mz_tol>2.0:
        raise SpectralKnowledgeError("Mass tolerance unusually broad: check measurement resolution")
    if not isinstance(min_matched_peaks,int) or not 6<=min_matched_peaks<=len(q):
        raise SpectralKnowledgeError("Require at least six reproducible matched ions")
    c=_positive_number(minimum_cosine,"spectral similarity")
    if c>1:raise SpectralKnowledgeError("Cosine score exceeds one")
    if acq=="GC_EI":
        if query.get("ionization")!="EI" or not query.get("derivatization"):
            raise SpectralKnowledgeError("GC query needs EI ionization and original derivatization")
        ri=_positive_number(query.get("retention_index"),"query GC RI")
        ri_method=query.get("retention_index_system")
        if not isinstance(ri_method,str) or not ri_method:
            raise SpectralKnowledgeError("GC query must have explicit RI method")
        ri_tol=_positive_number(retention_tolerance,"declared RI tolerance")
    else:
        if query.get("ionization")!="ESI" or query.get("polarity") not in ("positive","negative"):
            raise SpectralKnowledgeError("LC-MS/MS query ESI polarity is required")
        prec=_positive_number(query.get("precursor_mz"),"query precursor m/z")
        prec_tol=_positive_number(precursor_tolerance_da,"precursor tolerance")
    out=[]
    for (raw,) in db.execute("SELECT payload FROM spectral_reference WHERE acquisition=?",(acq,)):
        ref=json.loads(raw)
        if acq=="GC_EI":
            if (ref["derivatization"]!=query["derivatization"] or
                ref["retention_index_system"]!=ri_method or
                abs(ri-float(ref["retention_index"]))>ri_tol):
                continue
        else:
            if (ref["polarity"]!=query["polarity"] or
                abs(prec-float(ref["precursor_mz"]))>prec_tol):
                continue
        sim,nmatch=_cosine(q,ref["peaks"],mz_tol)
        if nmatch<min_matched_peaks or sim+1e-12<c:continue
        links=[{"gene":gene,"relation":relation,"organ":organ,
                "evidence_url":url,"source_association_verified":bool(verified)}
               for gene,relation,organ,url,verified in db.execute(
                   "SELECT gene,relation,organ,evidence_url,source_verified FROM mechanism_link WHERE chemical_id=?",
                   (ref["chemical_id"],))]
        out.append({"candidate_chemical_id":ref["chemical_id"],
                    "reference_source":ref["library_source"],
                    "reference_accession":ref["source_record_accession"],
                    "spectral_cosine":round(min(1.0,sim),6),
                    "n_matched_ions":nmatch,
                    "identification_grade":"PUTATIVE_SPECTRAL_MATCH_NOT_CONFIRMED",
                    "mechanistic_associations_if_chemical_ID_is_confirmed":links,
                    "organ_function_or_PK_change_demonstrated":False})
    out.sort(key=lambda v:(-v["spectral_cosine"],-v["n_matched_ions"],v["reference_accession"]))
    return {"status":"CANDIDATES_NOT_STRUCTURAL_CONFIRMATION" if out else "NO_QUALIFIED_SPECTRAL_MATCH",
            "candidates":out,"organ_function_inferred":False,
            "chemical_identity_confirmed":False,"in_vivo_PK_validated":False}


def triage_unassigned_feature_id(feature_id):
    """Original clinical matrix signal IDs are not mass-to-charge measurements."""
    if not isinstance(feature_id,str) or not feature_id.strip():
        raise SpectralKnowledgeError("Original feature identifier required")
    return {"feature_id":feature_id, "status":"UNASSIGNED_NEEDS_RAW_EI_SPECTRUM_AND_RETENTION_INDEX",
            "nominal_mass_da":None,"chemical_id":None,
            "candidate_ptm":None,"candidate_transporter":None,
            "candidate_enzyme":None,"mechanism_validated":False,
            "next_input":"source chromatographic deconvoluted EI peak list, experimental RI system, derivatization, original BinBase reference"}


def register_original_omic_feature(db, *, study, feature_id, label, assay, specimen,
                                   provenance_url, known_source_chemical_name=False):
    """Original study feature namespace; NO patient intensities or invented mass."""
    for name,value in (("study",study),("feature_id",feature_id),("label",label),
                       ("assay",assay),("specimen",specimen),("provenance_url",provenance_url)):
        if not isinstance(value,str) or not value.strip() or len(value)>250:
            raise SpectralKnowledgeError("Missing or invalid original feature metadata: "+name)
    if not provenance_url.startswith("https://"):
        raise SpectralKnowledgeError("Source original feature list must have an HTTPS provenance URL")
    if assay not in ("GC_EI","LC_MSMS","TRANSCRIPTOME","PROTEOME","PROTEOFORM"):
        raise SpectralKnowledgeError("Unsupported original omic feature assay")
    if not isinstance(known_source_chemical_name,bool):
        raise SpectralKnowledgeError("Source chemical name classification must be explicitly true/false")
    stage=("SOURCE_CHEMICAL_NAME_NOT_STANDARD_CONFIRMED" if known_source_chemical_name
           else "UNRESOLVED_BIN_ID_ONLY_NO_SPECTRUM")
    try:
        with db:
            db.execute("INSERT INTO omic_feature VALUES(?,?,?,?,?,?,?)",
                       (study,feature_id,assay,specimen,label,stage,provenance_url))
    except sqlite3.IntegrityError:
        raise SpectralKnowledgeError("Original study assay feature key duplicated") from None
    return stage


def summarize_registered_features(db, study):
    """Emit only aggregate safe counts, no individual rows or sample labels."""
    counts=dict(db.execute(
        "SELECT annotation_stage,COUNT(*) FROM omic_feature WHERE study=? GROUP BY annotation_stage",
        (study,)).fetchall())
    return {
      "study":study,
      "n_registered_source_features":sum(counts.values()),
      "annotation_stages":counts,
      "n_assigned_mass_spectra_in_this_source_feature_registry":0,
      "n_original_PK_measurements_in_registry":0,
      "n_transporter_or_enzyme_function_inferred":0,
    }
