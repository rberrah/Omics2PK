"""Read-only local audit of the original Degraeve 2024 Wiley supplement ZIP.

Original user/rights-provided bytes are NEVER uploaded or redistributed.
No individual clinical patient, PK record, ASV abundance vector or GraphPad
embedded source object is exported. Only known public aggregate metadata
and structure from the exact verified original archive is returned.

Run locally: python -m scripts.audit_local_degraeve_supplement /path/to/cpt3077-sup-0001-supinfo.zip
"""
from __future__ import annotations

from hashlib import sha256
import io
import json
from pathlib import Path
import re
from xml.etree import ElementTree as ET
from zipfile import ZipFile

DOI="10.1002/cpt.3077"
EXPECTED_SHA="03aab5deb292a03753073e493d7ae49896884be2cde22d5c2d1c9f3037488cc2"
EXPECTED_ARCHIVE_FILES={
    "2023-0506-s01.docx",
    "2023-0506-s02.pptx",
    "2023-0506-s03.docx",
}
WORD="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
PPT="{http://schemas.openxmlformats.org/presentationml/2006/main}"
MAX_COMPRESSED=10_000_000
MAX_UNCOMPRESSED=40_000_000


class SupplementSourceError(ValueError):
    pass


def _ensure_safe_zip(blob:bytes, expected:set[str]|None=None) -> ZipFile:
    if len(blob)>MAX_COMPRESSED or blob[:2]!=b"PK":
        raise SupplementSourceError("Original supplement must be a bounded ZIP")
    z=ZipFile(io.BytesIO(blob))
    members=z.infolist()
    names=[e.filename for e in members]
    if (len(names)>150 or len(names)!=len(set(names))
            or sum(e.file_size for e in members)>MAX_UNCOMPRESSED):
        raise SupplementSourceError("Supplement ZIP has unsafe size/count")
    if expected is not None and set(names)!=expected:
        raise SupplementSourceError("Original 3-file supplemental archive shape changed")
    if any(n.startswith("/") or n.startswith("\\") or
           ".." in n.replace("\\","/").split("/") for n in names):
        raise SupplementSourceError("ZIP contains unsafe traversal paths")
    return z


def _word_text(root:ET.Element)->str:
    return "".join(n.text or "" for n in root.iter(WORD+"t"))


def _docx_structure(blob:bytes):
    with _ensure_safe_zip(blob) as z:
        if "word/document.xml" not in z.namelist():
            raise SupplementSourceError("Missing genuine source DOCX document.xml")
        root=ET.fromstring(z.read("word/document.xml"))
        tables=root.findall(".//"+WORD+"tbl")
        extracted=[]
        for table in tables:
            rows=[]
            for r in table.findall("./"+WORD+"tr"):
                rows.append([_word_text(c) for c in r.findall("./"+WORD+"tc")])
            extracted.append(rows)
        text=_word_text(root)
        return text,extracted


def _powerpoint_structure(blob:bytes):
    with _ensure_safe_zip(blob) as z:
        names=z.namelist()
        slides=[n for n in names if re.fullmatch(r"ppt/slides/slide[0-9]+\.xml",n)]
        # GraphPad Prism9 OLE binary may contain graph points, but these are
        # NOT a verified n93 original patient ID ↔ PK ↔ 16S crosswalk.
        prism=[n for n in names if re.fullmatch(r"ppt/embeddings/oleObject[0-9]+\.bin",n)]
        embedded_prism=sum(b"Prism9.Document" in z.read(n) for n in prism)
        return {
            "n_supplementary_figure_slides":len(slides),
            "n_embedded_OLE_objects":len(prism),
            "n_ole_objects_recognized_as_GraphPad_Prism9":embedded_prism,
            "graphics_are_not_a_released_subject_pk_omics_matrix":True,
            "graphpad_original_data_not_exported_or_decoded":True,
        }


def audit_original_authorized_supplement(path: str|Path):
    file=Path(path)
    if not file.is_file() or file.stat().st_size>MAX_COMPRESSED:
        raise SupplementSourceError("Original locally authorized supplement file missing or too large")
    blob=file.read_bytes()
    exact_sha=sha256(blob).hexdigest()
    if exact_sha!=EXPECTED_SHA:
        raise SupplementSourceError("Original authorized Wiley supplement SHA256 changed; manually requalify")
    with _ensure_safe_zip(blob,EXPECTED_ARCHIVE_FILES) as z:
        method, method_tables=_docx_structure(z.read("2023-0506-s01.docx"))
        table_text,tables=_docx_structure(z.read("2023-0506-s03.docx"))
        figure=_powerpoint_structure(z.read("2023-0506-s02.pptx"))
        if method_tables:
            raise SupplementSourceError("Unexpected individual source tables in methods document")
        if len(tables)!=3 or [len(t) for t in tables]!=[15,30,5]:
            raise SupplementSourceError("Unexpected supplementary table counts; requalify source")
        if figure["n_supplementary_figure_slides"]!=8:
            raise SupplementSourceError("Unexpected original figure slide count")
        # This is a fixed aggregate figure already reported in the paper,
        # not a reconstructed per-patient source dataset.
        dose_rows=tables[2]
        dose_r2=[]
        for row in dose_rows[1:]:
            if len(row)<6:raise SupplementSourceError("Missing original dose model table columns")
            r2=row[5].strip().replace(",",".")
            if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?%",r2):
                raise SupplementSourceError("Dose supplemental R2 is not parseable")
            dose_r2.append(float(r2[:-1]))
        if dose_r2!=[18.5,26.3,34.9,38.1]:
            raise SupplementSourceError("Original supplementary DOSE model R2 deviates")
        needed=("In total 100 patients","Eventually, 93 patients",
                "3A5*3/*3","93 patients")
        if not all(term in method+table_text for term in needed):
            raise SupplementSourceError("Original clinical n100→93/CYP3A5 study evidence not found")
        if ("Table S3." not in table_text
                or "inter-individual variability (IIV) in tacrolimus (TAC) dose" not in table_text):
            raise SupplementSourceError("Cannot establish S3 is a DOSE model not an AUC model")
        if not all(term in dose_rows[-1][1] for term in ("Age","Haematocrit","ASV 1508","ASV 664")):
            raise SupplementSourceError("Original dose-table selected covariates differ")
    return {
       "status":"QUALIFIED_ORIGINAL_SUPPLEMENT_AGGREGATE_METHODS_NO_PATIENT_PK_16S_CROSSWALK",
       "article_doi":DOI,
       "archive_sha256":exact_sha,
       "archive_bytes":len(blob),
       "n_documents":3,
       "supplement_methods_tables":len(method_tables),
       "n_supplementary_word_tables":len(tables),
       "n_table_s1_rows":len(tables[0]),
       "n_table_s2_rows":len(tables[1]),
       "n_table_s3_rows":len(tables[2]),
       "reported_initial_recruited_patients":100,
       "excluded_unusable_microbiome":3,
       "excluded_protocol_ineligible":4,
       "n_analyzed_independent_humans":93,
       "clinical_PK_unit":"individual C0 measured in whole blood; AUC0-24 limited-sampling model estimated",
       "cyp3a5_study_population":"CYP3A5*3/*3 nonexpressers",
       "s3_endpoint":"TAC daily DOSE (mg), NOT AUC0-24",
       "s3_stepwise_dose_model_in_sample_total_R2_percent":dose_r2,
       "s3_final_covariates":["Age","Haematocrit","ASV 1508","ASV 664"],
       "main_article_auc_model_R2_percent":"52.1; separate main-article Table 2, not inside supplement S3",
       "n_actual_individual_C0_C1_C3_rows_in_supplement":0,
       "n_original_patient_ID_to_16S_SRA_accession_mapping_rows_in_supplement":0,
       "n_public_individual_ASV_abundance_matrices_in_supplement":0,
       "figure_archive":figure,
       "n_original_individual_patient_measurements_exported":0,
       "ready_for_original_omics_to_PopPK_fit":False,
       "next_required":"Authorized original n93 patient_id↔16S alias plus dose, formulation, observed C0/C1/C3, dose-clock, actual vs imputed AUC and subject-level ASV matrix",
    }


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("original_supplement_zip",type=Path)
    args=parser.parse_args()
    print(json.dumps(audit_original_authorized_supplement(args.original_supplement_zip),
                     sort_keys=True,indent=2,ensure_ascii=False))
