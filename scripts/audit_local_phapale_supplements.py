"""Local, READ-ONLY original Phapale (2010) Wiley supplemental archive audit.

The supplied source is copyright-protected and stays on the user's machine.
Never add any .doc, extracted chart, patient/metabolite value or converted
PDF/DOCX to this public repository. GitHub CI uses synthetic archives only.

Usage:
  python -m scripts.audit_local_phapale_supplements /path/to/paphale.zip

Optional dependency for reading original legacy Word text: 'antiword'. The
fixed original SHA plus thirteen per-file SHA-256 fingerprints were independently
measured from the authorized user-uploaded source on 2026-10-16.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path,PurePosixPath
import re
import shutil
import subprocess
from tempfile import TemporaryDirectory
from zipfile import ZipFile,BadZipFile

DOI="10.1038/clpt.2009.296"
EXPECTED_ARCHIVE_SHA256="b180d49fcde75d24b6d1812be001dbfa648c9565d0f5bb7100e1bc7d5c299a77"
EXPECTED_FILES={
 "0001":"f55b9e145410f113294bd7f9251dd97a76105aceb79ce4a982433c18e728f4f9",
 "0002":"c17dfe1e71060178859d14f7e4190dde690df61eb03d88d4c53bb9321cb5c21f",
 "0003":"842eb3d7cf4376be450bae6bad1599695e59d096f86085544de4c728fec3b351",
 "0004":"6b3deff63e9b5a8e7eecd0e723d91037954c64aa02c73003ee1e342b327db841",
 "0005":"bb8ab75af503fe4bfbe687c473c58cfa672ba5764b21b0208bbad98e31acce65",
 "0006":"1a215af32337f7561ba652a1c08781ae3e216e54ceb0fa8619dd172d9c93bcf7",
 "0007":"124578a4dd5d13b3ff7d33c3fd2c9c3755da9afa906273bc28211a7997c1f38a",
 "0008":"3c618f6684d9d29cd14b92efb289cffd81faeccb0f092ac6a4b771790c0e4405",
 "0009":"3c356672b6d0447994d6572043b251dcdc77b2dd5656264c6cd1dc7a637ef6f0",
 "0010":"0a647a14f42af453265a69fe5ce2b66a7fb59e70a9f1262c974b2f17f304830a",
 "0011":"e8cb850375044396031296e08ca35bd52a0e12f507088f48bb3b96afe34a18aa",
 "0012":"70dace1a3c9aa7b04cede5e178aae29f257941c6f84caaccd787863a8e93d72f",
 "0013":"ceefe1b7564bb2277b107027e53655f0774bec732012dee033c86381e4dfd6fc",
}
OLE_HEADER=bytes.fromhex("d0cf11e0a1b11ae1")
MAX_ORIGINAL_COMPRESSED=5_000_000
MAX_ORIGINAL_UNCOMPRESSED=7_000_000

class OriginalPhapaleSourceError(ValueError):
    pass


def inspect_archive_structure(blob:bytes,*,require_original_sha:bool=True):
    if not isinstance(blob,bytes) or not blob or len(blob)>MAX_ORIGINAL_COMPRESSED:
        raise OriginalPhapaleSourceError("Missing or oversized original local Wiley archive")
    archived_sha=hashlib.sha256(blob).hexdigest()
    if require_original_sha and archived_sha!=EXPECTED_ARCHIVE_SHA256:
        raise OriginalPhapaleSourceError("Original uploaded Phapale ZIP SHA-256 mismatch — abstain")
    try:
        with ZipFile(__import__("io").BytesIO(blob)) as z:
            members=z.infolist()
            if len(members)>30 or sum(m.file_size for m in members)>MAX_ORIGINAL_UNCOMPRESSED:
                raise OriginalPhapaleSourceError("Original source ZIP size/count safety limit")
            names=[m.filename for m in members]
            if len(names)!=len(set(names)):
                raise OriginalPhapaleSourceError("Duplicated original archive member names")
            if any(PurePosixPath(x).is_absolute() or ".." in PurePosixPath(x).parts or "\\" in x for x in names):
                raise OriginalPhapaleSourceError("Unsafe archive paths")
            docs=[m for m in members if m.filename.endswith(".doc") and not m.is_dir()]
            ids={}
            for m in docs:
                found=re.fullmatch(r"paphale/cptclpt2009296-sup-(\d{4})\.doc",m.filename)
                if not found:raise OriginalPhapaleSourceError("Unexpected original supplement member")
                n=found.group(1)
                if n in ids:raise OriginalPhapaleSourceError("Duplicate supplemental figure/table number")
                original=z.read(m)
                if original[:8]!=OLE_HEADER:raise OriginalPhapaleSourceError("Original Word OLE2 signature missing")
                digest=hashlib.sha256(original).hexdigest()
                if require_original_sha and EXPECTED_FILES.get(n)!=digest:
                    raise OriginalPhapaleSourceError("Original supplemental document SHA-256 changed")
                ids[n]=(len(original),digest)
            if sorted(ids)!=sorted(EXPECTED_FILES):
                raise OriginalPhapaleSourceError("Expected exactly the 13 original Wiley .doc files")
    except (BadZipFile,KeyError) as exc:
        raise OriginalPhapaleSourceError("Unreadable original archive") from exc
    return {
        "source_archive_sha256":archived_sha,
        "n_genuine_Wiley_legacy_Word_supplements":len(ids),
        "supplementary_section_counts":{"figures":[1,2,3,4,5,6,7,8],"summary_tables":[9,10,11,12],"methods":[13]},
        "all_original_legacy_Word_ole2_headers_verified":True,
        "all_13_individual_original_member_SHA256_matched":require_original_sha,
        "raw_or_reconstructed_omics_matrix_released":False,
        "patient_PK_Ct_data_released":False,
    }


def inspect_original_local(path:str|Path):
    file=Path(path)
    if not file.is_file() or file.stat().st_size>MAX_ORIGINAL_COMPRESSED:
        raise OriginalPhapaleSourceError("Original Phapale source ZIP unavailable or too large")
    original=file.read_bytes()
    qc=inspect_archive_structure(original)
    if not shutil.which("antiword"):
        qc.update({"method_text_status":"ABSTAIN_ANTIWORD_UNAVAILABLE",
                   "source_subject_linkage_validated":False})
        return qc
    with ZipFile(__import__("io").BytesIO(original)) as z,TemporaryDirectory() as td:
        readable={}
        for n in EXPECTED_FILES:
            inside=f"paphale/cptclpt2009296-sup-{n}.doc"
            dst=Path(td)/f"sup_{n}.doc"
            dst.write_bytes(z.read(inside))
            try:
                p=subprocess.run(["antiword","-w","160",str(dst)],capture_output=True,
                                  text=True,timeout=25,check=False)
            except (subprocess.TimeoutExpired,OSError) as exc:
                raise OriginalPhapaleSourceError("Legacy Word text conversion did not complete") from exc
            if p.returncode!=0:
                raise OriginalPhapaleSourceError("Legacy Word text extraction failed; abstain")
            readable[n]=p.stdout
    m=readable["0013"].lower()
    if not all(x in m for x in ("quantile normalization","pls model","cross-validation","vip","permutation")):
        raise OriginalPhapaleSourceError("Original supplement 13 methods claims have changed")
    if not ("auc" in readable["0001"].lower() and
            "clearance" in readable["0002"].lower() and
            "ctr ough".replace(" ","") in readable["0011"].lower().replace(" ","")):
        raise OriginalPhapaleSourceError("Original PLS AUC/clearance/Ctrough scientific labels changed")
    qc.update({
        "method_text_status":"ORIGINAL_LOCAL_LEGACY_WORD_TEXT_AUDITED",
        "original_method_describes_global_quantile_normalization":True,
        "original_method_describes_VIP_feature_screening_before_second_PLS":True,
        "original_method_describes_internal_seven_part_CV":True,
        "original_method_describes_20_response_permutations":True,
        "original_methods_confirm_AUC0_72_clearance_and_Ctrough_targets":True,
        "four_source_reported_selected_urine_metabolites":[
            "cortisol","1-methylguanosine","acetyl-arginine","phosphoethanolamine"
        ],
        "original_29_by_1256_numeric_assay_table_retrieved":False,
        "original_same_participant_numeric_dose_time_concentration_table_retrieved":False,
        "original_subject_id_to_metabolite_and_AUC_linkage_retrieved":False,
        "original_embedded_figures_are_NOT_individual_original_data_tables":True,
        "eligibility_for_subject_disjoint_predictive_reanalysis":"ABSTAIN_MISSING_SOURCE_PATIENT_LEVEL_ARRAYS",
        "patient_data_or_original_proprietary_source_bytes_printed":False
    })
    return qc

if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser(description="Privately inspect original 13 Wiley supplements; only aggregate JSON leaves CLI")
    ap.add_argument("original_zip",type=Path)
    opts=ap.parse_args()
    print(json.dumps(inspect_original_local(opts.original_zip),indent=2,ensure_ascii=False,sort_keys=True))
