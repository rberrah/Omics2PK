"""Original GNPS MGF: public nonpatient m/z/RT for selected source MS1 feature IDs.

Find original task spectra for selected feature IDs from nested
patient-heldout train-only predictor ranking; chemistry remains putative.
This file never reads/prints subjects or metabolite intensity observations.
"""
from __future__ import annotations

from hashlib import sha256
from zipfile import ZipFile
import json
import re
import sys

from scripts.validate_local_jarmusch14_gnps_biom_pk import (
    GNPS_TASK,REPACK_CANONICAL_MEMBER_SHA256,
)

FOCUSED_SOURCE_FEATURE_IDS=(
    "369","3018","4321","6135", # J9 urine full MS1 selected
    "6171","936",               # J1 plasma full MS1 selected
    "740","1889","2496",       # J1 urine putative bile path
    "454","6205",              # J1 plasma source annotated
)


def _source_selected_spectra(archive):
    with ZipFile(archive) as z:
        for member,digest in REPACK_CANONICAL_MEMBER_SHA256.items():
            if sha256(z.read(member)).hexdigest()!=digest:
                raise ValueError("Unrecognized public original GNPS source members")
        if GNPS_TASK not in z.read("params.xml").decode("utf-8"):
            raise ValueError("Original GNPS task mismatch")
        raw=z.read("spectra/specs_ms.mgf").decode("utf-8","replace")
    spectra={}
    current=None
    for line in raw.splitlines():
        val=line.strip()
        if val=="BEGIN IONS":
            if current is not None:raise ValueError("Malformed nested original spectra")
            current={}
        elif val=="END IONS":
            if current is None:raise ValueError("Original MGF spectrum parse mismatch")
            scan=current.get("SCANS","")
            if not scan:
                # A source title may contain original scan number:
                title=current.get("TITLE","")
                m=re.search(r"(?:scan[=: _-]*)(\d+)",title,re.I)
                scan=m.group(1) if m else ""
            if scan in FOCUSED_SOURCE_FEATURE_IDS:
                if scan in spectra:raise ValueError("Repeated original MS2 feature scan ID")
                mass=current.get("PEPMASS","").split()
                rt=current.get("RTINSECONDS")
                try:
                    mz=float(mass[0]) if mass else None
                    rtsec=float(rt) if rt is not None else None
                except ValueError:
                    mz=rtsec=None
                if mz is None or mz<=0:
                    raise ValueError("Invalid original precursor mass")
                spectra[scan]={
                    "original_public_GNPS_source_scan_id":scan,
                    "original_precursor_mz":round(mz,6),
                    "original_MS2_retention_time_seconds":(
                        round(rtsec,2) if rtsec is not None else None),
                    "original_precursor_charge":current.get("CHARGE"),
                    "original_source_identification_not_reference_standard":True,
                }
            current=None
        elif current is not None and "=" in val:
            key,value=val.split("=",1)
            current[key.strip().upper()]=value.strip()
    return {
        "status":"PUBLIC_GNPS_FEATURE_MZ_RT_IDENTITY_AUDIT_NO_PATIENT_DATA",
        "original_GNPS_task":GNPS_TASK,
        "n_requested_original_MS1_feature_ids":len(FOCUSED_SOURCE_FEATURE_IDS),
        "n_found_original_MGF_MS2_features":len(spectra),
        "feature_mz_rt":spectra,
        "not_found_source_scan_ids":[i for i in FOCUSED_SOURCE_FEATURE_IDS if i not in spectra],
        "no_patient_data_or_original_intensities_in_output":True,
        "chemical_identity_not_proven":True,
    }

if __name__=="__main__":
    if len(sys.argv)!=2:raise SystemExit("source ZIP path required")
    print(json.dumps(_source_selected_spectra(sys.argv[1]),indent=2,sort_keys=True))
