#!/usr/bin/env python3
"""Public BinDiscover API candidate metadata-only scan for source Rotroff bins.

The historical official documentation specifies:
  POST http://bindiscover.metabolomics.us:4999/binresource/
  JSON: {"bin_id": NUMBER}
This is a LIVE feasibility audit, not chemical identification. Source 6-digit
Rotroff bins might not share BinDiscover accession namespace/revision; outputs
never assert original-study bin identity without original method crosswalk.
No patient PK/omics observations, nor full mass spectra, are exported.
"""
from __future__ import annotations
import hashlib
import json
import re
import urllib.error
import urllib.request

BINS=("222169","221574","218765")
DOCUMENTATION="https://metabolomics-us.github.io/bindiscover/binbrowser.html"
OFFICIAL_ENDPOINT="http://bindiscover.metabolomics.us:4999/binresource/"
ALTERNATE_TLS_ENDPOINT="https://bindiscover.metabolomics.us/binresource/"
LIMIT=500_000

def summarize_source_reply(blob,source_bin):
    """Schema reconnaissance only: do not assume spectral provenance/identity."""
    if len(blob)>LIMIT or not blob:
        raise ValueError("Unqualified public BinDiscover payload")
    try:response=json.loads(blob)
    except (ValueError,UnicodeError) as error:
        raise ValueError("Not a documented BinDiscover JSON response") from error
    rows=response if isinstance(response,list) else [response]
    if len(rows)!=1 or not isinstance(rows[0],dict):
        raise ValueError("Unexpected BinDiscover response shape")
    row=rows[0]
    spectrum=row.get("spectrum")
    # Doc shows semicolon? actual binbrowser doc '85:0.007 86:..'
    n_peaks=None
    if isinstance(spectrum,str):
        n_peaks=sum(bool(re.fullmatch(r"\d+(?:\.\d+)?:\d+(?:\.\d+)?",token)) for token in spectrum.split())
    accession=row.get("bin_id",row.get("binId"))
    echoed=(str(accession)==source_bin) if accession is not None else None
    values={}
    for key in ("retentionIndex","kovats","quantMass","uniqueMass","purity","splash"):
        val=row.get(key)
        if val is None:continue
        if not isinstance(val,(int,float,str)) or len(str(val))>80:continue
        values[key]=val
    return {
      "source_original_study_feature_id":source_bin,
      "response_field_names":sorted(k for k in row.keys() if isinstance(k,str))[:30],
      "documented_bin_id_explicitly_echoed_and_matches":echoed,
      "n_spectral_ions_in_public_response":n_peaks,
      "reported_public_metadata":values,
      "source_response_sha256":hashlib.sha256(blob).hexdigest(),
      "spectrum_raw_values_exported":False,
      "original_Rotroff_bin_namespace_correspondence_proved":False,
      "original_Rotroff_ionization_derivatization_and_RI_method_equivalence_proved":False,
      "chemical_structure_identified":False,
      "transport_or_enzyme_activity_inferred":False,
    }

def probe_one(bin_id):
    if not re.fullmatch(r"\d{5,8}",bin_id):
        raise ValueError("Only original source numeric bin label accepted")
    attempts=[]
    for url in (OFFICIAL_ENDPOINT,ALTERNATE_TLS_ENDPOINT):
        data=json.dumps({"bin_id":int(bin_id)}).encode()
        req=urllib.request.Request(url,data=data,method="POST",
            headers={"Content-Type":"application/json","Accept":"application/json",
                     "User-Agent":"Omics2PK original Rotroff unknown BinDiscover public metadata audit"})
        try:
            with urllib.request.urlopen(req,timeout=8) as response:
                blob=response.read(LIMIT+1)
            parsed=summarize_source_reply(blob,bin_id)
            parsed["status"]="PUBLIC_BIN_METADATA_RETURNED_NOT_VERIFIED_AS_ORIGINAL_ROTR0FF_SOURCE"
            parsed["endpoint_kind"]="official_documented" if url==OFFICIAL_ENDPOINT else "tls_candidate_undocumented"
            return parsed
        except (urllib.error.URLError,TimeoutError,OSError,ValueError) as error:
            attempts.append({
                "endpoint_kind":"official_documented" if url==OFFICIAL_ENDPOINT else "tls_candidate_undocumented",
                "error_class":type(error).__name__,
                "http_code":getattr(error,"code",None)
            })
    return {"status":"ABSTAIN_NO_PUBLIC_CROSSWALKABLE_BIN_SPECTRUM",
            "source_original_study_feature_id":bin_id,
            "attempts":attempts,"chemical_identity_claimed":False}

def main():
    hits={bin_id:probe_one(bin_id) for bin_id in BINS}
    return {
      "status":"PUBLIC_BINDISCOVER_ID_API_FEASIBILITY_SCREEN_ONLY",
      "api_documentation":DOCUMENTATION,
      "n_original_SHAP_ranked_bins_queried":len(BINS),
      "n_returned_public_metadata_not_validated_as_Rotroff":
        sum(str(value.get("status","")).startswith("PUBLIC_BIN_METADATA_RETURNED") for value in hits.values()),
      "candidate_metadata":hits,
      "no_deconvoluted_spectra_copied":True,
      "no_original_study_to_spectrum_crosswalk_proved":True,
      "no_new_molecule_or_transporter_assigned":True,
    }
if __name__=="__main__":
    print(json.dumps(main(),indent=2,sort_keys=True))
