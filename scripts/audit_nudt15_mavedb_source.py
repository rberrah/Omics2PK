#!/usr/bin/env python3
"""Print aggregate-only provenance and schema of public NUDT15 paired MaveDB.

Original Suiter et al. 2020 PNAS source and MaveDB paired meta-score set.
This is a schema audit, NOT yet a model or patient PD result.
No individual variant, score, patient, or matrix is printed or committed.
"""
import csv
import io
import json
import hashlib
import re
from collections import Counter
from urllib.request import Request,urlopen

URN="urn:mavedb:00000055-0-1"
ENDPOINT="https://api.mavedb.org/api/v1/score-sets/"+URN+"/scores"
def source_audit():
    request=Request(ENDPOINT,headers={"User-Agent":"Omics2PK-source-audit-20261105/1.0",
                                      "Accept":"text/csv"})
    with urlopen(request,timeout=40) as r:
        b=r.read(5_000_001)
        ct=r.headers.get("Content-Type","")
    if len(b)>5_000_000 or len(b)<500:
        raise RuntimeError("MaveDB public paired data body size invalid")
    try:
        s=b.decode("utf-8-sig")
    except UnicodeDecodeError as e:
        raise RuntimeError("Original public source not UTF-8") from e
    reader=csv.DictReader(io.StringIO(s))
    headers=reader.fieldnames
    if headers is None or "score" not in " ".join(headers).lower():
        raise RuntimeError("Original source missing score headers")
    rows=list(reader)
    if not 2800<=len(rows)<=3400:
        raise RuntimeError("Unexpected MaveDB paired dataset row count")
    from omics2pk.public_nudt15_variant_pd_benchmark import _variant
    def sig(hgvs):
        # Aggregate only: remove all residue names and amino-acid
        # positions, no original variant identity leaves the job.
        if hgvs is None:return "absent"
        v=re.sub(r"[0-9]+","N",hgvs)
        v=re.sub(r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)","AAA",v)
        v=re.sub(r"[ACDEFGHIKLMNPQRSTVWY]","X",v)
        return v[:80]
    missed=Counter(sig(z["hgvs_pro"]) for z in rows
                   if z["abundance_score"] not in ("","NA") and
                      z["sensitivity_score"] not in ("","NA") and
                      _variant(z["hgvs_pro"]) is None)
    return {
        "nonparsed_paired_protein_variant_HGVS_SHAPES_NO_IDS":dict(missed),
        "status":"ORIGINAL_MAVEDB_URN_SCHEME_QUALIFIED_NO_NUMERICAL_BENCHMARK_YET",
        "accession":URN,
        "original_source":"Suiter et al 2020 PNAS 10.1073/pnas.1915680117",
        "http_content_type":ct,
        "source_SHA256":hashlib.sha256(b).hexdigest(),
        "n_source_rows":len(rows),
        "source_headers":headers,
        "first_10_rows_not_exported":True,
        "original_value_records_exported":False,
        "source_is_real_patient_PD":False,
        "measures":"variant protein stability and thiopurine sensitivity in engineered HEK293T",
    }
if __name__=="__main__":
    print(json.dumps(source_audit(),indent=2,sort_keys=True))
