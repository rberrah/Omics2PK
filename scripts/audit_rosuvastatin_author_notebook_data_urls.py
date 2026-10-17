"""Read SOURCE-CODE CELLS ONLY of original rosuvastatin first-author GitHub notebook.

Silveira et al 2021 DOI 10.3389/fphar.2021.752960 explicitly references
https://github.com/GustavoHBDuarte/Metabolomics_PK_proj
as the original Elastic Net model implementation. Its 21 MB notebook was
not readable through GitHub Contents; official git-pinned raw notebook
can reveal publicly shared *processed* patient×feature matrix links.

We NEVER print, rehost, parse or fit with notebook cell outputs, embedded
private medical data, subject IDs, entire source cell text or raw data bytes.
Only show counts/hostnames/URL-type and allowlisted public repository-path URLs.
Private Google Drive IDs/query credentials remain masked.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
import json
import re
from urllib.parse import urlsplit
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

PIN="dbf9182b0f72f8b290f623f7cc6560c5a5f99d73"
SOURCE=("https://raw.githubusercontent.com/GustavoHBDuarte/"
        "Metabolomics_PK_proj/"+PIN+"/Metabolomics_Elastic_net.ipynb")
CAP=28_000_000
REGEX=re.compile(r"https?://[^\s'\"<>\\)\]}]+",re.I)
FILE_EXT=re.compile(r"\.(csv|xlsx|xls|tsv|json|zip|parquet)(?:$|\?)",re.I)


def classify_notebook_payload(blob):
    if len(blob)>CAP:
        raise ValueError("Original notebook exceeds guarded size limit")
    doc=json.loads(blob)
    if not isinstance(doc,dict) or not isinstance(doc.get("cells"),list):
        raise ValueError("Not a Jupyter original notebook")
    counters=Counter()
    safe_references=set()
    data_loaders=Counter()
    contains_embedded_outputs=0
    for c in doc["cells"]:
        if not isinstance(c,dict):
            continue
        if c.get("outputs"):
            contains_embedded_outputs+=1
        raw=c.get("source",[])
        code="".join(raw) if isinstance(raw,list) else str(raw)
        if c.get("cell_type")!="code":
            continue
        counters["code_cells"]+=1
        for loader in ("pd.read_csv(","pd.read_excel(","read_csv(","read_excel(",
                       "wget ","curl ","pd.read_parquet(","requests.get("):
            data_loaders[loader]+=code.count(loader)
        for hit in REGEX.finditer(code):
            url=hit.group(0)
            parsed=urlsplit(url)
            host=(parsed.hostname or "").lower()
            counters["source_url_"+host]+=1
            if (host in ("raw.githubusercontent.com","github.com","www.ebi.ac.uk",
                         "ebi.ac.uk","www.omicsdi.org")
                    and FILE_EXT.search(parsed.path)):
                # Only public repository filenames; stripped of tokens/fragments.
                safe_references.add(parsed.scheme+"://"+parsed.netloc+parsed.path)
            elif host in ("drive.google.com","docs.google.com","www.dropbox.com",
                          "dropbox.com","colab.research.google.com"):
                counters["opaque_drive_or_colab_reference"]+=1
            elif host in ("raw.githubusercontent.com","github.com"):
                counters["other_github_url"]+=1
    return {
        "author_github_repo":"GustavoHBDuarte/Metabolomics_PK_proj",
        "original_notebook_commit":PIN,
        "original_notebook_sha256":sha256(blob).hexdigest(),
        "notebook_total_cells":len(doc["cells"]),
        "code_cells":counters["code_cells"],
        "n_cells_with_uninspected_embedded_outputs":contains_embedded_outputs,
        "code_data_loader_counts":{k:v for k,v in data_loaders.items() if v},
        "code_source_url_host_counts":{
            k.removeprefix("source_url_"):v for k,v in counters.items()
            if k.startswith("source_url_")},
        "n_masked_cloud_drive_references":counters["opaque_drive_or_colab_reference"],
        "public_source_data_url_paths":sorted(safe_references)[:35],
        "no_original_notebook_outputs_or_patient_data_exported":True,
        "disclaimer":"A URL in code is not a verified original subject×omics×PK matrix or a permitted reshare.",
    }


def run():
    report={
        "study":"rosuvastatin human n40 predose metabolomics and quantitative AUC/Cmax",
        "article_doi":"10.3389/fphar.2021.752960",
        "author_repo":"https://github.com/GustavoHBDuarte/Metabolomics_PK_proj",
        "actual_processed_omic_pk_matrix_verified":False,
    }
    try:
        with urlopen(Request(SOURCE,headers={"User-Agent":"Omics2PK-original-author-code-source-qualification/1.0"}),timeout=65) as r:
            raw=r.read(CAP+1)
        report["source"]="ORIGINAL_PINNED_AUTHOR_GITHUB_NOTEBOOK_READ_ONLY"
        report["original_source_structure"]=classify_notebook_payload(raw)
    except (HTTPError,URLError,OSError,TimeoutError,ValueError,UnicodeError,json.JSONDecodeError) as err:
        report["source"]="SOURCE_UNAVAILABLE_ABSTAIN"
        report["error_type"]=type(err).__name__
    return report


if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,sort_keys=True,indent=2))
