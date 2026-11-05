"""Download original publicly linked GNPS job output in ephemeral CI only.

Official published Jarmusch 2020 GNPS task, DOI 10.1111/cts.12785.
Check exact SHA of previously independently qualified public export before use.
Never publish the original archive or individual volunteer/sample metadata.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen, Request

SOURCE_URL=(
  "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?"
  "task=b3966445abad4658a3cdc63c8198a666&view=download_cytoscape_data"
)
ARCHIVE_SHA256="2047ea2032d73e3b48e2696fb1f652a7dedaa6f8cd5d40d0926d810bdd46c026"
MAX_BYTES=150_000_000
MIN_BYTES=50_000_000


def get_public_archive(path: Path):
    digest=hashlib.sha256()
    size=0
    req=Request(SOURCE_URL,headers={
        "User-Agent":"Omics2PK-academic-source-reproduction/1.0",
        "Accept":"application/zip,application/octet-stream,*/*"})
    try:
        with urlopen(req,timeout=160) as response, path.open("wb") as dest:
            while True:
                data=response.read(256*1024)
                if not data:break
                size+=len(data)
                if size>MAX_BYTES:
                    raise ValueError("GNPS original source larger than allowed original maximum")
                digest.update(data)
                dest.write(data)
        if not MIN_BYTES<size<MAX_BYTES or digest.hexdigest()!=ARCHIVE_SHA256:
            raise ValueError("Public GNPS download is not byte-identical to qualified 2019 source")
        return {
            "status":"ORIGINAL_PUBLIC_GNPS_ARTIFACT_SHA256_VERIFIED",
            "n_downloaded_bytes":size,
            "original_SHA256":ARCHIVE_SHA256,
            "no_sample_rows_exported":True,
        }
    except Exception:
        path.unlink(missing_ok=True)
        raise


if __name__=="__main__":
    if len(sys.argv)!=2:
        raise SystemExit("usage: python -m scripts.download_public_jarmusch_gnps /tmp/original.zip")
    print(json.dumps(get_public_archive(Path(sys.argv[1])),sort_keys=True))
