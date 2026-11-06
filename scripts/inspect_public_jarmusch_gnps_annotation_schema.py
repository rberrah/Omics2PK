"""Public GNPS task annotated-spectral chemical schema audit; no subject rows."""
from __future__ import annotations
import csv
import io
import json
from pathlib import Path
from zipfile import ZipFile

# File names of spectral lists / column names only (not row contents).
SOURCE_LIMIT=150_000_000
STEROID_SIGNATURES=("cortisol","cortisone","dehydroepiandrosterone",
    "hydroxycort","bile acid","cholate","cholic acid","chenodeoxychol",
    "ursodeoxychol","taurochol","glycochol","testosterone",
    "progesterone","pregnenolone","steroid","hydroxysteroid")
FOLD_REPORT_COLUMNS=("compound_name","libraryid","cluster index","cluster_index",
                     "mzerrorppm","mqscore","ionmode","file","smiles","inchi",
                     "spectrumid")


def inspect_gnps_public_annotations(archive):
    from scripts.download_public_jarmusch_gnps import ARCHIVE_SHA256
    import hashlib
    path=Path(archive)
    if path.stat().st_size > SOURCE_LIMIT:
        raise ValueError("Unqualified public source ZIP")
    new_sha=hashlib.sha256(path.read_bytes()).hexdigest()
    if len(new_sha)!=64:
        raise ValueError("Invalid original archive checksum")
    # DownloadResult can recompress the same source members with different ZIP bytes.
    # This is a SOURCE REQUALIFICATION audit; not a supervised model execution.

    with ZipFile(path) as z:
        entries=[]
        for item in z.infolist():
            name=item.filename
            obj={
                "path":name[:180],
                "original_file_bytes":item.file_size,
                "crc32":f"{item.CRC:08x}",
            }
            if name in ("params.xml","qiime2_output/qiime2_metadata.tsv","qiime2_output/qiime2_table.qza"):
                obj["uncompressed_member_sha256"]=hashlib.sha256(z.read(item)).hexdigest()

            # No sample metadata columns or records are emitted.
            if name.lower().endswith((".tsv",".csv",".txt")) and not any(
                word in name.lower() for word in
                ("metadata", "sample", "attribute", "quant")):
                with z.open(item) as handle:
                    line=handle.readline(20_000).decode("utf-8-sig","replace")
                cols=next(csv.reader([line],delimiter="\t" if "\t" in line else ","))
                obj["header_columns"]=cols[:80]
                obj["n_header_fields"]=len(cols)
            entries.append(obj)
        return {
            "source":"original 2019 GNPS task b3966445abad4658a3cdc63c8198a666",
            "new_official_download_SHA256":new_sha,
            "archive_old_SHA256_matches":new_sha==ARCHIVE_SHA256,
            "n_original_zip_members":len(entries),
            "member_schemas_with_no_subject_rows":entries,
            "steroid_bile_candidate_keywords_for_future_targeted_spectral_hypothesis":
                list(STEROID_SIGNATURES),
            "no_individual_source_values_exported":True,
        }

if __name__=="__main__":
    import sys
    print(json.dumps(inspect_gnps_public_annotations(sys.argv[1]),sort_keys=True,indent=2))
