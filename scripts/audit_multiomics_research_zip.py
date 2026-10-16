"""Audit actual MultiOmics ZIP source metadata versus exchange, aggregate only."""
import csv
import hashlib
import io
import json
import sys
import zipfile


def audit_zip(path):
    if __import__("pathlib").Path(path).stat().st_size > 140_000_000:
        raise ValueError("Research ZIP too large for bounded source verifier")
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if (len(names) > 500 or len(set(names)) != len(names)
                or sum(f.file_size for f in z.infolist()) > 400_000_000):
            raise ValueError("Unsafe or oversized ZIP envelope")
        if any(f.is_dir() or __import__("pathlib").PurePosixPath(f.filename).is_absolute()
               or ".." in __import__("pathlib").PurePosixPath(f.filename).parts for f in z.infolist()):
            raise ValueError("Unsafe research ZIP member name")
        if len(names) > 500 or len(set(names)) != len(names):
            raise ValueError("Invalid ZIP member count")
        required = ("provenance/manifest.json", "source/projet_reproductible.json",
                    "exchange/manifest.json", "exchange/samples.csv")
        if not all(k in names for k in required):
            raise ValueError("Missing original archive sections")
        outer = json.loads(z.read(required[0]))
        if outer.get("downstream_exchange",{}).get("status")!="included":
            raise ValueError("No compatible native exchange")
        for name, info in outer["files"].items():
            if name not in names:
                raise ValueError("Missing declared archive member")
            if hashlib.sha256(z.read(name)).hexdigest()!=info["sha256"]:
                raise ValueError("Archive file integrity mismatch")
        original = json.loads(z.read(required[1]))
        metadata = original["inputs"]["metadata"]["content"]
        source = list(csv.DictReader(io.StringIO(metadata)))
        exchange = list(csv.DictReader(io.StringIO(z.read(required[3]).decode("utf-8-sig"))))
        if len(source)!=len(exchange) or {r["subject_id"] for r in source}!={r["subject_id"] for r in exchange}:
            raise ValueError("Source/exchange subject or assay row mismatch")
        result={
          "status":"NATIVE_RESEARCH_ZIP_SOURCE_COMPARISON_NO_PK_PROOF",
          "n_source_omic_assays":len(source),
          "n_independent_subjects":len({r["subject_id"] for r in source}),
          "n_source_assays_with_timepoint":sum(bool(r.get("timepoint")) for r in source),
          "n_exported_assays_with_timepoint":sum(bool(r.get("timepoint")) for r in exchange),
          "n_exported_numeric_collection_times":sum(bool(r.get("collection_time")) for r in exchange),
          "original_categorical_timepoint_field_present":"timepoint" in source[0],
          "exchange_categorical_timepoint_field_present":"timepoint" in exchange[0],
          "source_metadata_preserved_in_full_research_zip":True,
          "all_original_archive_hashes_verified":True,
          "source_records_or_ids_exported":False,
        }
        return result


if __name__=="__main__":
    print(json.dumps(audit_zip(sys.argv[1]),indent=2,sort_keys=True))
