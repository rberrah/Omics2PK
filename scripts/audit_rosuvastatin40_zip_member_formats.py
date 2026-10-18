"""Bounded original EBI MTBLS3128 ZIP central-directory audit (no mzML download).

Inspects SOURCE ARCHIVE MEMBER EXTENSIONS via HTTP range requests only, when
remote FTP HTTPS mirrors honor byte ranges. This checks for processed CSV or
quantitative PK tables INSIDE ZIPs without downloading large raw experiments.

No person/sample/ZIP/member filename, original measured value or individual
identifier is printed/committed. An uninspected ZIP is never classified empty.
"""
from __future__ import annotations

from collections import Counter
from ftplib import FTP, error_perm
from pathlib import PurePosixPath
from statistics import median
from urllib.parse import quote
from urllib.request import Request, urlopen
import json
import re
import struct

from scripts.audit_rosuvastatin40_public_ftp_census import HOST, ROOT

TAIL_BYTES = 65557
MAX_CENTRAL_BYTES = 2_000_000
MAX_ARCHIVES = 200
MAX_ENTRIES = 20000
ZIP_CENTRAL_SIGNATURE = b"PK\x01\x02"
ZIP_END_SIGNATURE = b"PK\x05\x06"
TARGET_EXT = frozenset({".csv",".tsv",".xls",".xlsx",".parquet",".rds",".h5",".hdf5"})

class SourceZipError(ValueError):
    pass


def archive_member_extensions(central_directory: bytes, *, expected_count: int) -> dict:
    """Read ZIP central directory, never member content; report extensions only."""
    if expected_count <= 0 or expected_count > MAX_ENTRIES:
        raise SourceZipError("Unsafe source archive directory cardinality")
    counter = Counter()
    target = 0
    at = 0
    for _ in range(expected_count):
        if central_directory[at:at+4] != ZIP_CENTRAL_SIGNATURE or at + 46 > len(central_directory):
            raise SourceZipError("Malformed ZIP central directory or unsupported split archive")
        n, extra, comment = struct.unpack_from("<HHH", central_directory, at+28)
        end = at + 46 + n + extra + comment
        if end > len(central_directory) or n > 4096:
            raise SourceZipError("Unsafe ZIP filename or central entry length")
        name = central_directory[at+46:at+46+n].decode("utf-8","replace")
        suffix = PurePosixPath(name).suffix.lower() or "no-extension"
        if not name.endswith("/"):
            counter[suffix] += 1
            if suffix in TARGET_EXT: target += 1
        at = end
    return {"n_original_archive_members":sum(counter.values()),
            "original_member_extension_counts":dict(sorted(counter.items())),
            "n_candidate_tabular_data_members":target,
            "no_original_member_filenames_or_values_exported":True}


def zip_end_metadata(tail: bytes, *, file_size: int):
    pos = tail.rfind(ZIP_END_SIGNATURE)
    if pos < 0 or len(tail) - pos < 22:
        raise SourceZipError("ZIP end marker not in bounded source tail")
    (disk, start_disk, disk_entries, entries, central_bytes, central_offset, comment) = struct.unpack_from(
        "<HHHHIIH",tail,pos+4
    )
    if disk or start_disk or disk_entries != entries or entries in (0, 65535) or (
        central_bytes == 0xffffffff or central_offset == 0xffffffff
    ):
        raise SourceZipError("ZIP64 or split archive not supported by bounded reader")
    if central_bytes <= 0 or central_bytes > MAX_CENTRAL_BYTES:
        raise SourceZipError("Central directory too large")
    if central_offset + central_bytes > file_size or pos+22+comment > len(tail):
        raise SourceZipError("Invalid ZIP footer offsets")
    return central_offset, central_bytes, entries


def _range(url: str, spec: str, *, maximum: int) -> bytes:
    with urlopen(Request(url, headers={
        "User-Agent": "Omics2PK-EBI-original-ZIP-metadata-only/1.0",
        "Range": "bytes=" + spec,
        "Accept-Encoding": "identity"
    }),timeout=25) as connection:
        if connection.status != 206:
            # Stop BEFORE reading an unexpectedly full multi-GB raw archive.
            raise SourceZipError("Remote HTTPS server does not support bounded byte ranges")
        data = connection.read(maximum + 1)
    if len(data) > maximum:
        raise SourceZipError("Remote archive exceeded range response byte cap")
    return data


def probe_zip_central(url: str, size: int) -> dict:
    if not 50 <= size <= 30 * 1024**3:
        raise SourceZipError("Invalid original ZIP file size")
    tail_length = min(size, TAIL_BYTES)
    tail = _range(url, f"-{tail_length}", maximum=TAIL_BYTES)
    if len(tail) != tail_length:
        raise SourceZipError("Incomplete bounded archive footer response")
    offset, length, entries = zip_end_metadata(tail, file_size=size)
    central = _range(url, f"{offset}-{offset+length-1}", maximum=MAX_CENTRAL_BYTES)
    if len(central) != length:
        raise SourceZipError("Incomplete central directory range response")
    return archive_member_extensions(central, expected_count=entries)


def list_original_zip_sizes():
    ftp = FTP(timeout=35)
    archives = []
    folders = [ROOT]
    inaccessible_subfolders = 0
    try:
        # Match the source-verified direct FTP sequence used by the public census.
        ftp.connect(HOST, 21, timeout=35)
        ftp.login()
        ftp.voidcmd("TYPE I")
        while folders:
            folder=folders.pop()
            # An inaccessible nested vendor folder is not evidence of absence.
            try:
                ftp.cwd(folder)
                try:
                    entries = list(ftp.mlsd())
                except (error_perm, NotImplementedError):
                    # Original EBI server may reject MLSD; verified census uses NLST.
                    entries = []
                    for listed in ftp.nlst():
                        name = PurePosixPath(listed).name
                        if not name or name in (".", ".."):
                            continue
                        full = folder.rstrip("/") + "/" + name
                        try:
                            ftp.cwd(full)
                            ftp.cwd(folder)
                            kind = "dir"
                        except error_perm:
                            kind = "file"
                        size = 0
                        if kind == "file" and name.lower().endswith(".zip"):
                            try:
                                size = ftp.size(full) or 0
                            except error_perm:
                                pass
                        entries.append((name, {"type": kind, "size": str(size)}))
            except error_perm:
                if folder == ROOT:
                    raise
                inaccessible_subfolders += 1
                continue
            for name, f in entries:
                if name in (".","..") or "/" in name or "\\" in name:
                    continue
                source=folder.rstrip("/")+"/"+name
                if not source.startswith(ROOT+"/"):
                    continue
                typ=f.get("type","").lower()
                if typ=="dir":
                    folders.append(source)
                elif typ=="file" and name.lower().endswith(".zip"):
                    try: size=int(f.get("size",0))
                    except (ValueError,TypeError):size=0
                    archives.append((source,size))
                if len(archives)>20000:
                    raise SourceZipError("Original FTP inventory too large")
    finally:
        try:ftp.quit()
        except Exception:ftp.close()
    return archives, inaccessible_subfolders


def inspect_ebi_archive_sample(*, sample_size=8) -> dict:
    if sample_size is not None and not 0 <= sample_size <= MAX_ARCHIVES:
        raise SourceZipError("Invalid maximum number of inspected ZIP archives")
    archives, inaccessible_subfolders = list_original_zip_sizes()
    if sample_size is None:
        sample_size = min(MAX_ARCHIVES, len(archives))
    sizes=sorted(size for _,size in archives if size>0)
    selection=sorted((path,size) for path,size in archives if size>0)
    # Cover small and large archive sizes, not filename-derived patient labels.
    selection.sort(key=lambda x:x[1])
    if selection and sample_size:
        indices={round(i*(len(selection)-1)/max(1,sample_size-1)) for i in range(sample_size)}
        chosen=[selection[i] for i in sorted(indices)]
    else:
        chosen=[]
    extensions=Counter()
    n_members=0
    tabular=0
    success=0
    failed=Counter()
    for path,size in chosen:
        url="https://"+HOST+quote(path, safe="/")
        try:
            result=probe_zip_central(url,size)
            success+=1
            n_members+=result["n_original_archive_members"]
            tabular+=result["n_candidate_tabular_data_members"]
            extensions.update(result["original_member_extension_counts"])
        except (SourceZipError,OSError,ValueError) as error:
            failed[type(error).__name__]+=1
    return {
        "study":"MTBLS3128",
        "source":"Official EBI public FTP archive inventory + bounded HTTPS central directory probes",
        "n_original_zip_archives_listed":len(archives),
        "n_inaccessible_original_nested_directories":inaccessible_subfolders,
        "is_complete_original_public_zip_directory_census":inaccessible_subfolders == 0,
        "n_original_zip_sizes_reported":len(sizes),
        "min_zip_size_bytes":min(sizes) if sizes else None,
        "median_zip_size_bytes":int(median(sizes)) if sizes else None,
        "max_zip_size_bytes":max(sizes) if sizes else None,
        "n_archives_within_10MiB":sum(n <= 10*1024**2 for n in sizes),
        "n_zip_archives_selected_for_bounded_probe":len(chosen),
        "n_zip_archives_central_directories_actually_inspected":success,
        "n_member_files_in_inspected_archives":n_members,
        "member_extension_counts_in_inspected_archives":dict(sorted(extensions.items())),
        "n_candidate_processed_tables_in_inspected_archives":tabular,
        "failed_probe_categories":dict(failed),
        "not_all_archives_proven_if_subset":success < len(archives),
        "still_no_40_human_individual_AUC_Cmax_crosswalk":True,
        "no_source_filenames_or_subject_values_logged":True,
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="Attempt each official ZIP central-directory (metadata only; no member content)")
    args = parser.parse_args()
    try:
        print(json.dumps(inspect_ebi_archive_sample(sample_size=None if args.all else 8),indent=2,sort_keys=True))
    except (SourceZipError,OSError,ValueError,error_perm) as error:
        print(json.dumps({"status":"EBI_ARCHIVE_INVENTORY_UNAVAILABLE_ABSTAIN",
                          "error_type":type(error).__name__}))
