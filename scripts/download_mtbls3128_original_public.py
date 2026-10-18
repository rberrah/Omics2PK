"""Download PUBLIC original MTBLS3128 ISA metadata and optional LC-MS mzML locally.

Study: Silveira et al., Front Pharmacol (2021), doi:10.3389/fphar.2021.752960.
Requires local internet access to ftp.ebi.ac.uk. No authentication, no
upload to GitHub, and no original participant data is committed by this tool.

IMPORTANT: Public mzML/ISA files do not supply a verified individual
AUC/Cmax ↔ omics subject crosswalk or the authors' five processed CSVs.
"""
from __future__ import annotations

import argparse
from ftplib import FTP, error_perm
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import tempfile

HOST = "ftp.ebi.ac.uk"
STUDY = "MTBLS3128"
ROOT = f"/pub/databases/metabolights/studies/public/{STUDY}"
METADATA_PREFIXES = ("i_", "s_", "a_", "m_")
METADATA_SUFFIXES = (".txt", ".tsv", ".csv")
MAX_DIRECTORY_DEPTH = 6
MAX_DIRECTORY_ENTRIES = 20000


class DownloadGateError(ValueError):
    """Download contract error without exposing private source measurements."""


def _safe_relative(remote: str) -> str:
    if not remote.startswith(ROOT + "/"):
        raise DownloadGateError("Source path outside official MTBLS3128 study")
    rel = remote[len(ROOT) + 1:]
    parts = PurePosixPath(rel).parts
    if not parts or any(x in ("..", ".", "") for x in parts) or "\\" in rel:
        raise DownloadGateError("Unsafe source path")
    return "/".join(parts)


def _list_folder(ftp, directory: str):
    """Prefer MLSD; fallback handles public EBI servers without MLSD."""
    try:
        entries = list(ftp.mlsd(directory))
        for name, meta in entries:
            if name in (".", "..") or "/" in name or "\\" in name:
                continue
            yield name, str(meta.get("type", "")).lower()
    except (error_perm, AttributeError, NotImplementedError):
        current = ftp.pwd()
        try:
            ftp.cwd(directory)
            filenames = ftp.nlst()
            for path in filenames:
                name = PurePosixPath(path).name
                if not name or name in (".", ".."):
                    continue
                try:
                    ftp.cwd(directory.rstrip("/") + "/" + name)
                    ftp.cwd(directory)
                    typ = "dir"
                except error_perm:
                    typ = "file"
                yield name, typ
        finally:
            ftp.cwd(current)


def source_inventory(ftp, *, include_raw: bool):
    """Enumerate authentic FTP paths; filenames are kept local, never printed."""
    metadata = []
    mzml = []
    zip_files = []
    pending = [(ROOT, 0)]
    visited = set()
    count = 0
    while pending:
        directory, depth = pending.pop()
        if directory in visited:
            continue
        visited.add(directory)
        if depth > MAX_DIRECTORY_DEPTH:
            raise DownloadGateError("Exceeded study folder depth")
        for name, kind in _list_folder(ftp, directory):
            count += 1
            if count > MAX_DIRECTORY_ENTRIES:
                raise DownloadGateError("Too many entries in source folder")
            remote = directory.rstrip("/") + "/" + name
            rel = _safe_relative(remote)
            if kind == "dir":
                if name == "FILES" or directory != ROOT:
                    if include_raw:
                        pending.append((remote, depth + 1))
            elif kind in ("file", ""):
                if directory == ROOT and name.startswith(METADATA_PREFIXES) and name.lower().endswith(METADATA_SUFFIXES):
                    metadata.append(remote)
                elif include_raw and rel.startswith("FILES/"):
                    if rel.lower().endswith(".mzml"):
                        mzml.append(remote)
                    elif rel.lower().endswith(".zip"):
                        zip_files.append(remote)
    if len(metadata) < 4:
        raise DownloadGateError("Expected original ISA root metadata files unavailable")
    return sorted(metadata), sorted(mzml), sorted(zip_files)


def select_files(metadata, mzml, zip_files, *, mzml_max_files: int, zip_max_files: int):
    """Select by source format only, never by inferred patient/sample numbers."""
    if mzml_max_files < 0 or zip_max_files < 0:
        raise DownloadGateError("File limits must be nonnegative")
    if mzml_max_files > len(mzml) or zip_max_files > len(zip_files):
        raise DownloadGateError("Requested more original files than source inventory contains")
    return [*metadata, *mzml[:mzml_max_files], *zip_files[:zip_max_files]]


def _download_one(ftp, remote: str, root: Path, *, max_file_bytes: int, overwrite: bool):
    rel = _safe_relative(remote)
    dest = root.joinpath(*PurePosixPath(rel).parts)
    if dest.exists() and not overwrite:
        raise DownloadGateError("Local source file exists: use --overwrite to verify a fresh remote download")
    dest.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    digest = sha256()
    temp_name = None
    try:
        with tempfile.NamedTemporaryFile(dir=dest.parent, prefix=".partial-", delete=False) as handle:
            temp_name = handle.name

            def receive(chunk):
                nonlocal total
                total += len(chunk)
                if total > max_file_bytes:
                    raise DownloadGateError("Remote file exceeded explicit local maximum bytes")
                handle.write(chunk)
                digest.update(chunk)

            ftp.retrbinary("RETR " + remote, receive, blocksize=131072)
        os.replace(temp_name, dest)
    except Exception:
        if temp_name is not None:
            Path(temp_name).unlink(missing_ok=True)
        raise
    return {
        "relative_original_filename": rel,
        "downloaded_bytes": total,
        "original_download_sha256": digest.hexdigest(),
        "local_storage_only": True,
    }


def download_study(
    ftp,
    *,
    output_dir: Path,
    mzml_max_files: int = 0,
    zip_max_files: int = 0,
    max_total_bytes: int = 10 * 1024**3,
    max_file_bytes: int = 3 * 1024**3,
    overwrite: bool = False,
    list_only: bool = False,
):
    if min(max_total_bytes, max_file_bytes) <= 0:
        raise DownloadGateError("Explicit positive disk/download limits required")
    metadata, mzml, zip_files = source_inventory(
        ftp, include_raw=bool(mzml_max_files or zip_max_files or list_only)
    )
    selected = select_files(
        metadata, mzml, zip_files,
        mzml_max_files=mzml_max_files, zip_max_files=zip_max_files,
    )
    summary = {
        "study": STUDY,
        "original_article_doi": "10.3389/fphar.2021.752960",
        "original_public_isa_files": len(metadata),
        "original_public_mzml_files_listed": len(mzml) if (mzml_max_files or zip_max_files or list_only) else None,
        "original_public_zip_files_listed": len(zip_files) if (mzml_max_files or zip_max_files or list_only) else None,
        "selected_files": len(selected),
        "selected_mzml_files": mzml_max_files,
        "selected_zip_files": zip_max_files,
        "source": "official anonymous EBI FTP; local-only original retrieval",
        "no_proof_of_individual_AUC_Cmax_crosswalk": True,
        "no_individual_patient_rows_exported_to_stdout": True,
    }
    if list_only:
        return summary | {"status": "SOURCE_INVENTORY_ONLY_NO_DATA_DOWNLOADED"}
    if len(selected) == 0:
        raise DownloadGateError("No original public files selected")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    downloaded = 0
    for idx, remote in enumerate(selected, 1):
        item = _download_one(
            ftp, remote, output_dir,
            max_file_bytes=min(max_file_bytes, max_total_bytes - downloaded),
            overwrite=overwrite,
        )
        downloaded += item["downloaded_bytes"]
        if downloaded > max_total_bytes:
            raise DownloadGateError("Remote source exceeded global download quota")
        manifest.append(item)
        print(f"Original public source {idx}/{len(selected)}: verified local SHA256", flush=True)
    data = summary | {
        "status": "PUBLIC_ISA_AND_OPTIONAL_RAW_MS_DOWNLOADED_LOCALLY",
        "downloaded_file_count": len(manifest),
        "downloaded_bytes_total": downloaded,
        "original_source_manifest": manifest,
        "clinical_E3_results_qualified": False,
    }
    with (output_dir / "private_source_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
    return summary | {
        "status": data["status"],
        "downloaded_file_count": len(manifest),
        "downloaded_bytes_total": downloaded,
        "manifest_stored_locally_only": True,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, default=Path("private/MTBLS3128"))
    p.add_argument("--list-only", action="store_true")
    p.add_argument("--mzml-max-files", type=int, default=0,
                   help="Download first N source-sorted mzML; 0=none (default)")
    p.add_argument("--zip-max-files", type=int, default=0,
                   help="Download first N source-sorted original ZIP; 0=none (may be large vendor data)")
    p.add_argument("--max-total-gib", type=float, default=10)
    p.add_argument("--max-file-gib", type=float, default=3)
    p.add_argument("--overwrite", action="store_true",
                   help="Replace existing local files with independently re-downloaded source")
    args = p.parse_args(argv)
    try:
        with FTP(HOST, timeout=40) as ftp:
            ftp.login()
            ftp.voidcmd("TYPE I")
            result = download_study(
                ftp,
                output_dir=args.output_dir,
                mzml_max_files=args.mzml_max_files,
                zip_max_files=args.zip_max_files,
                max_total_bytes=int(args.max_total_gib * 1024**3),
                max_file_bytes=int(args.max_file_gib * 1024**3),
                overwrite=args.overwrite,
                list_only=args.list_only,
            )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (DownloadGateError, OSError, EOFError, error_perm, ValueError) as err:
        p.exit(2, f"ABSTAIN: {type(err).__name__} — no original patient-level data published\n")


if __name__ == "__main__":
    main()
