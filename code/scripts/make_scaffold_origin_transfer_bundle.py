#!/usr/bin/env python3
"""Create a checksummed archive of scaffold-origin audit outputs."""

import argparse
import csv
import hashlib
import tarfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--archive", required=True)
    return parser.parse_args()


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    files = [Path(value).resolve() for value in args.input]
    rows = []
    for source in files:
        if not source.is_file() or source.stat().st_size == 0:
            raise ValueError(f"Missing or empty transfer input: {source}")
        rows.append(
            {
                "source_file": str(source),
                "archive_path": str(source.relative_to(run_dir)),
                "size_bytes": source.stat().st_size,
                "sha256": sha256(source),
            }
        )

    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    with open(args.manifest, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source_file", "archive_path", "size_bytes", "sha256"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    with tarfile.open(args.archive, "w:gz") as archive:
        for row in rows:
            archive.add(row["source_file"], arcname=row["archive_path"])
        archive.add(args.manifest, arcname="transfer_manifest.tsv")


if __name__ == "__main__":
    main()
