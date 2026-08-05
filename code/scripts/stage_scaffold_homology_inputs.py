#!/usr/bin/env python3
"""Copy the frozen homology target set into a run-scoped input directory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scaffolds", required=True, type=Path)
    parser.add_argument("--genes", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--output-scaffolds", required=True, type=Path)
    parser.add_argument("--output-genes", required=True, type=Path)
    parser.add_argument("--output-provenance", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    pairs = [
        (args.scaffolds, args.output_scaffolds, "candidate_scaffolds"),
        (args.genes, args.output_genes, "candidate_genes"),
        (args.provenance, args.output_provenance, "candidate_provenance"),
    ]
    for source, target, _ in pairs:
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty frozen input: {source}")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["role", "source", "staged", "bytes", "sha256"],
            delimiter="\t",
        )
        writer.writeheader()
        for source, target, role in pairs:
            writer.writerow(
                {
                    "role": role,
                    "source": str(source.resolve()),
                    "staged": str(target),
                    "bytes": target.stat().st_size,
                    "sha256": sha256(target),
                }
            )


if __name__ == "__main__":
    main()
