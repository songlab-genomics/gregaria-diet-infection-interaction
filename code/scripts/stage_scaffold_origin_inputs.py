#!/usr/bin/env python3
"""Validate and copy frozen scaffold-audit inputs into a new run folder."""

import argparse
import csv
import hashlib
import shutil
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-list", required=True)
    parser.add_argument("--candidate-genes", required=True)
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--fcs-report", required=True)
    parser.add_argument("--output-list", required=True)
    parser.add_argument("--output-genes", required=True)
    parser.add_argument("--output-provenance", required=True)
    parser.add_argument("--output-fcs-report", required=True)
    parser.add_argument("--manifest", required=True)
    return parser.parse_args()


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_scaffolds(path):
    with open(path, encoding="utf-8") as handle:
        values = [
            line.strip()
            for line in handle
            if line.strip() and not line.startswith("#")
        ]
    if not values:
        raise ValueError("Candidate scaffold list is empty.")
    if len(values) != len(set(values)):
        raise ValueError("Candidate scaffold list contains duplicates.")
    invalid = [value for value in values if not value.startswith("NW_")]
    if invalid:
        raise ValueError("Candidate list contains non-NW_ accessions: " + ", ".join(invalid))
    return set(values)


def table_scaffolds(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return {
            row["scaffold"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("scaffold")
        }


def gene_scaffolds(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return {
            row["seqid"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("seqid")
        }


def main():
    args = parse_args()
    expected = read_scaffolds(args.candidate_list)
    if table_scaffolds(args.provenance) != expected:
        raise ValueError("Provenance scaffold set does not match the candidate list.")
    if gene_scaffolds(args.candidate_genes) != expected:
        raise ValueError("Candidate-gene scaffold set does not match the candidate list.")

    pairs = [
        (args.candidate_list, args.output_list, "candidate_scaffolds"),
        (args.candidate_genes, args.output_genes, "candidate_deg_genes"),
        (args.provenance, args.output_provenance, "candidate_provenance"),
        (args.fcs_report, args.output_fcs_report, "ncbi_fcs_report"),
    ]
    rows = []
    for source, destination, role in pairs:
        destination_path = Path(destination)
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        rows.append(
            {
                "input_role": role,
                "source_file": str(Path(source).resolve()),
                "staged_file": str(destination_path.resolve()),
                "size_bytes": destination_path.stat().st_size,
                "sha256": sha256(destination_path),
            }
        )

    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    with open(args.manifest, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "input_role",
                "source_file",
                "staged_file",
                "size_bytes",
                "sha256",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
