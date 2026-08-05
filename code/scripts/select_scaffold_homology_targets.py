#!/usr/bin/env python3
"""Freeze a contrast-specific unresolved scaffold/gene set for homology tests."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--candidate-genes", required=True, type=Path)
    parser.add_argument("--contrast", required=True)
    parser.add_argument(
        "--evidence-status",
        default="Unresolved by Kraken/FCS",
    )
    parser.add_argument("--output-scaffolds", required=True, type=Path)
    parser.add_argument("--output-genes", required=True, type=Path)
    parser.add_argument("--output-provenance", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def contrast_memberships(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def main() -> None:
    args = parse_args()
    evidence_fields, evidence_rows = read_tsv(args.evidence)
    gene_fields, gene_rows = read_tsv(args.candidate_genes)

    required_evidence = {"scaffold", "preliminary_evidence_status"}
    required_genes = {"gene_id", "seqid", "contrasts"}
    if not required_evidence.issubset(evidence_fields):
        raise ValueError(f"Evidence table is missing {required_evidence - set(evidence_fields)}")
    if not required_genes.issubset(gene_fields):
        raise ValueError(f"Candidate gene table is missing {required_genes - set(gene_fields)}")

    unresolved_scaffolds = {
        row["scaffold"]
        for row in evidence_rows
        if row["preliminary_evidence_status"] == args.evidence_status
    }
    selected_genes = [
        row
        for row in gene_rows
        if row["seqid"] in unresolved_scaffolds
        and args.contrast in contrast_memberships(row["contrasts"])
    ]
    selected_genes.sort(key=lambda row: (row["seqid"], row["gene_id"]))
    selected_scaffolds = sorted({row["seqid"] for row in selected_genes})

    if not selected_genes:
        raise ValueError("The requested evidence/contrast filters selected no genes.")

    for path in (
        args.output_scaffolds,
        args.output_genes,
        args.output_provenance,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)

    args.output_scaffolds.write_text(
        "".join(f"{scaffold}\n" for scaffold in selected_scaffolds)
    )
    with args.output_genes.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=gene_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected_genes)

    provenance_fields = [
        "contrast",
        "evidence_status",
        "selected_scaffolds",
        "selected_genes",
        "source_evidence",
        "source_evidence_sha256",
        "source_candidate_genes",
        "source_candidate_genes_sha256",
    ]
    provenance_row = {
        "contrast": args.contrast,
        "evidence_status": args.evidence_status,
        "selected_scaffolds": str(len(selected_scaffolds)),
        "selected_genes": str(len(selected_genes)),
        "source_evidence": str(args.evidence),
        "source_evidence_sha256": sha256(args.evidence),
        "source_candidate_genes": str(args.candidate_genes),
        "source_candidate_genes_sha256": sha256(args.candidate_genes),
    }
    with args.output_provenance.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=provenance_fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(provenance_row)

    print(
        f"Selected {len(selected_genes)} genes on {len(selected_scaffolds)} scaffolds "
        f"for {args.contrast}."
    )


if __name__ == "__main__":
    main()
