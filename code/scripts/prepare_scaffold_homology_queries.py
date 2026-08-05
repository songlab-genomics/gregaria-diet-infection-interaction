#!/usr/bin/env python3
"""Extract one genomic query per candidate gene and all encoded proteins."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-genes", required=True, type=Path)
    parser.add_argument("--candidate-scaffolds-fasta", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--protein-fasta", required=True, type=Path)
    parser.add_argument("--output-nucleotide", required=True, type=Path)
    parser.add_argument("--output-protein", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--output-missing", required=True, type=Path)
    return parser.parse_args()


def parse_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in raw.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attributes[key] = unquote(value)
    return attributes


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current = ""
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                sequences[current] = []
            elif current:
                sequences[current].append(line.strip())
    return {key: "".join(value).upper() for key, value in sequences.items()}


def wrap(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def stream_fasta(path: Path):
    header = ""
    sequence: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence)
                header = line[1:].rstrip()
                sequence = []
            else:
                sequence.append(line.strip())
        if header:
            yield header, "".join(sequence)


def main() -> None:
    args = parse_args()
    with args.candidate_genes.open(newline="") as handle:
        candidate_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not candidate_rows:
        raise ValueError("Candidate gene table is empty.")

    candidate_by_id = {row["gene_id"]: row for row in candidate_rows}
    candidate_ids = set(candidate_by_id)
    gene_features: dict[str, tuple[str, int, int, str]] = {}
    proteins_by_gene: dict[str, set[str]] = defaultdict(set)

    with args.annotation.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, feature, start, end, strand, raw_attrs = (
                fields[0],
                fields[2],
                fields[3],
                fields[4],
                fields[6],
                fields[8],
            )
            attrs = parse_attributes(raw_attrs)
            gene_id = attrs.get("gene") or attrs.get("Name") or attrs.get("gene_id")
            if gene_id not in candidate_ids:
                continue
            if feature == "gene":
                gene_features[gene_id] = (seqid, int(start), int(end), strand)
            elif feature == "CDS" and attrs.get("protein_id"):
                proteins_by_gene[gene_id].add(attrs["protein_id"])

    scaffolds = read_fasta(args.candidate_scaffolds_fasta)
    missing_feature = sorted(candidate_ids - set(gene_features))
    if missing_feature:
        raise ValueError(
            f"{len(missing_feature)} candidate genes lack gene features; "
            f"first: {missing_feature[:5]}"
        )

    for path in (
        args.output_nucleotide,
        args.output_protein,
        args.output_manifest,
        args.output_missing,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, str]] = []
    with args.output_nucleotide.open("w") as nucleotide_handle:
        for gene_id in sorted(candidate_ids):
            seqid, start, end, strand = gene_features[gene_id]
            if seqid not in scaffolds:
                raise ValueError(f"Candidate scaffold sequence missing: {seqid}")
            sequence = scaffolds[seqid][start - 1 : end]
            if strand == "-":
                sequence = reverse_complement(sequence)
            query_id = f"{gene_id}|{seqid}:{start}-{end}|{strand}"
            nucleotide_handle.write(f">{query_id}\n{wrap(sequence)}\n")
            row = candidate_by_id[gene_id]
            manifest_rows.append(
                {
                    "query_id": query_id,
                    "gene_id": gene_id,
                    "seqid": seqid,
                    "gene_biotype": row.get("gene_biotype", ""),
                    "description": row.get("description", ""),
                    "query_type": "gene_genomic_DNA",
                    "source_accession": seqid,
                    "start": str(start),
                    "end": str(end),
                    "strand": strand,
                }
            )

    wanted_proteins = {
        protein_id: gene_id
        for gene_id, protein_ids in proteins_by_gene.items()
        for protein_id in protein_ids
    }
    found_proteins: set[str] = set()
    with args.output_protein.open("w") as protein_handle:
        for header, sequence in stream_fasta(args.protein_fasta):
            accession = header.split()[0]
            if accession not in wanted_proteins:
                continue
            gene_id = wanted_proteins[accession]
            found_proteins.add(accession)
            protein_handle.write(f">{accession}\n{wrap(sequence)}\n")
            row = candidate_by_id[gene_id]
            manifest_rows.append(
                {
                    "query_id": accession,
                    "gene_id": gene_id,
                    "seqid": row["seqid"],
                    "gene_biotype": row.get("gene_biotype", ""),
                    "description": row.get("description", ""),
                    "query_type": "protein",
                    "source_accession": accession,
                    "start": "",
                    "end": "",
                    "strand": "",
                }
            )

    manifest_fields = [
        "query_id",
        "gene_id",
        "seqid",
        "gene_biotype",
        "description",
        "query_type",
        "source_accession",
        "start",
        "end",
        "strand",
    ]
    with args.output_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    missing_rows = []
    for gene_id in sorted(candidate_ids):
        expected = proteins_by_gene.get(gene_id, set())
        missing = sorted(expected - found_proteins)
        row = candidate_by_id[gene_id]
        missing_rows.append(
            {
                "gene_id": gene_id,
                "seqid": row["seqid"],
                "gene_biotype": row.get("gene_biotype", ""),
                "protein_accessions_expected": ";".join(sorted(expected)),
                "protein_accessions_missing": ";".join(missing),
                "n_proteins_found": str(len(expected & found_proteins)),
            }
        )
    with args.output_missing.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(missing_rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(missing_rows)

    protein_gene_count = len(
        {
            wanted_proteins[protein_id]
            for protein_id in found_proteins
        }
    )
    protein_coding_genes = {
        gene_id
        for gene_id, row in candidate_by_id.items()
        if row.get("gene_biotype") == "protein_coding"
    }
    minimum_expected = int(0.9 * len(protein_coding_genes))
    if protein_gene_count < minimum_expected:
        raise ValueError(
            "Protein FASTA/GFF mismatch: proteins were recovered for "
            f"{protein_gene_count} of {len(protein_coding_genes)} "
            "protein-coding candidate genes."
        )
    print(
        f"Prepared {len(candidate_ids)} nucleotide gene queries and "
        f"{len(found_proteins)} protein queries from {protein_gene_count} genes."
    )


if __name__ == "__main__":
    main()
