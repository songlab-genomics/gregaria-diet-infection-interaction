#!/usr/bin/env python3
"""Prefix and combine host and fungal FASTA/GTF references.

Prefixing sequence and feature identifiers prevents accidental collisions and
makes every competitive alignment explicitly traceable to host or fungus.
"""

import argparse
import csv
import re
from pathlib import Path


GTF_ID_KEYS = ("gene_id", "transcript_id", "exon_id", "protein_id")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host-fasta", required=True)
    parser.add_argument("--host-gtf", required=True)
    parser.add_argument("--fungus-fasta", required=True)
    parser.add_argument("--fungus-gtf", required=True)
    parser.add_argument("--host-prefix", default="HOST__")
    parser.add_argument("--fungus-prefix", default="MR__")
    parser.add_argument("--combined-fasta", required=True)
    parser.add_argument("--combined-gtf", required=True)
    parser.add_argument("--host-prefixed-gtf", required=True)
    parser.add_argument("--fungus-prefixed-gtf", required=True)
    parser.add_argument("--sequence-map", required=True)
    parser.add_argument("--gene-map", required=True)
    return parser.parse_args()


def rewrite_fasta(source, destination, prefix, organism, sequence_rows, mode):
    with open(source) as src, open(destination, mode) as dst:
        for line in src:
            if not line.startswith(">"):
                dst.write(line)
                continue
            header = line[1:].rstrip("\n")
            original = header.split(maxsplit=1)[0]
            description = header[len(original):]
            prefixed = prefix + original
            dst.write(f">{prefixed}{description}\n")
            sequence_rows.append(
                {
                    "organism": organism,
                    "original_seqid": original,
                    "prefixed_seqid": prefixed,
                }
            )


def prefix_attributes(attributes, prefix, organism, gene_rows):
    seen_gene = None
    rewritten = attributes
    for key in GTF_ID_KEYS:
        pattern = re.compile(rf'({re.escape(key)}\s+")([^"]+)(")')

        def replace(match):
            nonlocal seen_gene
            original = match.group(2)
            if key == "gene_id":
                seen_gene = original
            return match.group(1) + prefix + original + match.group(3)

        rewritten = pattern.sub(replace, rewritten)

    if seen_gene is not None:
        gene_rows[(organism, seen_gene)] = {
            "organism": organism,
            "original_gene_id": seen_gene,
            "prefixed_gene_id": prefix + seen_gene,
        }
    return rewritten


def rewrite_gtf(source, destination, prefix, organism, gene_rows):
    with open(source) as src, open(destination, "w") as dst:
        for line in src:
            if line.startswith("#") or not line.strip():
                dst.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected 9-column GTF line in {source}: {line[:120]}")
            fields[0] = prefix + fields[0]
            fields[8] = prefix_attributes(fields[8], prefix, organism, gene_rows)
            dst.write("\t".join(fields) + "\n")


def concatenate_text(inputs, output):
    with open(output, "w") as dst:
        for path in inputs:
            with open(path) as src:
                for line in src:
                    dst.write(line)


def main():
    args = parse_args()
    output_paths = [
        args.combined_fasta,
        args.combined_gtf,
        args.host_prefixed_gtf,
        args.fungus_prefixed_gtf,
        args.sequence_map,
        args.gene_map,
    ]
    for path in output_paths:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    sequence_rows = []
    gene_rows = {}

    rewrite_fasta(
        args.host_fasta,
        args.combined_fasta,
        args.host_prefix,
        "Schistocerca gregaria",
        sequence_rows,
        "w",
    )
    rewrite_fasta(
        args.fungus_fasta,
        args.combined_fasta,
        args.fungus_prefix,
        "Metarhizium robertsii reference proxy",
        sequence_rows,
        "a",
    )
    rewrite_gtf(
        args.host_gtf,
        args.host_prefixed_gtf,
        args.host_prefix,
        "Schistocerca gregaria",
        gene_rows,
    )
    rewrite_gtf(
        args.fungus_gtf,
        args.fungus_prefixed_gtf,
        args.fungus_prefix,
        "Metarhizium robertsii reference proxy",
        gene_rows,
    )
    concatenate_text(
        [args.host_prefixed_gtf, args.fungus_prefixed_gtf],
        args.combined_gtf,
    )

    with open(args.sequence_map, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["organism", "original_seqid", "prefixed_seqid"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(sequence_rows)

    with open(args.gene_map, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["organism", "original_gene_id", "prefixed_gene_id"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(sorted(gene_rows.values(), key=lambda row: row["prefixed_gene_id"]))


if __name__ == "__main__":
    main()
