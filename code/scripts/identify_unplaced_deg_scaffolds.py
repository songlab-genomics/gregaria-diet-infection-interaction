#!/usr/bin/env python3
"""Link fresh host-only DEGs to unplaced S. gregaria scaffolds in the GTF."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


GENE_ID_PATTERN = re.compile(r'gene_id "([^"]+)"')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deseq-results", required=True)
    parser.add_argument("--mapping-reconciliation", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output-genes", required=True)
    parser.add_argument("--output-scaffolds", required=True)
    return parser.parse_args()


def as_bool(value):
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def main():
    args = parse_args()
    significant_contrasts = defaultdict(set)
    with open(args.deseq_results, newline="") as handle:
        for row in csv.DictReader(handle):
            if as_bool(row.get("significant", False)):
                significant_contrasts[row["gene_id"]].add(row["contrast_id"])

    mapping_classes = defaultdict(set)
    with open(args.mapping_reconciliation, newline="") as handle:
        for row in csv.DictReader(handle):
            if row["gene_id"] in significant_contrasts:
                mapping_classes[row["gene_id"]].add(
                    row["mapping_support_class"]
                )

    gene_to_scaffolds = defaultdict(set)
    with open(args.gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "transcript", "exon"}:
                continue
            match = GENE_ID_PATTERN.search(fields[8])
            if not match:
                continue
            gene_id = match.group(1)
            if gene_id in significant_contrasts and fields[0].startswith("NW_"):
                gene_to_scaffolds[gene_id].add(fields[0])

    rows = []
    for gene_id in sorted(gene_to_scaffolds):
        rows.append(
            {
                "gene_id": gene_id,
                "unplaced_scaffolds": ",".join(sorted(gene_to_scaffolds[gene_id])),
                "significant_contrasts": ",".join(
                    sorted(significant_contrasts[gene_id])
                ),
                "n_significant_contrasts": len(significant_contrasts[gene_id]),
                "mapping_support_classes": ",".join(
                    sorted(mapping_classes[gene_id])
                ),
            }
        )

    Path(args.output_genes).parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "gene_id",
        "unplaced_scaffolds",
        "significant_contrasts",
        "n_significant_contrasts",
        "mapping_support_classes",
    ]
    with open(args.output_genes, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    scaffolds = sorted(
        {scaffold for values in gene_to_scaffolds.values() for scaffold in values}
    )
    Path(args.output_scaffolds).write_text(
        "".join(f"{scaffold}\n" for scaffold in scaffolds)
    )


if __name__ == "__main__":
    main()
