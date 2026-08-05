#!/usr/bin/env python3
"""Join host DEG mapping support, count retention, and scaffold taxonomy."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deg-reconciliation", required=True)
    parser.add_argument("--gene-count-comparison", required=True)
    parser.add_argument("--unplaced-genes", required=True)
    parser.add_argument("--scaffold-taxonomy", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_keyed(path, key, delimiter="\t"):
    with open(path, newline="") as handle:
        return {
            row[key]: row
            for row in csv.DictReader(handle, delimiter=delimiter)
        }


def review_priority(mapping_class, is_unplaced):
    if is_unplaced and mapping_class in {
        "host_only_only",
        "shared_opposite_direction",
    }:
        return "high: unplaced and mapping-sensitive"
    if is_unplaced:
        return "review: unplaced scaffold"
    if mapping_class in {"host_only_only", "shared_opposite_direction"}:
        return "review: mapping-sensitive"
    return "standard host DEG evidence"


def main():
    args = parse_args()
    count_rows = read_keyed(args.gene_count_comparison, "gene_id")
    unplaced_rows = read_keyed(args.unplaced_genes, "gene_id")
    taxonomy_rows = read_keyed(args.scaffold_taxonomy, "scaffold")

    output_rows = []
    with open(args.deg_reconciliation, newline="") as handle:
        for row in csv.DictReader(handle):
            gene_id = row["gene_id"]
            count_row = count_rows.get(gene_id, {})
            unplaced_row = unplaced_rows.get(gene_id, {})
            scaffolds = [
                value
                for value in unplaced_row.get("unplaced_scaffolds", "").split(",")
                if value
            ]
            taxa = sorted(
                {
                    taxonomy_rows.get(scaffold, {}).get(
                        "kraken_taxon_name", "unclassified"
                    )
                    for scaffold in scaffolds
                }
            )
            mapping_class = row["mapping_support_class"]
            is_unplaced = bool(scaffolds)
            output_rows.append(
                {
                    **row,
                    "competitive_host_retained_fraction": count_row.get(
                        "competitive_host_retained_fraction", ""
                    ),
                    "host_counts_not_retained": count_row.get(
                        "host_counts_not_retained", ""
                    ),
                    "genome_placement_class": (
                        "unplaced_NW_scaffold" if is_unplaced else "not_unplaced_NW"
                    ),
                    "unplaced_scaffolds": ",".join(scaffolds),
                    "scaffold_kraken_taxa": ",".join(taxa),
                    "origin_review_priority": review_priority(
                        mapping_class, is_unplaced
                    ),
                }
            )

    fields = [
        "contrast_id",
        "contrast_family",
        "gene_id",
        "host_only_log2FoldChange",
        "competitive_log2FoldChange",
        "host_only_padj",
        "competitive_padj",
        "host_only_significant",
        "competitive_significant",
        "mapping_support_class",
        "competitive_host_retained_fraction",
        "host_counts_not_retained",
        "genome_placement_class",
        "unplaced_scaffolds",
        "scaffold_kraken_taxa",
        "origin_review_priority",
    ]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
