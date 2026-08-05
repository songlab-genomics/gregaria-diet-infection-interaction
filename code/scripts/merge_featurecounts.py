#!/usr/bin/env python3
"""Merge per-sample featureCounts tables into one raw-count matrix."""

import argparse
import csv
import sys
from pathlib import Path


def allow_large_annotation_fields():
    """Allow long featureCounts coordinate/annotation columns."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample-file",
        action="append",
        required=True,
        help="Repeated SAMPLE=featureCounts.txt entry.",
    )
    parser.add_argument("--strip-prefix", default="")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_featurecounts(path, strip_prefix):
    values = {}
    with open(path) as handle:
        reader = csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")
        header = next(reader)
        if len(header) < 7:
            raise ValueError(f"Unexpected featureCounts header in {path}")
        for row in reader:
            if not row:
                continue
            gene_id = row[0]
            if strip_prefix and gene_id.startswith(strip_prefix):
                gene_id = gene_id[len(strip_prefix):]
            values[gene_id] = int(float(row[-1]))
    return values


def main():
    allow_large_annotation_fields()
    args = parse_args()
    sample_paths = []
    for entry in args.sample_file:
        if "=" not in entry:
            raise ValueError(f"Expected SAMPLE=PATH, received: {entry}")
        sample, path = entry.split("=", 1)
        sample_paths.append((sample, path))

    sample_counts = {
        sample: read_featurecounts(path, args.strip_prefix)
        for sample, path in sample_paths
    }
    genes = sorted({gene for counts in sample_counts.values() for gene in counts})
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id"] + [sample for sample, _ in sample_paths])
        for gene in genes:
            writer.writerow(
                [gene] + [sample_counts[sample].get(gene, 0) for sample, _ in sample_paths]
            )


if __name__ == "__main__":
    main()
