#!/usr/bin/env python3
"""Compare host-only and competitive featureCounts matrices."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host-only", required=True)
    parser.add_argument("--competitive-host", required=True)
    parser.add_argument("--competitive-fungus", required=True)
    parser.add_argument("--sample-output", required=True)
    parser.add_argument("--gene-output", required=True)
    return parser.parse_args()


def read_matrix(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "gene_id":
            raise ValueError(f"Expected gene_id as the first column in {path}")
        samples = reader.fieldnames[1:]
        counts = {}
        for row in reader:
            counts[row["gene_id"]] = {
                sample: int(float(row[sample])) for sample in samples
            }
    return samples, counts


def safe_fraction(numerator, denominator):
    return numerator / denominator if denominator else ""


def main():
    args = parse_args()
    host_samples, host_only = read_matrix(args.host_only)
    competitive_samples, competitive_host = read_matrix(args.competitive_host)
    fungus_samples, competitive_fungus = read_matrix(args.competitive_fungus)

    if host_samples != competitive_samples or host_samples != fungus_samples:
        raise ValueError("The three count matrices do not have identical sample columns.")

    sample_rows = []
    for sample in host_samples:
        host_only_total = sum(row[sample] for row in host_only.values())
        competitive_host_total = sum(
            row[sample] for row in competitive_host.values()
        )
        competitive_fungus_total = sum(
            row[sample] for row in competitive_fungus.values()
        )
        sample_rows.append(
            {
                "sample_id": sample,
                "host_only_assigned_counts": host_only_total,
                "competitive_host_assigned_counts": competitive_host_total,
                "competitive_fungus_assigned_counts": competitive_fungus_total,
                "competitive_host_retained_fraction": safe_fraction(
                    competitive_host_total, host_only_total
                ),
                "host_counts_not_retained": host_only_total - competitive_host_total,
            }
        )

    gene_rows = []
    all_host_genes = sorted(set(host_only) | set(competitive_host))
    for gene_id in all_host_genes:
        host_values = host_only.get(gene_id, {})
        competitive_values = competitive_host.get(gene_id, {})
        host_total = sum(host_values.get(sample, 0) for sample in host_samples)
        competitive_total = sum(
            competitive_values.get(sample, 0) for sample in host_samples
        )
        gene_rows.append(
            {
                "gene_id": gene_id,
                "host_only_total_counts": host_total,
                "competitive_host_total_counts": competitive_total,
                "competitive_host_retained_fraction": safe_fraction(
                    competitive_total, host_total
                ),
                "host_counts_not_retained": host_total - competitive_total,
            }
        )

    Path(args.sample_output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.sample_output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=sample_rows[0].keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(sample_rows)

    Path(args.gene_output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.gene_output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=gene_rows[0].keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(gene_rows)


if __name__ == "__main__":
    main()
