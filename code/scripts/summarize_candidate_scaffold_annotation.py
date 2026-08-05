#!/usr/bin/env python3
"""Summarize all host annotations and DEG membership on candidate scaffolds."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-list", required=True)
    parser.add_argument("--candidate-genes", required=True)
    parser.add_argument("--annotation", required=True, help="NCBI GTF or GFF3.")
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--output-genes", required=True)
    parser.add_argument("--output-scaffolds", required=True)
    return parser.parse_args()


def parse_attributes(text):
    attributes = {}
    for item in text.strip().strip(";").split(";"):
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
        attributes[key.strip()] = value.strip()
    return attributes


def gene_id_from_attributes(attributes):
    gene_id = attributes.get("gene_id")
    if gene_id:
        return gene_id.replace("gene-", "", 1)
    for item in attributes.get("Dbxref", "").split(","):
        if item.startswith("GeneID:"):
            return "LOC" + item.split(":", 1)[1]
    gene_id = attributes.get("gene") or attributes.get("Name")
    if gene_id:
        return gene_id.replace("gene-", "", 1)
    return attributes.get("ID", "").replace("gene-", "", 1)


def read_ids(path, field=None):
    with open(path, encoding="utf-8", newline="") as handle:
        if field is None:
            return {
                line.strip()
                for line in handle
                if line.strip() and not line.startswith("#")
            }
        return {
            row[field]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get(field)
        }


def read_metrics(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return {
            row["scaffold"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def main():
    args = parse_args()
    candidates = read_ids(args.candidate_list)
    candidate_deg_genes = read_ids(args.candidate_genes, "gene_id")
    metrics = read_metrics(args.metrics)

    rows = []
    genes_by_scaffold = defaultdict(set)
    degs_by_scaffold = defaultdict(set)
    with open(args.annotation, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in candidates:
                continue
            attributes = parse_attributes(fields[8])
            gene_id = gene_id_from_attributes(attributes)
            if not gene_id:
                continue
            is_candidate_deg = gene_id in candidate_deg_genes
            genes_by_scaffold[fields[0]].add(gene_id)
            if is_candidate_deg:
                degs_by_scaffold[fields[0]].add(gene_id)
            rows.append(
                {
                    "scaffold": fields[0],
                    "gene_id": gene_id,
                    "start": fields[3],
                    "end": fields[4],
                    "strand": fields[6],
                    "gene_biotype": attributes.get(
                        "gene_biotype", attributes.get("gene_type", "")
                    ),
                    "description": attributes.get(
                        "description", attributes.get("product", "")
                    ),
                    "candidate_deg": "yes" if is_candidate_deg else "no",
                }
            )

    gene_fields = [
        "scaffold",
        "gene_id",
        "start",
        "end",
        "strand",
        "gene_biotype",
        "description",
        "candidate_deg",
    ]
    Path(args.output_genes).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_genes, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=gene_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["scaffold"], int(row["start"]))))

    scaffold_fields = [
        "scaffold",
        "length_bp",
        "gc_percent",
        "n_percent",
        "annotated_genes",
        "candidate_deg_genes",
        "annotated_genes_per_mb",
        "candidate_deg_genes_per_mb",
    ]
    with open(args.output_scaffolds, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=scaffold_fields, delimiter="\t")
        writer.writeheader()
        for scaffold in sorted(candidates):
            metric = metrics.get(scaffold, {})
            length = int(float(metric.get("length_bp", 0) or 0))
            scale = length / 1_000_000 if length else 0
            annotated_count = len(genes_by_scaffold[scaffold])
            deg_count = len(degs_by_scaffold[scaffold])
            writer.writerow(
                {
                    "scaffold": scaffold,
                    "length_bp": length,
                    "gc_percent": metric.get("gc_percent", ""),
                    "n_percent": metric.get("n_percent", ""),
                    "annotated_genes": annotated_count,
                    "candidate_deg_genes": deg_count,
                    "annotated_genes_per_mb": annotated_count / scale if scale else "",
                    "candidate_deg_genes_per_mb": deg_count / scale if scale else "",
                }
            )


if __name__ == "__main__":
    main()
