#!/usr/bin/env python3
"""Freeze a DEG-derived scaffold set with NCBI assembly/FCS provenance."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", required=True, help="Unified-model DEG gene TSV.")
    parser.add_argument("--assembly-report", required=True)
    parser.add_argument("--fcs-report", required=True)
    parser.add_argument("--source-run", required=True)
    parser.add_argument("--source-model", required=True)
    parser.add_argument("--output-list", required=True)
    parser.add_argument("--output-genes", required=True)
    parser.add_argument("--output-provenance", required=True)
    return parser.parse_args()


def read_assembly_report(path):
    header = None
    records = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("# Sequence-Name"):
                header = line[2:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                raise ValueError("Assembly report column header was not found.")
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            for accession_field in ("RefSeq-Accn", "GenBank-Accn", "Sequence-Name"):
                accession = row.get(accession_field, "")
                if accession and accession != "na":
                    records[accession] = row
    return records


def read_fcs_report(path):
    metadata = {}
    flagged = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                payload = json.loads(line[2:])
                if isinstance(payload, list) and len(payload) > 1:
                    metadata = payload[1]
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8:
                flagged[fields[0]] = {
                    "fcs_action": fields[4],
                    "fcs_contam_type": fields[5],
                    "fcs_coverage": fields[6],
                    "fcs_contam_details": fields[7],
                }
    return metadata, flagged


def numeric(values, reducer):
    parsed = []
    for value in values:
        try:
            parsed.append(float(value))
        except (TypeError, ValueError):
            continue
    return str(reducer(parsed)) if parsed else ""


def consolidate_gene_rows(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["seqid"], row["gene_id"])].append(row)

    consolidated = []
    for _, memberships in sorted(grouped.items()):
        row = dict(memberships[0])
        contrasts = {
            contrast.strip()
            for membership in memberships
            for contrast in membership.get("contrasts", "").split(";")
            if contrast.strip()
        }
        support = {
            membership.get("mapping_support", "").strip()
            for membership in memberships
            if membership.get("mapping_support", "").strip()
        }
        row["contrasts"] = ";".join(sorted(contrasts))
        row["n_contrasts"] = str(len(contrasts))
        row["mapping_support"] = "; ".join(sorted(support))
        row["min_padj"] = numeric(
            [membership.get("min_padj") for membership in memberships], min
        )
        row["max_abs_log2FoldChange"] = numeric(
            [
                membership.get("max_abs_log2FoldChange")
                for membership in memberships
            ],
            max,
        )
        row["competitive_host_retained_fraction"] = numeric(
            [
                membership.get("competitive_host_retained_fraction")
                for membership in memberships
            ],
            min,
        )
        row["host_counts_not_retained"] = numeric(
            [
                membership.get("host_counts_not_retained")
                for membership in memberships
            ],
            max,
        )
        consolidated.append(row)
    return consolidated


def main():
    args = parse_args()
    with open(args.genes, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        gene_fields = reader.fieldnames or []
        membership_rows = [
            row
            for row in reader
            if row.get("seqid", "").startswith("NW_")
        ]

    if not membership_rows or "gene_id" not in gene_fields or "seqid" not in gene_fields:
        raise ValueError("No NW_ scaffold DEG genes found in the input table.")
    gene_rows = consolidate_gene_rows(membership_rows)

    genes_by_scaffold = defaultdict(set)
    contrasts_by_scaffold = defaultdict(set)
    for row in gene_rows:
        genes_by_scaffold[row["seqid"]].add(row["gene_id"])
        for contrast in row.get("contrasts", "").split(";"):
            contrast = contrast.strip()
            if contrast:
                contrasts_by_scaffold[row["seqid"]].add(contrast)

    assembly = read_assembly_report(args.assembly_report)
    fcs_metadata, fcs_flagged = read_fcs_report(args.fcs_report)
    scaffolds = sorted(genes_by_scaffold)

    for output_path in (
        args.output_list,
        args.output_genes,
        args.output_provenance,
    ):
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(args.output_list, "w", encoding="utf-8") as handle:
        handle.write("\n".join(scaffolds) + "\n")

    with open(args.output_genes, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=gene_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(gene_rows, key=lambda row: (row["seqid"], row["gene_id"])))

    db_metadata = fcs_metadata.get("db", {})
    run_metadata = fcs_metadata.get("run-info", {})
    provenance_fields = [
        "scaffold",
        "source_run",
        "source_model",
        "n_unique_deg_genes",
        "n_deg_contrasts",
        "deg_contrasts",
        "assembly_sequence_role",
        "assembly_unit",
        "assembly_length_bp",
        "assembly_bioproject",
        "assembly_biosample",
        "ncbi_fcs_run_date",
        "ncbi_fcs_db_build_date",
        "ncbi_fcs_asserted_division",
        "ncbi_fcs_flagged",
        "fcs_action",
        "fcs_contam_type",
        "fcs_coverage",
        "fcs_contam_details",
    ]
    with open(args.output_provenance, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=provenance_fields, delimiter="\t")
        writer.writeheader()
        for scaffold in scaffolds:
            assembly_row = assembly.get(scaffold, {})
            fcs_row = fcs_flagged.get(scaffold, {})
            writer.writerow(
                {
                    "scaffold": scaffold,
                    "source_run": args.source_run,
                    "source_model": args.source_model,
                    "n_unique_deg_genes": len(genes_by_scaffold[scaffold]),
                    "n_deg_contrasts": len(contrasts_by_scaffold[scaffold]),
                    "deg_contrasts": ";".join(sorted(contrasts_by_scaffold[scaffold])),
                    "assembly_sequence_role": assembly_row.get("Sequence-Role", ""),
                    "assembly_unit": assembly_row.get("Assembly-Unit", ""),
                    "assembly_length_bp": assembly_row.get("Sequence-Length", ""),
                    "assembly_bioproject": "PRJNA782545",
                    "assembly_biosample": "SAMN23383904",
                    "ncbi_fcs_run_date": fcs_metadata.get("run-date", ""),
                    "ncbi_fcs_db_build_date": db_metadata.get("build-date", ""),
                    "ncbi_fcs_asserted_division": run_metadata.get("asserted-div", ""),
                    "ncbi_fcs_flagged": "yes" if scaffold in fcs_flagged else "no",
                    "fcs_action": fcs_row.get("fcs_action", ""),
                    "fcs_contam_type": fcs_row.get("fcs_contam_type", ""),
                    "fcs_coverage": fcs_row.get("fcs_coverage", ""),
                    "fcs_contam_details": fcs_row.get("fcs_contam_details", ""),
                }
            )

    print(
        f"Exported {len(gene_rows)} unique DEG genes from "
        f"{len(membership_rows)} memberships on {len(scaffolds)} candidate scaffolds."
    )


if __name__ == "__main__":
    main()
