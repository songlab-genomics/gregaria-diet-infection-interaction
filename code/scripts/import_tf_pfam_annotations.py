#!/usr/bin/env python3
"""Lift protein-level TF/Pfam annotations to LOCID/GeneID identifiers.

The comparative-genomics GRN files use protein accessions such as XP_049826910.1,
whereas the time-course DEG and WGCNA reports use LOCID/GeneID identifiers.
This helper joins the protein table to data/list/allspecies_protein2geneid.tsv
and writes a compact annotation layer that can be overlaid on coexpression
modules without changing any expression counts.
"""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tf-csv", required=True, help="CSV with accession, query_name, and target_name columns")
    parser.add_argument("--protein2gene", required=True, help="Project TSV with protein_id, GeneID, Description, Species")
    parser.add_argument("--species", required=True, help="Species label to keep from protein2gene")
    parser.add_argument("--out-table", required=True, help="Output LOCID-level TF/Pfam TSV")
    parser.add_argument("--manifest", required=True, help="Output import summary TSV")
    return parser.parse_args()


def read_protein_map(path, species):
    mapping = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Species") != species:
                continue
            protein = row.get("protein_id", "").strip()
            if protein:
                mapping[protein] = row
    return mapping


def main():
    args = parse_args()
    protein_map = read_protein_map(args.protein2gene, args.species)
    rows = []

    with open(args.tf_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            protein = row.get("query_name", "").strip()
            mapped = protein_map.get(protein, {})
            rows.append({
                "GeneID": mapped.get("GeneID", ""),
                "protein_id": protein,
                "pfam_accession": row.get("accession", ""),
                "tf_domain": row.get("target_name", ""),
                "Project_description": mapped.get("Description", ""),
                "Species": mapped.get("Species", ""),
                "mapped_to_locid": "yes" if mapped else "no",
            })

    out_table = Path(args.out_table)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    with open(out_table, "w", newline="") as handle:
        fieldnames = [
            "GeneID", "protein_id", "pfam_accession", "tf_domain",
            "Project_description", "Species", "mapped_to_locid",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    mapped_rows = [row for row in rows if row["mapped_to_locid"] == "yes"]
    unique_locids = {row["GeneID"] for row in mapped_rows if row["GeneID"]}
    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["source_tf_csv", args.tf_csv])
        writer.writerow(["protein2gene", args.protein2gene])
        writer.writerow(["species", args.species])
        writer.writerow(["rows", len(rows)])
        writer.writerow(["mapped_to_locid", len(mapped_rows)])
        writer.writerow(["unique_mapped_locids", len(unique_locids)])


if __name__ == "__main__":
    main()
