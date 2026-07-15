#!/usr/bin/env python3
"""Convert eggNOG-mapper protein annotations into gene-level tables.

The DEG reports use LOCID/GeneID identifiers, while eggNOG-mapper annotates
protein accessions such as XP_049849364.1. This script joins eggNOG output to
the project protein-to-gene map and writes compact GO and KEGG/KO tables for
clusterProfiler-style enrichment.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--emapper", required=True, help="eggNOG .emapper.annotations file")
    parser.add_argument("--protein2gene", required=True, help="TSV with protein_id, GeneID, Description, Species")
    parser.add_argument("--species", required=True, help="Species label to retain from protein2gene")
    parser.add_argument("--gene-annotations", required=True, help="Output gene-level annotation TSV")
    parser.add_argument("--gene-go", required=True, help="Output long GO TSV")
    parser.add_argument("--gene-kegg", required=True, help="Output long KEGG/KO TSV")
    parser.add_argument("--manifest", required=True, help="Output manifest/statistics TSV")
    parser.add_argument("--run-id", required=True, help="Run identifier")
    parser.add_argument("--protein-fasta", required=True, help="Protein FASTA used by eggNOG")
    parser.add_argument("--eggnog-data-dir", required=True, help="eggNOG database directory")
    return parser.parse_args()


def split_terms(value):
    """Split comma/semicolon-delimited annotation terms and drop empty values."""
    if value is None:
        return []
    value = value.strip()
    if value in {"", "-"}:
        return []
    pieces = []
    for chunk in value.replace(";", ",").split(","):
        chunk = chunk.strip()
        if chunk and chunk != "-":
            pieces.append(chunk)
    return sorted(set(pieces))


def read_protein_map(path, species):
    """Return protein_id -> project GeneID/LOCID for the requested species."""
    mapping = {}
    descriptions = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("Species") != species:
                continue
            protein = row.get("protein_id", "").strip()
            gene = row.get("GeneID", "").strip()
            if protein and gene:
                mapping[protein] = gene
                descriptions[protein] = row.get("Description", "").strip()
    return mapping, descriptions


def read_emapper(path):
    """Read eggNOG annotations while respecting its commented header format."""
    header = None
    rows = []
    with open(path, newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("#query"):
                header = line.lstrip("#").split("\t")
                continue
            if not line or line.startswith("#"):
                continue
            if header is None:
                continue
            rows.append(dict(zip(header, line.split("\t"))))
    if header is None:
        raise SystemExit(f"No #query header was found in {path}")
    return header, rows


def first_present(row, names):
    for name in names:
        if name in row:
            return row.get(name, "")
    return ""


def write_table(path, fieldnames, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()
    protein_to_gene, protein_descriptions = read_protein_map(args.protein2gene, args.species)
    header, emapper_rows = read_emapper(args.emapper)

    gene_rows = []
    go_rows = []
    kegg_rows = []
    genes_with_go = set()
    genes_with_kegg = set()
    proteins_mapped = 0
    genes_seen = set()

    preferred_columns = [
        "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
        "COG_category", "Description", "Preferred_name", "preferred_name",
        "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
        "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
        "BiGG_Reaction", "PFAMs",
    ]
    annotation_columns = [col for col in preferred_columns if col in header]
    gene_fieldnames = ["GeneID", "protein_id", "Project_description"] + annotation_columns

    for row in emapper_rows:
        protein = row.get("query", "").strip()
        gene = protein_to_gene.get(protein, protein)
        if protein in protein_to_gene:
            proteins_mapped += 1
        genes_seen.add(gene)

        gene_row = {
            "GeneID": gene,
            "protein_id": protein,
            "Project_description": protein_descriptions.get(protein, ""),
        }
        for col in annotation_columns:
            gene_row[col] = row.get(col, "")
        gene_rows.append(gene_row)

        go_terms = [term for term in split_terms(row.get("GOs", "")) if term.startswith("GO:")]
        for go_id in go_terms:
            go_rows.append({"GeneID": gene, "protein_id": protein, "GO_ID": go_id})
            genes_with_go.add(gene)

        kegg_sources = {
            "KEGG_ko": first_present(row, ["KEGG_ko", "KEGG.KO"]),
            "KEGG_Pathway": first_present(row, ["KEGG_Pathway", "KEGG.pathway", "KEGG_Pathways"]),
            "KEGG_Module": first_present(row, ["KEGG_Module", "KEGG.module"]),
            "KEGG_Reaction": first_present(row, ["KEGG_Reaction", "KEGG.reaction"]),
        }
        for source, value in kegg_sources.items():
            for term in split_terms(value):
                kegg_rows.append({
                    "GeneID": gene,
                    "protein_id": protein,
                    "KEGG_source": source,
                    "KEGG_term": term,
                })
                genes_with_kegg.add(gene)

    # Remove exact duplicate long-table rows while preserving deterministic order.
    go_rows = [dict(row) for row in sorted({tuple(row.items()) for row in go_rows})]
    kegg_rows = [dict(row) for row in sorted({tuple(row.items()) for row in kegg_rows})]

    write_table(args.gene_annotations, gene_fieldnames, gene_rows)
    write_table(args.gene_go, ["GeneID", "protein_id", "GO_ID"], go_rows)
    write_table(args.gene_kegg, ["GeneID", "protein_id", "KEGG_source", "KEGG_term"], kegg_rows)

    stats = {
        "run_id": args.run_id,
        "species": args.species,
        "protein_fasta": args.protein_fasta,
        "eggnog_data_dir": args.eggnog_data_dir,
        "emapper_annotations": args.emapper,
        "protein2gene": args.protein2gene,
        "emapper_protein_rows": len(emapper_rows),
        "proteins_mapped_to_project_gene": proteins_mapped,
        "unique_gene_ids": len(genes_seen),
        "unique_gene_ids_with_go": len(genes_with_go),
        "unique_gene_ids_with_kegg": len(genes_with_kegg),
        "go_rows": len(go_rows),
        "kegg_rows": len(kegg_rows),
    }
    write_table(args.manifest, ["metric", "value"], [{"metric": k, "value": v} for k, v in stats.items()])


if __name__ == "__main__":
    main()
