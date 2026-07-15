#!/usr/bin/env python3
"""Build DesertLocustR-style GO/KEGG universe tables from eggNOG outputs.

The older DesertLocustR workflow converted eggNOG protein annotations into
LOCID-based enrichment tables. This script reproduces that table layout from
the current Snakemake eggNOG outputs without changing any read-count files.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-annotations", required=True, help="LOCID-level annotation TSV")
    parser.add_argument("--gene-go", required=True, help="Long LOCID-to-GO TSV")
    parser.add_argument("--gene-kegg", required=True, help="Long LOCID-to-KEGG TSV")
    parser.add_argument("--go-name-table", default="", help="Optional Blast2GO custom table for GO names")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--prefix", required=True, help="Output filename prefix")
    parser.add_argument("--species", required=True, help="Species label")
    parser.add_argument("--run-id", required=True, help="eggNOG run identifier")
    return parser.parse_args()


def read_tsv(path):
    if not path or not Path(path).exists():
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def unique_rows(rows, keys):
    seen = set()
    out = []
    for row in rows:
        marker = tuple(row.get(key, "") for key in keys)
        if marker in seen:
            continue
        seen.add(marker)
        out.append(row)
    return out


def clean_value(value):
    value = (value or "").strip()
    return "" if value in {"", "-", "NA", "nan"} else value


def clean_go(go_id):
    go_id = clean_value(go_id)
    return go_id if go_id.startswith("GO:") else ""


def clean_ko(term):
    term = clean_value(term).replace("ko:", "")
    return term if term.startswith("K") and len(term) == 6 else ""


def clean_pathway(term):
    term = clean_value(term).replace("path:", "")
    if term.startswith("ko") and len(term) == 7:
        return "map" + term[-5:]
    return term if term.startswith(("map", "ko")) and len(term) == 8 else ""


def parse_go_name_table(path):
    """Read optional Blast2GO custom export: category, GO ID, and GO term."""
    go_info = {}
    if not path or not Path(path).exists():
        return go_info
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            raw = row.get("GO in Extended Format (GO Category GO ID GO Term)", "")
            parts = raw.split()
            if len(parts) < 2:
                continue
            ontology = parts[0]
            go_id = parts[1]
            name = " ".join(parts[2:]) if len(parts) > 2 else go_id
            if go_id.startswith("GO:"):
                go_info[go_id] = {"name": name, "ontology": ontology}
    return go_info


def collapse_terms(rows, group_key, term_key):
    grouped = defaultdict(set)
    for row in rows:
        group = clean_value(row.get(group_key))
        term = clean_value(row.get(term_key))
        if group and term:
            grouped[group].add(term)
    return [
        {group_key: group, term_key: ",".join(sorted(terms))}
        for group, terms in sorted(grouped.items())
    ]


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    annotations = read_tsv(args.gene_annotations)
    go_rows_raw = read_tsv(args.gene_go)
    kegg_rows_raw = read_tsv(args.gene_kegg)
    go_info = parse_go_name_table(args.go_name_table)

    annotation_by_gene = {}
    annotation_by_protein = {}
    for row in annotations:
        gene = clean_value(row.get("GeneID"))
        protein = clean_value(row.get("protein_id"))
        desc = clean_value(row.get("Project_description")) or clean_value(row.get("Description"))
        preferred = clean_value(row.get("Preferred_name")) or clean_value(row.get("preferred_name"))
        if gene and gene not in annotation_by_gene:
            annotation_by_gene[gene] = {"Description": desc, "Preferred_name": preferred}
        if protein and protein not in annotation_by_protein:
            annotation_by_protein[protein] = {"GeneID": gene, "Description": desc, "Preferred_name": preferred}

    go_rows = []
    for row in go_rows_raw:
        gene = clean_value(row.get("GeneID"))
        protein = clean_value(row.get("protein_id"))
        go_id = clean_go(row.get("GO_ID"))
        if gene and go_id:
            info = go_info.get(go_id, {})
            go_rows.append({
                "GeneID": gene,
                "protein_id": protein,
                "GO_ID": go_id,
                "GO_Name": info.get("name", go_id),
                "Ontology": info.get("ontology", "unknown"),
            })
    go_rows = unique_rows(go_rows, ["GeneID", "protein_id", "GO_ID"])

    term2gene = unique_rows(
        [{"term": row["GO_ID"], "gene": row["GeneID"]} for row in go_rows],
        ["term", "gene"],
    )
    term2name = unique_rows(
        [{"term": row["GO_ID"], "name": row["GO_Name"], "ontology": row["Ontology"]} for row in go_rows],
        ["term"],
    )
    term2gene.sort(key=lambda row: (row["term"], row["gene"]))
    term2name.sort(key=lambda row: row["term"])

    # DesertLocustR-compatible names preserve the original function's columns.
    goterm2loc = [
        {"GOs": row["term"], "Names": row["name"], "Ontology": row["ontology"]}
        for row in term2name
    ]
    go_all_terms = [
        {"GOs": row["term"], "X.query": row["gene"]}
        for row in term2gene
    ]
    go_by_ontology = {
        ontology: [
            {"GOs": row["GO_ID"], "X.query": row["GeneID"]}
            for row in go_rows
            if row["Ontology"] == ontology
        ]
        for ontology in ["P", "F", "C", "BP", "MF", "CC"]
    }

    kegg_ko = []
    kegg_pathway = []
    for row in kegg_rows_raw:
        gene = clean_value(row.get("GeneID"))
        source = clean_value(row.get("KEGG_source"))
        term = clean_value(row.get("KEGG_term"))
        if not gene:
            continue
        if source == "KEGG_ko":
            ko = clean_ko(term)
            if ko:
                kegg_ko.append({"KEGG_ko": ko, "X.query": gene})
        elif source == "KEGG_Pathway":
            pathway = clean_pathway(term)
            if pathway:
                kegg_pathway.append({"term": pathway, "gene": gene})
    kegg_ko = unique_rows(kegg_ko, ["KEGG_ko", "X.query"])
    kegg_pathway = unique_rows(kegg_pathway, ["term", "gene"])
    kegg_ko.sort(key=lambda row: (row["KEGG_ko"], row["X.query"]))
    kegg_pathway.sort(key=lambda row: (row["term"], row["gene"]))

    protein_go = collapse_terms(
        [{"protein": row["protein_id"], "GOs": row["GO_ID"]} for row in go_rows if row["protein_id"]],
        "protein",
        "GOs",
    )
    gene_go = collapse_terms(
        [{"GeneID": row["GeneID"], "GOs": row["GO_ID"]} for row in go_rows],
        "GeneID",
        "GOs",
    )
    gene_go_locid = []
    for row in gene_go:
        annot = annotation_by_gene.get(row["GeneID"], {})
        gene_go_locid.append({
            "GeneID": row["GeneID"],
            "GOs": row["GOs"],
            "Description": annot.get("Description", ""),
            "Preferred_name": annot.get("Preferred_name", ""),
        })

    protkey = []
    for protein, annot in sorted(annotation_by_protein.items()):
        protkey.append({
            "protein": protein,
            "gene": annot.get("GeneID", ""),
            "description": annot.get("Description", ""),
            "preferred_name": annot.get("Preferred_name", ""),
        })

    prefix = args.prefix
    write_tsv(out_dir / f"{prefix}.GO_Universe.tsv", ["protein", "GOs"], protein_go)
    write_tsv(out_dir / f"{prefix}.GO_Universe_LOCID.tsv", ["GeneID", "GOs", "Description", "Preferred_name"], gene_go_locid)
    write_tsv(out_dir / f"{prefix}.protkey.tsv", ["protein", "gene", "description", "preferred_name"], protkey)
    write_tsv(out_dir / f"{prefix}.GO_TERM2GENE.tsv", ["term", "gene"], term2gene)
    write_tsv(out_dir / f"{prefix}.GO_TERM2NAME.tsv", ["term", "name", "ontology"], term2name)
    write_tsv(out_dir / f"{prefix}.GOTERM2LOC.tsv", ["GOs", "Names", "Ontology"], goterm2loc)
    write_tsv(out_dir / f"{prefix}.GO_ALL_TERMS.tsv", ["GOs", "X.query"], go_all_terms)
    write_tsv(out_dir / f"{prefix}.GO_BP_TERMS.tsv", ["GOs", "X.query"], go_by_ontology["P"] + go_by_ontology["BP"])
    write_tsv(out_dir / f"{prefix}.GO_MF_TERMS.tsv", ["GOs", "X.query"], go_by_ontology["F"] + go_by_ontology["MF"])
    write_tsv(out_dir / f"{prefix}.GO_CC_TERMS.tsv", ["GOs", "X.query"], go_by_ontology["C"] + go_by_ontology["CC"])
    write_tsv(out_dir / f"{prefix}.KEGGTERM2LOC.tsv", ["KEGG_ko", "X.query"], kegg_ko)
    write_tsv(out_dir / f"{prefix}.KEGG_KO_TERM2GENE.tsv", ["term", "gene"], [{"term": row["KEGG_ko"], "gene": row["X.query"]} for row in kegg_ko])
    write_tsv(out_dir / f"{prefix}.KEGG_PATHWAY_TERM2GENE.tsv", ["term", "gene"], kegg_pathway)

    file_index = [
        {"file": f"{prefix}.GO_Universe.tsv", "purpose": "Bellini-style protein-to-GO universe"},
        {"file": f"{prefix}.GO_Universe_LOCID.tsv", "purpose": "LOCID-to-GO universe for DEG genes"},
        {"file": f"{prefix}.protkey.tsv", "purpose": "protein-to-LOCID and description key"},
        {"file": f"{prefix}.GO_TERM2GENE.tsv", "purpose": "clusterProfiler GO TERM2GENE"},
        {"file": f"{prefix}.GO_TERM2NAME.tsv", "purpose": "clusterProfiler GO TERM2NAME"},
        {"file": f"{prefix}.GOTERM2LOC.tsv", "purpose": "DesertLocustR GO term-name table"},
        {"file": f"{prefix}.GO_BP_TERMS.tsv", "purpose": "DesertLocustR biological-process TERM2GENE"},
        {"file": f"{prefix}.GO_MF_TERMS.tsv", "purpose": "DesertLocustR molecular-function TERM2GENE"},
        {"file": f"{prefix}.GO_CC_TERMS.tsv", "purpose": "DesertLocustR cellular-component TERM2GENE"},
        {"file": f"{prefix}.KEGGTERM2LOC.tsv", "purpose": "DesertLocustR KO-to-LOCID table"},
    ]
    write_tsv(out_dir / f"{prefix}.desertlocustr_file_index.tsv", ["file", "purpose"], file_index)

    manifest = [
        {"metric": "run_id", "value": args.run_id},
        {"metric": "species", "value": args.species},
        {"metric": "gene_annotation_rows", "value": str(len(annotations))},
        {"metric": "genes_with_go", "value": str(len({row["GeneID"] for row in go_rows}))},
        {"metric": "go_term2gene_rows", "value": str(len(term2gene))},
        {"metric": "go_terms_with_names", "value": str(sum(1 for row in term2name if row["name"] != row["term"]))},
        {"metric": "go_name_table", "value": args.go_name_table if args.go_name_table else "not provided"},
        {"metric": "kegg_ko_rows", "value": str(len(kegg_ko))},
        {"metric": "kegg_pathway_rows", "value": str(len(kegg_pathway))},
    ]
    write_tsv(out_dir / f"{prefix}.desertlocustr_manifest.tsv", ["metric", "value"], manifest)


if __name__ == "__main__":
    main()
