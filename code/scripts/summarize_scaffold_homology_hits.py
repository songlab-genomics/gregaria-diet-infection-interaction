#!/usr/bin/env python3
"""Summarize BLASTn/DIAMOND top-hit taxonomy by query, gene, and scaffold."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


HIT_FIELDS = [
    "qseqid",
    "qlen",
    "sseqid",
    "staxids",
    "pident",
    "length",
    "qcovhsp",
    "evalue",
    "bitscore",
    "stitle",
]

SUPPORT_ORDER = [
    "Locust/orthopteran support",
    "Broader arthropod support",
    "Other metazoan support",
    "Mixed/ambiguous",
    "Non-metazoan candidate",
    "No informative hit",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--candidate-genes", required=True, type=Path)
    parser.add_argument("--blastn-hits", required=True, type=Path)
    parser.add_argument("--diamond-hits", required=True, type=Path)
    parser.add_argument("--nodes", required=True, type=Path)
    parser.add_argument("--names", required=True, type=Path)
    parser.add_argument(
        "--exclude-taxid",
        action="append",
        type=int,
        default=[],
        help="Remove this taxid before selecting the top-score tier.",
    )
    parser.add_argument("--top-fraction", type=float, default=0.05)
    parser.add_argument("--output-top-hits", required=True, type=Path)
    parser.add_argument("--output-query-summary", required=True, type=Path)
    parser.add_argument("--output-gene-summary", required=True, type=Path)
    parser.add_argument("--output-scaffold-summary", required=True, type=Path)
    parser.add_argument("--output-non-metazoan", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_taxonomy(
    nodes_path: Path, names_path: Path
) -> tuple[dict[int, int], dict[int, str], dict[int, str]]:
    parents: dict[int, int] = {}
    ranks: dict[int, str] = {}
    with nodes_path.open() as handle:
        for line in handle:
            parts = [part.strip() for part in line.split("|")]
            if len(parts) < 3:
                continue
            taxid = int(parts[0])
            parents[taxid] = int(parts[1])
            ranks[taxid] = parts[2]

    names: dict[int, str] = {}
    with names_path.open() as handle:
        for line in handle:
            parts = [part.strip() for part in line.split("|")]
            if len(parts) < 4 or parts[3] != "scientific name":
                continue
            names[int(parts[0])] = parts[1]
    return parents, ranks, names


def lineage(
    taxid: int,
    parents: dict[int, int],
    names: dict[int, str],
) -> list[tuple[int, str]]:
    result: list[tuple[int, str]] = []
    seen: set[int] = set()
    current = taxid
    while current and current not in seen:
        seen.add(current)
        result.append((current, names.get(current, f"taxid_{current}")))
        parent = parents.get(current)
        if parent is None or parent == current:
            break
        current = parent
    result.reverse()
    return result


def broad_group(lineage_names: set[str]) -> str:
    if "Schistocerca" in lineage_names:
        return "Schistocerca"
    if "Acrididae" in lineage_names:
        return "Other Acrididae"
    if "Orthoptera" in lineage_names:
        return "Other Orthoptera"
    if "Arthropoda" in lineage_names:
        return "Other Arthropoda"
    if "Metazoa" in lineage_names:
        return "Other Metazoa"
    if "Fungi" in lineage_names:
        return "Fungi"
    if "Viridiplantae" in lineage_names or "Streptophyta" in lineage_names:
        return "Plants"
    if "Bacteria" in lineage_names:
        return "Bacteria"
    if "Archaea" in lineage_names:
        return "Archaea"
    if "Viruses" in lineage_names:
        return "Viruses"
    if "Eukaryota" in lineage_names:
        return "Other Eukaryota"
    return "Other/unclassified"


def support_from_groups(groups: set[str]) -> str:
    metazoan = groups & {
        "Schistocerca",
        "Other Acrididae",
        "Other Orthoptera",
        "Other Arthropoda",
        "Other Metazoa",
    }
    non_metazoan = groups & {
        "Fungi",
        "Plants",
        "Bacteria",
        "Archaea",
        "Viruses",
        "Other Eukaryota",
    }
    if non_metazoan and metazoan:
        return "Mixed/ambiguous"
    if groups & {"Schistocerca", "Other Acrididae", "Other Orthoptera"}:
        return "Locust/orthopteran support"
    if "Other Arthropoda" in groups:
        return "Broader arthropod support"
    if "Other Metazoa" in groups:
        return "Other metazoan support"
    if non_metazoan:
        return "Non-metazoan candidate"
    return "No informative hit"


def parse_hits(path: Path, method: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, fieldnames=HIT_FIELDS, delimiter="\t")
        for row in reader:
            if not row["qseqid"]:
                continue
            row["method"] = method
            rows.append(row)
    return rows


def parse_taxids(value: str) -> list[int]:
    taxids: list[int] = []
    for token in value.replace(",", ";").split(";"):
        token = token.strip()
        if token.isdigit() and int(token) > 0:
            taxids.append(int(token))
    return taxids


def combined_support(query_supports: set[str]) -> str:
    informative = query_supports - {"No informative hit"}
    if not informative:
        return "No informative hit"
    metazoan_support = informative & {
        "Locust/orthopteran support",
        "Broader arthropod support",
        "Other metazoan support",
    }
    if "Mixed/ambiguous" in informative:
        return "Mixed/ambiguous"
    if metazoan_support and "Non-metazoan candidate" in informative:
        return "Mixed/ambiguous"
    if "Locust/orthopteran support" in informative:
        return "Locust/orthopteran support"
    if "Broader arthropod support" in informative:
        return "Broader arthropod support"
    if "Other metazoan support" in informative:
        return "Other metazoan support"
    if informative == {"Non-metazoan candidate"}:
        return "Non-metazoan candidate"
    return "Mixed/ambiguous"


def main() -> None:
    args = parse_args()
    _, manifest_rows = read_tsv(args.manifest)
    _, gene_rows = read_tsv(args.candidate_genes)
    manifest_by_query = {row["query_id"]: row for row in manifest_rows}
    gene_by_id = {row["gene_id"]: row for row in gene_rows}
    parents, _, names = load_taxonomy(args.nodes, args.names)

    raw_hits = parse_hits(args.blastn_hits, "BLASTn core_nt")
    raw_hits.extend(parse_hits(args.diamond_hits, "DIAMOND nr_cluster_seq"))
    excluded_taxids = set(args.exclude_taxid)
    hits_by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
    for hit in raw_hits:
        hit_taxids = parse_taxids(hit["staxids"])
        retained_taxids = [
            taxid for taxid in hit_taxids if taxid not in excluded_taxids
        ]
        if hit_taxids and not retained_taxids:
            continue
        if retained_taxids != hit_taxids:
            hit = {**hit, "staxids": ";".join(map(str, retained_taxids))}
        if hit["qseqid"] in manifest_by_query:
            hits_by_query[hit["qseqid"]].append(hit)

    top_hit_rows: list[dict[str, object]] = []
    query_rows: list[dict[str, object]] = []
    for query_id, manifest in sorted(manifest_by_query.items()):
        hits = hits_by_query.get(query_id, [])
        top_hits: list[dict[str, str]] = []
        if hits:
            best_score = max(float(hit["bitscore"]) for hit in hits)
            cutoff = best_score * (1.0 - args.top_fraction)
            top_hits = [hit for hit in hits if float(hit["bitscore"]) >= cutoff]

        groups: set[str] = set()
        taxids: set[int] = set()
        taxon_names: set[str] = set()
        for hit in top_hits:
            hit_taxids = parse_taxids(hit["staxids"])
            if not hit_taxids:
                top_hit_rows.append(
                    {
                        **hit,
                        "gene_id": manifest["gene_id"],
                        "seqid": manifest["seqid"],
                        "taxid": "",
                        "scientific_name": "",
                        "broad_group": "Other/unclassified",
                        "lineage": "",
                    }
                )
                groups.add("Other/unclassified")
                continue
            for taxid in hit_taxids:
                tax_lineage = lineage(taxid, parents, names)
                lineage_names = {name for _, name in tax_lineage}
                group = broad_group(lineage_names)
                groups.add(group)
                taxids.add(taxid)
                taxon_names.add(names.get(taxid, f"taxid_{taxid}"))
                top_hit_rows.append(
                    {
                        **hit,
                        "gene_id": manifest["gene_id"],
                        "seqid": manifest["seqid"],
                        "taxid": taxid,
                        "scientific_name": names.get(taxid, f"taxid_{taxid}"),
                        "broad_group": group,
                        "lineage": ";".join(name for _, name in tax_lineage),
                    }
                )

        best_evalue = min((float(hit["evalue"]) for hit in top_hits), default=None)
        best_bitscore = max((float(hit["bitscore"]) for hit in top_hits), default=None)
        max_identity = max((float(hit["pident"]) for hit in top_hits), default=None)
        max_qcov = max((float(hit["qcovhsp"]) for hit in top_hits), default=None)
        query_rows.append(
            {
                **manifest,
                "top_hit_support": support_from_groups(groups),
                "top_hit_groups": ";".join(sorted(groups)),
                "top_hit_taxids": ";".join(map(str, sorted(taxids))),
                "top_hit_taxon_names": ";".join(sorted(taxon_names)),
                "top_tier_hit_records": len(top_hits),
                "best_evalue": "" if best_evalue is None else best_evalue,
                "best_bitscore": "" if best_bitscore is None else best_bitscore,
                "max_percent_identity": "" if max_identity is None else max_identity,
                "max_query_coverage_percent": "" if max_qcov is None else max_qcov,
            }
        )

    query_fields = list(query_rows[0])
    write_tsv(args.output_query_summary, query_rows, query_fields)
    top_fields = [
        "method",
        "qseqid",
        "gene_id",
        "seqid",
        "qlen",
        "sseqid",
        "taxid",
        "scientific_name",
        "broad_group",
        "lineage",
        "pident",
        "length",
        "qcovhsp",
        "evalue",
        "bitscore",
        "stitle",
    ]
    write_tsv(args.output_top_hits, top_hit_rows, top_fields)

    query_by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in query_rows:
        query_by_gene[str(row["gene_id"])].append(row)

    gene_summary_rows: list[dict[str, object]] = []
    for gene_id, source in sorted(gene_by_id.items()):
        rows = query_by_gene.get(gene_id, [])
        supports = {str(row["top_hit_support"]) for row in rows}
        nucleotide_supports = {
            str(row["top_hit_support"])
            for row in rows
            if row["query_type"] == "gene_genomic_DNA"
        }
        protein_supports = {
            str(row["top_hit_support"])
            for row in rows
            if row["query_type"] == "protein"
        }
        all_groups = {
            group
            for row in rows
            for group in str(row["top_hit_groups"]).split(";")
            if group
        }
        gene_summary_rows.append(
            {
                **source,
                "homology_origin_class": combined_support(supports),
                "nucleotide_support": ";".join(sorted(nucleotide_supports)),
                "protein_support": ";".join(sorted(protein_supports)),
                "top_hit_groups": ";".join(sorted(all_groups)),
                "n_nucleotide_queries": sum(
                    row["query_type"] == "gene_genomic_DNA" for row in rows
                ),
                "n_protein_queries": sum(row["query_type"] == "protein" for row in rows),
            }
        )

    gene_fields = list(gene_summary_rows[0])
    write_tsv(args.output_gene_summary, gene_summary_rows, gene_fields)
    non_metazoan_rows = [
        row
        for row in gene_summary_rows
        if row["homology_origin_class"] == "Non-metazoan candidate"
    ]
    write_tsv(args.output_non_metazoan, non_metazoan_rows, gene_fields)

    scaffold_to_genes: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in gene_summary_rows:
        scaffold_to_genes[str(row["seqid"])].append(row)

    scaffold_rows: list[dict[str, object]] = []
    for scaffold, rows in sorted(scaffold_to_genes.items()):
        counts = Counter(str(row["homology_origin_class"]) for row in rows)
        scaffold_row: dict[str, object] = {
            "scaffold": scaffold,
            "candidate_genes": len(rows),
        }
        for status in SUPPORT_ORDER:
            scaffold_row[status] = counts.get(status, 0)
        scaffold_row["non_metazoan_fraction"] = (
            counts.get("Non-metazoan candidate", 0) / len(rows)
        )
        scaffold_rows.append(scaffold_row)
    scaffold_fields = list(scaffold_rows[0])
    write_tsv(args.output_scaffold_summary, scaffold_rows, scaffold_fields)

    print(
        f"Summarized {len(query_rows)} queries, {len(gene_summary_rows)} genes, "
        f"and {len(scaffold_rows)} scaffolds; "
        f"{len(non_metazoan_rows)} genes are non-metazoan candidates."
    )


if __name__ == "__main__":
    main()
