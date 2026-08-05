#!/usr/bin/env python3
"""Combine assembly, annotation, FCS, and Kraken evidence conservatively."""

import argparse
import csv
from pathlib import Path


NON_HOST_GROUPS = {"Fungi", "Bacteria", "Archaea", "Viruses"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--annotation-summary", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--cross-species", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--follow-up", required=True)
    return parser.parse_args()


def read_keyed(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return {
            row["scaffold"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def read_cross_species(path):
    all_rows = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("comparison_species") == "not_run":
                continue
            scaffold = row["scaffold"]
            all_rows.setdefault(scaffold, []).append(row)

    rows_by_scaffold = {}
    for scaffold, scaffold_rows in all_rows.items():
        def high_conf_coverage(row):
            return to_float(
                row.get(
                    "high_conf_query_coverage_percent",
                    row.get("query_coverage_percent"),
                )
            )

        def high_conf_identity(row):
            return to_float(
                row.get(
                    "high_conf_weighted_identity_percent",
                    row.get("weighted_identity_percent"),
                )
            )

        best = max(
            scaffold_rows,
            key=lambda row: (high_conf_coverage(row), high_conf_identity(row)),
        )
        strong_species = [
            row["comparison_species"]
            for row in scaffold_rows
            if high_conf_coverage(row) >= 50 and high_conf_identity(row) >= 70
        ]
        partial_species = [
            row["comparison_species"]
            for row in scaffold_rows
            if high_conf_coverage(row) >= 10 and high_conf_identity(row) >= 70
        ]
        rows_by_scaffold[scaffold] = {
            "cross_species_best_species": best.get("comparison_species", ""),
            "cross_species_max_coverage_percent": high_conf_coverage(best),
            "cross_species_weighted_identity_percent": high_conf_identity(best),
            "cross_species_best_target_sequence": best.get(
                "best_target_sequence", ""
            ),
            "cross_species_strong_support_count": len(strong_species),
            "cross_species_strong_support_species": ";".join(strong_species),
            "cross_species_partial_support_count": len(partial_species),
            "cross_species_partial_support_species": ";".join(partial_species),
        }
    return rows_by_scaffold


def interpret(row):
    group = row.get("kraken_broad_group", "")
    fcs_flagged = row.get("ncbi_fcs_flagged", "").lower() == "yes"
    strong_support_count = int(
        to_float(row.get("cross_species_strong_support_count"))
    )
    partial_support_count = int(
        to_float(row.get("cross_species_partial_support_count"))
    )
    if strong_support_count >= 1 and (fcs_flagged or group in NON_HOST_GROUPS):
        return (
            "Conflicting host/non-host evidence; high-priority review",
            "Inspect alignment specificity, then confirm with BLASTn/DIAMOND and genomic DNA coverage.",
        )
    if fcs_flagged:
        return (
            "NCBI FCS flagged for review",
            "Inspect the FCS interval/action, then confirm with nucleotide and protein homology.",
        )
    if group in NON_HOST_GROUPS:
        return (
            "Non-host taxonomic signal; review",
            "Confirm with BLASTn/DIAMOND and original genomic DNA coverage before filtering.",
        )
    if strong_support_count >= 2:
        return (
            "Strong cross-Schistocerca support for host origin",
            "Confirm with original genomic DNA coverage/Hi-C; retain as host unless stronger contradictory evidence appears.",
        )
    if strong_support_count == 1 or partial_support_count >= 2:
        return (
            "Partial cross-Schistocerca support; likely host-associated",
            "Retain pending genomic DNA coverage/Hi-C confirmation and inspect repeat content.",
        )
    if group == "Arthropoda":
        return (
            "Arthropod-compatible taxonomic signal",
            "Test conservation in other Schistocerca genomes and original genomic DNA coverage.",
        )
    if group == "Metazoa_non_arthropod":
        return (
            "Metazoan but not arthropod-specific signal; unresolved",
            "Use cross-Schistocerca alignment plus BLASTn/DIAMOND and genomic DNA coverage.",
        )
    if group == "Other_eukaryote":
        return (
            "Eukaryotic but not host-specific signal; unresolved",
            "Use cross-Schistocerca alignment plus BLASTn/DIAMOND and genomic DNA coverage.",
        )
    return (
        "Unresolved by Kraken/FCS",
        "Use cross-Schistocerca alignment, BLASTn/DIAMOND, and original genomic DNA coverage.",
    )


def main():
    args = parse_args()
    provenance = read_keyed(args.provenance)
    annotation = read_keyed(args.annotation_summary)
    taxonomy = read_keyed(args.taxonomy)
    cross_species = read_cross_species(args.cross_species)
    scaffolds = sorted(
        set(provenance) | set(annotation) | set(taxonomy) | set(cross_species)
    )
    rows = []
    follow_up_rows = []

    for scaffold in scaffolds:
        row = {
            **provenance.get(scaffold, {}),
            **annotation.get(scaffold, {}),
            **taxonomy.get(scaffold, {}),
            **cross_species.get(scaffold, {}),
        }
        row["scaffold"] = scaffold
        status, next_test = interpret(row)
        row["preliminary_evidence_status"] = status
        row["recommended_next_test"] = next_test
        row["automatic_filter_decision"] = "retain_pending_independent_validation"
        rows.append(row)
        follow_up_rows.append(
            {
                "scaffold": scaffold,
                "preliminary_evidence_status": status,
                "recommended_next_test": next_test,
                "priority_reason": (
                    f"{row.get('candidate_deg_genes', row.get('n_unique_deg_genes', ''))} "
                    "candidate DEG genes; "
                    f"Kraken group={row.get('kraken_broad_group', 'unresolved')}"
                ),
            }
        )

    preferred_fields = [
        "scaffold",
        "source_run",
        "source_model",
        "assembly_sequence_role",
        "assembly_unit",
        "length_bp",
        "gc_percent",
        "n_percent",
        "annotated_genes",
        "candidate_deg_genes",
        "annotated_genes_per_mb",
        "candidate_deg_genes_per_mb",
        "n_deg_contrasts",
        "deg_contrasts",
        "kraken_status",
        "kraken_taxid",
        "kraken_taxon_name",
        "kraken_broad_group",
        "kraken_lineage",
        "cross_species_best_species",
        "cross_species_max_coverage_percent",
        "cross_species_weighted_identity_percent",
        "cross_species_best_target_sequence",
        "cross_species_strong_support_count",
        "cross_species_strong_support_species",
        "cross_species_partial_support_count",
        "cross_species_partial_support_species",
        "ncbi_fcs_run_date",
        "ncbi_fcs_db_build_date",
        "ncbi_fcs_asserted_division",
        "ncbi_fcs_flagged",
        "fcs_action",
        "fcs_contam_type",
        "fcs_coverage",
        "fcs_contam_details",
        "preliminary_evidence_status",
        "recommended_next_test",
        "automatic_filter_decision",
    ]
    extra_fields = sorted(
        set().union(*(row.keys() for row in rows)) - set(preferred_fields)
    )
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=preferred_fields + extra_fields,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    with open(args.follow_up, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "scaffold",
                "preliminary_evidence_status",
                "recommended_next_test",
                "priority_reason",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(follow_up_rows)


if __name__ == "__main__":
    main()
