#!/usr/bin/env python3
"""Join Kraken sequence classifications to unplaced-scaffold properties."""

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--kraken-output", required=True)
    parser.add_argument("--kraken-report", required=True)
    parser.add_argument("--nodes", default="", help="Optional NCBI taxonomy nodes.dmp.")
    parser.add_argument("--names", default="", help="Optional NCBI taxonomy names.dmp.")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_report_names(path):
    names = {}
    with open(path, errors="ignore") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            names[fields[4].strip()] = fields[5].strip()
    return names


def read_sequence_calls(path):
    calls = {}
    with open(path, errors="ignore") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            calls[fields[1]] = {
                "kraken_status": "Classified" if fields[0] == "C" else "Unclassified",
                "kraken_taxid": fields[2],
            }
    return calls


def read_taxonomy(nodes_path, names_path):
    parents = {}
    names = {}
    if not nodes_path or not names_path:
        return parents, names
    with open(nodes_path, errors="ignore") as handle:
        for line in handle:
            fields = [value.strip() for value in line.split("|")]
            if len(fields) >= 2:
                parents[fields[0]] = fields[1]
    with open(names_path, errors="ignore") as handle:
        for line in handle:
            fields = [value.strip() for value in line.split("|")]
            if len(fields) >= 4 and fields[3] == "scientific name":
                names[fields[0]] = fields[1]
    return parents, names


def lineage_for_taxid(taxid, parents, scientific_names):
    if taxid in {"", "0"} or not parents:
        return []
    lineage = []
    seen = set()
    current = taxid
    while current and current not in seen:
        seen.add(current)
        lineage.append(scientific_names.get(current, current))
        parent = parents.get(current)
        if not parent or parent == current:
            break
        current = parent
    return list(reversed(lineage))


def broad_group(lineage):
    names = set(lineage)
    if "Arthropoda" in names:
        return "Arthropoda"
    if "Metazoa" in names:
        return "Metazoa_non_arthropod"
    if "Fungi" in names:
        return "Fungi"
    if "Bacteria" in names:
        return "Bacteria"
    if "Archaea" in names:
        return "Archaea"
    if "Viruses" in names:
        return "Viruses"
    if "Eukaryota" in names:
        return "Other_eukaryote"
    return "Unresolved"


def main():
    args = parse_args()
    names = read_report_names(args.kraken_report)
    calls = read_sequence_calls(args.kraken_output)
    parents, scientific_names = read_taxonomy(args.nodes, args.names)
    rows = []
    with open(args.metrics, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            call = calls.get(
                row["scaffold"],
                {"kraken_status": "Unclassified", "kraken_taxid": "0"},
            )
            row.update(call)
            row["kraken_taxon_name"] = names.get(call["kraken_taxid"], "unclassified")
            lineage = lineage_for_taxid(
                call["kraken_taxid"], parents, scientific_names
            )
            row["kraken_lineage"] = ";".join(lineage)
            row["kraken_broad_group"] = (
                broad_group(lineage)
                if call["kraken_status"] == "Classified"
                else "Unclassified"
            )
            rows.append(row)

    fieldnames = [
        "scaffold",
        "length_bp",
        "gc_percent",
        "n_percent",
        "kraken_status",
        "kraken_taxid",
        "kraken_taxon_name",
        "kraken_lineage",
        "kraken_broad_group",
    ]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
