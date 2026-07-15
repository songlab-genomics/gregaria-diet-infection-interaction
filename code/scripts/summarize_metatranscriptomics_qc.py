#!/usr/bin/env python3
"""Summarize host RNA-seq QC and non-host metatranscriptomic signals.

Inputs are a Snakemake manifest plus per-sample QC/classification files. The
output tables are designed for downstream workflowR summaries across species,
tissues, and individuals.
"""

import argparse
import csv
import json
import re
from pathlib import Path


GROUP_TAXIDS = {
    "Bacteria": "2",
    "Archaea": "2157",
    "Eukaryota": "2759",
    "Fungi": "4751",
    "Viruses": "10239",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Sample manifest TSV from Snakemake")
    parser.add_argument("--sample-qc", required=True, help="Output per-sample QC TSV")
    parser.add_argument("--taxa-long", required=True, help="Output long taxonomy table TSV")
    parser.add_argument("--organism-summary", required=True, help="Output domain/group summary TSV")
    parser.add_argument("--top-taxa", required=True, help="Output top taxa per sample TSV")
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


def to_number(value):
    value = str(value).strip().replace("%", "")
    if value in {"", "NA", "nan", "None"}:
        return None
    try:
        if re.search(r"[.eE]", value):
            return float(value)
        return int(value)
    except ValueError:
        return None


def pct(numerator, denominator):
    if numerator is None or denominator in {None, 0}:
        return ""
    return 100.0 * float(numerator) / float(denominator)


def rpm(count, denominator):
    if count is None or denominator in {None, 0}:
        return ""
    return 1_000_000.0 * float(count) / float(denominator)


def read_fastp(path):
    if not path or not Path(path).exists():
        return {}
    try:
        with open(path) as handle:
            data = json.load(handle)
    except Exception:
        return {}
    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    return {
        "raw_reads_fastp": before.get("total_reads", ""),
        "raw_bases_fastp": before.get("total_bases", ""),
        "raw_q30_rate_fastp": before.get("q30_rate", ""),
        "trimmed_reads_fastp": after.get("total_reads", ""),
        "trimmed_bases_fastp": after.get("total_bases", ""),
        "trimmed_q30_rate_fastp": after.get("q30_rate", ""),
        "fastp_read_pairs_after_filtering": (
            after.get("total_reads", 0) / 2 if isinstance(after.get("total_reads"), (int, float)) else ""
        ),
    }


def read_star_log(path):
    out = {}
    if not path or not Path(path).exists():
        return out
    with open(path, errors="ignore") as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = [x.strip() for x in line.split("|", 1)]
            out[key] = value
    wanted = {
        "Number of input reads": "star_input_reads",
        "Uniquely mapped reads number": "star_unique_mapped_reads",
        "Uniquely mapped reads %": "star_unique_mapped_pct",
        "Number of reads mapped to multiple loci": "star_multimapped_reads",
        "% of reads mapped to multiple loci": "star_multimapped_pct",
        "Number of reads mapped to too many loci": "star_too_many_loci_reads",
        "% of reads mapped to too many loci": "star_too_many_loci_pct",
        "% of reads unmapped: too many mismatches": "star_unmapped_mismatch_pct",
        "% of reads unmapped: too short": "star_unmapped_tooshort_pct",
        "% of reads unmapped: other": "star_unmapped_other_pct",
    }
    return {dest: out.get(src, "") for src, dest in wanted.items()}


def read_featurecounts_summary(path):
    rows = read_tsv(path)
    if not rows:
        return {}
    counts = {}
    for row in rows:
        status = row.get("Status", "")
        values = [to_number(v) for k, v in row.items() if k != "Status"]
        values = [v for v in values if v is not None]
        counts[status] = sum(values)
    total = sum(counts.values()) if counts else None
    assigned = counts.get("Assigned")
    out = {
        "featurecounts_total_fragments": total if total is not None else "",
        "featurecounts_assigned_fragments": assigned if assigned is not None else "",
        "featurecounts_assigned_pct": pct(assigned, total),
    }
    for key, value in counts.items():
        out[f"featurecounts_{key}"] = value
    return out


def read_two_column_counts(path):
    counts = {}
    if not path or not Path(path).exists():
        return counts
    with open(path, errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            value = to_number(parts[-1])
            if value is not None:
                counts[parts[0]] = value
    return counts


def read_rrna_stats(counts_path, rrna_path):
    counts = read_two_column_counts(counts_path)
    if not counts:
        return {"rrna_count": "", "rrna_pct_of_assigned_gene_counts": ""}
    rrna_loci = set()
    if rrna_path and Path(rrna_path).exists():
        rrna_loci = {line.strip() for line in open(rrna_path) if line.strip()}
    total = sum(counts.values())
    rrna_count = sum(count for gene, count in counts.items() if gene in rrna_loci)
    return {
        "rrna_count": rrna_count,
        "rrna_pct_of_assigned_gene_counts": pct(rrna_count, total),
        "rrna_loci_in_reference_list": len(rrna_loci),
    }


def read_flagstat(path):
    out = {}
    if not path or not Path(path).exists():
        return out
    patterns = {
        "samtools_total_reads": r"^(\d+) \+ \d+ in total",
        "samtools_secondary": r"^(\d+) \+ \d+ secondary",
        "samtools_supplementary": r"^(\d+) \+ \d+ supplementary",
        "samtools_duplicates_flagged": r"^(\d+) \+ \d+ duplicates",
        "samtools_mapped_reads": r"^(\d+) \+ \d+ mapped",
        "samtools_properly_paired": r"^(\d+) \+ \d+ properly paired",
    }
    text = Path(path).read_text(errors="ignore").splitlines()
    for line in text:
        for key, pattern in patterns.items():
            m = re.search(pattern, line)
            if m:
                out[key] = int(m.group(1))
    out["samtools_mapped_pct"] = pct(out.get("samtools_mapped_reads"), out.get("samtools_total_reads"))
    out["samtools_duplicate_flagged_pct"] = pct(out.get("samtools_duplicates_flagged"), out.get("samtools_total_reads"))
    return out


def read_markdup_metrics(path):
    out = {}
    if not path or not Path(path).exists():
        return out
    for line in Path(path).read_text(errors="ignore").splitlines():
        if ":" not in line:
            continue
        key, value = [x.strip() for x in line.split(":", 1)]
        clean = key.lower().replace(" ", "_").replace("-", "_")
        parsed = to_number(value)
        out[f"markdup_{clean}"] = parsed if parsed is not None else value
    total_dup = out.get("markdup_duplicate_total")
    examined = out.get("markdup_examined")
    out["markdup_duplicate_pct"] = pct(total_dup, examined)
    return out


def read_unmapped_stats(path):
    rows = read_tsv(path)
    if not rows:
        return {}
    row = rows[0]
    return {
        "host_unmapped_read_pairs": to_number(row.get("Host_unmapped_read_pairs")) or "",
        "host_unmapped_reads_estimate": (to_number(row.get("Host_unmapped_read_pairs")) or 0) * 2,
    }


def read_kraken_report(path):
    rows = []
    if not path or not Path(path).exists():
        return rows
    with open(path, errors="ignore") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rows.append({
                "percent": to_number(parts[0]),
                "clade_reads": to_number(parts[1]) or 0,
                "taxon_reads": to_number(parts[2]) or 0,
                "rank": parts[3].strip(),
                "taxid": parts[4].strip(),
                "taxon_name": parts[5].strip(),
            })
    return rows


def metadata_lookup(rows):
    lookup = {}
    for row in rows:
        species = row.get("Species", "")
        sample = row.get("SampleName") or row.get("Sample") or row.get("sample") or row.get("locust") or ""
        tissue = row.get("Tissue", "")
        if sample:
            lookup[(species, sample)] = row
            if tissue:
                lookup[(species, tissue, sample)] = row
    return lookup


def enrich_with_metadata(base, meta):
    for key in ["TimePoint", "PhaseTransition", "Treatment", "Condition", "Tissue", "Species"]:
        if key in meta and meta[key] not in {"", None}:
            base[key] = meta[key]
    return base


def main():
    args = parse_args()
    manifest = read_tsv(args.manifest)
    metadata_rows = []
    for row in manifest:
        metadata_rows.extend(read_tsv(row.get("metadata_file", "")))
    meta_by_sample = metadata_lookup(metadata_rows)

    sample_rows = []
    taxa_long = []
    organism_rows = []
    top_rows = []

    for row in manifest:
        species = row["species"]
        tissue = row["tissue"]
        sample = row["sample"]
        meta = meta_by_sample.get((species, tissue, sample), meta_by_sample.get((species, sample), {}))

        qc = {
            "Species": species,
            "Tissue": tissue,
            "SampleName": sample,
            "metadata_file": row.get("metadata_file", ""),
        }
        enrich_with_metadata(qc, meta)
        qc.update(read_fastp(row.get("fastp_json")))
        qc.update(read_star_log(row.get("star_log")))
        qc.update(read_featurecounts_summary(row.get("featurecounts_summary")))
        qc.update(read_rrna_stats(row.get("featurecounts_counts"), row.get("rrna_list")))
        qc.update(read_flagstat(row.get("flagstat")))
        qc.update(read_markdup_metrics(row.get("markdup_metrics")))
        qc.update(read_unmapped_stats(row.get("host_unmapped_stats")))

        denominator_pairs = to_number(qc.get("star_input_reads")) or to_number(qc.get("fastp_read_pairs_after_filtering"))
        nonhost_pairs = to_number(qc.get("host_unmapped_read_pairs"))
        qc["host_unmapped_pct_of_star_input"] = pct(nonhost_pairs, denominator_pairs)

        kraken_rows = read_kraken_report(row.get("kraken2_report"))
        group_by_taxid = {k["taxid"]: k for k in kraken_rows}
        for group, taxid in GROUP_TAXIDS.items():
            hit = group_by_taxid.get(taxid, {})
            count = hit.get("clade_reads", 0)
            organism_rows.append({
                "Species": species,
                "Tissue": tissue,
                "SampleName": sample,
                "Organism_group": group,
                "Kraken2_clade_reads": count,
                "RPM_total_input_pairs": rpm(count, denominator_pairs),
                "RPM_host_unmapped_pairs": rpm(count, nonhost_pairs),
            })
            qc[f"kraken2_{group.lower()}_clade_reads"] = count
            qc[f"kraken2_{group.lower()}_rpm_total_input_pairs"] = rpm(count, denominator_pairs)
            qc[f"kraken2_{group.lower()}_rpm_host_unmapped_pairs"] = rpm(count, nonhost_pairs)

        for tax in kraken_rows:
            count = tax.get("clade_reads", 0)
            if count <= 0:
                continue
            taxa_long.append({
                "Species": species,
                "Tissue": tissue,
                "SampleName": sample,
                "Classifier": "kraken2",
                "Rank": tax.get("rank", ""),
                "TaxID": tax.get("taxid", ""),
                "Taxon": tax.get("taxon_name", ""),
                "Clade_reads": count,
                "Taxon_reads": tax.get("taxon_reads", 0),
                "RPM_total_input_pairs": rpm(count, denominator_pairs),
                "RPM_host_unmapped_pairs": rpm(count, nonhost_pairs),
            })

        ranked = sorted(
            [tax for tax in kraken_rows if tax.get("rank") in {"S", "S1", "G", "G1"} and tax.get("clade_reads", 0) > 0],
            key=lambda tax: tax.get("clade_reads", 0),
            reverse=True,
        )[:10]
        for rank_i, tax in enumerate(ranked, start=1):
            top_rows.append({
                "Species": species,
                "Tissue": tissue,
                "SampleName": sample,
                "Rank_order": rank_i,
                "Classifier": "kraken2",
                "Taxonomic_rank": tax.get("rank", ""),
                "TaxID": tax.get("taxid", ""),
                "Taxon": tax.get("taxon_name", ""),
                "Clade_reads": tax.get("clade_reads", 0),
                "RPM_total_input_pairs": rpm(tax.get("clade_reads", 0), denominator_pairs),
                "RPM_host_unmapped_pairs": rpm(tax.get("clade_reads", 0), nonhost_pairs),
            })

        sample_rows.append(qc)

    sample_fields = sorted({key for row in sample_rows for key in row.keys()})
    lead = ["Species", "Tissue", "SampleName", "PhaseTransition", "TimePoint", "Treatment", "Condition"]
    sample_fields = [x for x in lead if x in sample_fields] + [x for x in sample_fields if x not in lead]
    write_tsv(args.sample_qc, sample_fields, sample_rows)
    write_tsv(args.taxa_long, [
        "Species", "Tissue", "SampleName", "Classifier", "Rank", "TaxID", "Taxon",
        "Clade_reads", "Taxon_reads", "RPM_total_input_pairs", "RPM_host_unmapped_pairs",
    ], taxa_long)
    write_tsv(args.organism_summary, [
        "Species", "Tissue", "SampleName", "Organism_group", "Kraken2_clade_reads",
        "RPM_total_input_pairs", "RPM_host_unmapped_pairs",
    ], organism_rows)
    write_tsv(args.top_taxa, [
        "Species", "Tissue", "SampleName", "Rank_order", "Classifier", "Taxonomic_rank",
        "TaxID", "Taxon", "Clade_reads", "RPM_total_input_pairs", "RPM_host_unmapped_pairs",
    ], top_rows)


if __name__ == "__main__":
    main()
