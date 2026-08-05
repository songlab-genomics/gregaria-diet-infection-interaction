#!/usr/bin/env python3
"""Summarize host-only, competitive, prior DEG, and sequence-placement evidence."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


CONTRAST_LABELS = {
    "control_diet33_vs50": "Control 33 vs 50",
    "control_diet83_vs33": "Control 83 vs 33",
    "control_diet83_vs50": "Control 83 vs 50",
    "infected_diet33_vs50": "Infected 33 vs 50",
    "infected_diet83_vs33": "Infected 83 vs 33",
    "infected_diet83_vs50": "Infected 83 vs 50",
    "infection_diet33": "Diet 33: infected vs control",
    "infection_diet50": "Diet 50: infected vs control",
    "infection_diet83": "Diet 83: infected vs control",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host-results", required=True)
    parser.add_argument("--competitive-results", required=True)
    parser.add_argument("--previous-catalogue", required=True)
    parser.add_argument("--host-counts", required=True)
    parser.add_argument("--competitive-counts", required=True)
    parser.add_argument("--previous-counts", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--rrna-list", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--assembly-report", required=True)
    parser.add_argument("--sample-count-comparison", required=True)
    parser.add_argument("--host-mapping-summary", required=True)
    parser.add_argument("--competitive-mapping-summary", required=True)
    parser.add_argument("--gene-count-comparison", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def read_result(path: str) -> pd.DataFrame:
    table = pd.read_csv(path, low_memory=False)
    table["significant"] = (
        table["significant"].astype(str).str.lower().isin({"true", "t", "1", "yes"})
    )
    return table


def read_rrna_ids(path: str) -> set[str]:
    with open(path) as handle:
        return {line.strip() for line in handle if line.strip()}


def parse_attributes(value: str) -> dict[str, str]:
    attributes = {}
    for field in value.rstrip().split(";"):
        if "=" not in field:
            continue
        key, item = field.split("=", 1)
        attributes[key] = unquote(item)
    return attributes


def read_gene_annotation(path: str) -> pd.DataFrame:
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            gene_id = (
                attributes.get("gene")
                or attributes.get("Name")
                or attributes.get("ID", "").removeprefix("gene-")
            )
            if not gene_id:
                continue
            rows.append(
                {
                    "gene_id": gene_id,
                    "seqid": fields[0],
                    "gene_start": int(fields[3]),
                    "gene_end": int(fields[4]),
                    "strand": fields[6],
                    "description": attributes.get("description", "no GFF description"),
                    "gene_biotype": attributes.get("gene_biotype", "unknown"),
                }
            )
    return pd.DataFrame(rows).drop_duplicates("gene_id")


def read_assembly_report(path: str) -> pd.DataFrame:
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            role = fields[1]
            location_type = fields[3]
            if role == "assembled-molecule" and location_type == "Chromosome":
                sequence_class = "Nuclear chromosome"
            elif role == "assembled-molecule" and location_type == "Mitochondrion":
                sequence_class = "Mitochondrial"
            elif role == "unlocalized-scaffold":
                sequence_class = "Unlocalized scaffold"
            elif role == "unplaced-scaffold":
                sequence_class = "Unplaced scaffold"
            else:
                sequence_class = role
            rows.append(
                {
                    "sequence_name": fields[0],
                    "sequence_role": role,
                    "assigned_molecule": fields[2],
                    "location_type": location_type,
                    "genbank_accession": fields[4],
                    "seqid": fields[6],
                    "assembly_unit": fields[7],
                    "sequence_length": int(fields[8]),
                    "sequence_class": sequence_class,
                }
            )
    return pd.DataFrame(rows)


def direction_from_lfc(values: pd.Series) -> pd.Series:
    return np.where(values > 0, "Upregulated", "Downregulated")


def compare_degs(
    host: pd.DataFrame,
    competitive: pd.DataFrame,
    previous: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for branch, fresh in [
        ("Fresh host-only", host),
        ("Fresh competitive-host", competitive),
    ]:
        for contrast_id, previous_label in CONTRAST_LABELS.items():
            fresh_sig = fresh[
                (fresh["contrast_id"] == contrast_id) & fresh["significant"]
            ].copy()
            old_sig = previous[previous["comparison"] == previous_label].copy()
            fresh_direction = dict(
                zip(
                    fresh_sig["gene_id"],
                    direction_from_lfc(fresh_sig["log2FoldChange"]),
                )
            )
            old_direction = dict(zip(old_sig["gene_id"], old_sig["signed_direction"]))
            fresh_set = set(fresh_direction)
            old_set = set(old_direction)
            shared = fresh_set & old_set
            union = fresh_set | old_set
            rows.append(
                {
                    "mapping_branch": branch,
                    "contrast_id": contrast_id,
                    "contrast_label": previous_label,
                    "previous_degs": len(old_set),
                    "fresh_degs": len(fresh_set),
                    "fresh_minus_previous": len(fresh_set) - len(old_set),
                    "shared_same_direction": sum(
                        fresh_direction[gene] == old_direction[gene] for gene in shared
                    ),
                    "shared_opposite_direction": sum(
                        fresh_direction[gene] != old_direction[gene] for gene in shared
                    ),
                    "previous_only": len(old_set - fresh_set),
                    "fresh_only": len(fresh_set - old_set),
                    "jaccard": len(shared) / len(union) if union else math.nan,
                }
            )
    return pd.DataFrame(rows)


def summarize_fresh_degs(host: pd.DataFrame, competitive: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for branch, table in [("Host-only", host), ("Competitive-host", competitive)]:
        significant = table[table["significant"]].copy()
        significant["direction"] = direction_from_lfc(significant["log2FoldChange"])
        summary = (
            significant.groupby(["contrast_id", "direction"])
            .size()
            .rename("n_degs")
            .reset_index()
        )
        summary.insert(0, "mapping_branch", branch)
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def compare_count_matrices(
    old_path: str,
    host_path: str,
    competitive_path: str,
    rrna_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    old = pd.read_csv(old_path)
    old = old.rename(columns={old.columns[0]: "gene_id"}).set_index("gene_id")
    old.columns = (
        old.columns.str.removeprefix("mehreen_").str.removesuffix("_MERGE")
    )

    host = pd.read_csv(host_path, sep="\t").set_index("gene_id")
    competitive = pd.read_csv(competitive_path, sep="\t").set_index("gene_id")
    host = host.loc[~host.index.isin(rrna_ids)]
    competitive = competitive.loc[~competitive.index.isin(rrna_ids)]

    common_genes = old.index.intersection(host.index).intersection(competitive.index)
    common_samples = old.columns.intersection(host.columns).intersection(competitive.columns)
    rows = []
    for sample in common_samples:
        old_values = old.loc[common_genes, sample].to_numpy(dtype=float)
        host_values = host.loc[common_genes, sample].to_numpy(dtype=float)
        competitive_values = competitive.loc[common_genes, sample].to_numpy(dtype=float)
        rows.append(
            {
                "sample_id": sample,
                "common_genes": len(common_genes),
                "previous_total_common_genes": old_values.sum(),
                "host_only_total_common_genes": host_values.sum(),
                "competitive_total_common_genes": competitive_values.sum(),
                "host_to_previous_total_ratio": (
                    host_values.sum() / old_values.sum()
                    if old_values.sum()
                    else math.nan
                ),
                "competitive_to_previous_total_ratio": (
                    competitive_values.sum() / old_values.sum()
                    if old_values.sum()
                    else math.nan
                ),
                "previous_vs_host_log1p_pearson": np.corrcoef(
                    np.log1p(old_values), np.log1p(host_values)
                )[0, 1],
                "previous_vs_competitive_log1p_pearson": np.corrcoef(
                    np.log1p(old_values), np.log1p(competitive_values)
                )[0, 1],
                "host_vs_competitive_log1p_pearson": np.corrcoef(
                    np.log1p(host_values), np.log1p(competitive_values)
                )[0, 1],
            }
        )

    matrix_summary = pd.DataFrame(
        [
            {
                "matrix": "Previous",
                "genes": old.shape[0],
                "samples": old.shape[1],
                "total_counts": old.to_numpy(dtype=np.int64).sum(),
            },
            {
                "matrix": "Fresh host-only after rRNA removal",
                "genes": host.shape[0],
                "samples": host.shape[1],
                "total_counts": host.to_numpy(dtype=np.int64).sum(),
            },
            {
                "matrix": "Fresh competitive-host after rRNA removal",
                "genes": competitive.shape[0],
                "samples": competitive.shape[1],
                "total_counts": competitive.to_numpy(dtype=np.int64).sum(),
            },
        ]
    )
    return pd.DataFrame(rows), matrix_summary


def placement_summaries(
    host: pd.DataFrame,
    competitive: pd.DataFrame,
    annotation: pd.DataFrame,
    gene_count_comparison: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    annotated = annotation.set_index("gene_id")
    class_rows = []
    unplaced_rows = []

    for branch, table in [("Host-only", host), ("Competitive-host", competitive)]:
        for contrast_id, contrast in table.groupby("contrast_id", sort=False):
            contrast = contrast.join(
                annotated[["sequence_class"]],
                on="gene_id",
                how="left",
            )
            contrast["sequence_class"] = contrast["sequence_class"].fillna(
                "Missing placement annotation"
            )
            tested_total = len(contrast)
            deg_total = int(contrast["significant"].sum())
            for sequence_class, class_table in contrast.groupby("sequence_class"):
                tested_class = len(class_table)
                deg_class = int(class_table["significant"].sum())
                nonclass_deg = deg_total - deg_class
                nonclass_tested_nondeg = (
                    tested_total - tested_class - nonclass_deg
                )
                table_2x2 = [
                    [deg_class, tested_class - deg_class],
                    [nonclass_deg, nonclass_tested_nondeg],
                ]
                odds_ratio, pvalue = fisher_exact(table_2x2)
                background_fraction = tested_class / tested_total
                deg_fraction = deg_class / deg_total if deg_total else math.nan
                class_rows.append(
                    {
                        "mapping_branch": branch,
                        "contrast_id": contrast_id,
                        "sequence_class": sequence_class,
                        "tested_genes_in_class": tested_class,
                        "degs_in_class": deg_class,
                        "all_tested_genes": tested_total,
                        "all_degs": deg_total,
                        "background_fraction": background_fraction,
                        "deg_fraction": deg_fraction,
                        "enrichment_ratio": (
                            deg_fraction / background_fraction
                            if deg_total and background_fraction
                            else math.nan
                        ),
                        "fisher_odds_ratio": odds_ratio,
                        "fisher_pvalue": pvalue,
                    }
                )

            significant_unplaced = contrast[
                contrast["significant"]
                & (contrast["sequence_class"] == "Unplaced scaffold")
            ]
            for row in significant_unplaced.itertuples(index=False):
                unplaced_rows.append(
                    {
                        "mapping_branch": branch,
                        "contrast_id": contrast_id,
                        "gene_id": row.gene_id,
                        "log2FoldChange": row.log2FoldChange,
                        "padj": row.padj,
                        "direction": (
                            "Upregulated"
                            if row.log2FoldChange > 0
                            else "Downregulated"
                        ),
                    }
                )

    placement = pd.DataFrame(class_rows)
    unplaced = pd.DataFrame(unplaced_rows)
    if unplaced.empty:
        return placement, unplaced, pd.DataFrame(), pd.DataFrame()

    unplaced = (
        unplaced.merge(annotation, on="gene_id", how="left")
        .merge(gene_count_comparison, on="gene_id", how="left")
        .sort_values(["contrast_id", "gene_id", "mapping_branch"])
    )

    support = (
        unplaced.pivot_table(
            index=["contrast_id", "gene_id"],
            columns="mapping_branch",
            values="direction",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(columns=None)
    )
    for branch in ["Host-only", "Competitive-host"]:
        if branch not in support:
            support[branch] = np.nan
    support["mapping_support"] = np.select(
        [
            support["Host-only"].notna() & support["Competitive-host"].notna(),
            support["Host-only"].notna(),
            support["Competitive-host"].notna(),
        ],
        [
            "Significant in both mappings",
            "Host-only significant only",
            "Competitive-host significant only",
        ],
        default="Unknown",
    )

    unplaced = unplaced.merge(
        support[["contrast_id", "gene_id", "mapping_support"]],
        on=["contrast_id", "gene_id"],
        how="left",
    )
    unplaced = unplaced.drop_duplicates(
        ["mapping_branch", "contrast_id", "gene_id"]
    )

    scaffold_summary = (
        unplaced.groupby(["mapping_branch", "contrast_id"])
        .agg(
            unplaced_deg_calls=("gene_id", "nunique"),
            unplaced_scaffolds=("seqid", "nunique"),
        )
        .reset_index()
    )

    gene_summary = (
        unplaced.groupby(
            [
                "gene_id",
                "seqid",
                "sequence_length",
                "description",
                "gene_biotype",
                "mapping_support",
            ],
            dropna=False,
        )
        .agg(
            n_contrasts=("contrast_id", "nunique"),
            contrasts=("contrast_id", lambda x: "; ".join(sorted(set(x)))),
            min_padj=("padj", "min"),
            max_abs_log2FoldChange=("log2FoldChange", lambda x: np.max(np.abs(x))),
            competitive_host_retained_fraction=(
                "competitive_host_retained_fraction",
                "first",
            ),
            host_counts_not_retained=("host_counts_not_retained", "first"),
        )
        .reset_index()
        .sort_values(["n_contrasts", "min_padj"], ascending=[False, True])
    )
    return placement, unplaced, scaffold_summary, gene_summary


def mapping_metric_summary(
    sample_counts: pd.DataFrame,
    host_mapping: pd.DataFrame,
    competitive_mapping: pd.DataFrame,
) -> pd.DataFrame:
    comparable = host_mapping[
        ["sample_id", "input_reads", "uniquely_mapped_percent"]
    ].merge(
        competitive_mapping[
            [
                "sample_id",
                "unique_read_alignments",
                "unique_host_unplaced",
                "unique_fungus",
            ]
        ],
        on="sample_id",
        how="inner",
    )
    # The competitive BAM table counts both mates; divide alignments by two.
    comparable["competitive_unique_percent_of_input"] = (
        100
        * comparable["unique_read_alignments"]
        / (2 * comparable["input_reads"])
    )
    comparable["competitive_minus_host_unique_percentage_points"] = (
        comparable["competitive_unique_percent_of_input"]
        - comparable["uniquely_mapped_percent"]
    )
    comparable["competitive_unplaced_percent_of_input"] = (
        100 * comparable["unique_host_unplaced"] / (2 * comparable["input_reads"])
    )
    comparable["competitive_fungus_percent_of_input"] = (
        100 * comparable["unique_fungus"] / (2 * comparable["input_reads"])
    )

    metrics = {
        "Host counts retained after competitive mapping": sample_counts[
            "competitive_host_retained_fraction"
        ],
        "Host counts not retained": sample_counts["host_counts_not_retained"],
        "Competitive fungal assigned counts": sample_counts[
            "competitive_fungus_assigned_counts"
        ],
        "Host-only uniquely mapped percent": host_mapping[
            "uniquely_mapped_percent"
        ],
        "Competitive uniquely mapped percent of input": comparable[
            "competitive_unique_percent_of_input"
        ],
        "Competitive minus host-only unique mapping percentage points": comparable[
            "competitive_minus_host_unique_percentage_points"
        ],
        "Competitive host-unplaced unique alignments percent of input": comparable[
            "competitive_unplaced_percent_of_input"
        ],
        "Competitive fungal unique alignments percent of input": comparable[
            "competitive_fungus_percent_of_input"
        ],
        "Competitive unique alignments on host chromosomes percent": competitive_mapping[
            "pct_unique_host_chromosome"
        ],
        "Competitive unique alignments on host unplaced scaffolds percent": competitive_mapping[
            "pct_unique_host_unplaced"
        ],
        "Competitive unique alignments on fungus percent": competitive_mapping[
            "pct_unique_fungus"
        ],
    }
    rows = []
    for metric, values in metrics.items():
        values = pd.to_numeric(values, errors="coerce")
        rows.append(
            {
                "metric": metric,
                "minimum": values.min(),
                "median": values.median(),
                "mean": values.mean(),
                "maximum": values.max(),
            }
        )
    return pd.DataFrame(rows)


def summarize_dominant_unplaced_expression(
    unplaced: pd.DataFrame,
    host_counts_path: str,
    metadata_path: str,
    rrna_ids: set[str],
) -> pd.DataFrame:
    host_unplaced = unplaced[unplaced["mapping_branch"] == "Host-only"]
    if host_unplaced.empty:
        return pd.DataFrame()
    dominant_contrast = (
        host_unplaced.groupby("contrast_id")["gene_id"].nunique().idxmax()
    )
    genes = set(
        host_unplaced.loc[
            host_unplaced["contrast_id"] == dominant_contrast,
            "gene_id",
        ]
    )
    counts = pd.read_csv(host_counts_path, sep="\t").set_index("gene_id")
    counts = counts.loc[~counts.index.isin(rrna_ids)]
    cluster_counts = counts.loc[counts.index.intersection(genes)].sum(axis=0)
    total_counts = counts.sum(axis=0)
    burden = pd.DataFrame(
        {
            "sample_id": counts.columns.astype(str),
            "dominant_unplaced_contrast": dominant_contrast,
            "unplaced_deg_genes": len(genes),
            "unplaced_cluster_counts": cluster_counts.to_numpy(),
            "total_host_non_rrna_counts": total_counts.to_numpy(),
            "unplaced_cluster_cpm": (
                1e6 * cluster_counts.to_numpy() / total_counts.to_numpy()
            ),
        }
    )
    metadata = pd.read_csv(
        metadata_path,
        sep="\t",
        dtype={"sample_id": str, "diet": str},
    )
    return (
        metadata.merge(burden, on="sample_id", how="right")
        .sort_values("unplaced_cluster_cpm", ascending=False)
        .reset_index(drop=True)
    )


def make_deg_count_plot(comparison: pd.DataFrame, output_path: Path) -> None:
    if plt is None:
        return
    plot_table = comparison[
        [
            "mapping_branch",
            "contrast_id",
            "contrast_label",
            "previous_degs",
            "fresh_degs",
        ]
    ].copy()
    previous = (
        plot_table[["contrast_id", "contrast_label", "previous_degs"]]
        .drop_duplicates()
        .rename(columns={"previous_degs": "n_degs"})
    )
    previous["analysis"] = "Previous"
    fresh = plot_table.rename(
        columns={"mapping_branch": "analysis", "fresh_degs": "n_degs"}
    )[["contrast_id", "contrast_label", "analysis", "n_degs"]]
    plot_table = pd.concat([previous, fresh], ignore_index=True)
    order = list(CONTRAST_LABELS)
    plot_table["contrast_id"] = pd.Categorical(
        plot_table["contrast_id"], categories=order, ordered=True
    )
    plot_table = plot_table.sort_values(["contrast_id", "analysis"])

    analyses = ["Previous", "Fresh host-only", "Fresh competitive-host"]
    colors = ["#B8B1A8", "#C85A54", "#315C8A"]
    x = np.arange(len(order))
    width = 0.25
    fig, axis = plt.subplots(figsize=(14, 7))
    for index, (analysis, color) in enumerate(zip(analyses, colors)):
        values = (
            plot_table[plot_table["analysis"] == analysis]
            .set_index("contrast_id")
            .reindex(order)["n_degs"]
            .fillna(0)
        )
        axis.bar(x + (index - 1) * width, values, width, label=analysis, color=color)
    axis.set_xticks(x)
    axis.set_xticklabels(
        [CONTRAST_LABELS[item] for item in order],
        rotation=35,
        ha="right",
    )
    axis.set_ylabel("Significant DEGs")
    axis.set_title("Previous and fresh DEG burdens using matched DESeq2 models")
    axis.legend(frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    host = read_result(args.host_results)
    competitive = read_result(args.competitive_results)
    previous = pd.read_csv(args.previous_catalogue, low_memory=False)
    rrna_ids = read_rrna_ids(args.rrna_list)

    comparison = compare_degs(host, competitive, previous)
    comparison.to_csv(
        output_dir / "previous_vs_fresh_deg_comparison.tsv",
        sep="\t",
        index=False,
    )
    summarize_fresh_degs(host, competitive).to_csv(
        output_dir / "fresh_deg_counts_by_direction.tsv",
        sep="\t",
        index=False,
    )

    count_metrics, matrix_summary = compare_count_matrices(
        args.previous_counts,
        args.host_counts,
        args.competitive_counts,
        rrna_ids,
    )
    count_metrics.to_csv(
        output_dir / "previous_vs_fresh_sample_count_metrics.tsv",
        sep="\t",
        index=False,
    )
    matrix_summary.to_csv(
        output_dir / "count_matrix_summary.tsv",
        sep="\t",
        index=False,
    )

    assembly = read_assembly_report(args.assembly_report)
    annotation = read_gene_annotation(args.gff).merge(
        assembly,
        on="seqid",
        how="left",
    )
    annotation["sequence_class"] = annotation["sequence_class"].fillna(
        "Missing placement annotation"
    )
    annotation.to_csv(
        output_dir / "host_gene_sequence_placement.tsv",
        sep="\t",
        index=False,
    )

    gene_count_comparison = pd.read_csv(args.gene_count_comparison, sep="\t")
    placement, unplaced, scaffold_summary, gene_summary = placement_summaries(
        host,
        competitive,
        annotation,
        gene_count_comparison,
    )
    placement.to_csv(
        output_dir / "deg_sequence_class_enrichment.tsv",
        sep="\t",
        index=False,
    )
    unplaced.to_csv(
        output_dir / "unplaced_scaffold_deg_catalogue.tsv",
        sep="\t",
        index=False,
    )
    if not unplaced.empty:
        scaffold_detail = (
            unplaced.assign(
                upregulated=unplaced["direction"].eq("Upregulated").astype(int),
                downregulated=unplaced["direction"].eq("Downregulated").astype(int),
            )
            .groupby(
                [
                    "mapping_branch",
                    "contrast_id",
                    "seqid",
                    "sequence_length",
                ],
                dropna=False,
            )
            .agg(
                unplaced_deg_genes=("gene_id", "nunique"),
                upregulated=("upregulated", "sum"),
                downregulated=("downregulated", "sum"),
                median_competitive_host_retained_fraction=(
                    "competitive_host_retained_fraction",
                    "median",
                ),
            )
            .reset_index()
            .sort_values(
                ["mapping_branch", "contrast_id", "unplaced_deg_genes"],
                ascending=[True, True, False],
            )
        )
    else:
        scaffold_detail = pd.DataFrame()
    scaffold_detail.to_csv(
        output_dir / "unplaced_scaffold_deg_by_scaffold.tsv",
        sep="\t",
        index=False,
    )
    scaffold_summary.to_csv(
        output_dir / "unplaced_scaffold_deg_summary.tsv",
        sep="\t",
        index=False,
    )
    gene_summary.to_csv(
        output_dir / "unplaced_scaffold_unique_deg_genes.tsv",
        sep="\t",
        index=False,
    )
    summarize_dominant_unplaced_expression(
        unplaced,
        args.host_counts,
        args.metadata,
        rrna_ids,
    ).to_csv(
        output_dir / "dominant_unplaced_deg_sample_burden.tsv",
        sep="\t",
        index=False,
    )

    sample_counts = pd.read_csv(args.sample_count_comparison, sep="\t")
    host_mapping = pd.read_csv(args.host_mapping_summary, sep="\t")
    competitive_mapping = pd.read_csv(args.competitive_mapping_summary, sep="\t")
    mapping_metric_summary(
        sample_counts,
        host_mapping,
        competitive_mapping,
    ).to_csv(
        output_dir / "mapping_metric_summary.tsv",
        sep="\t",
        index=False,
    )

    make_deg_count_plot(
        comparison,
        output_dir / "previous_vs_fresh_deg_counts.png",
    )

    with open(output_dir / "README.txt", "w") as handle:
        handle.write(
            "Host mapping audit\n"
            "==================\n"
            "All samples were retained. Annotated rRNA loci were removed.\n"
            "DEGs use padj < 0.05 and absolute log2 fold change >= 1.\n"
            "Fresh legacy-matched models reproduce the workflowR subset strategy.\n"
            "Sequence roles come from the official NCBI assembly report.\n"
            "Unplaced loci are retained and reported, not interpreted as contaminants by default.\n"
        )


if __name__ == "__main__":
    main()
