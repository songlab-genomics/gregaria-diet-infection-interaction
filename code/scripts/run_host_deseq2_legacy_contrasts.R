#!/usr/bin/env Rscript

# Reproduce the all-samples contrast strategy used by the workflowR site.
suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) {
    stop("Missing required argument: ", flag)
  }
  args[[index + 1]]
}

optional_arg_value <- function(flag, default) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) {
    return(default)
  }
  args[[index + 1]]
}

counts_path <- arg_value("--counts")
metadata_path <- arg_value("--metadata")
rrna_path <- arg_value("--rrna-list")
mapping_branch <- arg_value("--mapping-branch")
output_dir <- arg_value("--output-dir")
alpha <- as.numeric(arg_value("--alpha"))
lfc_threshold <- as.numeric(arg_value("--lfc-threshold"))
min_total_count <- as.integer(arg_value("--min-total-count"))
outlier_mode <- optional_arg_value("--outlier-mode", "legacy-default")
if (!outlier_mode %in% c("legacy-default", "retain-all-counts")) {
  stop("--outlier-mode must be legacy-default or retain-all-counts.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts_table <- read.delim(
  counts_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!"gene_id" %in% names(counts_table)) {
  stop("Count matrix must contain a gene_id column.")
}

counts <- as.matrix(counts_table[, setdiff(names(counts_table), "gene_id")])
rownames(counts) <- counts_table$gene_id
storage.mode(counts) <- "integer"

metadata <- read.delim(
  metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_metadata <- c("sample_id", "treatment", "diet")
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0) {
  stop("Metadata is missing columns: ", paste(missing_metadata, collapse = ", "))
}

metadata$sample_id <- as.character(metadata$sample_id)
metadata$treatment <- factor(metadata$treatment, levels = c("Control", "Infected"))
metadata$diet <- factor(as.character(metadata$diet), levels = c("33", "50", "83"))
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
if (anyNA(metadata$sample_id)) {
  stop("Count matrix sample columns do not all occur in metadata.")
}
rownames(metadata) <- metadata$sample_id

rrna_ids <- unique(trimws(readLines(rrna_path, warn = FALSE)))
rrna_ids <- rrna_ids[nzchar(rrna_ids)]
counts <- counts[!rownames(counts) %in% rrna_ids, , drop = FALSE]

format_result <- function(result, contrast_id, contrast_family, numerator,
                          denominator, diet = NA_character_) {
  table <- as.data.frame(result)
  table$gene_id <- rownames(table)
  table$mapping_branch <- mapping_branch
  table$contrast_id <- contrast_id
  table$contrast_family <- contrast_family
  table$numerator <- numerator
  table$denominator <- denominator
  table$diet <- diet
  table$significant <- (
    !is.na(table$padj) &
      table$padj < alpha &
      abs(table$log2FoldChange) >= lfc_threshold
  )
  table$direction <- ifelse(
    !table$significant,
    "Not DEG",
    ifelse(table$log2FoldChange > 0, "Upregulated", "Downregulated")
  )
  table
}

run_infection_within_diet <- function(diet_level) {
  sample_keep <- metadata$diet == diet_level
  metadata_sub <- droplevels(metadata[sample_keep, , drop = FALSE])
  counts_sub <- counts[, sample_keep, drop = FALSE]
  gene_keep <- rowSums(counts_sub) >= min_total_count

  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub[gene_keep, , drop = FALSE],
    colData = metadata_sub,
    design = ~ treatment
  )
  if (outlier_mode == "retain-all-counts") {
    dds <- DESeq(dds, minReplicatesForReplace = Inf, quiet = TRUE)
    result <- results(
      dds,
      contrast = c("treatment", "Infected", "Control"),
      cooksCutoff = FALSE
    )
  } else {
    dds <- DESeq(dds, quiet = TRUE)
    result <- results(dds, contrast = c("treatment", "Infected", "Control"))
  }

  format_result(
    result,
    paste0("infection_diet", diet_level),
    "infection_within_diet",
    "Infected",
    "Control",
    diet_level
  )
}

run_diet_pairwise <- function(treatment_level) {
  sample_keep <- metadata$treatment == treatment_level
  metadata_sub <- droplevels(metadata[sample_keep, , drop = FALSE])
  counts_sub <- counts[, sample_keep, drop = FALSE]
  gene_keep <- rowSums(counts_sub) >= min_total_count

  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub[gene_keep, , drop = FALSE],
    colData = metadata_sub,
    design = ~ diet
  )
  if (outlier_mode == "retain-all-counts") {
    dds <- DESeq(dds, minReplicatesForReplace = Inf, quiet = TRUE)
  } else {
    dds <- DESeq(dds, quiet = TRUE)
  }

  contrasts <- list(
    diet33_vs50 = c("diet", "33", "50"),
    diet83_vs33 = c("diet", "83", "33"),
    diet83_vs50 = c("diet", "83", "50")
  )

  do.call(rbind, lapply(names(contrasts), function(contrast_name) {
    contrast <- contrasts[[contrast_name]]
    result <- if (outlier_mode == "retain-all-counts") {
      results(dds, contrast = contrast, cooksCutoff = FALSE)
    } else {
      results(dds, contrast = contrast)
    }
    format_result(
      result,
      paste0(tolower(treatment_level), "_", contrast_name),
      paste0("diet_among_", tolower(treatment_level)),
      contrast[[2]],
      contrast[[3]]
    )
  }))
}

all_results <- do.call(
  rbind,
  c(
    lapply(c("33", "50", "83"), run_infection_within_diet),
    list(run_diet_pairwise("Control")),
    list(run_diet_pairwise("Infected"))
  )
)
rownames(all_results) <- NULL

write.csv(
  all_results,
  file.path(output_dir, "legacy_matched_all_contrasts.csv"),
  row.names = FALSE
)

significant <- all_results[all_results$significant, , drop = FALSE]
summary_table <- as.data.frame(
  xtabs(
    ~ contrast_id + direction,
    data = significant,
    drop.unused.levels = TRUE
  )
)
names(summary_table)[names(summary_table) == "Freq"] <- "n_degs"
summary_table <- summary_table[summary_table$n_degs > 0, , drop = FALSE]

write.table(
  summary_table,
  file.path(output_dir, "legacy_matched_contrast_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

audit <- data.frame(
  mapping_branch = mapping_branch,
  samples = ncol(counts),
  genes_after_rrna_removal = nrow(counts),
  rrna_ids_removed = sum(rrna_ids %in% counts_table$gene_id),
  alpha = alpha,
  lfc_threshold = lfc_threshold,
  min_total_count_per_subset = min_total_count,
  model_strategy = "separate DESeq2 model for each diet or treatment subset",
  outlier_mode = outlier_mode,
  min_replicates_for_replace = ifelse(
    outlier_mode == "retain-all-counts",
    "Inf",
    "DESeq2 default"
  ),
  cooks_cutoff = ifelse(
    outlier_mode == "retain-all-counts",
    "FALSE",
    "DESeq2 default"
  ),
  independent_filtering = "DESeq2 default",
  stringsAsFactors = FALSE
)
write.table(
  audit,
  file.path(output_dir, "legacy_matched_analysis_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
