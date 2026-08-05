#!/usr/bin/env Rscript

# Test M. robertsii expression among diets only when fungal RNA coverage and
# biological replication are sufficient for a defensible comparison.

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) {
    stop("Missing argument: ", flag)
  }
  args[[position + 1]]
}

count_file <- arg_value("--counts")
metadata_file <- arg_value("--metadata")
output_dir <- arg_value("--output-dir")
min_assigned <- as.numeric(arg_value("--min-assigned"))
min_genes <- as.numeric(arg_value("--min-genes"))
min_replicates <- as.integer(arg_value("--min-replicates"))
alpha <- as.numeric(arg_value("--alpha"))
lfc_threshold <- as.numeric(arg_value("--lfc-threshold"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_table <- readr::read_tsv(count_file, show_col_types = FALSE)
gene_ids <- count_table$gene_id
count_matrix <- as.matrix(count_table[, setdiff(names(count_table), "gene_id")])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- gene_ids

metadata <- readr::read_tsv(metadata_file, show_col_types = FALSE) |>
  mutate(
    sample_id = as.character(.data$sample_id),
    treatment = as.character(.data$treatment),
    diet = as.character(.data$diet)
  )
common_samples <- intersect(metadata$sample_id, colnames(count_matrix))
metadata <- metadata |>
  filter(.data$sample_id %in% common_samples) |>
  arrange(match(.data$sample_id, colnames(count_matrix)))
count_matrix <- count_matrix[, metadata$sample_id, drop = FALSE]

sample_qc <- tibble(
  sample_id = colnames(count_matrix),
  fungal_assigned_counts = colSums(count_matrix),
  fungal_genes_count_ge_10 = colSums(count_matrix >= 10)
) |>
  left_join(metadata, by = "sample_id") |>
  mutate(
    eligible_for_fungal_de = (
      .data$treatment == "Infected" &
        .data$fungal_assigned_counts >= min_assigned &
        .data$fungal_genes_count_ge_10 >= min_genes
    ),
    exclusion_reason = case_when(
      .data$treatment != "Infected" ~ "Control sample: used for background only",
      .data$fungal_assigned_counts < min_assigned ~ "Insufficient fungal assigned counts",
      .data$fungal_genes_count_ge_10 < min_genes ~ "Insufficient fungal genes detected",
      TRUE ~ "Eligible"
    )
  )
readr::write_tsv(
  sample_qc,
  file.path(output_dir, "fungal_deseq2_sample_eligibility.tsv")
)

eligible_samples <- sample_qc |>
  filter(.data$eligible_for_fungal_de) |>
  pull(.data$sample_id)

empty_results <- tibble(
  gene_id = character(),
  baseMean = numeric(),
  log2FoldChange = numeric(),
  lfcSE = numeric(),
  stat = numeric(),
  pvalue = numeric(),
  padj = numeric(),
  contrast = character(),
  numerator = character(),
  denominator = character(),
  higher_in = character(),
  significant = logical()
)
result_rows <- empty_results
summary_rows <- tibble(
  contrast = character(),
  numerator = character(),
  denominator = character(),
  n_numerator = integer(),
  n_denominator = integer(),
  status = character(),
  upregulated = integer(),
  downregulated = integer()
)

placeholder_plot <- function(title, subtitle) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = subtitle, size = 5) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = title) +
    theme_void(base_size = 14)
}

eligible_metadata <- metadata |>
  filter(.data$sample_id %in% eligible_samples) |>
  mutate(Diet = factor(.data$diet))

supported_diets <- eligible_metadata |>
  count(.data$diet, name = "eligible_replicates") |>
  filter(.data$eligible_replicates >= min_replicates) |>
  pull(.data$diet)
analysis_metadata <- eligible_metadata |>
  filter(.data$diet %in% supported_diets) |>
  mutate(Diet = droplevels(factor(.data$diet)))
analysis_counts <- count_matrix[, analysis_metadata$sample_id, drop = FALSE]

if (nrow(analysis_metadata) >= 2 && n_distinct(analysis_metadata$Diet) >= 2) {
  keep_genes <- rowSums(analysis_counts >= 10) >= min_replicates
  analysis_counts <- analysis_counts[keep_genes, , drop = FALSE]

  if (nrow(analysis_counts) > 0) {
    dds <- DESeqDataSetFromMatrix(
      countData = analysis_counts,
      colData = as.data.frame(analysis_metadata),
      design = ~ Diet
    )
    dds <- tryCatch(
      DESeq(
        dds,
        sfType = "poscounts",
        minReplicatesForReplace = Inf,
        quiet = TRUE
      ),
      error = function(error) {
        message(
          "Standard DESeq2 dispersion fit failed; using gene-wise estimates: ",
          conditionMessage(error)
        )
        fallback <- estimateSizeFactors(dds, type = "poscounts")
        fallback <- estimateDispersionsGeneEst(fallback, quiet = TRUE)
        dispersions(fallback) <- mcols(fallback)$dispGeneEst
        nbinomWaldTest(fallback, quiet = TRUE)
      }
    )

    contrast_definitions <- list(
      c("33", "50"),
      c("83", "50"),
      c("83", "33")
    )
    for (definition in contrast_definitions) {
      numerator <- definition[[1]]
      denominator <- definition[[2]]
      n_numerator <- sum(eligible_metadata$diet == numerator)
      n_denominator <- sum(eligible_metadata$diet == denominator)
      contrast_label <- paste0("Diet ", numerator, " vs ", denominator)

      if (n_numerator < min_replicates || n_denominator < min_replicates) {
        summary_rows <- bind_rows(
          summary_rows,
          tibble(
            contrast = contrast_label,
            numerator = numerator,
            denominator = denominator,
            n_numerator = n_numerator,
            n_denominator = n_denominator,
            status = "Not tested: insufficient eligible replicates",
            upregulated = NA_integer_,
            downregulated = NA_integer_
          )
        )
        next
      }

      result <- results(
        dds,
        contrast = c("Diet", numerator, denominator),
        alpha = alpha
      ) |>
        as.data.frame() |>
        tibble::rownames_to_column("gene_id") |>
        as_tibble() |>
        mutate(
          contrast = contrast_label,
          numerator = numerator,
          denominator = denominator,
          higher_in = if_else(
            .data$log2FoldChange >= 0,
            numerator,
            denominator
          ),
          significant = (
            !is.na(.data$padj) &
              .data$padj < alpha &
              abs(.data$log2FoldChange) >= lfc_threshold
          )
        )
      result_rows <- bind_rows(result_rows, result)
      summary_rows <- bind_rows(
        summary_rows,
        tibble(
          contrast = contrast_label,
          numerator = numerator,
          denominator = denominator,
          n_numerator = n_numerator,
          n_denominator = n_denominator,
          status = "Tested",
          upregulated = sum(
            result$significant & result$log2FoldChange > 0,
            na.rm = TRUE
          ),
          downregulated = sum(
            result$significant & result$log2FoldChange < 0,
            na.rm = TRUE
          )
        )
      )
    }

    pca_result <- tryCatch(
      {
        transformed <- varianceStabilizingTransformation(dds, blind = FALSE)
        data <- plotPCA(
          transformed,
          intgroup = "Diet",
          returnData = TRUE
        )
        list(
          data = data,
          percent = round(100 * attr(data, "percentVar"), 1)
        )
      },
      error = function(error) {
        message(
          "VST PCA failed; using log2 normalized fungal counts: ",
          conditionMessage(error)
        )
        normalized <- counts(dds, normalized = TRUE)
        pca <- prcomp(t(log2(normalized + 1)), center = TRUE, scale. = FALSE)
        percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
        data <- tibble(
          name = rownames(pca$x),
          PC1 = pca$x[, 1],
          PC2 = if (ncol(pca$x) >= 2) pca$x[, 2] else 0
        ) |>
          left_join(
            analysis_metadata |>
              select("sample_id", "Diet"),
            by = c("name" = "sample_id")
          )
        list(data = data, percent = percent)
      }
    )
    pca_data <- pca_result$data
    percent_variance <- pca_result$percent
    pca_plot <- ggplot(
      pca_data,
      aes(x = .data$PC1, y = .data$PC2, color = .data$Diet)
    ) +
      geom_point(size = 4) +
      geom_text(aes(label = .data$name), nudge_y = 0.4, size = 3) +
      labs(
        title = "M. robertsii transcriptome PCA",
        subtitle = "Only infected libraries passing fungal coverage thresholds",
        x = paste0("PC1: ", percent_variance[[1]], "%"),
        y = paste0("PC2: ", percent_variance[[2]], "%"),
        color = "Diet"
      ) +
      theme_bw(base_size = 13)
  } else {
    pca_plot <- placeholder_plot(
      "M. robertsii transcriptome PCA",
      "No fungal genes passed the minimum count filter."
    )
  }
} else {
  pca_plot <- placeholder_plot(
    "M. robertsii transcriptome PCA",
    "Insufficient fungal-positive diets or samples for formal analysis."
  )
}

if (nrow(summary_rows) > 0 && any(summary_rows$status == "Tested")) {
  bar_rows <- summary_rows |>
    filter(.data$status == "Tested") |>
    select("contrast", "upregulated", "downregulated") |>
    tidyr::pivot_longer(
      cols = c("upregulated", "downregulated"),
      names_to = "direction",
      values_to = "n_degs"
    )
  bar_plot <- ggplot(
    bar_rows,
    aes(x = .data$contrast, y = .data$n_degs, fill = .data$direction)
  ) +
    geom_col(position = "dodge") +
    scale_fill_manual(
      values = c(upregulated = "#D73027", downregulated = "#2878B5"),
      labels = c(upregulated = "Higher in numerator", downregulated = "Lower in numerator")
    ) +
    labs(
      title = "M. robertsii differential expression among host diets",
      x = NULL,
      y = "Significant fungal DEGs",
      fill = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
} else {
  bar_plot <- placeholder_plot(
    "M. robertsii differential expression among host diets",
    "No contrast had enough eligible fungal-positive replicates."
  )
}

readr::write_csv(
  result_rows,
  file.path(output_dir, "fungal_deseq2_all_contrasts.csv")
)
readr::write_tsv(
  summary_rows,
  file.path(output_dir, "fungal_deseq2_contrast_summary.tsv")
)
ggsave(
  file.path(output_dir, "fungal_transcriptome_pca.png"),
  pca_plot,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(output_dir, "fungal_deg_counts_by_diet.png"),
  bar_plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)
