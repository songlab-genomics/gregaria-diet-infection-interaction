#!/usr/bin/env Rscript

# Apply one identical all-samples DESeq2 analysis to either host mapping branch.
# Samples are never excluded or replaced; annotated rRNA loci are removed first.

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(tibble)
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
rrna_file <- arg_value("--rrna-list")
mapping_branch <- arg_value("--mapping-branch")
output_dir <- arg_value("--output-dir")
alpha <- as.numeric(arg_value("--alpha"))
lfc_threshold <- as.numeric(arg_value("--lfc-threshold"))
min_count <- as.integer(arg_value("--min-count"))
min_samples <- as.integer(arg_value("--min-samples"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_table <- readr::read_tsv(count_file, show_col_types = FALSE)
if (!"gene_id" %in% names(count_table)) {
  stop("Count matrix must contain a gene_id column: ", count_file)
}
gene_ids <- as.character(count_table$gene_id)
count_matrix <- as.matrix(count_table[, setdiff(names(count_table), "gene_id")])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- gene_ids

metadata <- readr::read_tsv(metadata_file, show_col_types = FALSE) |>
  transmute(
    sample_id = as.character(.data$sample_id),
    diet = as.character(.data$diet),
    treatment = as.character(.data$treatment)
  )

if (!setequal(metadata$sample_id, colnames(count_matrix))) {
  missing_metadata <- setdiff(colnames(count_matrix), metadata$sample_id)
  missing_counts <- setdiff(metadata$sample_id, colnames(count_matrix))
  stop(
    "Metadata/count sample mismatch. Missing metadata: ",
    paste(missing_metadata, collapse = ","),
    "; missing counts: ",
    paste(missing_counts, collapse = ",")
  )
}
metadata <- metadata |>
  arrange(match(.data$sample_id, colnames(count_matrix)))
count_matrix <- count_matrix[, metadata$sample_id, drop = FALSE]

rrna_ids <- unique(trimws(readLines(rrna_file, warn = FALSE)))
rrna_ids <- rrna_ids[nzchar(rrna_ids) & !startsWith(rrna_ids, "#")]
is_rrna <- rownames(count_matrix) %in% rrna_ids
counts_no_rrna <- count_matrix[!is_rrna, , drop = FALSE]
keep_expression <- rowSums(counts_no_rrna >= min_count) >= min_samples
analysis_counts <- counts_no_rrna[keep_expression, , drop = FALSE]

if (nrow(analysis_counts) == 0) {
  stop("No host genes passed the rRNA and minimum-expression filters.")
}

readr::write_tsv(
  as_tibble(counts_no_rrna, rownames = "gene_id"),
  file.path(output_dir, "host_counts_no_rrna.tsv")
)
readr::write_tsv(
  tibble(
    metric = c(
      "mapping_branch",
      "samples_retained",
      "input_genes",
      "rrna_ids_in_exclusion_list",
      "rrna_genes_removed_from_matrix",
      "genes_after_rrna_removal",
      "minimum_count",
      "minimum_samples",
      "genes_entering_deseq2",
      "sample_outlier_policy",
      "cooks_filtering"
    ),
    value = as.character(c(
      mapping_branch,
      ncol(count_matrix),
      nrow(count_matrix),
      length(rrna_ids),
      sum(is_rrna),
      nrow(counts_no_rrna),
      min_count,
      min_samples,
      nrow(analysis_counts),
      "All samples retained; automatic replacement disabled",
      "Disabled for planned contrasts"
    ))
  ),
  file.path(output_dir, "host_deseq2_filter_audit.tsv")
)

metadata <- metadata |>
  mutate(
    diet = factor(.data$diet, levels = c("50", "33", "83")),
    treatment = factor(.data$treatment, levels = c("Control", "Infected")),
    group = factor(
      paste0("diet", .data$diet, "_", tolower(.data$treatment)),
      levels = c(
        "diet50_control",
        "diet33_control",
        "diet83_control",
        "diet50_infected",
        "diet33_infected",
        "diet83_infected"
      )
    )
  )

if (anyNA(metadata$group)) {
  stop("Unexpected diet or treatment values in metadata.")
}

dds <- DESeqDataSetFromMatrix(
  countData = analysis_counts,
  colData = as.data.frame(metadata),
  design = ~ 0 + group
)
dds <- DESeq(
  dds,
  minReplicatesForReplace = Inf,
  quiet = TRUE
)

pair_contrast <- function(id, family, numerator, denominator, numerator_label, denominator_label) {
  list(
    id = id,
    family = family,
    numerator = numerator_label,
    denominator = denominator_label,
    weights = setNames(c(1, -1), c(numerator, denominator))
  )
}

contrast_definitions <- list(
  pair_contrast(
    "infection_diet33", "Infection within diet",
    "diet33_infected", "diet33_control", "Infected diet 33", "Control diet 33"
  ),
  pair_contrast(
    "infection_diet50", "Infection within diet",
    "diet50_infected", "diet50_control", "Infected diet 50", "Control diet 50"
  ),
  pair_contrast(
    "infection_diet83", "Infection within diet",
    "diet83_infected", "diet83_control", "Infected diet 83", "Control diet 83"
  ),
  pair_contrast(
    "control_diet33_vs50", "Diet among controls",
    "diet33_control", "diet50_control", "Control diet 33", "Control diet 50"
  ),
  pair_contrast(
    "control_diet83_vs50", "Diet among controls",
    "diet83_control", "diet50_control", "Control diet 83", "Control diet 50"
  ),
  pair_contrast(
    "control_diet83_vs33", "Diet among controls",
    "diet83_control", "diet33_control", "Control diet 83", "Control diet 33"
  ),
  pair_contrast(
    "infected_diet33_vs50", "Diet among infected",
    "diet33_infected", "diet50_infected", "Infected diet 33", "Infected diet 50"
  ),
  pair_contrast(
    "infected_diet83_vs50", "Diet among infected",
    "diet83_infected", "diet50_infected", "Infected diet 83", "Infected diet 50"
  ),
  pair_contrast(
    "infected_diet83_vs33", "Diet among infected",
    "diet83_infected", "diet33_infected", "Infected diet 83", "Infected diet 33"
  ),
  list(
    id = "infection_average_all_diets",
    family = "Average infection response",
    numerator = "Average infected",
    denominator = "Average control",
    weights = c(
      diet33_infected = 1 / 3,
      diet50_infected = 1 / 3,
      diet83_infected = 1 / 3,
      diet33_control = -1 / 3,
      diet50_control = -1 / 3,
      diet83_control = -1 / 3
    )
  ),
  list(
    id = "interaction_diet33_vs50",
    family = "Diet by infection interaction",
    numerator = "Infection response diet 33",
    denominator = "Infection response diet 50",
    weights = c(
      diet33_infected = 1, diet33_control = -1,
      diet50_infected = -1, diet50_control = 1
    )
  ),
  list(
    id = "interaction_diet83_vs50",
    family = "Diet by infection interaction",
    numerator = "Infection response diet 83",
    denominator = "Infection response diet 50",
    weights = c(
      diet83_infected = 1, diet83_control = -1,
      diet50_infected = -1, diet50_control = 1
    )
  ),
  list(
    id = "interaction_diet83_vs33",
    family = "Diet by infection interaction",
    numerator = "Infection response diet 83",
    denominator = "Infection response diet 33",
    weights = c(
      diet83_infected = 1, diet83_control = -1,
      diet33_infected = -1, diet33_control = 1
    )
  )
)

coefficient_names <- resultsNames(dds)
result_rows <- vector("list", length(contrast_definitions))
summary_rows <- vector("list", length(contrast_definitions))

for (index in seq_along(contrast_definitions)) {
  definition <- contrast_definitions[[index]]
  contrast_vector <- setNames(rep(0, length(coefficient_names)), coefficient_names)
  requested_names <- paste0("group", names(definition$weights))
  if (!all(requested_names %in% coefficient_names)) {
    stop(
      "Missing DESeq2 group coefficients for contrast ",
      definition$id,
      ": ",
      paste(setdiff(requested_names, coefficient_names), collapse = ",")
    )
  }
  contrast_vector[requested_names] <- unname(definition$weights)

  result <- results(
    dds,
    contrast = contrast_vector,
    alpha = alpha,
    cooksCutoff = FALSE
  ) |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    as_tibble() |>
    mutate(
      mapping_branch = mapping_branch,
      contrast_id = definition$id,
      contrast_family = definition$family,
      numerator = definition$numerator,
      denominator = definition$denominator,
      higher_in = case_when(
        is.na(.data$log2FoldChange) ~ NA_character_,
        .data$log2FoldChange >= 0 ~ definition$numerator,
        TRUE ~ definition$denominator
      ),
      significant = (
        !is.na(.data$padj) &
          .data$padj < alpha &
          abs(.data$log2FoldChange) >= lfc_threshold
      )
    )
  result_rows[[index]] <- result
  summary_rows[[index]] <- tibble(
    mapping_branch = mapping_branch,
    contrast_id = definition$id,
    contrast_family = definition$family,
    numerator = definition$numerator,
    denominator = definition$denominator,
    samples_retained = ncol(dds),
    upregulated = sum(
      result$significant & result$log2FoldChange > 0,
      na.rm = TRUE
    ),
    downregulated = sum(
      result$significant & result$log2FoldChange < 0,
      na.rm = TRUE
    )
  )
}

all_results <- bind_rows(result_rows)
contrast_summary <- bind_rows(summary_rows)
readr::write_csv(
  all_results,
  file.path(output_dir, "host_deseq2_all_contrasts.csv")
)
readr::write_tsv(
  contrast_summary,
  file.path(output_dir, "host_deseq2_contrast_summary.tsv")
)

vst <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_data <- plotPCA(
  vst,
  intgroup = c("diet", "treatment"),
  returnData = TRUE
)
percent_variance <- round(100 * attr(pca_data, "percentVar"), 1)
diet_colors <- c("33" = "#D9822B", "50" = "#2E8B57", "83" = "#7B4FA3")
treatment_shapes <- c("Control" = 21, "Infected" = 24)
pca_plot <- ggplot(
  pca_data,
  aes(
    x = .data$PC1,
    y = .data$PC2,
    color = .data$diet,
    shape = .data$treatment,
    fill = .data$diet
  )
) +
  geom_point(size = 4, stroke = 0.8) +
  scale_color_manual(values = diet_colors) +
  scale_fill_manual(values = diet_colors) +
  scale_shape_manual(values = treatment_shapes) +
  labs(
    title = paste("Host transcriptome PCA:", mapping_branch),
    subtitle = "All samples retained; annotated rRNA loci removed",
    x = paste0("PC1: ", percent_variance[[1]], "% variance"),
    y = paste0("PC2: ", percent_variance[[2]], "% variance"),
    color = "Diet",
    fill = "Diet",
    shape = "Treatment"
  ) +
  theme_bw(base_size = 13)

bar_data <- contrast_summary |>
  select(
    "contrast_id",
    "contrast_family",
    "upregulated",
    "downregulated"
  ) |>
  pivot_longer(
    cols = c("upregulated", "downregulated"),
    names_to = "direction",
    values_to = "n_degs"
  )
bar_plot <- ggplot(
  bar_data,
  aes(x = .data$contrast_id, y = .data$n_degs, fill = .data$direction)
) +
  geom_col(position = "dodge") +
  facet_wrap(vars(.data$contrast_family), scales = "free_x") +
  scale_fill_manual(
    values = c(upregulated = "#D73027", downregulated = "#2878B5"),
    labels = c(upregulated = "Higher in numerator", downregulated = "Lower in numerator")
  ) +
  labs(
    title = paste("Host DEG counts:", mapping_branch),
    subtitle = paste0("padj < ", alpha, "; |log2FC| >= ", lfc_threshold),
    x = NULL,
    y = "Significant DEGs",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "host_transcriptome_pca.png"),
  pca_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(output_dir, "host_deg_counts.png"),
  bar_plot,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)
