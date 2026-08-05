#!/usr/bin/env Rscript

# Build colleague-style family composition plots from Kraken2 or Bracken files.

suppressPackageStartupMessages({
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

manifest_file <- arg_value("--manifest")
output_dir <- arg_value("--output-dir")
top_n <- as.integer(arg_value("--top-families"))
focus_families <- strsplit(
  arg_value("--focus-families"),
  ",",
  fixed = TRUE
)[[1]]
focus_families <- trimws(focus_families[nzchar(trimws(focus_families))])

manifest <- readr::read_tsv(manifest_file, show_col_types = FALSE)

read_family_file <- function(path, format) {
  if (format == "bracken") {
    table <- readr::read_tsv(path, show_col_types = FALSE)
    count_columns <- intersect(
      c("new_est_reads", "estimated_reads", "kraken_assigned_reads"),
      names(table)
    )
    if (length(count_columns) == 0) {
      stop("No supported read-count column found in Bracken file: ", path)
    }
    count_column <- count_columns[[1]]
    return(
      table |>
        transmute(
          Family = trimws(.data$name),
          Estimated_reads = as.numeric(.data[[count_column]])
        )
    )
  }

  table <- readr::read_tsv(
    path,
    col_names = c(
      "percent", "clade_reads", "taxon_reads", "rank", "taxid", "name"
    ),
    show_col_types = FALSE,
    trim_ws = FALSE
  )
  table |>
    filter(trimws(.data$rank) == "F") |>
    transmute(
      Family = trimws(.data$name),
      Estimated_reads = as.numeric(.data$clade_reads)
    )
}

all_rows <- lapply(seq_len(nrow(manifest)), function(index) {
  metadata <- manifest[index, ]
  family_rows <- read_family_file(
    metadata$abundance_file,
    metadata$abundance_format
  )
  if (nrow(family_rows) == 0) {
    family_rows <- tibble(
      Family = "No family-level assignment",
      Estimated_reads = 0
    )
  }
  family_rows |>
    mutate(
      sample_id = as.character(metadata$sample_id),
      source = as.character(metadata$source),
      Treatment = as.character(metadata$treatment),
      Diet = as.character(metadata$diet)
    )
}) |>
  bind_rows()

source_labels <- c(
  all_trimmed = "All trimmed reads",
  host_unmapped = "Reads unmapped to S. gregaria",
  competitive_unmapped = "Reads unmapped to locust and Metarhizium"
)

for (current_source in unique(as.character(manifest$source))) {
  source_rows <- all_rows |>
    filter(.data$source == current_source)

  leading_families <- source_rows |>
    group_by(.data$Family) |>
    summarise(total = sum(.data$Estimated_reads, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(.data$total)) |>
    slice_head(n = top_n) |>
    pull(.data$Family)
  leading_families <- union(
    leading_families,
    intersect(
      c(focus_families, "No family-level assignment"),
      unique(source_rows$Family)
    )
  )

  plot_rows <- source_rows |>
    mutate(
      Family_plot = if_else(
        .data$Family %in% leading_families,
        .data$Family,
        "Other"
      )
    ) |>
    group_by(
      .data$sample_id, .data$source, .data$Treatment, .data$Diet,
      .data$Family_plot
    ) |>
    summarise(
      Estimated_reads = sum(.data$Estimated_reads, na.rm = TRUE),
      .groups = "drop"
    ) |>
    group_by(.data$sample_id) |>
    mutate(
      Relative_abundance = (
        100 * .data$Estimated_reads /
          max(sum(.data$Estimated_reads), 1)
      )
    ) |>
    ungroup()

  ordered_samples <- plot_rows |>
    distinct(.data$sample_id) |>
    mutate(sample_number = suppressWarnings(as.numeric(.data$sample_id))) |>
    arrange(.data$sample_number, .data$sample_id) |>
    pull(.data$sample_id)
  plot_rows$sample_label <- factor(
    paste0("mehreen_", plot_rows$sample_id, "_MERGE"),
    levels = paste0("mehreen_", ordered_samples, "_MERGE")
  )

  observed_families <- unique(plot_rows$Family_plot)
  family_levels <- c(
    sort(setdiff(observed_families, "Other")),
    intersect("Other", observed_families)
  )
  plot_rows$Family_plot <- factor(plot_rows$Family_plot, levels = family_levels)
  family_colors <- setNames(
    grDevices::hcl.colors(length(family_levels), palette = "Dynamic"),
    family_levels
  )
  if ("Other" %in% names(family_colors)) {
    family_colors[["Other"]] <- "#B8B8B8"
  }

  common_theme <- theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1, size = 9),
      legend.position = "right",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )

  source_title <- unname(source_labels[current_source])
  if (is.na(source_title)) {
    source_title <- current_source
  }

  reads_plot <- ggplot(
    plot_rows,
    aes(x = .data$sample_label, y = .data$Estimated_reads, fill = .data$Family_plot)
  ) +
    geom_col(width = 0.9) +
    facet_grid(~ Treatment + Diet, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = family_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    labs(
      title = paste("Family composition:", source_title),
      subtitle = paste(
        "Top", top_n,
        "families plus focal Clavicipitaceae; remaining families pooled as Other"
      ),
      x = NULL,
      y = "Estimated family-level reads"
    ) +
    common_theme

  percent_plot <- ggplot(
    plot_rows,
    aes(
      x = .data$sample_label,
      y = .data$Relative_abundance,
      fill = .data$Family_plot
    )
  ) +
    geom_col(width = 0.9) +
    facet_grid(~ Treatment + Diet, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = family_colors, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 25),
      labels = function(x) paste0(x, "%")
    ) +
    labs(
      title = paste("Relative family composition:", source_title),
      subtitle = "Percent of reads assigned at family level",
      x = NULL,
      y = "Relative abundance"
    ) +
    common_theme

  source_dir <- file.path(output_dir, current_source)
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(source_dir, "Family_Composition_Reads.png"),
    reads_plot,
    width = 18,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(source_dir, "Family_Composition_Reads.svg"),
    reads_plot,
    width = 18,
    height = 9,
    bg = "white"
  )
  ggsave(
    file.path(source_dir, "Family_Composition_Percent.png"),
    percent_plot,
    width = 18,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(source_dir, "Family_Composition_Percent.svg"),
    percent_plot,
    width = 18,
    height = 9,
    bg = "white"
  )
  readr::write_tsv(
    plot_rows,
    file.path(source_dir, "family_composition_plot_data.tsv")
  )
}
