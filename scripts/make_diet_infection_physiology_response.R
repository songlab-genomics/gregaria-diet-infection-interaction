#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
repo_root <- if (!is.na(script_file)) normalizePath(file.path(dirname(script_file), "..")) else getwd()
if (!dir.exists(file.path(repo_root, "analysis"))) {
  repo_root <- normalizePath(getwd())
}

results_root <- file.path(repo_root, "output", "rmd_runs")

latest_dir <- function(pattern) {
  dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(pattern, basename(dirs))]
  if (length(dirs) == 0) {
    stop("No output folder matching ", pattern, " under ", results_root, call. = FALSE)
  }
  dirs[order(file.info(dirs)$mtime, decreasing = TRUE)][[1]]
}

first_nonempty <- function(x, fallback = NA_character_) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) fallback else x[[1]]
}

classify_physiology_theme <- function(annotation_text) {
  text <- stringr::str_to_lower(dplyr::coalesce(annotation_text, ""))
  dplyr::case_when(
    stringr::str_detect(
      text,
      "antimicrobial|defensin|attacin|cecropin|diptericin|gloverin|metchnikowin|thaumatin|termicin|lysozyme"
    ) ~ "Antimicrobial peptides",
    stringr::str_detect(
      text,
      "toll|imd|relish|myd88|cactus|dorsal|nf-kappa|nf kappa|jak|stat"
    ) ~ "Immune-related",
    stringr::str_detect(
      text,
      "immune|peptidoglycan|pgrp|gram-negative|gram negative|beta-1,3-glucan|glucan binding|phenoloxidase|prophenoloxidase|melanization|lectin|hemocyanin|complement|serpin"
    ) ~ "Immune-related",
    stringr::str_detect(
      text,
      "lipid|fatty acid|fatty-acid|acyl-coa|acyl coa|acetyl-coa carboxylase|elongase|desaturase|glycerolipid|phospholipid|triacylglycerol|diacylglycerol|sterol|lipase"
    ) ~ "Lipid synthesis / metabolism",
    stringr::str_detect(
      text,
      "protein synthesis|translation|ribosom|elongation factor|initiation factor|aminoacyl|trna synthetase|vitellogenin|hexamerin|arylphorin|storage protein"
    ) ~ "Protein synthesis / storage",
    stringr::str_detect(
      text,
      "detox|cytochrome p450|\\bcyp[0-9a-z]|glutathione|gst|carboxylesterase|esterase|udp-glucuronosyltransferase|\\bugt\\b|abc transporter|oxidoreductase|peroxidase|catalase|superoxide|heat shock|stress"
    ) ~ "Detoxification",
    TRUE ~ "Other"
  )
}

wrap_gene_label <- function(description, gene_id, width = 58) {
  label <- ifelse(is.na(description) | description == "",
                  gene_id,
                  paste0(description, " [", gene_id, "]"))
  stringr::str_wrap(label, width = width)
}

theme_meeting <- function(base_size = 15) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(color = "#374151"),
      strip.background = element_rect(fill = "#F8F6F1", color = "#111827", linewidth = 0.8),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      axis.text = element_text(color = "#111827"),
      axis.title = element_text(face = "bold")
    )
}

infection_dir <- latest_dir("^all-samples-infection-within-diet-deg_")
catalogue_dir <- latest_dir("^all-deg-catalogue_")

infection_file <- file.path(infection_dir, "all_infection_within_diet_results.csv")
catalogue_file <- file.path(catalogue_dir, "all_samples_deg_catalogue.csv")

infection_results <- read.csv(infection_file, check.names = FALSE)
catalogue <- read.csv(catalogue_file, check.names = FALSE, na.strings = c("", "NA"))

gene_annotation <- catalogue |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(
    Description = first_nonempty(Description, "no GFF description"),
    gene_biotype = first_nonempty(gene_biotype, "unknown"),
    GOs = first_nonempty(GOs, ""),
    KEGG_ko = first_nonempty(KEGG_ko, ""),
    KEGG_Pathway = first_nonempty(KEGG_Pathway, ""),
    PFAMs = first_nonempty(PFAMs, ""),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    annotation_text = stringr::str_squish(paste(
      Description, GOs, KEGG_ko, KEGG_Pathway, PFAMs, sep = " | "
    )),
    detailed_theme = classify_physiology_theme(annotation_text),
    meeting_theme = dplyr::case_when(
      detailed_theme %in% c("Antimicrobial peptides", "Immune-related") ~
        "Immune pathway",
      detailed_theme == "Lipid synthesis / metabolism" ~
        "Lipid metabolism",
      detailed_theme == "Protein synthesis / storage" ~
        "Protein storage and synthesis",
      detailed_theme == "Detoxification" ~
        "Detoxification",
      TRUE ~ NA_character_
    )
  )

theme_levels <- c(
  "Immune pathway",
  "Protein storage and synthesis",
  "Lipid metabolism",
  "Detoxification"
)
theme_labels <- c(
  "Immune pathway" = "Immune\npathway",
  "Protein storage and synthesis" = "Protein\nstorage\nand synthesis",
  "Lipid metabolism" = "Lipid\nmetabolism",
  "Detoxification" = "Detoxification"
)
direction_levels <- c("Higher in infected", "Higher in control")
direction_colors <- c(
  "Higher in infected" = "#B91C1C",
  "Higher in control" = "#1D4ED8"
)
diet_colors <- c("33" = "#E98A2C", "50" = "#2C9C5B", "83" = "#7A4BB3")

infection_annotated <- infection_results |>
  dplyr::left_join(gene_annotation, by = "gene_id") |>
  dplyr::mutate(
    diet = factor(as.character(diet), levels = c("33", "50", "83")),
    significant = as.logical(significant),
    direction = dplyr::if_else(log2FoldChange >= 0,
                               "Higher in infected",
                               "Higher in control"),
    direction = factor(direction, levels = direction_levels),
    meeting_theme = factor(meeting_theme, levels = theme_levels)
  )

candidate_results <- infection_annotated |>
  dplyr::filter(significant, !is.na(meeting_theme))

if (nrow(candidate_results) == 0) {
  stop("No significant infection-within-diet DEGs matched the physiology themes.", call. = FALSE)
}

run_id <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_dir <- file.path(results_root, paste0("diet-infection-physiology-response_", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate_counts <- candidate_results |>
  dplyr::count(meeting_theme, diet, direction, name = "n") |>
  tidyr::complete(
    meeting_theme = factor(theme_levels, levels = theme_levels),
    diet = factor(c("33", "50", "83"), levels = c("33", "50", "83")),
    direction = factor(direction_levels, levels = direction_levels),
    fill = list(n = 0)
  ) |>
  dplyr::arrange(meeting_theme, diet, direction) |>
  dplyr::mutate(
    theme_label = factor(theme_labels[as.character(meeting_theme)],
                         levels = unname(theme_labels[theme_levels]))
  )

candidate_counts |>
  readr::write_csv(file.path(out_dir, "diet_infection_candidate_response_counts.csv"))

diet_50_83_delta <- candidate_counts |>
  dplyr::filter(diet %in% c("50", "83")) |>
  dplyr::mutate(diet = as.character(diet)) |>
  tidyr::pivot_wider(names_from = diet, values_from = n, values_fill = 0) |>
  dplyr::mutate(
    diet50_minus_diet83 = `50` - `83`,
    theme_label = factor(theme_labels[as.character(meeting_theme)],
                         levels = unname(theme_labels[theme_levels])),
    more_in = dplyr::case_when(
      diet50_minus_diet83 > 0 ~ "More in diet 50",
      diet50_minus_diet83 < 0 ~ "More in diet 83",
      TRUE ~ "Same count"
    )
  )

diet_50_83_delta |>
  readr::write_csv(file.path(out_dir, "diet50_vs_diet83_candidate_response_delta.csv"))

gene_level_table <- infection_annotated |>
  dplyr::filter(!is.na(meeting_theme)) |>
  dplyr::select(
    gene_id, Description, gene_biotype, meeting_theme, diet, direction,
    significant, log2FoldChange, padj, baseMean
  ) |>
  dplyr::arrange(meeting_theme, gene_id, diet)

gene_level_table |>
  readr::write_csv(file.path(out_dir, "diet_infection_candidate_gene_table_all_tests.csv"))

candidate_gene_summary <- candidate_results |>
  dplyr::group_by(gene_id, Description, gene_biotype, meeting_theme) |>
  dplyr::summarise(
    called_in_diets = paste(sort(unique(as.character(diet))), collapse = ";"),
    n_called_diets = dplyr::n_distinct(diet),
    max_abs_log2FC = max(abs(log2FoldChange), na.rm = TRUE),
    min_padj = min(padj, na.rm = TRUE),
    response_directions = paste(sort(unique(as.character(direction))), collapse = "; "),
    .groups = "drop"
  ) |>
  dplyr::arrange(meeting_theme, min_padj, dplyr::desc(max_abs_log2FC))

candidate_gene_summary |>
  readr::write_csv(file.path(out_dir, "diet_infection_candidate_gene_summary_significant.csv"))

p_counts <- candidate_counts |>
  dplyr::mutate(label_y = n + max(candidate_counts$n) * 0.025) |>
  ggplot(aes(x = diet, y = n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.62, color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(n > 0, n, "")),
    position = position_dodge(width = 0.74),
    vjust = -0.35,
    size = 4.2
  ) +
  facet_wrap(~ theme_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = direction_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(
    title = "Diet-specific infection-response DEGs",
    subtitle = "Counts are significant infected-vs-control DEGs within each diet",
    x = "Diet",
    y = "Number of DEGs",
    fill = "Direction"
  ) +
  theme_meeting(15)

p_delta <- diet_50_83_delta |>
  ggplot(aes(x = diet50_minus_diet83, y = theme_label, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280", linewidth = 0.8) +
  geom_segment(aes(x = 0, xend = diet50_minus_diet83, yend = theme_label),
               linewidth = 1.3, alpha = 0.72) +
  geom_point(size = 5.4) +
  geom_text(
    aes(label = sprintf("%+d", diet50_minus_diet83)),
    nudge_x = ifelse(diet_50_83_delta$diet50_minus_diet83 >= 0, 4, -4),
    size = 4.8,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~ direction, nrow = 1) +
  scale_color_manual(values = direction_colors) +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  labs(
    title = "Diet 50 versus diet 83",
    subtitle = "Positive values mean more infection-response DEGs in diet 50; negative values mean more in diet 83",
    x = "DEG count difference: diet 50 - diet 83",
    y = NULL,
    color = "Direction"
  ) +
  theme_meeting(15) +
  theme(legend.position = "none")

top_gene_ids <- candidate_gene_summary |>
  dplyr::group_by(meeting_theme) |>
  dplyr::slice_min(order_by = min_padj, n = 8, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::pull(gene_id)

heatmap_data <- infection_annotated |>
  dplyr::filter(gene_id %in% top_gene_ids) |>
  dplyr::mutate(
    gene_label = wrap_gene_label(Description, gene_id),
    meeting_theme = factor(meeting_theme, levels = theme_levels),
    theme_label = factor(theme_labels[as.character(meeting_theme)],
                         levels = unname(theme_labels[theme_levels]))
  ) |>
  dplyr::left_join(
    candidate_gene_summary |>
      dplyr::select(gene_id, min_padj, max_abs_log2FC),
    by = "gene_id"
  ) |>
  dplyr::group_by(meeting_theme) |>
  dplyr::mutate(
    gene_label = factor(
      gene_label,
      levels = unique(gene_label[order(min_padj, -max_abs_log2FC)])
    )
  ) |>
  dplyr::ungroup()

p_heat <- heatmap_data |>
  ggplot(aes(x = diet, y = gene_label, fill = log2FoldChange)) +
  geom_tile(aes(alpha = significant),
            width = 0.55,
            height = 0.9,
            color = "white",
            linewidth = 0.25) +
  geom_point(
    data = dplyr::filter(heatmap_data, significant),
    aes(x = diet, y = gene_label),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.9,
    stroke = 0.65,
    fill = "black",
    color = "black"
  ) +
  facet_grid(theme_label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#1D4ED8",
    mid = "#F8F6F1",
    high = "#B91C1C",
    midpoint = 0,
    oob = scales::squish,
    name = "Infected vs control\nlog2FC"
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.28), guide = "none") +
  labs(
    title = "Candidate genes driving diet-specific infection responses",
    subtitle = "Black dots mark diet-specific significant DEGs; tile color shows infection log2FC in every diet",
    x = "Diet",
    y = NULL
  ) +
  theme_meeting(13) +
  theme(
    axis.text.y = element_text(size = 10.3, lineheight = 0.92,
                               margin = margin(r = 8)),
    legend.position = "right",
    strip.text.y = element_text(angle = 0, size = 12),
    panel.spacing.y = unit(0.55, "lines"),
    plot.margin = margin(5.5, 12, 5.5, 24)
  )

composite <- (p_counts | p_delta) / p_heat +
  plot_layout(heights = c(0.82, 1.9), widths = c(1.08, 1), guides = "collect") +
  plot_annotation(
    title = "Diet modifies the fat body transcriptional response to infection",
    subtitle = "Focused on immune pathway, protein storage and synthesis, lipid metabolism, and detoxification candidate DEGs",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 26),
      plot.subtitle = element_text(size = 16, color = "#374151"),
      plot.tag = element_text(face = "bold", size = 22)
    )
  )

save_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")),
                  plot, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".svg")),
                  plot, width = width, height = height, bg = "white")
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")),
                  plot, width = width, height = height, bg = "white")
}

save_plot(composite, "diet_infection_physiology_response_composite", 24, 18)
save_plot(p_counts, "diet_infection_physiology_response_counts", 17, 5.2)
save_plot(p_delta, "diet50_vs_diet83_physiology_response_delta", 13, 5.2)
save_plot(p_heat, "diet_infection_physiology_candidate_gene_heatmap", 12, 13)

message("Wrote diet-infection physiology response outputs to: ", out_dir)
message("Input infection DEG table: ", infection_file)
message("Input annotation/catalogue table: ", catalogue_file)
message("")
message("Diet 50 minus diet 83 candidate DEG count differences:")
print(diet_50_83_delta |> dplyr::select(meeting_theme, direction, `50`, `83`, diet50_minus_diet83))
