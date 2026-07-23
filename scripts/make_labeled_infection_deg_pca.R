#!/usr/bin/env Rscript

# Presentation export: labeled PCA panels for infected-vs-control DEGs within each diet.
# This script reuses the saved all-samples infection-within-diet PCA score tables
# so it does not rerun DESeq2; it only replots the existing PCA coordinates.

required <- c("DESeq2", "SummarizedExperiment", "ggplot2", "ggrepel", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install required packages before running: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_dir) == "scripts") {
  project_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
}

results_root <- file.path(project_dir, "output", "rmd_runs")
run_dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
source_dirs <- run_dirs[grepl("^all-samples-infection-within-diet-deg_", basename(run_dirs))]
if (length(source_dirs) == 0) {
  stop("No all-samples infection-within-diet run folder found under ", results_root,
       call. = FALSE)
}
source_dir <- source_dirs[which.max(file.info(source_dirs)$mtime)]

run_stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_dir <- file.path(results_root, paste0("presentation-labeled-infection-deg-pca_", run_stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

treatment_colors <- c("Control" = "#7B8794", "Infected" = "#C85A54")
treatment_fills <- c("Control" = "#D8DEE2", "Infected" = "#F0C7BF")
treatment_shapes <- c("Control" = 16, "Infected" = 17)

summary_file <- file.path(source_dir, "infection_within_diet_summary.csv")
deg_summary <- read.csv(summary_file, check.names = FALSE)

metadata <- read.delim(file.path(project_dir, "data", "metadata",
                                 "mehreen_metadata_fixed.txt"),
                       check.names = FALSE)
metadata$label <- as.character(metadata$label)
metadata$Diet <- factor(metadata$Diet, levels = c("33", "50", "83"))
metadata$Treatment <- relevel(factor(metadata$Treatment), ref = "Control")

counts_raw <- read.csv(file.path(project_dir, "data", "raw_read_counts",
                                 "master_counts_no_rRNA.csv"),
                       check.names = FALSE)
counts <- as.matrix(counts_raw[, -1, drop = FALSE])
rownames(counts) <- counts_raw[[1]]
storage.mode(counts) <- "integer"

count_labels_all <- sub("_MERGE$", "", sub("^mehreen_", "", colnames(counts)))
keep_samples <- count_labels_all %in% metadata$label
counts <- counts[, keep_samples, drop = FALSE]
count_labels <- count_labels_all[keep_samples]
metadata <- metadata[match(count_labels, metadata$label), , drop = FALSE]
stopifnot(identical(metadata$label, count_labels))
rownames(metadata) <- metadata$label
colnames(counts) <- count_labels

keep_vst <- rowSums(counts) >= 10
dds_vst <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts[keep_vst, , drop = FALSE],
  colData = metadata,
  design = ~ Diet + Treatment
)
vsd <- DESeq2::vst(dds_vst, blind = TRUE)
vst_mat <- SummarizedExperiment::assay(vsd)

make_hull <- function(df) {
  if (nrow(df) < 3) return(df)
  df[grDevices::chull(df$PC1, df$PC2), , drop = FALSE]
}

save_png_svg <- function(plot, file_base, width, height, dpi = 350) {
  png_path <- file.path(out_dir, paste0(file_base, ".png"))
  svg_path <- file.path(out_dir, paste0(file_base, ".svg"))
  ggplot2::ggsave(png_path, plot, width = width, height = height,
                  dpi = dpi, bg = "white", limitsize = FALSE)
  grDevices::svg(svg_path, width = width, height = height, onefile = FALSE)
  print(plot)
  grDevices::dev.off()
  invisible(c(png = png_path, svg = svg_path))
}

read_diet_scores <- function(diet_level) {
  deg_genes <- readLines(
    file.path(source_dir, paste0("significant_genes_infection_diet_",
                                 diet_level, ".txt")),
    warn = FALSE
  )
  deg_genes <- intersect(deg_genes, rownames(vst_mat))
  metadata_sub <- droplevels(metadata[metadata$Diet == diet_level, , drop = FALSE])
  mat <- vst_mat[deg_genes, metadata_sub$label, drop = FALSE]
  mat <- mat[apply(mat, 1, var) > 0, , drop = FALSE]
  if (nrow(mat) < 3 || ncol(mat) < 4) {
    stop("Diet ", diet_level, " does not have enough DEG expression data for PCA.",
         call. = FALSE)
  }

  pca_local <- stats::prcomp(t(mat), scale. = FALSE)
  percent <- round(100 * (pca_local$sdev^2 / sum(pca_local$sdev^2)), 1)
  scores <- data.frame(
    label = rownames(pca_local$x),
    PC1 = pca_local$x[, "PC1"],
    PC2 = pca_local$x[, "PC2"],
    Diet = metadata_sub$Diet,
    Treatment = metadata_sub$Treatment,
    check.names = FALSE
  )
  centroids <- aggregate(cbind(PC1, PC2) ~ Treatment, scores, mean)
  names(centroids)[names(centroids) == "PC1"] <- "centroid_PC1"
  names(centroids)[names(centroids) == "PC2"] <- "centroid_PC2"

  # PCA axis signs are arbitrary. Flip PC1 when needed so control appears left.
  control_x <- centroids$centroid_PC1[centroids$Treatment == "Control"][1]
  infected_x <- centroids$centroid_PC1[centroids$Treatment == "Infected"][1]
  if (is.finite(control_x) && is.finite(infected_x) && control_x > infected_x) {
    scores$PC1 <- -scores$PC1
    centroids$centroid_PC1 <- -centroids$centroid_PC1
  }

  scores$Treatment <- factor(scores$Treatment, levels = c("Control", "Infected"))
  centroids$Treatment <- factor(centroids$Treatment, levels = c("Control", "Infected"))
  list(scores = scores, centroids = centroids, percent = percent,
       n_genes = nrow(mat))
}

make_labeled_pca <- function(diet_level, show_legend = TRUE) {
  pca_parts <- read_diet_scores(diet_level)
  scores <- pca_parts$scores
  centroids <- pca_parts$centroids
  hulls <- do.call(rbind, lapply(split(scores, scores$Treatment), make_hull))
  if (!("PC1" %in% names(scores)) || !("PC2" %in% names(scores))) {
    stop("PCA score table must contain PC1 and PC2 columns.", call. = FALSE)
  }

  n_sig <- deg_summary$n_significant_genes[as.character(deg_summary$diet) == diet_level][1]

  write.csv(scores,
            file.path(out_dir, paste0("diet_", diet_level,
                                      "_infection_deg_pca_labeled_scores.csv")),
            row.names = FALSE)
  write.csv(centroids,
            file.path(out_dir, paste0("diet_", diet_level,
                                      "_infection_deg_pca_labeled_centroids.csv")),
            row.names = FALSE)

  ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_polygon(
      data = hulls,
      ggplot2::aes(fill = Treatment, group = Treatment),
      alpha = 0.20,
      color = NA
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = Treatment, shape = Treatment),
      size = 4.3,
      alpha = 0.96
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 4.6,
      color = "#111827",
      min.segment.length = 0,
      segment.color = "#B8B8B8",
      segment.size = 0.25,
      box.padding = 0.38,
      point.padding = 0.28,
      max.overlaps = Inf,
      seed = 24
    ) +
    ggplot2::geom_point(
      data = centroids,
      ggplot2::aes(x = centroid_PC1, y = centroid_PC2, color = Treatment),
      shape = 21,
      fill = "white",
      stroke = 1.2,
      size = 6.2,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_label_repel(
      data = centroids,
      ggplot2::aes(x = centroid_PC1, y = centroid_PC2, label = Treatment),
      color = "#1F2937",
      fill = "white",
      label.size = 0.25,
      size = 5.0,
      min.segment.length = 0,
      box.padding = 0.8,
      point.padding = 0.45,
      seed = 25,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(values = treatment_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = treatment_fills, drop = FALSE) +
    ggplot2::scale_shape_manual(values = treatment_shapes, drop = FALSE) +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = ggplot2::element_text(size = 20),
      legend.text = ggplot2::element_text(size = 18),
      plot.title = ggplot2::element_text(size = 30, face = "plain"),
      plot.subtitle = ggplot2::element_text(size = 23, color = "#111827"),
      axis.title = ggplot2::element_text(size = 24),
      axis.text = ggplot2::element_text(size = 18, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      plot.margin = ggplot2::margin(12, 18, 12, 12)
    ) +
    ggplot2::labs(
      title = paste("Diet", diet_level, "PCA using infection DEGs only"),
      subtitle = paste(format(n_sig, big.mark = ","), "significant genes"),
      x = paste0("PC1: ", pca_parts$percent[1], "% variance"),
      y = paste0("PC2: ", pca_parts$percent[2], "% variance"),
      color = "Treatment",
      fill = "Treatment",
      shape = "Treatment"
    )
}

plots <- lapply(c("33", "50", "83"), make_labeled_pca)
names(plots) <- c("33", "50", "83")

for (diet_level in names(plots)) {
  save_png_svg(
    plots[[diet_level]],
    paste0("diet_", diet_level, "_infection_deg_pca_labeled"),
    width = 13,
    height = 9
  )
}

combined_portrait <- plots[["33"]] / plots[["50"]] / plots[["83"]]
combined_portrait <- combined_portrait +
  patchwork::plot_layout(heights = c(1, 1, 1), guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

combined_landscape <- plots[["33"]] | plots[["50"]] | plots[["83"]]
combined_landscape <- combined_landscape +
  patchwork::plot_layout(widths = c(1, 1, 1), guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

save_png_svg(combined_portrait, "all_diets_infection_deg_pca_labeled",
             width = 13, height = 24)
save_png_svg(combined_landscape, "all_diets_infection_deg_pca_labeled_landscape",
             width = 30, height = 9.5)

manifest <- data.frame(
  item = c("source_run", "output_dir",
           paste0("diet_", names(plots), "_png"),
           paste0("diet_", names(plots), "_svg"),
           "combined_portrait_png", "combined_portrait_svg",
           "combined_landscape_png", "combined_landscape_svg"),
  path = c(source_dir, out_dir,
           file.path(out_dir, paste0("diet_", names(plots), "_infection_deg_pca_labeled.png")),
           file.path(out_dir, paste0("diet_", names(plots), "_infection_deg_pca_labeled.svg")),
           file.path(out_dir, "all_diets_infection_deg_pca_labeled.png"),
           file.path(out_dir, "all_diets_infection_deg_pca_labeled.svg"),
           file.path(out_dir, "all_diets_infection_deg_pca_labeled_landscape.png"),
           file.path(out_dir, "all_diets_infection_deg_pca_labeled_landscape.svg")),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_dir, "presentation_labeled_infection_deg_pca_manifest.csv"),
          row.names = FALSE)

message("Wrote labeled PCA exports to: ", out_dir)
