#!/usr/bin/env Rscript

# Presentation workflow: diagnostic steps to test whether unplaced-scaffold DEGs
# are host biology, Metarhizium-derived signal, or ambiguous host-pathogen mapping.

required <- c("ggplot2", "grid", "svglite")
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
publication_dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
publication_dirs <- publication_dirs[
  grepl("^publication-figures_", basename(publication_dirs)) &
    file.exists(file.path(publication_dirs, "figureS2_deg_coordinates_by_context.csv"))
]
if (length(publication_dirs) == 0) {
  stop("No publication-figures folder with Figure S2 coordinates was found.",
       call. = FALSE)
}
source_dir <- publication_dirs[which.max(file.info(publication_dirs)$mtime)]

coords <- read.csv(file.path(source_dir, "figureS2_deg_coordinates_by_context.csv"),
                   check.names = FALSE)
unplaced <- coords[
  coords$chromosome_class == "Unplaced scaffold" | grepl("^NW_", coords$seqid),
  ,
  drop = FALSE
]
counts <- aggregate(
  list(n_degs = unplaced$gene_id, n_scaffolds = unplaced$seqid),
  list(DEG_context = unplaced$DEG_context),
  function(x) length(unique(x))
)

focus_context <- "Diet DEGs in infected"
focus <- counts[counts$DEG_context == focus_context, , drop = FALSE]
if (nrow(focus) == 0) {
  focus_n_degs <- 0
  focus_n_scaffolds <- 0
} else {
  focus_n_degs <- focus$n_degs[[1]]
  focus_n_scaffolds <- focus$n_scaffolds[[1]]
}

run_stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_dir <- file.path(results_root, paste0("unplaced-scaffold-origin-workflow_", run_stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wrap_text <- function(x, width = 24) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1))
}

nodes <- data.frame(
  id = c(
    "observation", "extract", "host", "pathogen", "mapping",
    "sensitivity", "host_out", "fungal_out", "ambiguous_out", "report",
    "qpcr"
  ),
  x = c(1.55, 4.25, 7.15, 7.15, 7.15, 10.25, 15.75, 15.75, 15.75, 12.9,
        4.25),
  y = c(6.35, 6.35, 7.05, 5.55, 4.05, 5.55, 6.95, 5.55, 4.15, 5.55,
        2.45),
  w = c(2.45, 2.55, 2.75, 2.75, 2.75, 2.55, 2.75, 2.75, 2.75, 2.65,
        3.15),
  h = c(1.15, 1.15, 1.05, 1.05, 1.05, 1.25, 0.90, 0.90, 0.90, 1.05,
        1.05),
  fill = c(
    "#FFF3C4", "#F8FAFC", "#E7F5EA", "#FDE8E8", "#EAF2FF",
    "#F3F4F6", "#E7F5EA", "#FDE8E8", "#EAF2FF", "#FFF7ED",
    "#F9FAFB"
  ),
  border = c(
    "#B7791F", "#374151", "#2F8F5B", "#B91C1C", "#2563EB",
    "#111827", "#2F8F5B", "#B91C1C", "#2563EB", "#C2410C",
    "#6B7280"
  ),
  label = c(
    paste0(
      "Unplaced signal\nDiet DEGs in infected\n",
      format(focus_n_degs, big.mark = ","), " DEGs on ",
      focus_n_scaffolds, " scaffolds"
    ),
    "Candidate set\nDEG genes + NW scaffolds\n+ raw-read alignments",
    "1. Host assembly audit\nGFF/NCBI taxonomy\nbiotype, GC, repeats, coverage",
    "2. Fungal similarity test\nBLAST/DIAMOND vs\nM. robertsii + fungi",
    "3. Competitive remapping\nS. gregaria + M. robertsii\nunique vs ambiguous reads",
    "Classify each gene/scaffold\nhost-like\nfungal-like\nambiguous",
    "Host-like\nretain in host DEG biology",
    "Fungal-like\nanalyze as pathogen load",
    "Ambiguous\nexclude or sensitivity-only",
    "Rerun / compare DESeq2\nDoes the diet-infected\nsignal remain?",
    "Ideal parallel validation\nqPCR / targeted fungal-load assay\nnot possible: no tissue remains"
  ),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c("observation", "extract", "extract", "extract",
           "host", "pathogen", "mapping", "sensitivity",
           "report", "report", "report"),
  to = c("extract", "host", "pathogen", "mapping",
         "sensitivity", "sensitivity", "sensitivity", "report",
         "host_out", "fungal_out", "ambiguous_out"),
  stringsAsFactors = FALSE
)

node_lookup <- nodes
rownames(node_lookup) <- node_lookup$id
edge_xy <- transform(
  edges,
  x = node_lookup[from, "x"] + node_lookup[from, "w"] / 2,
  y = node_lookup[from, "y"],
  xend = node_lookup[to, "x"] - node_lookup[to, "w"] / 2,
  yend = node_lookup[to, "y"]
)

qpcr_edge_xy <- data.frame(
  x = node_lookup["extract", "x"],
  y = node_lookup["extract", "y"] - node_lookup["extract", "h"] / 2,
  xend = node_lookup["qpcr", "x"],
  yend = node_lookup["qpcr", "y"] + node_lookup["qpcr", "h"] / 2
)

annotation_text <- paste(
  "Goal: decide whether the unplaced-scaffold signal is locust biology,",
  "Metarhizium signal, or ambiguous mapping before interpreting diet effects."
)

workflow_plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = edge_xy,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.75,
    color = "#4B5563",
    arrow = grid::arrow(length = grid::unit(0.18, "inches"), type = "closed")
  ) +
  ggplot2::geom_segment(
    data = qpcr_edge_xy,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.65,
    color = "#6B7280",
    linetype = "dashed",
    arrow = grid::arrow(length = grid::unit(0.14, "inches"), type = "closed")
  ) +
  ggplot2::geom_rect(
    data = nodes,
    ggplot2::aes(xmin = x - w / 2, xmax = x + w / 2,
                 ymin = y - h / 2, ymax = y + h / 2,
                 fill = fill, color = border),
    linewidth = 1.0
  ) +
  ggplot2::geom_text(
    data = nodes,
    ggplot2::aes(x = x, y = y, label = label),
    size = 3.45,
    lineheight = 0.90,
    color = "#111827",
    fontface = "bold"
  ) +
  ggplot2::annotate(
    "text",
    x = 8.7,
    y = 8.55,
    label = "Resolving the origin of unplaced-scaffold DEGs",
    size = 7.6,
    fontface = "bold",
    color = "#111827"
  ) +
  ggplot2::annotate(
    "text",
    x = 8.7,
    y = 8.12,
    label = wrap_text(annotation_text, width = 100),
    size = 3.85,
    color = "#374151"
  ) +
  ggplot2::annotate(
    "text",
    x = 8.7,
    y = 1.25,
    label = "Practical readout: a per-gene/per-scaffold origin table plus DESeq2 sensitivity plots with host-like, fungal-like, and ambiguous genes separated.",
    size = 3.35,
    color = "#4B5563"
  ) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_color_identity() +
  ggplot2::coord_cartesian(xlim = c(0.0, 17.3), ylim = c(0.95, 8.85),
                           expand = FALSE, clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(20, 24, 18, 24)
  )

png_path <- file.path(out_dir, "unplaced_scaffold_origin_workflow.png")
svg_path <- file.path(out_dir, "unplaced_scaffold_origin_workflow.svg")
pdf_path <- file.path(out_dir, "unplaced_scaffold_origin_workflow.pdf")

ggplot2::ggsave(png_path, workflow_plot, width = 18, height = 10,
                dpi = 350, bg = "white", limitsize = FALSE)
ggplot2::ggsave(svg_path, workflow_plot, width = 18, height = 10,
                bg = "white", limitsize = FALSE)
ggplot2::ggsave(pdf_path, workflow_plot, width = 18, height = 10,
                bg = "white", limitsize = FALSE)

write.csv(counts,
          file.path(out_dir, "unplaced_scaffold_counts_by_context.csv"),
          row.names = FALSE)
write.csv(data.frame(source_dir = source_dir, output_dir = out_dir),
          file.path(out_dir, "unplaced_scaffold_origin_workflow_manifest.csv"),
          row.names = FALSE)

message("Wrote workflow figure to: ", out_dir)
