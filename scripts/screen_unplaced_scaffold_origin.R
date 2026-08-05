#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(forcats)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
repo_root <- if (!is.na(script_file)) normalizePath(file.path(dirname(script_file), "..")) else getwd()
if (!dir.exists(file.path(repo_root, "analysis"))) {
  repo_root <- normalizePath(getwd())
}

results_root <- file.path(repo_root, "output", "rmd_runs")
gff_file <- file.path(repo_root, "data", "reference",
                      "GCF_023897955.1_iqSchGreg1.2_genomic.gff")

latest_dir <- function(pattern) {
  dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(pattern, basename(dirs))]
  if (length(dirs) == 0) {
    stop("No output folder matching ", pattern, " under ", results_root, call. = FALSE)
  }
  dirs[order(file.info(dirs)$mtime, decreasing = TRUE)][[1]]
}

first_nonempty <- function(x, fallback = "") {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (length(x) == 0) fallback else x[[1]]
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- unique(x[!is.na(x) & x != "" & x != "NA"])
  if (length(x) == 0) "" else paste(sort(x), collapse = "; ")
}

extract_attr <- function(attributes, key) {
  out <- stringr::str_match(attributes, paste0("(?:^|;)", key, "=([^;]+)"))[, 2]
  utils::URLdecode(out)
}

read_gff_minimal <- function(path) {
  gff <- read.delim(path, comment.char = "#", header = FALSE, sep = "\t",
                    quote = "", stringsAsFactors = FALSE)
  names(gff)[1:9] <- c("seqid", "source", "type", "start", "end", "score",
                       "strand", "phase", "attributes")
  gff
}

classify_origin_evidence <- function(annotation_text, has_domain_annotation) {
  text <- stringr::str_to_lower(dplyr::coalesce(annotation_text, ""))
  nonhost_keyword <- stringr::str_detect(
    text,
    paste(
      "metarhizium|fungal|fungus|fungi|mycosis|yeast|aspergillus|fusarium",
      "beauveria|cordyceps|bacteria|bacterial|prokaryot|microbial|virus",
      "viral|phage|mycoplasma|bacillus|pseudomonas",
      sep = "|"
    )
  )
  repeat_keyword <- stringr::str_detect(
    text,
    paste(
      "transposon|retrotransposon|reverse transcriptase|integrase|gag-pol",
      "polyprotein|mobile element|tigger|mariner|gypsy|copia|\\bline\\b|\\bsine\\b",
      sep = "|"
    )
  )
  low_info <- stringr::str_detect(
    text,
    "uncharacterized|hypothetical|unknown protein|no gff description|predicted protein"
  )
  dplyr::case_when(
    nonhost_keyword ~ "Possible non-host keyword",
    repeat_keyword ~ "Repeat/mobile-element-like",
    low_info & !has_domain_annotation ~ "Low-information RefSeq gene",
    low_info & has_domain_annotation ~ "Uncharacterized but has domains",
    TRUE ~ "RefSeq annotated host-like"
  )
}

publication_dir <- latest_dir("^publication-figures_")
catalogue_dir <- latest_dir("^all-deg-catalogue_")

coordinate_file <- file.path(publication_dir, "figureS2_deg_coordinates_by_context.csv")
catalogue_file <- file.path(catalogue_dir, "all_samples_deg_catalogue.csv")

if (!file.exists(coordinate_file)) {
  stop("Missing coordinate file: ", coordinate_file, call. = FALSE)
}
if (!file.exists(catalogue_file)) {
  stop("Missing DEG catalogue file: ", catalogue_file, call. = FALSE)
}
if (!file.exists(gff_file)) {
  stop("Missing GFF file: ", gff_file, call. = FALSE)
}

run_id <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_dir <- file.path(results_root, paste0("unplaced-scaffold-origin-screen_", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading DEG coordinates: ", coordinate_file)
coordinates <- read.csv(coordinate_file, check.names = FALSE)
unplaced_memberships <- coordinates |>
  dplyr::filter(
    chromosome_class == "Unplaced scaffold" |
      stringr::str_detect(seqid, "^NW_")
  ) |>
  dplyr::distinct(DEG_context, gene_id, seqid, start, end, midpoint, strand,
                  seq_start, seq_end, chromosome, genome, chromosome_class)

message("Reading DEG catalogue: ", catalogue_file)
catalogue <- read.csv(catalogue_file, check.names = FALSE, na.strings = c("", "NA"))

catalogue_annotation <- catalogue |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(
    Description = first_nonempty(Description, "no GFF description"),
    gene_biotype = first_nonempty(gene_biotype, "unknown"),
    GOs = collapse_unique(GOs),
    KEGG_ko = collapse_unique(KEGG_ko),
    KEGG_Pathway = collapse_unique(KEGG_Pathway),
    PFAMs = collapse_unique(PFAMs),
    exact_deg_contrasts = collapse_unique(paste(comparison_type, comparison,
                                                signed_direction, sep = " | ")),
    .groups = "drop"
  )

message("Reading RefSeq GFF gene and scaffold metadata.")
gff <- read_gff_minimal(gff_file)

region_metadata <- gff |>
  dplyr::filter(type == "region") |>
  dplyr::transmute(
    seqid,
    region_name = extract_attr(attributes, "Name"),
    region_chromosome = extract_attr(attributes, "chromosome"),
    region_genome = extract_attr(attributes, "genome"),
    taxon = stringr::str_match(attributes, "Dbxref=taxon:([^;]+)")[, 2],
    isolate = extract_attr(attributes, "isolate"),
    mol_type = extract_attr(attributes, "mol_type"),
    tissue_type = extract_attr(attributes, "tissue-type"),
    region_length_bp = end - start + 1
  ) |>
  dplyr::distinct(seqid, .keep_all = TRUE)

gff_gene_metadata <- gff |>
  dplyr::filter(type == "gene") |>
  dplyr::transmute(
    seqid,
    source,
    gene_id = dplyr::coalesce(
      extract_attr(attributes, "gene"),
      extract_attr(attributes, "Name")
    ),
    gff_gene_id = stringr::str_match(attributes, "GeneID:([^,;]+)")[, 2],
    gff_description = extract_attr(attributes, "description"),
    gff_biotype = extract_attr(attributes, "gene_biotype"),
    gff_gbkey = extract_attr(attributes, "gbkey")
  ) |>
  dplyr::filter(!is.na(gene_id), gene_id != "") |>
  dplyr::mutate(
    gff_description = dplyr::coalesce(gff_description, "no GFF description"),
    gff_biotype = dplyr::coalesce(gff_biotype, "unknown")
  ) |>
  dplyr::distinct(gene_id, .keep_all = TRUE)

rm(gff)
gc(verbose = FALSE)

unplaced_screen <- unplaced_memberships |>
  dplyr::left_join(region_metadata, by = "seqid") |>
  dplyr::left_join(gff_gene_metadata, by = c("gene_id", "seqid")) |>
  dplyr::left_join(catalogue_annotation, by = "gene_id") |>
  dplyr::mutate(
    Description = dplyr::coalesce(Description, gff_description, "no GFF description"),
    gene_biotype = dplyr::coalesce(gene_biotype, gff_biotype, "unknown"),
    annotation_text = stringr::str_squish(paste(
      Description, GOs, KEGG_ko, KEGG_Pathway, PFAMs, exact_deg_contrasts,
      sep = " | "
    )),
    has_domain_annotation = (dplyr::coalesce(GOs, "") != "") |
      (dplyr::coalesce(KEGG_ko, "") != "") |
      (dplyr::coalesce(KEGG_Pathway, "") != "") |
      (dplyr::coalesce(PFAMs, "") != ""),
    origin_evidence_class = classify_origin_evidence(annotation_text,
                                                     has_domain_annotation),
    origin_screen_priority = dplyr::case_when(
      origin_evidence_class == "Possible non-host keyword" ~
        "Highest priority for sequence-level taxonomy",
      origin_evidence_class == "Repeat/mobile-element-like" ~
        "Ambiguous repeat-rich locus",
      origin_evidence_class %in% c("Low-information RefSeq gene",
                                   "Uncharacterized but has domains") ~
        "Needs sequence-level taxonomy",
      TRUE ~ "No annotation red flag"
    ),
    region_claim = dplyr::case_when(
      taxon == "7010" & region_genome == "genomic" ~
        "RefSeq S. gregaria unplaced genomic region",
      taxon == "7010" ~
        "RefSeq S. gregaria region",
      is.na(taxon) ~
        "No taxon parsed from region",
      TRUE ~ paste("Region taxon", taxon)
    )
  )

unplaced_gene_table <- unplaced_screen |>
  dplyr::group_by(gene_id, seqid) |>
  dplyr::summarise(
    region_claim = first_nonempty(region_claim),
    region_length_bp = dplyr::first(region_length_bp),
    taxon = first_nonempty(taxon),
    isolate = first_nonempty(isolate),
    mol_type = first_nonempty(mol_type),
    source = first_nonempty(source),
    gff_gene_id = first_nonempty(gff_gene_id),
    gene_biotype = first_nonempty(gene_biotype, "unknown"),
    Description = first_nonempty(Description, "no GFF description"),
    origin_evidence_class = first_nonempty(origin_evidence_class),
    origin_screen_priority = first_nonempty(origin_screen_priority),
    DEG_contexts = collapse_unique(DEG_context),
    exact_deg_contrasts = first_nonempty(exact_deg_contrasts),
    GOs = first_nonempty(GOs),
    KEGG_ko = first_nonempty(KEGG_ko),
    KEGG_Pathway = first_nonempty(KEGG_Pathway),
    PFAMs = first_nonempty(PFAMs),
    start = min(start, na.rm = TRUE),
    end = max(end, na.rm = TRUE),
    strand = first_nonempty(strand, "."),
    .groups = "drop"
  ) |>
  dplyr::arrange(origin_screen_priority, origin_evidence_class, seqid, gene_id)

unplaced_scaffold_summary <- unplaced_screen |>
  dplyr::group_by(seqid, region_claim, taxon, isolate, mol_type, region_length_bp) |>
  dplyr::summarise(
    n_deg_memberships = dplyr::n(),
    n_distinct_deg_genes = dplyr::n_distinct(gene_id),
    DEG_contexts = collapse_unique(DEG_context),
    evidence_classes = collapse_unique(origin_evidence_class),
    n_possible_nonhost = dplyr::n_distinct(gene_id[origin_evidence_class == "Possible non-host keyword"]),
    n_repeat_mobile = dplyr::n_distinct(gene_id[origin_evidence_class == "Repeat/mobile-element-like"]),
    n_low_information = dplyr::n_distinct(gene_id[origin_evidence_class %in%
                                                    c("Low-information RefSeq gene",
                                                      "Uncharacterized but has domains")]),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_distinct_deg_genes), seqid)

evidence_by_context <- unplaced_screen |>
  dplyr::distinct(DEG_context, gene_id, origin_evidence_class) |>
  dplyr::count(DEG_context, origin_evidence_class, name = "n_genes") |>
  dplyr::group_by(DEG_context) |>
  dplyr::mutate(
    context_total = sum(n_genes),
    percent = n_genes / context_total
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(DEG_context, dplyr::desc(n_genes))

region_claim_summary <- unplaced_gene_table |>
  dplyr::count(region_claim, taxon, isolate, mol_type, name = "n_unplaced_deg_genes") |>
  dplyr::arrange(dplyr::desc(n_unplaced_deg_genes))

write_csv(unplaced_screen,
          file.path(out_dir, "unplaced_deg_origin_screen_memberships.csv"))
write_csv(unplaced_gene_table,
          file.path(out_dir, "unplaced_deg_origin_screen_gene_table.csv"))
write_csv(unplaced_scaffold_summary,
          file.path(out_dir, "unplaced_scaffold_origin_screen_summary.csv"))
write_csv(evidence_by_context,
          file.path(out_dir, "unplaced_origin_evidence_by_context.csv"))
write_csv(region_claim_summary,
          file.path(out_dir, "unplaced_region_claim_summary.csv"))

unplaced_gene_bed <- unplaced_gene_table |>
  dplyr::transmute(
    seqid,
    bed_start = pmax(start - 1L, 0L),
    bed_end = end,
    name = paste(gene_id, origin_evidence_class, sep = "|"),
    score = 0,
    strand
  )

readr::write_tsv(
  unplaced_gene_bed,
  file.path(out_dir, "unplaced_deg_gene_intervals.bed"),
  col_names = FALSE
)

readr::write_lines(
  sort(unique(unplaced_gene_table$seqid)),
  file.path(out_dir, "unplaced_deg_scaffolds.txt")
)

top_scaffolds <- unplaced_scaffold_summary |>
  dplyr::slice_max(n_distinct_deg_genes, n = 20, with_ties = FALSE) |>
  dplyr::pull(seqid)

readr::write_lines(
  top_scaffolds,
  file.path(out_dir, "top_unplaced_deg_scaffolds.txt")
)

evidence_colors <- c(
  "RefSeq annotated host-like" = "#2C7A51",
  "Uncharacterized but has domains" = "#A3A948",
  "Low-information RefSeq gene" = "#D9A21B",
  "Repeat/mobile-element-like" = "#8C6BB1",
  "Possible non-host keyword" = "#C51B1D"
)

p_evidence <- evidence_by_context |>
  dplyr::mutate(
    DEG_context = forcats::fct_reorder(DEG_context, context_total),
    origin_evidence_class = factor(origin_evidence_class,
                                   levels = names(evidence_colors))
  ) |>
  ggplot(aes(x = DEG_context, y = n_genes, fill = origin_evidence_class)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  geom_text(aes(label = n_genes),
            position = position_stack(vjust = 0.5),
            color = "white",
            fontface = "bold",
            size = 4,
            check_overlap = TRUE) +
  coord_flip() +
  scale_fill_manual(values = evidence_colors, drop = TRUE) +
  labs(
    title = "Annotation screen for DEGs on unplaced scaffolds",
    subtitle = "Counts are distinct DEG genes per broad DEG context; no microbial/fungal or repeat keyword flags were detected",
    caption = "This is a triage screen based on RefSeq/GFF and existing eggNOG-style annotations, not taxonomic proof.",
    x = NULL,
    y = "Distinct DEG genes on unplaced scaffolds",
    fill = "Evidence"
  ) +
  theme_classic(base_size = 14) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.caption = element_text(color = "#4B5563"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

p_scaffolds <- unplaced_screen |>
  dplyr::filter(seqid %in% top_scaffolds) |>
  dplyr::distinct(seqid, gene_id, DEG_context) |>
  dplyr::count(seqid, DEG_context, name = "n_genes") |>
  dplyr::left_join(
    unplaced_scaffold_summary |> dplyr::select(seqid, n_distinct_deg_genes),
    by = "seqid"
  ) |>
  dplyr::mutate(seqid = forcats::fct_reorder(seqid, n_distinct_deg_genes)) |>
  ggplot(aes(x = seqid, y = n_genes, fill = DEG_context)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  coord_flip() +
  labs(
    title = "Top unplaced scaffolds contributing DEG genes",
    subtitle = "Stacked by broad DEG context",
    x = NULL,
    y = "Distinct DEG genes",
    fill = "DEG context"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".png")),
         plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".svg")),
         plot, width = width, height = height, bg = "white")
}

save_plot(p_evidence, "unplaced_origin_evidence_by_context", 14, 7.2)
save_plot(p_scaffolds, "top_unplaced_scaffolds_by_deg_context", 12, 8)

screen_manifest <- data.frame(
  key = c("output_dir", "coordinate_file", "catalogue_file", "gff_file",
          "n_unplaced_deg_memberships", "n_unplaced_deg_genes",
          "n_unplaced_scaffolds"),
  value = c(out_dir, coordinate_file, catalogue_file, gff_file,
            nrow(unplaced_screen),
            dplyr::n_distinct(unplaced_screen$gene_id),
            dplyr::n_distinct(unplaced_screen$seqid))
)
write_csv(screen_manifest,
          file.path(out_dir, "unplaced_scaffold_origin_screen_manifest.csv"))

message("Wrote unplaced scaffold origin screen to: ", out_dir)
message("Summary:")
print(screen_manifest)
message("")
message("Evidence classes by context:")
print(evidence_by_context)
