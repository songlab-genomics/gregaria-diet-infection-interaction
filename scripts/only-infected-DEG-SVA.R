library(DESeq2)
library(tidyverse)
library(sva)
library(uwot)
library(EnhancedVolcano)
library(pheatmap)
library(SARTools)

# ----------------------
# PARAMETERS
# ----------------------

workDir <- "/scratch/ebaker48/transcriptome/mehreen"
workDir_DEseq2 <- file.path(workDir,
                            "03-GCF_023897955.1_iqSchGreg1.2-DESeq2")

target_file <- file.path(workDir, "mehreen_metadata_no_outliers.txt")
counts_file <- file.path(workDir_DEseq2, "master_counts.csv")

outDir <- file.path(workDir, "rRNA_removed/DEG_InfectedOnly_DietEffects_SVA_NO_OUTLIERS")
dir.create(outDir, showWarnings = FALSE)

# ----------------------
# LOAD DATA
# ----------------------

target <- read.delim(target_file, header = TRUE)
counts <- read.csv(counts_file, row.names = 1)

target <- target %>% filter(Treatment == "Infected")

target$sample <- gsub("_counts\\.txt$", "", target$files)
rownames(target) <- target$sample

common <- intersect(colnames(counts), target$sample)

counts <- counts[, common]

target <- target %>%
  filter(sample %in% common)

target <- target[match(colnames(counts), target$sample), ]

target$Diet <- factor(target$Diet, levels = c("33","50","83"))

# ----------------------
# DESEQ2
# ----------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = target,
  design = ~ Diet
)

dds <- dds[rowSums(counts(dds)) > 10, ]

# ----------------------
# SVA
# ----------------------

vsd <- vst(dds, blind = TRUE)
expr <- assay(vsd)

mod  <- model.matrix(~ Diet, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

n.sv <- min(num.sv(expr, mod, method = "leek"), 2)

if(n.sv > 0){
  
  svobj <- sva(expr, mod, mod0, n.sv = n.sv)
  
  for(i in 1:n.sv){
    colData(dds)[[paste0("SV", i)]] <- svobj$sv[, i]
  }
  
}

if(n.sv > 0){
  design(dds) <- as.formula(
    paste("~", paste(c(paste0("SV", 1:n.sv),
                       "Diet"),
                     collapse = " + "))
  )
} else {
  design(dds) <- ~ Diet
}

dds <- DESeq(dds, fitType = "local")

# ----------------------
# REMOVE rRNA GENES
# ----------------------

library(rtracklayer)

gff_file <- file.path(
  workDir,
  "reference/GCF_023897955.1_iqSchGreg1.2_genomic.gff"
)

gff <- import(gff_file)
gff_df <- as.data.frame(gff)

# Restrict to gene features
gene_annot <- gff_df %>%
  filter(type == "gene")

# Identify rRNA genes
rrna_genes <- gene_annot %>%
  filter(
    grepl("rRNA|ribosomal RNA|18S|28S|5.8S|5S",
          gene_biotype,
          ignore.case = TRUE) |
    grepl("rRNA|ribosomal RNA|18S|28S|5.8S|5S",
          product,
          ignore.case = TRUE) |
    grepl("rRNA|ribosomal RNA|18S|28S|5.8S|5S",
          Name,
          ignore.case = TRUE)
  ) %>%
  pull(gene) %>%         # <-- FIXED (was gene_id)
  unique()

cat("rRNA genes found:", length(rrna_genes), "\n")

# Check overlap with count matrix IDs
cat("Matching genes in count matrix:",
    sum(rownames(dds) %in% rrna_genes), "\n")

# Remove rRNA genes
dds <- dds[
  !rownames(dds) %in% rrna_genes,
]

cat("Genes remaining:", nrow(dds), "\n")

# ----------------------
# FULL HEATMAP (TOP VARIABLE GENES)
# ----------------------

make_full_heatmap <- function(dds, vsd, outDir) {
  
  mat <- assay(vsd)
  
  # 🔹 Use top variable genes across ALL infected samples
  vars <- matrixStats::rowVars(mat)
  top_genes <- names(sort(vars, decreasing = TRUE))[1:100]
  
  heat_mat <- mat[top_genes, , drop = FALSE]
  
  annotation <- as.data.frame(colData(dds)[, "Diet", drop = FALSE])
  
  pheatmap::pheatmap(
    heat_mat,
    scale = "row",
    annotation_col = annotation,
    show_rownames = FALSE,
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    filename = file.path(outDir, "Heatmap_ALL_Infected_no_rRNA.pdf")
  )
}


# ----------------------
# VST FOR PLOTS
# ----------------------

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
make_full_heatmap(dds, vsd, outDir)

# ----------------------
# PCA (ALL SAMPLES)
# ----------------------

pca_data <- plotPCA(vsd, intgroup = "Diet", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Diet)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA - Infected Only",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  )

ggsave(file.path(outDir, "PCA_Infected_no_rRNA.png"), p_pca, width = 6, height = 5)

# ----------------------
# UMAP
# ----------------------

set.seed(123)
umap_res <- umap(t(mat))

umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  Diet = colData(dds)$Diet
)

p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Diet)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "UMAP - Infected Only")

ggsave(file.path(outDir, "UMAP_Infected_no_rRNA.png"), p_umap, width = 6, height = 5)

# ----------------------
# SAFE VOLCANO FUNCTION
# ----------------------

make_volcano <- function(res_df, name, od) { 
  p <- EnhancedVolcano( res_df, 
                        lab = res_df$gene,
                        x = "log2FoldChange",
                        y = "padj",
                        title = name,
                        pCutoff = 0.05,
                        FCcutoff = 1 ) 
  
  pdf(file.path(od, paste0("Volcano_", name, "_no_rRNA.pdf")),
      width = 8, 
      height = 6) 
  
  print(p) 
  
  dev.off() 
  
}

# ----------------------
# MA PLOT (TOP 500)
# ----------------------

make_ma <- function(res_df, name, od) {
  
  res_df <- res_df %>% filter(!is.na(padj))
  
  res_df$group <- "NotSig"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Up"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Down"
  
  top500 <- res_df %>%
    arrange(desc(baseMean)) %>%
    slice(1:500)
  
  p <- ggplot(top500,
              aes(x = log10(baseMean + 1),
                  y = log2FoldChange)) +
    geom_point(aes(color = group), size = 1.4, alpha = 0.8) +
    scale_color_manual(values = c(
      "Up" = "red",
      "Down" = "blue",
      "NotSig" = "grey70"
    )) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    labs(title = name)
  
  pdf(file.path(od, paste0("MA_", name, "_no_rRNA.pdf")), width = 8, height = 6)
  print(p)
  dev.off()
}

# ----------------------
# HEATMAP (TOP VARIABLE GENES)
# ----------------------

make_heatmap <- function(dds_subset, vsd_subset, res_df, name, od) {
  
  mat <- assay(vsd_subset)
  
  sig_genes <- res_df %>%
    filter(!is.na(padj), padj < 0.05) %>%
    arrange(padj) %>%
    slice(1:50) %>%
    pull(gene)
  
  if (length(sig_genes) < 10) {
    vars <- matrixStats::rowVars(mat)
    sig_genes <- names(sort(vars, decreasing = TRUE))[1:50]
  }
  
  heat_mat <- mat[sig_genes, , drop = FALSE]
  
  annotation <- as.data.frame(colData(dds_subset)[, "Diet", drop = FALSE])
  
  pheatmap::pheatmap(
    heat_mat,
    scale = "row",
    annotation_col = annotation,
    show_rownames = FALSE,
    filename = file.path(od, paste0("Heatmap_", name, "_no_rRNA.pdf"))
  )
}

# ----------------------
# DEG OUTPUT FUNCTION
# ----------------------

save_deg_outputs <- function(res, name, outDir) {
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  write.csv(res_df,
            file.path(outDir, paste0("FULL_", name, "_no_rRNA.csv")),
            row.names = FALSE)
  
  sig <- res_df %>%
    filter(!is.na(padj), padj < 0.05)
  
  write.csv(sig,
            file.path(outDir, paste0("SIG_", name, "_no_rRNA.csv")),
            row.names = FALSE)
  
  write.table(sig$gene,
              file.path(outDir, paste0("SIG_", name, ".txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ----------------------
# CLUSTER GRAPHS FUNCTION
# ----------------------

make_pca <- function(vsd_subset, name, od) {
  
  pca_data <- plotPCA(vsd_subset, intgroup = "Diet", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  p <- ggplot(pca_data, aes(PC1, PC2, color = Diet)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = paste0("PCA - ", name),
      x = paste0("PC1: ", percentVar[1], "%"),
      y = paste0("PC2: ", percentVar[2], "%")
    )
  
  ggsave(file.path(od, paste0("PCA_", name, "_no_rRNA.png")), p, width = 6, height = 5)
}

make_umap <- function(vsd_subset, name, od) {
  
  mat <- assay(vsd_subset)
  n_samples <- ncol(mat)
  
  # 🔹 Ensure valid n_neighbors
  n_neighbors <- min(10, n_samples - 1)
  
  # Optional: skip if too few samples
  if (n_samples < 3) {
    message("Skipping UMAP for ", name, " (too few samples)")
    return(NULL)
  }
  
  set.seed(123)
  umap_res <- umap(
    t(mat),
    n_neighbors = n_neighbors
  )
  
  umap_df <- data.frame(
    UMAP1 = umap_res[,1],
    UMAP2 = umap_res[,2],
    Diet = colData(vsd_subset)$Diet
  )
  
  p <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Diet)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = paste0("UMAP - ", name))
  
  ggsave(file.path(od, paste0("UMAP_", name, "_no_rRNA.png")), p, width = 6, height = 5)
}

# ----------------------
# CONTRAST PIPELINE
# ----------------------

run_contrast <- function(dds, contrast, outDir, vsd) {
  
  d1 <- contrast[2]
  d2 <- contrast[3]
  name <- paste0(d1, "_vs_", d2)
  
  od <- file.path(outDir, name)
  dir.create(od, showWarnings = FALSE)
  
  # Subset to only the two diets in this contrast
  keep <- colData(dds)$Diet %in% c(d1, d2)
  
  dds_subset <- dds[, keep]
  vsd_subset <- vsd[, keep]
  
  message("Running: ", name)
  
  res <- results(dds, contrast = contrast)
  
  save_deg_outputs(res, name, od)
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  make_volcano(res_df, name, od)
  make_ma(res_df, name, od)
  make_heatmap(dds_subset, vsd_subset, res_df, name, od)
  make_pca(vsd_subset, name, od)
  make_umap(vsd_subset, name, od)
  
  return(res)
}

# ----------------------
# RUN ALL CONTRASTS
# ----------------------

res_33_50 <- run_contrast(dds, c("Diet","33","50"), outDir, vsd)
res_33_83 <- run_contrast(dds, c("Diet","33","83"), outDir, vsd)
res_50_83 <- run_contrast(dds, c("Diet","50","83"), outDir, vsd)
