library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(ggplot2)
library(uwot)

## PARAMETERS
workDir <- "/scratch/ebaker48/transcriptome/mehreen"
workDir_DEseq2 <- "/scratch/ebaker48/transcriptome/mehreen/03-GCF_023897955.1_iqSchGreg1.2-DESeq2"

target_file <- file.path(workDir, "mehreen_metadata_fixed.txt")
counts_file <- file.path(workDir_DEseq2, "master_counts.csv")

outDir <- file.path(workDir, "Diet_Pairwise_DEG_PCA")
dir.create(outDir, showWarnings = FALSE)

# ----------------------
# LOAD DATA
# ----------------------

target <- read.delim(target_file, header = TRUE)
counts <- read.csv(counts_file, row.names = 1)

target$Diet <- factor(target$Diet, levels = c("33","50","83"))

rownames(target) <- target$label


# ----------------------
# PAIRWISE COMPARISONS
# ----------------------

comparisons <- list(
  c("33", "50"),
  c("33", "83"),
  c("50", "83")
)

for (cmp in comparisons) {

  d1 <- cmp[1]
  d2 <- cmp[2]

  message("Processing: ", d1, " vs ", d2)

  # subset
  sel <- target$Diet %in% c(d1, d2)
  target_sub <- target[sel, ]
  counts_sub <- counts[, sel]

  target_sub$Diet <- droplevels(target_sub$Diet)

  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = target_sub,
    design = ~ Diet
  )

  dds <- dds[rowSums(counts(dds)) > 10, ]

  dds <- DESeq(dds, fitType = "local")

  # ----------------------
  # DEG RESULTS
  # ----------------------

  res <- results(dds, contrast = c("Diet", d1, d2))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)

  # ----------------------
# SAVE SIGNIFICANT DEGs TO TXT
# ----------------------

sig_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  pull(gene)

# SAVE SIGNIFICANT DEGs TO TXT (FIXED)
write.table(sig_genes,
            file.path(outDir,
                      paste0("Significant_DEGs_", d1, "_vs_", d2, ".txt")),
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


# ----------------------
# SAVE SIGNIFICANT DEGs TO CSV (filtered)
# ----------------------

deg_df <- res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::filter(padj < 0.05)

write.csv(deg_df,
          file.path(outDir,
                    paste0("Significant_DEGs_", d1, "_vs_", d2, ".csv")),
          row.names = FALSE)

  write.csv(res_df,
            file.path(outDir, paste0("DEG_", d1, "_vs_", d2, ".csv")),
            row.names = FALSE)

  # volcano
  p1 <- EnhancedVolcano(
    res_df,
    lab = res_df$gene,
    x = "log2FoldChange",
    y = "padj",
    title = paste(d1, "vs", d2)
  )

  ggsave(file.path(outDir, paste0("Volcano_", d1, "_vs_", d2, ".png")),
         p1, width = 10, height = 8, dpi = 300)

  # ----------------------
  # PCA (ONLY ADDITION)
  # ----------------------

  vsd <- vst(dds, blind = TRUE)

  pcaData <- plotPCA(vsd, intgroup = "Diet", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  p2 <- ggplot(pcaData, aes(PC1, PC2, color = Diet)) +
    geom_point(size = 4) +
    theme_classic() +
    labs(
      title = paste("PCA:", d1, "vs", d2),
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance")
    )

  ggsave(file.path(outDir, paste0("PCA_", d1, "_vs_", d2, ".png")),
         p2, width = 8, height = 6, dpi = 300)


# ----------------------
  # UMAP (ADDED)
  # ----------------------

  mat <- t(assay(vsd))

  set.seed(123)
  umap_res <- umap(mat,
                   n_neighbors = min(5, nrow(mat) - 1),
                   min_dist = 0.3)

  umap_df <- data.frame(
    UMAP1 = umap_res[,1],
    UMAP2 = umap_res[,2],
    Diet = target_sub$Diet
  )

  umap_plot <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Diet)) +
    geom_point(size = 4) +
    theme_classic() +
    labs(title = paste("UMAP:", d1, "vs", d2))

  ggsave(file.path(outDir, paste0("UMAP_", d1, "_vs_", d2, ".png")),
         umap_plot, width = 8, height = 6, dpi = 300)


# ----------------------
# MA PLOT
# ----------------------

  # ----------------------
  # CUSTOM MA PLOT (TOP 500 + SIGNIFICANCE COLORS)
  # ----------------------
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # remove NA
  res_df <- res_df %>% dplyr::filter(!is.na(padj))
  
  # define significance
  res_df$group <- "NotSig"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Up"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Down"
  
  # keep top 500 most expressed genes (baseMean = "A" in MA plot)
  top500 <- res_df %>%
    dplyr::arrange(desc(baseMean)) %>%
    dplyr::slice(1:500)
  
  # MA plot
  p_ma <- ggplot(top500, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    
    geom_point(aes(color = group), size = 1.5, alpha = 0.8) +
    
    scale_color_manual(values = c(
      "Up" = "red",
      "Down" = "blue",
      "NotSig" = "grey70"
    )) +
    
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    theme_classic() +
    
    labs(
      title = paste("MA Plot (Top 500 genes):", d1, "vs", d2),
      x = "log10(Base Mean + 1)",
      y = "log2 Fold Change",
      color = "Status"
    )
  
  ggsave(
    filename = file.path(outDir, paste0("MA_top500_", d1, "_vs_", d2, ".png")),
    plot = p_ma,
    width = 8,
    height = 6,
    dpi = 300
  )

# ----------------------
# SAMPLE DISTANCE HEATMAP
# ----------------------

library(ComplexHeatmap)
library(circlize)

# ----------------------
# SAMPLE DISTANCE HEATMAP (ComplexHeatmap)
# ----------------------

mat <- assay(vsd)

sample_dist <- as.matrix(as.dist(1 - cor(mat)))

# annotation
ha <- HeatmapAnnotation(
  Diet = target_sub$Diet,
  col = list(Diet = c("33" = "#1b9e77",
                      "50" = "#d95f02",
                      "83" = "#7570b3"))
)

# color scale
col_fun <- colorRamp2(
  c(min(sample_dist), mean(sample_dist), max(sample_dist)),
  c("#2166ac", "white", "#b2182b")
)

ht <- Heatmap(
  sample_dist,
  name = "1 - cor",
  col = col_fun,
  top_annotation = ha,
  row_names_side = "left",
  column_names_side = "top",
  clustering_distance_rows = as.dist(sample_dist),
  clustering_distance_columns = as.dist(sample_dist),
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  show_row_dend = TRUE,
  show_column_dend = TRUE
)

pdf(file.path(outDir,
              paste0("ComplexHeatmap_sample_distance_", d1, "_vs_", d2, ".pdf")),
    width = 7, height = 6)

draw(ht)

# ----------------------
# GENE EXPRESSION HEATMAP (TOP VARIABLE GENES)
# ----------------------

library(ComplexHeatmap)
library(circlize)

mat <- assay(vsd)

# ---- select top variable genes ----
gene_var <- apply(mat, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:50]

heat_mat <- mat[top_genes, ]

# z-score scaling per gene
heat_mat <- t(scale(t(heat_mat)))

# annotation
ha <- HeatmapAnnotation(
  Diet = target_sub$Diet,
  col = list(Diet = c("33" = "#1b9e77",
                      "50" = "#d95f02",
                      "83" = "#7570b3"))
)

# color scale
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166ac", "white", "#b2182b")
)

ht <- Heatmap(
  heat_mat,
  name = "Z-score",
  col = col_fun,
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "average",
  column_title = paste("Top Variable Genes:", d1, "vs", d2)
)

pdf(file.path(outDir,
              paste0("GeneHeatmap_Top50Var_", d1, "_vs_", d2, ".pdf")),
    width = 8, height = 10)

draw(ht)

dev.off()

}
