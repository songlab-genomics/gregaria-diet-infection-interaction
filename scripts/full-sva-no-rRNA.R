library(DESeq2)
library(tidyverse)
library(sva)
library(uwot)
library(EnhancedVolcano)
library(pheatmap)
library(SARTools)

workDir <- "/scratch/ebaker48/transcriptome/mehreen"
rawDir <- file.path(workDir, "03-GCF_023897955.1_iqSchGreg1.2-DESeq2")
outDir <- file.path(workDir, "DEseq2_Mehreen_Interaction_Model_No_rRNA_NO_OUTLIERS")

dir.create(outDir, showWarnings = FALSE)

target <- loadTargetFile(
  targetFile = file.path(workDir, "mehreen_metadata_no_outliers.txt"),
  varInt = "Treatment",
  condRef = "Control",
  batch = NULL
)

#######################
#   HELPER FUNCTIONS
#######################

# ----------------------
# MA PLOT FUNCTION (MULTI)
# ----------------------

make_ma_plots <- function(res_df, name, outDir) {
  
  res_df <- res_df %>% filter(!is.na(padj))
  
  # classify genes
  res_df$group <- "NotSig"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Up"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Down"
  
  plot_ma <- function(df, suffix) {
    
    p <- ggplot(df,
                aes(x = log10(baseMean + 1),
                    y = log2FoldChange)) +
      geom_point(aes(color = group), size = 1.2, alpha = 0.7) +
      scale_color_manual(values = c(
        "Up" = "red",
        "Down" = "blue",
        "NotSig" = "grey70"
      )) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_classic() +
      labs(title = paste0(name, " - ", suffix),
           x = "log10(baseMean)",
           y = "log2FoldChange")
    
    ggsave(
      file.path(outDir, paste0("MA_", name, "_", suffix, ".png")),
      p, width = 7, height = 5, dpi = 300
    )
  }
  
  # ALL genes
  plot_ma(res_df, "ALL")
  
  # TOP 500 by baseMean
  top500 <- res_df %>%
    arrange(desc(baseMean)) %>%
    slice(1:500)
  
  plot_ma(top500, "TOP500")
  
  # TOP 1000 by baseMean
  top1000 <- res_df %>%
    arrange(desc(baseMean)) %>%
    slice(1:1000)
  
  plot_ma(top1000, "TOP1000")
}

# ----------------------
# HEATMAP TOP 500 DEGs
# ----------------------

make_deg_heatmap <- function(res_df, expr, target, name, outDir) {
  
  sig <- res_df %>%
    filter(!is.na(padj), padj < 0.05) %>%
    arrange(padj) %>%
    slice(1:500)
  
  genes <- sig$gene
  genes <- intersect(genes, rownames(expr))
  
  if(length(genes) < 10){
    message("Not enough DEGs for heatmap: ", name)
    return(NULL)
  }
  
  mat <- expr[genes, , drop = FALSE]
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), ]
  mat_scaled <- t(scale(t(mat)))
  
  annotation_col <- data.frame(
    Diet = target$Diet,
    Treatment = target$Treatment
  )
  rownames(annotation_col) <- colnames(mat_scaled)
  
  pheatmap(
    mat_scaled,
    show_rownames = FALSE,
    annotation_col = annotation_col,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = paste0("Top 500 DEGs - ", name),
    filename = file.path(outDir, paste0("Heatmap_TOP500_", name, ".pdf"))
  )
}

###################################

counts <- loadCountData(target = target, rawDir = rawDir)

target$Treatment <- factor(target$Treatment)
target$Diet <- factor(target$Diet, levels = c("50","33","83"))  # Diet 50 reference

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = target,
  design = ~ 1
)

dds <- dds[rowSums(counts(dds)) > 10, ]

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
  pull(ID) %>%         # <-- FIXED (was gene_id)
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

vsd <- vst(dds, blind = TRUE)
expr <- assay(vsd)

mod  <- model.matrix(~ Diet * Treatment, colData(dds))
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
                       "Diet * Treatment"),
                     collapse = " + "))
  )
} else {
  design(dds) <- ~ Diet * Treatment
}

dds <- DESeq(dds, fitType = "local")

resultsNames(dds)

res_infect <- results(dds, name = "Treatment_Infected_vs_Control")
res_diet33 <- results(dds, name = "Diet_33_vs_50")
res_diet83 <- results(dds, name = "Diet_83_vs_50")

int33 <- results(dds, name = "Diet33.TreatmentInfected")
int83 <- results(dds, name = "Diet83.TreatmentInfected")

sig_int <- unique(c(
  rownames(int33)[which(int33$padj < 0.05)],
  rownames(int83)[which(int83$padj < 0.05)]
))

write.table(sig_int,
            file.path(outDir, "Interaction_Genes_DietxInfection.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
inf_50 <- res_infect
inf_33 <- results(dds, contrast = list("Treatment_Infected_vs_Control",
                                       "Diet33.TreatmentInfected"))
inf_83 <- results(dds, contrast = list("Treatment_Infected_vs_Control",
                                       "Diet83.TreatmentInfected"))
write.csv(as.data.frame(res_infect),
          file.path(outDir, "Infection_effect_Diet50.csv"))

write.csv(as.data.frame(int33),
          file.path(outDir, "Interaction_Diet33.csv"))

write.csv(as.data.frame(int83),
          file.path(outDir, "Interaction_Diet83.csv"))

res33 <- as.data.frame(int33)
res33$gene <- rownames(res33)

p_int33 <- EnhancedVolcano(
  res33,
  lab = res33$gene,
  x = "log2FoldChange",
  y = "padj",
  title = "Interaction: Diet 33 vs 50 (Infection response shift)",
  pCutoff = 0.05,
  FCcutoff = 1
)

ggsave(file.path(outDir, "Volcano_Interaction_Diet33.png"),
       p_int33, width = 16, height = 12, dpi = 300)

res83 <- as.data.frame(int83)
res83$gene <- rownames(res83)

p_int83 <- EnhancedVolcano(
  res83,
  lab = res83$gene,
  x = "log2FoldChange",
  y = "padj",
  title = "Interaction: Diet 83 vs 50 (Infection response shift)",
  pCutoff = 0.05,
  FCcutoff = 1
)

ggsave(file.path(outDir, "Volcano_Interaction_Diet83.png"),
       p_int83, width = 16, height = 12, dpi = 300)


all_int_genes <- unique(c(
  rownames(int33)[which(int33$padj < 0.05)],
  rownames(int83)[which(int83$padj < 0.05)]
))

all_int_genes <- intersect(all_int_genes, rownames(expr))

mat <- expr[all_int_genes, , drop = FALSE]

mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), ]
mat_scaled <- t(scale(t(mat)))

annotation_col <- data.frame(
  Diet = target$Diet,
  Treatment = target$Treatment
)
rownames(annotation_col) <- colnames(mat_scaled)

ann_colors <- list(
  Diet = c("33" = "#B31B21", "50" = "#1465AC", "83" = "#219a25"),
  Treatment = c("Control" = "grey40", "Infected" = "black")
)

pheatmap(
  mat_scaled,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Interaction DEGs (Diet × Infection)",
  filename = file.path(outDir, "Heatmap_Interaction_DEGs.pdf")
)

# ----------------------
# PCA
# ----------------------

pca <- prcomp(t(expr), center = TRUE, scale. = FALSE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Diet = target$Diet,
  Treatment = target$Treatment
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Diet, shape = Treatment)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title = "PCA (VST expression)")

ggsave(
  filename = file.path(outDir, "PCA_VST.png"),
  plot = p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

# ----------------------
# UMAP
# ----------------------

set.seed(123)

umap_res <- uwot::umap(
  t(expr),
  n_neighbors = min(15, ncol(expr) - 1),
  min_dist = 0.3,
  metric = "euclidean"
)

umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  Diet = target$Diet,
  Treatment = target$Treatment
)

p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Diet, shape = Treatment)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title = "UMAP (VST expression)")

ggsave(
  filename = file.path(outDir, "UMAP_VST.png"),
  plot = p_umap,
  width = 8,
  height = 6,
  dpi = 300
)

### MA Plots + Heatmaps
make_ma_plots(res33, "Interaction_Diet33", outDir)
make_deg_heatmap(res33, expr, target, "Interaction_Diet33", outDir)

make_ma_plots(res83, "Interaction_Diet83", outDir)
make_deg_heatmap(res83, expr, target, "Interaction_Diet83", outDir)

res_inf_df <- as.data.frame(res_infect)
res_inf_df$gene <- rownames(res_inf_df)

make_ma_plots(res_inf_df, "Infection", outDir)
make_deg_heatmap(res_inf_df, expr, target, "Infection", outDir)


# ----------------------
# MASTER DESEQ2 TABLE
# ----------------------

make_df <- function(res_obj, prefix){
  df <- as.data.frame(res_obj)
  df$GeneID <- rownames(df)
  
  df <- df %>%
    dplyr::select(
      GeneID,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    )
  
  colnames(df)[-1] <- paste0(prefix, "_", colnames(df)[-1])
  return(df)
}

# main effects
df_infect <- make_df(res_infect, "INFECT")
df_diet33 <- make_df(res_diet33, "DIET33")
df_diet83 <- make_df(res_diet83, "DIET83")

# interactions
df_int33 <- make_df(int33, "INT33")
df_int83 <- make_df(int83, "INT83")

master_table <- Reduce(function(x, y){
  full_join(x, y, by = "GeneID")
}, list(df_infect, df_diet33, df_diet83, df_int33, df_int83))


master_table$GeneType <- NA
master_table$Description <- NA
master_table$Species <- "Schistocerca_gregaria"


sig_cols <- grep("padj", names(master_table), value = TRUE)

master_table$NumSig <- rowSums(
  master_table[, sig_cols] < 0.05,
  na.rm = TRUE
)

write.csv(
  master_table,
  file.path(outDir, "DESeq2_MASTER_TABLE_AllContrasts.csv"),
  row.names = FALSE
)