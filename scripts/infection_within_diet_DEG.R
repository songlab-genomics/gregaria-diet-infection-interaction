library(DESeq2)
library(tidyverse)
library(sva)
library(uwot)
library(EnhancedVolcano)
library(pheatmap)
library(SARTools)

## PARAMETERS
workDir <- "/scratch/ebaker48/transcriptome/mehreen"
workDir_DEseq2 <- "/scratch/ebaker48/transcriptome/mehreen/03-GCF_023897955.1_iqSchGreg1.2-DESeq2"

target_file <- file.path(workDir, "mehreen_metadata_no_outliers.txt")
counts_file <- file.path(workDir_DEseq2, "master_counts.csv")

outDir <- file.path(workDir, "rRNA_removed/Treatment_within_Diet_DEG_NO_OUTLIERS")
dir.create(outDir, showWarnings = FALSE)

# ----------------------
# LOAD DATA
# ----------------------

target <- read.delim(target_file, header = TRUE)
counts <- read.csv(counts_file, row.names = 1)

target$Diet <- factor(target$Diet, levels = c("33","50","83"))
target$Treatment <- factor(target$Treatment, levels = c("Control","Infected"))

# ----------------------
# LOOP THROUGH DIETS
# ----------------------

for (d in levels(target$Diet)) {
  
  message("Processing Diet: ", d)
  
  od <- file.path(outDir, paste0("Diet_", d))
  dir.create(od, showWarnings = FALSE)
  
  # subset
  sel <- target$Diet == d
  target_sub <- target[sel, ]

 # Construct sample IDs matching count matrix columns
target_sub$sample_id <- paste0(
  "mehreen_",
  target_sub$label,
  "_MERGE"
)

  sample_ids <- target_sub$sample_id

  # Keep only matching samples
  sample_ids <- intersect(sample_ids, colnames(counts))

  counts_sub <- counts[, sample_ids, drop = FALSE]

  # Reorder metadata to match counts
  target_sub <- target_sub[
    match(colnames(counts_sub),
          target_sub$sample_id),
  ]

rownames(target_sub) <- target_sub$sample_id
  
  target_sub$Treatment <- droplevels(target_sub$Treatment)


  cat("\n====================\n")
  cat("Diet:", d, "\n")
  cat("====================\n")

  cat("target rows:", nrow(target_sub), "\n")
  cat("count cols :", ncol(counts_sub), "\n\n")

  cat("Metadata sample names:\n")
  print(target_sub$sample_id)

  cat("\nCount matrix column names:\n")
  print(colnames(counts_sub))

  cat("\nIntersection size:\n")
  print(length(intersect(target_sub$sample_id,
                        colnames(counts))))

  cat("\nSamples missing from counts:\n")
  print(setdiff(target_sub$sample_id,
                colnames(counts)))

  cat("\nSamples missing from metadata:\n")
  print(setdiff(colnames(counts),
        target_sub$sample_id))
  
  # ----------------------
  # DESEQ2
  # ----------------------

  cat("FINAL CHECK\n")
  cat("counts columns:", ncol(counts_sub), "\n")
  cat("metadata rows:", nrow(target_sub), "\n")

  print(colnames(counts_sub))
  print(rownames(target_sub))
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = target_sub,
    design = ~ Treatment
  )
  
  dds <- dds[rowSums(counts(dds)) > 10, ]


# ----------------------
# SVA
# ----------------------

  vsd <- vst(dds, blind = TRUE)
  expr <- assay(vsd)

  mod  <- model.matrix(~ Treatment, colData(dds))
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
                        "Treatment"),
                      collapse = " + "))
    )
  } else {
    design(dds) <- ~ Treatment
  }


   dds <- DESeq(dds, fitType = "local")
  
  # ----------------------
  # RESULTS
  # ----------------------
  
  res <- results(dds, contrast = c("Treatment", "Infected", "Control"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # ----------------------
  # SAVE TABLES
  # ----------------------
  
  write.csv(res_df,
            file.path(od, paste0("FULL_Diet_", d, "_no_rRNA.csv")),
            row.names = FALSE)
  
  sig <- res_df %>%
    filter(!is.na(padj), padj < 0.05)
  
  write.csv(sig,
            file.path(od, paste0("SIG_Diet_", d, "_no_rRNA.csv")),
            row.names = FALSE)
  
  write.table(sig$gene,
              file.path(od, paste0("SIG_Diet_", d, ".txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # ----------------------
  # VST (ONCE PER DIET)
  # ----------------------
  
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)
  
  # ----------------------
  # VOLCANO (SAFE)
  # ----------------------
  
  pdf(file.path(od, paste0("Volcano_Diet_", d, "_no_rRNA.pdf")), 8, 6)
  print(
    EnhancedVolcano(
      res_df,
      lab = res_df$gene,
      x = "log2FoldChange",
      y = "padj",
      title = paste("Diet", d, ": Infected vs Control"),
      pCutoff = 0.05,
      FCcutoff = 1
    )
  )
  dev.off()
  
  # ----------------------
  # MA PLOT (TOP 500)
  # ----------------------
  
  res_df <- res_df %>% filter(!is.na(padj))
  
  res_df$group <- "NotSig"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Up"
  res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Down"
  
  top500 <- res_df %>%
    arrange(desc(baseMean)) %>%
    slice(1:500)
  
  p_ma <- ggplot(top500,
                 aes(log10(baseMean + 1), log2FoldChange)) +
    geom_point(aes(color = group), size = 1.5) +
    scale_color_manual(values = c("Up"="red","Down"="blue","NotSig"="grey70")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    labs(title = paste("MA Plot Diet", d))
  
  ggsave(file.path(od, paste0("MA_Diet_", d, "_no_rRNA.png")),
         p_ma, width = 7, height = 6)
  
  # ----------------------
  # PCA
  # ----------------------
  
  pcaData <- plotPCA(vsd, intgroup = "Treatment", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  p_pca <- ggplot(pcaData, aes(PC1, PC2, color = Treatment)) +
    geom_point(size = 4) +
    theme_classic() +
    labs(
      title = paste("PCA Diet", d),
      x = paste0("PC1: ", percentVar[1], "%"),
      y = paste0("PC2: ", percentVar[2], "%")
    )
  
  ggsave(file.path(od, paste0("PCA_Diet_", d, "_no_rRNA.png")),
         p_pca, width = 7, height = 6)
  
  # ----------------------
  # UMAP
  # ----------------------
  
  set.seed(123)
  umap_res <- umap(t(mat), n_neighbors = min(5, ncol(mat)-1))
  
  umap_df <- data.frame(
    UMAP1 = umap_res[,1],
    UMAP2 = umap_res[,2],
    Treatment = target_sub$Treatment
  )
  
  p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Treatment)) +
    geom_point(size = 4) +
    theme_classic() +
    labs(title = paste("UMAP Diet", d))
  
  ggsave(file.path(od, paste0("UMAP_Diet_", d, "_no_rRNA.png")),
         p_umap, width = 7, height = 6)
  
  # ----------------------
  # HEATMAP (FIXED)
  # ----------------------
  
  top_var <- head(order(rowVars(mat), decreasing = TRUE), 50)
  heat_mat <- mat[top_var, ]
  
  annotation <- as.data.frame(colData(dds)[, "Treatment", drop = FALSE])
  
  pheatmap::pheatmap(
    heat_mat,
    scale = "row",
    annotation_col = annotation,
    show_rownames = FALSE,
    filename = file.path(od, paste0("Heatmap_Diet_", d, "_no_rRNA.pdf"))
  )
}