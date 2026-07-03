library(DESeq2)
library(tximport)
library(txdbmaker)
library(tidyverse)
library(data.table)
library(DT)
library(plotly)
library(ggthemes)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(apeglm)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(SARTools)
library(pheatmap)
library(clusterProfiler)
library(sva)
library(cowplot)
library(ashr)
library(ggforce)
library(ggConvexHull)
library(uwot)
library(limma)

## PARAMETERS
homeDir <- "/scratch/ebaker48/transcriptome/mehreen"
workDir <- "/scratch/ebaker48/transcriptome/mehreen/DEG"
rawDir <- file.path(workDir, "03-GCF_023897955.1_iqSchGreg1.2-DESeq2")
projectName <- "Mehreen_Gregaria_Transcriptome_Interaction"

Dirname <- paste("DEseq2_", projectName, sep = "")
dir.create(file.path(workDir, Dirname), showWarnings = FALSE)
workDir_DEseq2 <- file.path(workDir, Dirname)
setwd(workDir_DEseq2)

targetFile <- file.path(workDir, "mehreen_metadata_fixed.txt")
colors <- c("#B31B21", "#1465AC", "#219a25")
shape_values <- c("Control" = 16, "Infected" = 17)

# Load target and counts
target <- loadTargetFile(targetFile=targetFile,
                         varInt="Treatment",
                         condRef="Control",
                         batch=NULL)

counts <- loadCountData(target=target,
                        rawDir=rawDir,
                        featuresToRemove=NULL)

# Ensure factors
target$Treatment <- factor(target$Treatment)
target$Diet <- factor(target$Diet, levels=c("33","50","83"))

# Description plots
majSequences <- descriptionPlots(counts=counts,
                                 group=target$Diet,
                                 col=colors)

# ----------------------
# DESEQ2 MODEL (needed for MA plot)
# ----------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = target,
  design = ~ Treatment + Diet
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds, fitType = "local")

# choose a contrast (EDIT THIS IF NEEDED)
res <- results(dds, contrast = c("Diet", "33", "50"))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# ----------------------
# CUSTOM MA PLOT (TOP 500 + COLOR BY SIGNIFICANCE)
# ----------------------

res_df <- res_df %>% dplyr::filter(!is.na(padj))

res_df$group <- "NotSig"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Up"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Down"

# top 500 expressed genes
top500 <- res_df %>%
  dplyr::arrange(desc(baseMean)) %>%
  dplyr::slice(1:500)

p_ma <- ggplot(top500,
               aes(x = log10(baseMean + 1),
                   y = log2FoldChange)) +
  
  geom_point(aes(color = group),
             size = 1.5,
             alpha = 0.8) +
  
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "NotSig" = "grey70"
  )) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic() +
  
  labs(
    title = "MA Plot (Top 500 genes)",
    x = "log10(Base Mean + 1)",
    y = "log2 Fold Change",
    color = "Status"
  )

ggsave(
  filename = file.path(workDir_DEseq2,
                       "MAplot_top500_custom.png"),
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 300
)