# ============================================
# SYDNEY TRANSCRIPTOME DESEQ2 ANALYSIS WITH SVA,
# ============================================

# -------------------------------
# Load libraries
# -------------------------------
library(DESeq2)
library(tximport)
library(txdbmaker)
library(knitr)
library(rmdformats)
library(tidyverse)
library(data.table)
library(DT)
library(plotly)
library(ggthemes)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
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
library(limma)

# -------------------------------
# PARAMETERS
# -------------------------------
homeDir <- "/scratch/ebaker48/transcriptome/mehreen"
workDir <- homeDir
rawDir <- file.path(workDir, "03-GCF_023897955.1_iqSchGreg1.2-DESeq2") 
projectName <- "Mehreen_Transcriptome_Infection_Only"
author <- "Emily Baker"
targetFile <- file.path(workDir, "mehreen_metadata_fixed.txt")
colors <- c("#B31B21", "#1465AC", "#219a25")
shape_values <- c("Control" = 16, "Infected" = 17)
species <- "03-GCF_023897955.1_iqSchGreg1.2-DESeq2"  # Example species
varInt <- "Treatment"
condRef <- "Control"
batch <- NULL
fitType <- "parametric"
cooksCutoff <- TRUE
independentFiltering <- TRUE
typeTrans <- "rlog"
locfunc <- "median"
tresh_logfold <- 1
tresh_padj <- 0.05
alpha_DEseq2 <- 0.05
pAdjustMethod_DEseq2 <- "BH"
featuresToRemove <- c(NULL)

# -------------------------------
# Setup working directories
# -------------------------------
setwd(workDir)
Dirname <- paste0("DEseq2_", projectName)
dir.create(Dirname, showWarnings = FALSE)
setwd(Dirname)
workDir_DEseq2 <- getwd()

# -------------------------------
# Load data
# -------------------------------
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# -------------------------------
# Basic descriptive plots
# -------------------------------
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# Path and name of targetfile containing conditions and file names


# checking parameters
setwd(workDir_DEseq2)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# -------------------------------
# loading counts
# -------------------------------
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# -------------------------------
# CREATE DESEQ OBJECT
# -------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = target,
  design    = as.formula(paste("~", varInt))
)

# Prefilter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Normalize
dds <- estimateSizeFactors(dds)

# -------------------------------
# RUN SVA HERE
# -------------------------------
norm_counts <- counts(dds, normalized = TRUE)

mod  <- model.matrix(as.formula(paste("~", varInt)), data = colData(dds))
mod0 <- model.matrix(~ 1, data = colData(dds))

svobj <- sva(norm_counts, mod, mod0)

# Add surrogate variables
for (i in 1:svobj$n.sv) {
  colData(dds)[, paste0("SV", i)] <- svobj$sv[, i]
}

# Update design to include SVs BEFORE variable of interest
sv_names <- paste0("SV", 1:svobj$n.sv)

design(dds) <- as.formula(
  paste("~", paste(c(sv_names, varInt), collapse = " + "))
)

# -------------------------------
# RUN DESEQ2 WITH SVA INCLUDED
# -------------------------------
dds <- DESeq(dds)

rld <- rlog(dds, blind = FALSE)
mat <- assay(rld)
mat_t <- t(mat)

# UMAP
set.seed(123)
umap_res <- uwot::umap(mat_t)

umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  Treatment = colData(dds)$Treatment,
  Phase = colData(dds)$Phase
)

umap_plot <- ggplot(umap_df,
                    aes(x = UMAP1, y = UMAP2,
                        color = Phase,
                        shape = Treatment)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title = "UMAP (rlog + SVA)",
       color = "Diet",
       shape = "Treatment")

ggsave("UMAP_rlog_SVA.png",
       plot = umap_plot,
       width = 6,
       height = 5,
       dpi = 300)

res <- results(dds, alpha = alpha_DEseq2)

# -------------------------------
# PCA + Clustering
# -------------------------------

pcaData <- plotPCA(rld, intgroup = varInt, returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData,
                   aes(PC1, PC2,
                       color = Treatment,
                       shape = Phase)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  labs(title = "PCA (rlog + SVA)")

# Save PCA
ggsave("PCA_rlog_SVA.png",
       plot = pca_plot,
       width = 6,
       height = 5,
       dpi = 300)


# MA Plot
plotMA(res, ylim = c(-5,5))

# Volcano Plot
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1)


# Dispersion Plot
plotDispEsts(dds)

# Export Results Table 
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered),
          file="DESeq2_results_with_SVA.csv")


# Heatmap
sigGenes <- rownames(res)[which(res$padj < 0.05)]

mat <- assay(rld)[sigGenes, ]
mat <- mat - rowMeans(mat)

pheatmap(mat,
         annotation_col = as.data.frame(colData(dds)[, varInt, drop=FALSE]))

