library(tidyverse)

# ---------------------------
# PATHS
# ---------------------------

outDir <- "/scratch/ebaker48/transcriptome/mehreen/rRNA_removed"
annotDir <- "/scratch/ebaker48/transcriptome/mehreen/reference"

file_50 <- file.path(outDir, "Treatment_within_Diet_DEG_SVA/Diet_50/SIG_Diet_50_no_rRNA.csv")
file_83 <- file.path(outDir, "Treatment_within_Diet_DEG_SVA/Diet_83/SIG_Diet_83_no_rRNA.csv")

annotation_file <- file.path(annotDir, "gene_annotations.csv")

# ---------------------------
# LOAD DATA
# ---------------------------

df50 <- read.csv(file_50, stringsAsFactors = FALSE) %>%
  rename(GeneID = gene)

df83 <- read.csv(file_83, stringsAsFactors = FALSE) %>%
  rename(GeneID = gene)

annotation_df <- read.csv(annotation_file, stringsAsFactors = FALSE)

# ---------------------------
# CLEAN IDS
# ---------------------------

df50$GeneID <- trimws(df50$GeneID)
df83$GeneID <- trimws(df83$GeneID)
annotation_df$GeneID <- trimws(annotation_df$GeneID)

# ---------------------------
# ADD DIRECTION
# ---------------------------

df50 <- df50 %>%
  mutate(
    Direction_50 = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Neutral"
    )
  ) %>%
  select(GeneID, log2FC_50 = log2FoldChange, padj_50 = padj, Direction_50)

df83 <- df83 %>%
  mutate(
    Direction_83 = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Neutral"
    )
  ) %>%
  select(GeneID, log2FC_83 = log2FoldChange, padj_83 = padj, Direction_83)

# ---------------------------
# MERGE
# ---------------------------

combined <- full_join(df50, df83, by = "GeneID")

# ---------------------------
# ADD ANNOTATION
# ---------------------------

annotation_df <- annotation_df %>%
  select(GeneID, GeneType = Biotype, Description) %>%
  distinct()

final_table <- combined %>%
  left_join(annotation_df, by = "GeneID")

# ---------------------------
# OPTIONAL FILTER
# ---------------------------

final_table <- final_table %>%
  filter(Direction_50 != "Neutral" | Direction_83 != "Neutral")

# ---------------------------
# ADD PATTERN COLUMN
# ---------------------------

final_table <- final_table %>%
  mutate(
    Pattern = paste(Direction_50, Direction_83, sep = "_")
  )

# ---------------------------
# ORDER COLUMNS
# ---------------------------

final_table <- final_table %>%
  select(GeneID, Description, Pattern, everything())

# ---------------------------
# WRITE CSV
# ---------------------------

write.csv(
  final_table,
  file.path(outDir, "Within_50_83_Combined_Annotated.csv"),
  row.names = FALSE
)