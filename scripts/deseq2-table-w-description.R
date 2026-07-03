library(tidyverse)
library(DT)

# ---------------------------
# FILE PATHS
# ---------------------------

outDir <- "/scratch/ebaker48/transcriptome/mehreen/rRNA_removed"
annotDir <- "/scratch/ebaker48/transcriptome/mehreen/reference"

master_file <- file.path(
  outDir,
  "DEG_InfectedOnly_DietEffects_SVA/33_vs_83/SIG_33_vs_83_no_rRNA.csv"
)

annotation_file <- file.path(
  annotDir,
  "gene_annotations.csv"
)

# ---------------------------
# LOAD FILES
# ---------------------------

master <- read.csv(
  master_file,
  stringsAsFactors = FALSE
)

master <- read.csv(
  master_file,
  stringsAsFactors = FALSE
) %>%
  rename(GeneID = gene)

annotation_df <- read.csv(
  annotation_file,
  stringsAsFactors = FALSE
)

# ---------------------------
# CLEAN IDS
# ---------------------------

master$GeneID <- trimws(master$GeneID)
annotation_df$GeneID <- trimws(annotation_df$GeneID)

# ---------------------------
# KEEP ONLY ANNOTATION COLUMNS YOU WANT
# ---------------------------

annotation_df <- annotation_df %>%
  select(
    GeneID,
    GeneType = Biotype,
    Description
  ) %>%
  distinct()

# ---------------------------
# MERGE ANNOTATION
# ---------------------------

final_table <- master %>%
  left_join(annotation_df, by="GeneID")

# ---------------------------
# FILTER |log2FC| > 1 IN ANY CONTRAST
# ---------------------------

fc_cols <- grep(
  "log2FoldChange",
  names(final_table),
  value=TRUE
)

final_table <- final_table %>%
  filter(
    if_any(
      all_of(fc_cols),
      ~ abs(.x) > 1
    )
  )

# ---------------------------
# REMOVE UNUSED COLUMNS
# ---------------------------

final_table <- final_table %>%
  select(
    -any_of(c(
      "GeneType.x",
      "Species",
      "NumSig"
    ))
  )

# ---------------------------
# REORDER COLUMNS (SAFE)
# ---------------------------

cols <- colnames(final_table)

if ("Description" %in% cols) {
  final_table <- final_table %>%
    select(GeneID, Description, everything())
} else {
  final_table <- final_table %>%
    select(GeneID, everything())
}

# ---------------------------
# WRITE OUTPUT
# ---------------------------

write.csv(
  final_table,
  file.path(
    outDir,
    "Infected_33_vs_83_Annotated_Table.csv"
  ),
  row.names=FALSE
)

# ---------------------------
# VIEW INTERACTIVE TABLE
# ---------------------------

datatable(
  final_table,
  options=list(
    pageLength=25,
    scrollX=TRUE
  ),
  rownames=FALSE
)