library(tidyverse)
library(DT)

# ---------------------------
# FILE PATHS
# ---------------------------

outDir <- "/scratch/ebaker48/transcriptome/mehreen/DEseq2_Mehreen_Interaction_Model"

master_file <- file.path(
  outDir,
  "DESeq2_MASTER_TABLE_AllContrasts.csv"
)

annotation_file <- file.path(
  outDir,
  "Gene_Annotations.csv"
)

# ---------------------------
# LOAD FILES
# ---------------------------

master <- read.csv(
  master_file,
  stringsAsFactors = FALSE
)

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
      "GeneType.y",
      "Species",
      "NumSig"
    ))
  )

# ---------------------------
# REORDER COLUMNS
# ---------------------------

other_cols <- setdiff(
  names(final_table),
  c("GeneID","Description")
)

final_table <- final_table %>%
  select(
    GeneID,
    Description,
    all_of(other_cols)
  )

# ---------------------------
# WRITE OUTPUT
# ---------------------------

write.csv(
  final_table,
  file.path(
    outDir,
    "DESeq2_FINAL_Annotated_Table.csv"
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