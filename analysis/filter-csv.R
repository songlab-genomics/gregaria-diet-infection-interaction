library(dplyr)

# Load full master table
df <- read.csv(
  "docs/assets/rRNA_removed/DESeq2_FINAL_Annotated_Table.csv",
  stringsAsFactors = FALSE
)

# Keep genes significant in infection AND at least one biased diet
filtered_genes <- df %>%
  filter(
    INFECT_padj < 0.05 &
      (
        DIET33_padj < 0.05 |
          DIET83_padj < 0.05
      )
  )

# Write filtered table with ALL original columns retained
write.csv(
  filtered_genes,
  "analysis/filtered_table_NO_rRNA.csv",
  row.names = FALSE
)

cat("Genes retained:", nrow(filtered_genes), "\n")
