library(tidyverse)

# ---------------------------
# FILE PATHS
# ---------------------------

file1 <- "analysis/Control_50_83_Combined_Annotated.csv"
file2 <- "analysis/Infected_50_83_Combined_Annotated.csv"

out_file <- "non_overlapping_rows.csv"

# ---------------------------
# LOAD DATA
# ---------------------------

df1 <- read.csv(file1, stringsAsFactors = FALSE)
df2 <- read.csv(file2, stringsAsFactors = FALSE)

# ---------------------------
# CLEAN IDS
# ---------------------------

df1$GeneID <- trimws(df1$GeneID)
df2$GeneID <- trimws(df2$GeneID)

# ---------------------------
# TAG SOURCE
# ---------------------------

df1$Source <- "File1"
df2$Source <- "File2"

# ---------------------------
# GET NON-SHARED ROWS
# ---------------------------

only_df1 <- left_join(df1, df2, by = "GeneID")
only_df2 <- left_join(df2, df1, by = "GeneID")

diff_df <- bind_rows(only_df1, only_df2)

# ---------------------------
# WRITE OUTPUT
# ---------------------------

write.csv(diff_df, out_file, row.names = FALSE)
