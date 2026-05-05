library(tidyverse)

# ---------------------------
# INPUT FILES
# ---------------------------
file1 <- "analysis/Control_50_83_Combined_Annotated.csv"
file2 <- "analysis/Infected_50_83_Combined_Annotated.csv"
file3 <- "analysis/Within_50_83_Combined_Annotated.csv"

# ---------------------------
# FUNCTION TO PROCESS A FILE
# ---------------------------
process_file <- function(file, gene_col_name) {
  read_csv(file, show_col_types = FALSE) %>%

    # Remove "uncharacterized"
    filter(!str_detect(Description, regex("uncharacterized", ignore_case = TRUE))) %>%

    # Clean Description
    mutate(Description_clean = trimws(tolower(Description))) %>%

    # Rename GeneID column
    rename(!!gene_col_name := GeneID) %>%

    # Collapse GeneIDs per Description
    group_by(Description_clean) %>%
    summarise(
      !!gene_col_name := paste(unique(.data[[gene_col_name]]), collapse = ";"),
      .groups = "drop"
    )
}

# ---------------------------
# PROCESS EACH FILE
# ---------------------------
df1 <- process_file(file1, "GeneID_File1")
df2 <- process_file(file2, "GeneID_File2")
df3 <- process_file(file3, "GeneID_File3")

# ---------------------------
# JOIN ALL THREE
# ---------------------------
output <- df1 %>%
  inner_join(df2, by = "Description_clean") %>%
  inner_join(df3, by = "Description_clean")

# ---------------------------
# RESTORE ORIGINAL Description (from file1)
# ---------------------------
desc_lookup <- read_csv(file1, show_col_types = FALSE) %>%
  filter(!str_detect(Description, regex("uncharacterized", ignore_case = TRUE))) %>%
  mutate(Description_clean = trimws(tolower(Description))) %>%
  select(Description_clean, Description) %>%
  distinct()

output <- output %>%
  left_join(desc_lookup, by = "Description_clean") %>%
  select(Description, GeneID_File1, GeneID_File2, GeneID_File3)

# ---------------------------
# WRITE OUTPUT
# ---------------------------
write_csv(output, "collapsed_GeneIDs_3files.csv")

cat("Done! Output written to collapsed_GeneIDs_3files.csv\n")
