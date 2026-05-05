library(tidyverse)

# ---------------------------
# INPUT FILES
# ---------------------------
file1 <- "analysis/Control_50_83_Combined_Annotated.csv"
file2 <- "analysis/Infected_50_83_Combined_Annotated.csv"

# ---------------------------
# READ DATA
# ---------------------------
df1 <- read_csv(file1, show_col_types = FALSE)
df2 <- read_csv(file2, show_col_types = FALSE)

# ---------------------------
# REMOVE "uncharacterized"
# ---------------------------
df1 <- df1 %>%
  filter(!str_detect(Description, regex("uncharacterized", ignore_case = TRUE)))

df2 <- df2 %>%
  filter(!str_detect(Description, regex("uncharacterized", ignore_case = TRUE)))

# ---------------------------
# CLEAN + STANDARDIZE
# ---------------------------
df1 <- df1 %>%
  mutate(Description_clean = trimws(tolower(Description))) %>%
  rename(GeneID_File1 = GeneID)

df2 <- df2 %>%
  mutate(Description_clean = trimws(tolower(Description))) %>%
  rename(GeneID_File2 = GeneID)

# ---------------------------
# COLLAPSE EACH FILE FIRST
# ---------------------------
df1_collapsed <- df1 %>%
  group_by(Description_clean) %>%
  summarise(
    GeneID_File1 = paste(unique(GeneID_File1), collapse = ";"),
    .groups = "drop"
  )

df2_collapsed <- df2 %>%
  group_by(Description_clean) %>%
  summarise(
    GeneID_File2 = paste(unique(GeneID_File2), collapse = ";"),
    .groups = "drop"
  )

# ---------------------------
# JOIN COLLAPSED TABLES
# ---------------------------
output <- inner_join(df1_collapsed, df2_collapsed, by = "Description_clean")

# ---------------------------
# RESTORE ORIGINAL Description
# ---------------------------
desc_lookup <- df1 %>%
  select(Description_clean, Description) %>%
  distinct()

output <- output %>%
  left_join(desc_lookup, by = "Description_clean") %>%
  select(Description, GeneID_File1, GeneID_File2)

# ---------------------------
# WRITE OUTPUT
# ---------------------------
write_csv(output, "collapsed_GeneIDs_filtered.csv")

cat("Done! Output written to collapsed_GeneIDs_filtered.csv\n")
