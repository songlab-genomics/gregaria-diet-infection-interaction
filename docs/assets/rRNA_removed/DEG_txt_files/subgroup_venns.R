library(tidyverse)

file1 <- "/scratch/ebaker48/transcriptome/mehreen/rRNA_removed/DEG_txt_files/infection_only_no_rRNA_SVA/SIG_33_vs_50.txt"
file2 <- "/scratch/ebaker48/transcriptome/mehreen/rRNA_removed/DEG_txt_files/infection_only_no_rRNA_SVA/SIG_33_vs_83.txt"
file3 <- "/scratch/ebaker48/transcriptome/mehreen/rRNA_removed/DEG_txt_files/infection_only_no_rRNA_SVA/SIG_50_vs_83.txt"

read_genes <- function(f) {
  x <- readLines(f)
  x <- trimws(x)
  x[x != ""]
}

clean_genes <- function(x) {
  x <- trimws(x)
  x <- x[!is.na(x)]
  x <- x[x != ""]
  unique(x)
}

set1 <- unique(read_genes(file1))
set2 <- unique(read_genes(file2))
set3 <- unique(read_genes(file3))

set1 <- clean_genes(set1)
set2 <- clean_genes(set2)
set3 <- clean_genes(set3)

length(set1)
length(set2)
length(set3)

venn_list <- list(
  '33_vs_50' = set1,
  '33_vs_83' = set2,
  '50_vs_83' = set3
)

graphics.off()

library(VennDiagram)
library(grid)

venn_list <- list(
  `33_vs_50` = set1,
  `33_vs_83` = set2,
  `50_vs_83` = set3
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("#4E79A7", "#E15759", "#59A14F"),
  alpha = 0.5,
  cex = 1.3,
  cat.cex = 1.2,
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)

