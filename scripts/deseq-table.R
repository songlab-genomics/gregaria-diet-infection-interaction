library(tidyverse)
library(DT)

base_dir <- "/scratch/ebaker48/transcriptome/mehreen"

parent_dirs <- c(
  Treatment_within_Diet_DEG = file.path(base_dir,"Treatment_within_Diet_DEG"),
  Diet_Pairwise_DEG_PCA = file.path(base_dir,"Diet_Pairwise_DEG_PCA"),
  DEG_ControlOnly_DietEffects = file.path(base_dir,"DEG_ControlOnly_DietEffects"),
  DEG_InfectedOnly_DietEffects = file.path(base_dir,"DEG_InfectedOnly_DietEffects")
)

all_results <- list()

# ----------------------------
# LOAD EVERY CSV AS A SUBTEST
# ----------------------------

for(parent_name in names(parent_dirs)){
  
  files <- list.files(
    parent_dirs[[parent_name]],
    pattern="FULL.*\\.csv$|DESeq2.*\\.csv$",
    recursive=TRUE,
    full.names=TRUE
  )
  
  for(f in files){
    
    df <- read.csv(f, stringsAsFactors=FALSE)
    
    # identify gene column
    if("gene" %in% names(df)){
      genes <- df$gene
    } else if("GeneID" %in% names(df)){
      genes <- df$GeneID
    } else if("X" %in% names(df)){
      genes <- df$X
    } else {
      next
    }
    
    df$Gene <- genes
    
    df$Status <- "NS"
    
    df$Status[
      !is.na(df$padj) &
        df$padj < 0.05 &
        df$log2FoldChange > 1
    ] <- "Up"
    
    df$Status[
      !is.na(df$padj) &
        df$padj < 0.05 &
        df$log2FoldChange < -1
    ] <- "Down"
    
    rel_path <- gsub(parent_dirs[[parent_name]], "", f)
    rel_path <- gsub("^/", "", rel_path)
    cname <- gsub("\\.csv$","", rel_path)
    cname <- gsub("/","__", cname)
    
    cname <- paste(parent_name, cname, sep="__")
    
    all_results[[cname]] <- df[,c("Gene","Status")]
    
    
  }
}

# ----------------------------
# MASTER GENE LIST
# ----------------------------

all_genes <- unique(
  unlist(lapply(all_results, function(x) x$Gene))
)

summary_table <- data.frame(
  GeneID = all_genes,
  stringsAsFactors=FALSE
)

# ----------------------------
# FILL EACH COLUMN WITH UP/DOWN/NS/ABSENT
# ----------------------------

for(test in names(all_results)){
  
  sub <- all_results[[test]]
  
  status_vec <- rep("Absent", nrow(summary_table))
  
  m <- match(summary_table$GeneID, sub$Gene)
  
  status_vec[!is.na(m)] <- sub$Status[m[!is.na(m)]]
  
  summary_table[[test]] <- status_vec
}

# ----------------------------
# TEST COLUMNS
# ----------------------------

test_cols <- setdiff(names(summary_table), c("GeneID", "NumSig"))

# ensure clean data frame
sig_df <- summary_table[, test_cols, drop = FALSE]
sig_df[] <- lapply(sig_df, as.character)

# ----------------------------
# NUMSIG (robust + safe)
# ----------------------------

summary_table$NumSig <- apply(sig_df, 1, function(x) {
  sum(x %in% c("Up","Down"))
})

# ----------------------------
# OPTIONAL: ADD STRONG SUMMARY METRICS
# ----------------------------

summary_table$UpCount <- apply(sig_df, 1, function(x) {
  sum(x == "Up")
})

summary_table$DownCount <- apply(sig_df, 1, function(x) {
  sum(x == "Down")
})

summary_table$DirectionalConsistency <- apply(sig_df, 1, function(x) {
  
  x <- x[x %in% c("Up","Down")]
  
  if(length(x) == 0) return(NA)
  
  if(all(x == "Up")) return("All_Up")
  if(all(x == "Down")) return("All_Down")
  
  return("Mixed")
})

# ----------------------------
# ORDERING (IMPORTANT FOR VIEWING)
# ----------------------------

summary_table <- summary_table %>%
  arrange(desc(NumSig), desc(UpCount + DownCount))

datatable(
  summary_table,
  options = list(
    pageLength = 50,
    scrollX = TRUE
  ),
  rownames = FALSE
) %>%
  formatStyle(
    columns = test_cols,
    backgroundColor = styleEqual(
      c("Up","Down","NS","Absent"),
      c("lightcoral","lightblue","grey85","white")
    )
  )
