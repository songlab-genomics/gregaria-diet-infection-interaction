# prep_data.R

#' Lift Annotations from Eggnog & Create GO/KEGG Tables
#'
#' A function that replaces protein IDs from Schistocerca with LOCID from RefSeq Annotation.
#' It also creates GO/KEGG annotation tables for use in pathway enrichment.
#'
#' @param eggnog_file .annotaitons file output from the EggNog Mapper tool.
#' @param gff_file .gff file downloaded from NCBI FTP server
#' @return a named list of objects needed for downstream pathway enrichment
#'
#'  @examples
#'  lift_annotations("schisctocerca_greg.annotations", "GCF_023897955.1_iqSchGreg1.2_genomic.gff")
#'
#'  @export
lift_annotaitons <- function(eggnog_file, gff_file) {
  eggnog_annots <- read.delim(eggnog_file,
                              sep = "\t", skip = 4, header = T)
  eggnog_annots <- eggnog_annots[1:(nrow(eggnog_annots)-3),]
  gff.df <- as.data.frame( import(gff_file))

  protein_2_gene <- gff.df[, c("Name","gene")] #changed from 26 -> 27 as RefSeq GFFs have more columns pushing the gene column up
  uniq_protein_2_gene <- unique(protein_2_gene)
  sum(grepl("^XP", uniq_protein_2_gene$Name))
  protein_2_gene_df <- subset(uniq_protein_2_gene, grepl("^XP", uniq_protein_2_gene$Name))
  eggnog_annots$Name <- eggnog_annots$X.query
  eggnog_annots <- left_join(eggnog_annots, protein_2_gene_df, by = "Name") #join locID back to eggnog_annots as gene vector from DESeq2 has loc ids from RefSeq
  eggnog_annots$X.query <- eggnog_annots$gene

  GO_terms <- eggnog_annots[, c("X.query", "GOs")] #switched to isolate by names
  colnames(GO_terms) <- c("X.query", "GOs")
  GO_terms <- data.table(GO_terms)
  GO_terms <- GO_terms[, list(GOs = unlist(strsplit(GOs , ","))), by = X.query]
  go_final <- GO_terms[,c(2,1)]
  term2name = go_final
  term2name$Names <- mapIds(GO.db,
                            keys = term2name$GOs,
                            column = "TERM",
                            keytype = "GOID",
                            multiVals = "first"
  )
  term2name$Ontology <- mapIds(GO.db,
                               keys = term2name$GOs,
                               column = "ONTOLOGY",
                               keytype = "GOID",
                               multiVals = "first"
  )
  term2name <- as.data.frame(term2name)
  go_bp <- term2name[term2name$Ontology == "BP", c("GOs", "X.query")]
  go_cc <- term2name[term2name$Ontology == "CC", c("GOs", "X.query")]
  go_mf <- term2name[term2name$Ontology == "MF", c("GOs", "X.query")]
  term2name <- term2name[c(1,3)]

  KO_terms <- eggnog_annots[, c("X.query", "KEGG_ko")]
  KO_terms$KEGG_ko <- gsub( "ko:", "", as.character(KO_terms$KEGG_ko))
  KO_terms <- data.table(KO_terms)
  KO_terms <- KO_terms[, list(KEGG_ko = unlist(strsplit(KEGG_ko , ","))), by = X.query]
  kegg_final <- KO_terms[,c(2,1)]
  return(list(GOTERM2LOC = term2name, GO_BP_TERMS = go_bp, GO_CC_TERMS = go_cc,
              GO_MF_TERMS = go_mf, KEGGTERM2LOC = kegg_final))
}
