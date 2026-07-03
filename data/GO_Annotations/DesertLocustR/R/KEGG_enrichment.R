# KEGG_enrichment.R

#' KEGG_enrichment
#'
#' A small function that performs GO enrichment.
#' Built upon the clusterProfiler enrichKEGG function
#'
#'  @param gene_vector character vector of genes to perform enrichment on.
#'  @param KEGGTERM2LOC Mapping of KEGG terms to gene ID. This is stored in the object list returned by prep_data() included in this package.
#'  @param pval ajusted p-value to threshold enriched terms on, same as pvalueCutoff from clusterProfiler
#'  @param qval FDR to threshold enriched terms on, same as qvalueCutoff from clusterProfiler

#'  @return clusterProfiler enrichment result
#'
#'  @examples
#'  KEGG_enrichment("dge_genes", prep_data.list$KEGGTERm2LOC,0.05, 0.2)
#'
#'  @export
KEGG_enrichment <- function(gene_vector, KEGGTERM2LOC, pval, qval){
  gene_mapping.df <- data.frame(X.query = gene_vector)
  dge_with_kegg_ids <- left_join(gene_mapping.df, KEGGTERM2LOC, by = "X.query")
  KEGG_ids <- dge_with_kegg_ids$KEGG_ko[grepl("^K", dge_with_kegg_ids$KEGG_ko)]
  KEGGenrichment_results <- enrichKEGG(KEGG_ids,
                                       organism = "ko",
                                       pvalueCutoff = pval,
                                       qvalueCutoff = qval,
                                       pAdjustMethod = "BH")
  return(KEGGenrichment_results)
}
