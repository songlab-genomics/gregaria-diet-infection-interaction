# GO_enrichment.R

#' GO_enrichment
#'
#' A small function that performs GO enrichment.
#' Built upon the clusterProfiler enricher function
#'
#'  @param gene_vector character vector of genes to perform enrichment on.
#'  @param sub_ontology GO sub-onology to perform enrichment on. This is stored in the object list returned by prep_data() included in this package.
#'  @param term2name the dataframe linking gene names to GO terms. Required by enricher from clusterProfiler. This is stored in the object list returned by prep_data() included in this package
#'  @param pval ajusted p-value to threshold enriched terms on, same as pvalueCutoff from clusterProfiler
#'  @param qval FDR to threshold enriched terms on, same as qvalueCutoff from clusterProfiler

#'  @return clusterProfiler enrichment result
#'
#'  @examples
#'  GO_enrichment("dge_genes", prep_data.list$GO_BP_TERMS, prep_data.list$GOTERM2LOC, 0.05, 0.2)
#'
#'  @export
GO_enrichment <- function(gene_vector, sub_ontology, term2name, pval, qval){
  GOenrichment_results <- enricher(gene_vector, TERM2GENE = sub_ontology,
                                   TERM2NAME = term2name, pvalueCutoff = pval,
                                   pAdjustMethod = "BH", qvalueCutoff = qval)
  return(GOenrichment_results)
}
