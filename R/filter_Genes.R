#' Filter lowly abundant genes.
#'
#' \code{filter_genes} filters genes, according to the overall number of counts (across all samples) compatible with the gene.
#' The filtering also applies to groups of genes with reads/fragments compatible with >1 gene;
#' in this case, the number of counts considered is across all genes in the group.
#' 
#' The function inputs a 'BANDITS_data' object, and returns again a 'BANDITS_data' object after filtering genes and groups of genes.
#'
#' @param BANDITS_data a 'BANDITS_data' object, created with the \code{\link{create_data}} function.
#' @param min_counts_per_gene the minimum number of counts compatible with a gene (across all samples).
#' 
#' @return A \code{\linkS4class{BANDITS_data}} object.
#' @examples
#' # load the pre-computed data:
#' data("input_data", package = "BANDITS")
#' input_data
#' 
#' # Filter lowly abundant genes:
#' input_data = filter_genes(input_data, min_counts_per_gene = 20)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#'  
#' @seealso \code{\link{filter_transcripts}}, \code{\link{create_data}}, \code{\linkS4class{BANDITS_data}}
#' 
#' @export
filter_genes  = function(BANDITS_data, min_counts_per_gene = 10){
  if( !is(BANDITS_data, "BANDITS_data") ){
    message("'BANDITS_data' must be a 'BANDITS_data' object created via the create_data function")
    return(NULL)
  }
  
  # I filter out lowly expressed genes: at least 1 count per sample and at least 11 counts per condition:
  tot_counts = vapply( counts(BANDITS_data), sum, FUN.VALUE = numeric(1) )
  # sapply( counts(BANDITS_data), sum)
  SEL = tot_counts >= min_counts_per_gene
  
  # filter genes/groups from BANDITS_data:
  if(mean(SEL) < 1){ # if mean(SEL) == 1, no genes/groups were filtered.
    BANDITS_data_new = new("BANDITS_data",
                           genes       = genes(BANDITS_data)[SEL], 
                           transcripts = transcripts(BANDITS_data)[SEL],
                           effLen      = effLen(BANDITS_data)[SEL],
                           classes     = classes(BANDITS_data)[SEL],
                           counts      = counts(BANDITS_data)[SEL], 
                           uniqueId    = uniqueId(BANDITS_data)[SEL],
                           all_genes   = all_genes(BANDITS_data)[ all_genes(BANDITS_data) %in% unlist(genes(BANDITS_data)[SEL]) ] )
    
    message(paste0("Initial number of genes: ", length(all_genes(BANDITS_data)), "; number of selected genes: ", length(all_genes(BANDITS_data_new)) ) )
    
    return(BANDITS_data_new)
  }
  
  message(paste0("Initial number of genes: ", length(all_genes(BANDITS_data)), "; number of selected genes: ", length(all_genes(BANDITS_data)) ) )
  
  BANDITS_data
}
