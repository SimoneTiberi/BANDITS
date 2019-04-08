#' Filter lowly abundant genes.
#'
#' \code{filter_genes} filters genes, according to the overall number of counts (across all samples) compatible with the gene.
#' The filtering also applies to groups of genes with reads/fragments compatible with >1 gene;
#' in this case, the number of counts considered is across all genes in the group.
#' 
#' The function inputs a 'BANDITS_data' object, and returns agaain a 'BANDITS_data' object after filtering genes and groups of genes.
#'
#' @param BANDITS_data a 'BANDITS_data' object, created with the \code{\link{create_data}} function.
#' @param min_counts_per_gene the minimum number of counts compatible with a gene (across all samples).
#' 
#' @return A \code{\linkS4class{BANDITS_data}} object.
#' @examples
#' ## Preliminary information
#' 
#' # specify the directory of the internal data:
#' data_dir = system.file("extdata", package = "BANDITS")
#' data_dir
#' 
#' # load gene_to_transcript matching:
#' data("gene_tr_id", package = "BANDITS")
#' # gene_tr_id contains transcripts ids on the first column
#' # and the corresponding gene ids on the second column:
#' head(gene_tr_id)
#' 
#' # Specify the directory of the transcript level estimated counts.
#' sample_names = paste0("sample", seq_len(4))
#' quant_files = file.path(data_dir, sample_names, "quant.sf")
#' file.exists(quant_files)
#' 
#' # Load the transcript level estimated counts via tximport:
#' library(tximport)
#' txi = tximport(files = quant_files, type = "salmon", txOut = TRUE)
#' counts = txi$counts
#' head(counts)
#' 
#' 
#' 
#' ## Optional (recommended): transcript pre-filtering
#' 
#' transcripts_to_keep = filter_transcripts(gene_to_transcript = gene_tr_id,
#'                                          transcript_counts = counts,
#'                                          min_transcript_proportion = 0.01,
#'                                          min_transcript_counts = 10,
#'                                          min_gene_counts = 20)
#' head(transcripts_to_keep)
#' 
#' 
#' 
#' ## Load the data:
#' 
#' # compute the Median estimated effective length for each transcript:
#' eff_len = eff_len_compute(x_eff_len = txi$length)
#' 
#' # specify the path to the equivalence classes:
#' equiv_classes_files = file.path(data_dir, sample_names, "aux_info", "eq_classes.txt")
#' file.exists(equiv_classes_files)
#' 
#' # create data and filter internally lowly abundant transcripts:
#' #input_data = create_data(gene_to_transcript = gene_tr_id,
#' #                           path_to_eq_classes = equiv_classes_files, eff_len = eff_len, 
#' #                           n_cores = 2,
#' #                           transcripts_to_keep = transcripts_to_keep)
#' 
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
  tot_counts = vapply( BANDITS_data@counts, sum, FUN.VALUE = numeric(1) )
  # sapply( BANDITS_data@counts, sum)
  SEL = tot_counts >= min_counts_per_gene
  
  n_initial = length(BANDITS_data@all_genes) # only genes with > 1 transcript can be analyzed for DTU
  
  # filter genes/groups from BANDITS_data:
  if(mean(SEL) < 1){ # if mean(SEL) == 1, no genes/groups were filtered.
    BANDITS_data@genes       = BANDITS_data@genes[SEL]
    BANDITS_data@transcripts = BANDITS_data@transcripts[SEL]
    BANDITS_data@effLen      = BANDITS_data@effLen[SEL]
    BANDITS_data@classes     = BANDITS_data@classes[SEL]
    BANDITS_data@counts      = BANDITS_data@counts[SEL]
    BANDITS_data@uniqueId    = BANDITS_data@uniqueId[SEL]
  }
  BANDITS_data@all_genes = BANDITS_data@all_genes[ BANDITS_data@all_genes %in% unlist(BANDITS_data@genes) ]
  
  message(paste0("Initial number of genes: ", n_initial, "; number of selected genes: ", length(BANDITS_data@all_genes) ) )
  
  return(BANDITS_data)
}
