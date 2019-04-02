#' Results of the DTU test, generated with \code{\link{test_DTU}}
#' 
#' @rdname results
#' @name results
#' @aliases results
#'
#' @param results a \code{\linkS4class{BANDITS_test}} object containing the results of the DTU test,
#' generated with \code{\link{test_DTU}}.
#' 
#' @examples
#' # Computed as shown in the vignettes, see: browseVignettes("BANDITS")
#' 
#' #set.seed(61217)
#' #results = test_DTU(BANDITS_data = BANDITS_data,
#' #             prior_precision = precision$prior,
#' #             samples_design = samples_design,
#' #             R = 10^4, burn_in = 2*10^3, n_cores = 2,
#' #             gene_to_transcript = gene_tr_id)
#' 
#' # load the pre-computed results:
#' data("results", package = "BANDITS")
#' results
#' 
#' # Visualize the most significant Genes, sorted by gene level significance.
#' head(top_genes(results))
#' 
#' # Alternatively, gene-level results can also be sorted according to DTU_measure, 
#' # which is a measure of the strength of the change between the 
#' # average relative abundances of the two groups.
#' head(top_genes(results, sort_by = "DTU_measure"))
#' 
#' # Visualize the most significant transcripts, sorted by transcript level significance.
#' head(top_transcripts(results, sort_by = "transcript"))
#' 
#' # Visualize the convergence output for the most significant genes, 
#' # sorted by gene level significance.
#' head(convergence(results))
#' 
#' # We can further use the \code{gene} function to gather all output for a specific gene:
#' # gene level, transcript level and convergence results.
#' top_gene = top_genes(results, n = 1)
#' gene(results, top_gene$Gene_id)
#' 
#' # Similarly we can use the \code{transcript} function to gather all output 
#' # for a specific transcript.
#' top_transcript = top_transcripts(results, n = 1)
#' transcript(results, top_transcript$Transcript_id)
#' 
#' #Finally, we can plot the estimated average transcript relative expression 
#' # in the two groups for a specific gene via \code{plot_proportions}.
#' library(ggplot2)
#' plot_proportions(results, top_gene$Gene_id)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{test_DTU}}, \code{\linkS4class{BANDITS_test}}
NULL