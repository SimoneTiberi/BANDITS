#' Estimates for the log-precision parameter
#' 
#' @rdname DTU_results
#' @name DTU_results
#' @aliases DTU_results
#' 
#' @description A  \code{\linkS4class{BANDITS_test}} object containing the results of the DTU test,
#' obtained by running \code{\link{test_DTU}}.
#' 
#' @examples
#' # load the pre-computed results:
#' data("DTU_results", package = "BANDITS")
#' x
#' 
#' # Visualize the most significant Genes, sorted by gene level significance.
#' head(top_genes(x))
#' 
#' # Alternatively, gene-level results can also be sorted according to DTU_measure, 
#' # which is a measure of the strength of the change between the 
#' # average relative abundances of the two groups.
#' head(top_genes(x, sort_by = "DTU_measure"))
#' 
#' # Visualize the most significant transcripts, sorted by transcript level significance.
#' head(top_transcripts(x, sort_by = "transcript"))
#' 
#' # Visualize the convergence output for the most significant genes, 
#' # sorted by gene level significance.
#' head(convergence(x))
#' 
#' # We can further use the \code{gene} function to gather all output for a specific gene:
#' # gene level, transcript level and convergence results.
#' top_gene = top_genes(x, n = 1)
#' gene(x, top_gene$Gene_id)
#' 
#' # Similarly we can use the \code{transcript} function to gather all output 
#' # for a specific transcript.
#' top_transcript = top_transcripts(x, n = 1)
#' transcript(x, top_transcript$Transcript_id)
#' 
#' #Finally, we can plot the estimated average transcript relative expression 
#' # in the two groups for a specific gene via \code{plot_proportions}.
#' library(ggplot2)
#' plot_proportions(x, top_gene$Gene_id)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{test_DTU}}, \code{\linkS4class{BANDITS_test}}
NULL