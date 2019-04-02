#' Estimates for the log-precision parameter, generated with \code{\link{prior_precision}}
#' 
#' @rdname precision
#' @name precision
#' @aliases precision
#' 
#' @param precision a \code{list} with two elements:
#' \itemize{
#'   \item 'prior', a numeric \code{vector} containing the mean and standard deviation prior estimates
#'   for the log-precision parameter of the Dirichlet-multinomial distribution;
#'   \item 'genewise_log_precision', a numeric \code{vector} containing the gene-wise precision estimates
#'   of the log-precision parameter of the Dirichlet-multinomial distribution.
#' }
#
#' @examples
#' # Computed as shown in the vignettes, see: browseVignettes("BANDITS")
#'
#' #set.seed(61217)
#' #precision = prior_precision(gene_to_transcript = gene_tr_id, transcript_counts = counts,
#' #                       min_transcript_proportion = 0.01, min_transcript_counts = 10,
#' #                       min_gene_counts = 20, n_cores = 2)
#' 
#' # load the precision estimates
#' data(precision, package = "BANDITS")
#' 
#' # mean and standard deviation estimates of the log-precision parameter:
#' precision$prior
#' 
#' # gene-wise estimates of the log-precision parameter:
#' precision$genewise_log_precision
#' 
#' # Plot the histogram of the genewise log-precision estimates.
#' # The black solid line represents the normally distributed prior distribution 
#' # for the log-precision parameter.
#' plot_precision(precision)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{prior_precision}}
NULL