#' Estimates for the log-precision parameter, generated with \code{\link{prior_precision}}
#' 
#' @rdname precision
#' @name precision
#' @aliases precision
#' 
#' @param precision a \code{list}, generated via \code{\link{prior_precision}}, with two elements:
#' \itemize{
#'   \item 'prior', a numeric \code{vector} containing the mean and standard deviation prior estimates
#'   for the log-precision parameter of the Dirichlet-multinomial distribution;
#'   \item 'genewise_log_precision', a numeric \code{vector} containing the gene-wise precision estimates
#'   of the log-precision parameter of the Dirichlet-multinomial distribution.
#' }
#
#' @examples
#' # Object 'precision' is generated via 'prior_precision' as shown in the vignettes: 
#' # see browseVignettes("BANDITS").
#' 
#' # load the pre-computed precision estimates:
#' data(precision, package = "BANDITS")
#' 
#' precision$prior
#' head(precision$genewise_log_precision)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{prior_precision}}, \code{\link{plot_precision}}
NULL