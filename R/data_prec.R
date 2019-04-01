#' Estimates for the log-precision parameter
#' 
#' @rdname prec
#' @name prec
#' @aliases prec
#' 
#' @description A list containing the estimates for the log-precision parameter of the Dirichlet-multinomial distribution.
#' 
#' @param prec a \code{list} with two elements:
#' \itemize{
#'   \item 'prior', a numeric \code{vector} containing the mean and standard deviation prior estimates
#'   for the log-precision parameter;
#'   \item 'genewise_log_precision', a numeric \code{vector} containing the gene-wise precision estimates
#'   of the log-precision parameter.
#' }
#
#' @examples
#' # load the precision estimates
#' data(prec, package = "BANDITS")
#' 
#' # mean and standard deviation estimates of the log-precision parameter:
#' prec$prior
#' 
#' # gene-wise estimates of the log-precision parameter:
#' prec$genewise_log_precision
#' 
#' # Plot the histogram of the genewise log-precision estimates.
#' # The black solid line represents the normally distributed prior distribution 
#' # for the log-precision parameter.
#' plot_precision(prec)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{prior_precision}}
NULL