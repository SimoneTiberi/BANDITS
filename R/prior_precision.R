#' Infer an informative prior for the precision
#'
#' \code{prior_precision} uses \code{DRIMSeq}'s pipeline to infer an informative prior for the
#' precision parameter of the Dirichlet-Multinomial distribution.
#' The function computes the genewise estimates for the precision via \code{DRIMSeq::dmPrecision},
#' and calculates the mean and standard deviation of the log-precision estimates.
#' 
#' @param gene_to_transcript a matrix or data.frame with a list of gene-to-transcript correspondances.
#' The first column represents the gene id, while the second one contains the transcript id.
#' @param transcript_counts a matrix or data.frame, with 1 column per sample and 1 row per transcript, 
#' containing the estimated abundances for each transcript in each sample.
#' @param min_transcript_proportion the minimum relative abundance (i.e., proportion) of a transcript in a gene.
#' @param min_transcript_counts the minimum overall abundance of a transcript (adding counts from all samples).
#' @param min_gene_counts the minimum overall abundance of a gene (adding counts from all samples).
#' @param n_cores the number of cores to parallelize the tasks on.
#' 
#' @return A list with 2 objects containing:
#' \itemize{
#' \item prior: a vector containing the mean and standard deviation of the log-precision, used to formulate an informative prior in \code{\link{test_DTU}};
#' \item genewise_log_precision: a numeric vector with the individual genewise estimates for the log-precision.
#' }
#' 
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
#' ## Optional (recommended): infer an informative prior for the precision parameter
#' 
#' # Use the same filtering criteria as in \code{\link{filter_transcripts}}; 
#' # if transcript pre-filtering is not performed, set \code{min_transcript_proportion},
#' # \code{min_transcript_counts} and \code{min_gene_counts} to 0.
#' 
#' #set.seed(61217)
#' #precision = prior_precision(gene_to_transcript = gene_tr_id, transcript_counts = counts,
#' #                       min_transcript_proportion = 0.01, min_transcript_counts = 10,
#' #                       min_gene_counts = 20, n_cores = 2)
#' 
#' # load the pre-computes precision estimates:
#' data(precision, package = "BANDITS")
#' 
#' # Plot the histogram of the genewise log-precision estimates.
#' # The black solid line represents the normally distributed prior distribution 
#' # for the log-precision parameter.
#' plot_precision(precision)
#' 
#' @author Simone Tiberi
#'  
#' @seealso \code{\link{test_DTU}}, \code{\link{plot_precision}}
#' 
#' @export
prior_precision = function(gene_to_transcript, transcript_counts,
                           min_transcript_proportion = 0.01, min_transcript_counts = 1, 
                           min_gene_counts = 10, n_cores = 1){
  # check that gene_to_transcript is a matrix or data.frame object
  if( !is.data.frame(gene_to_transcript) & !is.matrix(gene_to_transcript)  ){
    message("'gene_to_transcript' must be a matrix or data.frame")
    return(NULL)
  }
  
  if( ncol(gene_to_transcript) != 2 ){
    message("'gene_to_transcript' must be a 2 column matrix or data.frame")
    return(NULL)
  }
  
  # check that transcript_counts is a matrix or data.frame object
  if( !is.data.frame(transcript_counts) & !is.matrix(transcript_counts) ){
    message("'transcript_counts' must be a matrix or data.frame")
    return(NULL)
  }
  
  if( !all( rownames(transcript_counts) %in% gene_to_transcript[,2])  ){
    message("All transcript names in 'rownames(transcript_counts)' must be in 'gene_to_transcript[,2]'")
    return(NULL)
  }
  
  if( !all(dim(transcript_counts) > 0.5) ){
    message("'transcript_counts' must have at least 1 row (transcripts) and 1 column (samples)")
    return(NULL)
  }
  
  N = ncol(transcript_counts)
  
  matches = match( rownames(transcript_counts), gene_to_transcript[,2] )
  gene_id = as.character( gene_to_transcript[matches,1] )
  
  colnames(transcript_counts) = paste("sample", seq_len(N))
  counts_df = data.frame(transcript_counts, gene_id = gene_id, feature_id = rownames(transcript_counts))
  
  # maybe design and samples not needed for precision estimate...check the function internally!
  samples = data.frame(sample_id = colnames(counts_df)[seq_len(N)],
                       group = "A")
  
  d = dmDSdata(counts = counts_df, samples = samples)
  
  d = dmFilter(d, min_gene_expr = min_gene_counts, min_feature_expr = min_transcript_counts,
               min_feature_prop = min_transcript_proportion)
  
  design = model.matrix(~ 1, data = samples(d))
  
  message("Estimating gene-wise precision parameters")
  suppressMessages({
    d = dmPrecision(d, genewise_precision = TRUE, 
                    design = design, BPPARAM = MulticoreParam(workers = n_cores))
  })
  message("Estimation completed")
  
  log_prec = log( genewise_precision(d)$genewise_precision )
  mean_log_precision = mean(log_prec, na.rm = TRUE)
  sd_log_precision   = sd(log_prec, na.rm = TRUE)
  
  names(log_prec) = genewise_precision(d)$gene_id
  list( prior = c(mean_log_precision, sd_log_precision), genewise_log_precision = log_prec )
}

#' Plot the log-precision estimates
#'
#' \code{\link{{plot_precision}} plots a histogram of the estimates for the log-precision parameter 
#' of the Dirichlet-Multinomial distribution, obtained via \code{\link{prior_precision}}.
#' The solid line represents the normal prior for the log-precision parameter.
#' 
#' @param prior the prior of the log-precision parameter, computed via \code{\link{prior_precision}}.
#' 
#' @return A plot.
#' 
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
#' ## Optional (recommended): infer an informative prior for the precision parameter
#' 
#' # Use the same filtering criteria as in \code{\link{filter_transcripts}}; 
#' # if transcript pre-filtering is not performed, set \code{min_transcript_proportion},
#' # \code{min_transcript_counts} and \code{min_gene_counts} to 0.
#' 
#' #set.seed(61217)
#' #precision = prior_precision(gene_to_transcript = gene_tr_id, transcript_counts = counts,
#' #                       min_transcript_proportion = 0.01, min_transcript_counts = 10,
#' #                       min_gene_counts = 20, n_cores = 2)
#' 
#' # load the pre-computed precision estimates:
#' data(precision, package = "BANDITS")
#' 
#' # Plot the histogram of the genewise log-precision estimates.
#' # The black solid line represents the normally distributed prior distribution 
#' # for the log-precision parameter.
#' plot_precision(precision)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{test_DTU}}, \code{\link{prior_precision}}
#' 
#' @export
plot_precision = function(prior){
  if( !is.list(prior) ){
    message("'prior' must be a list of length 2, created via 'prior_precision'")
    return(NULL)
  }

  if( length(prior) != 2 ){
    message("'prior' must be a list of length 2, created via 'prior_precision'")
    return(NULL)
  }
  
  hist( prior[[2]], main = "Log-prior precision estimates", 
        xlab = "log-prior", freq = FALSE)
  curve(dnorm(x, prior[[1]][1], prior[[1]][2]), add = TRUE, lwd = 3)
}
