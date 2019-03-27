#' Compute the median effective length of transcripts.
#'
#' \code{eff_len_compute} inputs the estimated effective length of transcripts from every sample, and
#' computes the median effective length of each transcript across samples.
#'
#' @param x_eff_len is a list: each element of the list refers to a specific sample and is a matrix or data.frame
#'  with the estimated effective length under the column 'EffectiveLength' and the transcript name under the column 'Name'.
#' 
#' @return A vector containing the effective length of transcripts; the vector names indicate the transcript ids.
#' @examples
#' # specify the directory of the internal data:
#' data_dir = system.file("extdata", package = "BANDITS")
#' data_dir
#' 
#' # Specify the directory of the transcript level estimated counts.
#' sample_names = paste0("sample", seq_len(4))
#' quant_files = file.path(data_dir, sample_names, "quant.sf")
#' file.exists(quant_files)
#' 
#' # Load the transcript level estimated counts via tximport:
#' library(tximport)
#' txi = tximport(files = quant_files, type = "salmon", txOut = TRUE)
#' 
#' # compute the Median estimated effective length for each transcript:
#' eff_len = eff_len_compute(x_eff_len = txi$length)
#' 
#' @author Simone Tiberi
#'  
#' @seealso \code{\link{filter_transcripts}}, \code{\link{create_data}}
#' @export
eff_len_compute = function(x_eff_len){
#  tr_id_eff_len = unique( c( sapply( x_eff_len, function(x) x$Name ) ) )
#  y_eff_len = sapply(x_eff_len, function(y){ y$EffectiveLength[match(tr_id_eff_len, y$Name)] } )

  EffLen = apply(x_eff_len, 1, median, na.rm = TRUE)
#  names(EffLen) = tr_id_eff_len
  return( EffLen )
}
