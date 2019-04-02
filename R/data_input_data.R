#' A \code{\linkS4class{BANDITS_data}} object, generated with \code{\link{create_data}}
#' 
#' @rdname input_data
#' @name input_data
#' @aliases input_data
#' 
#' @param input_data a \code{\linkS4class{BANDITS_data}} object.
#
#' @examples
#' # Computed as shown in the vignettes, see: browseVignettes("BANDITS")
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
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{create_data}}, \code{\linkS4class{BANDITS_data}}
NULL