#' BANDITS_data class
#'
#' \code{BANDITS_data} contains all the information required to perform differential transcript usage (DTU).
#' \code{BANDITS_data} associates each gene (genes), to its transcript ids (transcripts),
#' effective transcript lengths (effLen), equivalence classes (classes) and respective counts (counts).
#' The same structure is also used for groups of genes with reads/fragments compatible with >1 gene (with uniqueId == FALSE); 
#' in this case the 'genes' field contains all the genes ids in the group.
#' Created via \code{\link{create_data}}.
#' 
#' @return
#' \itemize{
#' \item \code{show(object)}: returns the number of genes and transcripts in the \code{BANDITS_data} object.
#' }
#' 
#' @slot genes \code{list} of gene names: each element is a vector of 1 or 
#' more gene names indicating the genes to be analyzed together.
#' @slot transcripts \code{list} of transcript names: each element is a vector of 1 or
#'  more transcript names indicating the transcripts matching the gene names in the corresponding element of @genes object.
#' @slot effLen \code{list} of transcript effective lengths: each element is a vector of 1 or
#'  more numbers, indicating the effective length of the transcripts in the corresponding element of @transcripts object.
#' @slot classes \code{list} of matrices: the (i,j) element of each matrix is 1 if the i-th transcript
#' is present in the j-th equivalence class, 0 otherwise.
#' @slot counts \code{list} of matrices: the (i,j) element indicates the reads/fragments compatible with 
#' the i-th equivalence class in sample j.
#' @slot uniqueId \code{logical}, it indicates if the element contains one gene to be analyzed alone (TRUE),
#'  or more genes to be analyzed jointly (FALSE).
#' @slot all_genes \code{vector}, it lists all the genes to be analyzed (with at least 2 transcripts).
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
#' 
#' # create data and filter internally lowly abundant transcripts:
#' BANDITS_data = create_data(gene_to_transcript = gene_tr_id,
#'                            path_to_eq_classes = equiv_classes_files, eff_len = eff_len, 
#'                            n_cores = 2,
#'                            transcripts_to_keep = transcripts_to_keep)
#' 
#' # If transcripts pre-filtering is not wanted, do not specify \code{transcripts\_to\_keep} parameter.
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#'
#' @seealso \code{\link{create_data}}, \code{\link{filter_transcripts}},  \code{\link{eff_len_compute}}
#' 
#' @aliases BANDITS_data
#' 
#' @export
setClass("BANDITS_data", 
         slots = representation(genes = "list", transcripts = "list", effLen = "list", classes = "list", counts = "list", uniqueId = "logical",
                                all_genes = "vector" ))

#' @rdname BANDITS_data-class
#' @param object a 'BANDITS_data' object.
#' @export
setMethod("show", "BANDITS_data", function(object){
  message(paste0("A 'BANDITS_data' object of length ", length(object@uniqueId), "."))
  message(paste0("Number of samples: ", ncol(object@counts[[1]]), "."))
  message(paste0("Number of genes: ", length(object@all_genes), "."))
})
