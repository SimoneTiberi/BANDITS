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
#' # load the pre-computed data:
#' data("input_data", package = "BANDITS")
#' show(input_data)
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
  message(paste0("A 'BANDITS_data' object with ", ncol(counts(object)[[1]]), " samples and ", length(all_genes(object)), " genes."))
})

###############################################################################
### Set validity of the object
###############################################################################
setValidity("BANDITS_data", function(object){
  # Has to return TRUE for a valid object!
  n_groups = length(genes(object))
  if( (n_groups != length(transcripts(object))) | (n_groups != length(effLen(object))) | (n_groups != length(classes(object))) | (n_groups != length(counts(object))) | (n_groups != length(uniqueId(object))) ){
    return("The length of @genes, @transcripts, @effLen, @classes, @counts and @uniqueId slots must be the same")
  }
  
  return(TRUE)
})


###############################################################################
### accessing methods, private
###############################################################################
setGeneric("genes", function(x) 
  standardGeneric("genes") )
setMethod("genes", "BANDITS_data", function(x) x@genes)

setGeneric("transcripts", function(x) 
  standardGeneric("transcripts") )
setMethod("transcripts", "BANDITS_data", function(x) x@transcripts)

setGeneric("effLen", function(x) 
  standardGeneric("effLen") )
setMethod("effLen", "BANDITS_data", function(x) x@effLen)

setGeneric("classes", function(x) 
  standardGeneric("classes") )
setMethod("classes", "BANDITS_data", function(x) x@classes)

setGeneric("counts", function(x) 
  standardGeneric("counts") )
setMethod("counts", "BANDITS_data", function(x) x@counts)

setGeneric("uniqueId", function(x) 
  standardGeneric("uniqueId") )
setMethod("uniqueId", "BANDITS_data", function(x) x@uniqueId)

setGeneric("all_genes", function(x) 
  standardGeneric("all_genes") )
setMethod("all_genes", "BANDITS_data", function(x) x@all_genes)
