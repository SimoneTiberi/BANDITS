#' Gene-transcript matching
#' 
#' @rdname gene_tr_id
#' @name gene_tr_id
#' @aliases gene_tr_id
#' 
#' @param gene_tr_id a \code{data.frame} containing the matching between gene (1st column) and transcript identifiers (2nd column).
#' The gtf file used was downloaded from the ARMOR github repository 
#' \href{https://github.com/csoneson/ARMOR/blob/master/example_data/reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.93.1.1.10M.gtf}{here}.
#' 
#' @examples
#' # Compute the 'gene_tr_id' object from the gtf file as shown in the vignettes,
#' # see: browseVignettes("BANDITS").
#' # suppressMessages(library(GenomicFeatures))
#' # tx = makeTxDbFromGFF("Homo_sapiens.GRCh38.93.1.1.10M.gtf")
#' # ss = unlist(transcriptsBy(tx, by="gene"))
#' # gene_tr_id_gtf = data.frame(gene_id = names(ss), transcript_id = ss$tx_name )
#' # gene_tr_id_gtf = gene_tr_id_gtf[ rowSums( is.na(gene_tr_id_gtf)) == 0, ] # remove eventual NA's
#' # gene_tr_id_gtf = unique(gene_tr_id_gtf) # remove eventual duplicated rows
#' 
#' # load the Gene-Transcript data.frame and visualize its top
#' data(gene_tr_id, package = "BANDITS")
#' head(gene_tr_id)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
NULL