#' Filter lowly abundant transcripts.
#'
#' \code{filter_transcripts} filters transcripts, before loading the data, according to estimated transcript level counts.
#' The function outputs a vector containing the list of transcripts which respect the filtering criteria across all samples 
#' (i.e., min_transcript_proportion, min_transcript_counts and min_gene_counts).
#' 
#' Transcript pre-filtering is highly suggested: it both improves the performance of the method 
#' and decreases its computational cost.
#'
#' @param gene_to_transcript a matrix or data.frame with a list of gene-to-transcript correspondances.
#' The first column represents the gene id, while the second one contains the transcript id.
#' @param transcript_counts a matrix or data.frame, with 1 column per sample and 1 row per transcript, 
#' containing the estimated abundances for each transcript in each sample.
#' @param min_transcript_proportion the minimum relative abundance (i.e., proportion) of a transcript in a gene.
#' @param min_transcript_counts the minimum overall abundance of a transcript (adding counts from all samples).
#' @param min_gene_counts the minimum overall abundance of a gene (adding counts from all samples).
#' 
#' @return A vector containing the list of transcripts which respect the filtering criteria.
#' @examples
#' ## Preliminary information
#' 
#' # specify the directory of the internal data:
#' data_dir = system.file("extdata", package = "BANDITS")
#' data_dir
#' 
#' # load gene_to_transcript matching:
#' data("GeneTr_id", package = "BANDITS")
#' # GeneTr_id contains transcripts ids on the first column and the corresponding gene ids on the second column:
#' head(GeneTr_id)
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
#' transcripts_to_keep = filter_transcripts(gene_to_transcript = GeneTr_id,
#'                                          transcript_counts = counts, min_transcript_proportion = 0.01,
#'                                          min_transcript_counts = 10, min_gene_counts = 20)
#' head(transcripts_to_keep)
#' 
#' @author Simone Tiberi
#'
#' @seealso \code{\link{filter_genes}}, \code{\link{create_data}}, \code{\linkS4class{BANDITS_data}}
#'
#' @export
filter_transcripts = function(gene_to_transcript, transcript_counts,
                              min_transcript_proportion = 0.01, min_transcript_counts = 1, 
                              min_gene_counts = 10){
  n_tr_initial = nrow(transcript_counts)
  Gene_id = gene_to_transcript[,1]
  Tr_id   = gene_to_transcript[,2]
  
  # Compute the relative expression of transcripts (the pi's), from Salmon estimated transcript level transcript_counts:
  # I associate each transcript in "transcript_counts" to the respective gene:
  transcript_counts_gene_id = Gene_id[ match(rownames(transcript_counts), Tr_id) ]
  
  transcript_split_by_gene = split(data.frame(transcript_counts), transcript_counts_gene_id)
  transcript_split_by_gene = lapply(transcript_split_by_gene, data.frame)
  
  # (ev.) PARALLELIZE:
  transcript_proportions = lapply(transcript_split_by_gene, compute_proportions)
  
  # (ev.) PARALLELIZE:
  # Min 
  sel_transcripts_proportions = unlist( lapply(transcript_proportions, filter_tr_proportions_OneN, min_transcript_proportion = min_transcript_proportion) )
  sel_transcripts_totalCounts = rownames(transcript_counts)[ rowSums( transcript_counts) >= min_transcript_counts ]

  All_sel_tr  = unique(sel_transcripts_proportions, sel_transcripts_totalCounts)
  sel_in_both = {All_sel_tr %in% sel_transcripts_proportions} & {All_sel_tr %in% sel_transcripts_totalCounts}
  
  transcripts_to_keep = All_sel_tr[sel_in_both]
  
  # gene-filtering (if over-all estimated gene count is < N):
  if(min_gene_counts > 0){
    # remove the transcripts which were filtered first!
    transcript_counts = transcript_counts[rownames(transcript_counts) %in% transcripts_to_keep,]
    
    transcript_counts_gene_id = Gene_id[ match(rownames(transcript_counts), Tr_id) ]
    transcript_split_by_gene = split(data.frame(transcript_counts), as.character(transcript_counts_gene_id))
    transcript_split_by_gene = lapply(transcript_split_by_gene, data.frame)
    
    gene_sums = sapply(transcript_split_by_gene, sum)
    genes_sel = names(gene_sums)[ gene_sums > min_gene_counts ]
    
    transcripts_to_keep = transcripts_to_keep[ transcripts_to_keep %in% Tr_id[Gene_id %in% genes_sel] ]
  }
  
  message(paste0("After filtering, ", round(100 * length(transcripts_to_keep)/n_tr_initial, 2), "% of transcripts are kept"))
  
  return(transcripts_to_keep)
}

# function to compute the relative abundance of transcripts
compute_proportions = function(x){
  if(nrow(x) > 1){ # if > 1 transcripts
    res = apply(x, 2, function(y) y/sum(y))
  }else{
    res = x/x
  }
  res[is.na(res)] = 0
  res
}

filter_tr_proportions_OneN = function(x, min_transcript_proportion){
  avg_exp = rowMeans(x)             # avg relative abundance ovarall
  
  filter_out = avg_exp < min_transcript_proportion # filter transcripts if avg abundance < 0.01 overall
  rownames(x)[filter_out == FALSE] # return tr ids, except the ones filtered out.
}
