#' BANDITS_test class
#' 
#' \code{BANDITS_test} contains the results of the differential transcript usage (DTU) test.
#' \code{BANDITS_test} is organized in three data.frames containing: gene-level results, 
#' transcript-level results and 
#' convergence diagnostics of the Markov chain Monte Carlo (MCMC) posterior chains.
#' Created via \code{\link{test_DTU}}.
#' To test for convergence, we use Heidelberger and Welch's convergence diagnostic,
#' implemented in \code{coda::heidel.diag}, to test for the stationarity of the 
#' chain for the full log-posterior density;
#' we use a 0.01 threshold on the p.value to reject the null hypotehsis of stationarity.
#' 
#' @return
#' \itemize{
#' \item \code{show(object)}: returns the \code{head} of the gene and transcript level results and of the convergence diagnostic.
#' \item \code{top_genes(x, n = Inf,  sort_by = "p.value")}: returns the gene-level results of the DTU test for the top 'n' significant genes.
#' By default n = Inf and all results will be reported.
#' sort_by = "gene" for sorting results according to gene-level significance; sort_by = "DTU_measure" for sorting results according to the 'DTU_measure'.
#' \item \code{top_transcripts(x, n = Inf, sort_by = "gene")}: returns the transcript-level results of the DTU test for the top 'n' significant genes.
#' By default n = Inf and all results will be reported.
#' sort_by = "gene" for sorting results according to gene-level significance; sort_by = "transcript" for sorting results according to transcript-level significance.
#' \item \code{convergence(x)}: returns the convergence diagnostic of the posterior MCMC chains for every gene.
#' \item \code{gene(x, gene_id)}: returns a list with all results for the gene(s) specified in 'gene_id': gene results, corresponding transcript results and convergence diagnostic.
#' \item \code{transcript(x, transcript_id)}: returns a list with all results for the trancript specified in 'transcript_id': transcript results, corresponding gene results and convergence diagnostic.
#' \item \code{plot_proportions(x, gene_id, group_names = NULL)}: plots the posterior means of the average transcripts
#'  relative expression (i.e., the proportions) of each condition, for the gene specified in 'gene_id'.
#'  group_names is a carachter vector with the same length as the number of groups; if not provided, the first letters of the alphabet are used.
#' }
#' 
#' @slot Gene_results a \code{data.frame} containing the gene-level results of the DTU test, structured in the following columns:
#' \itemize{
#' \item Gene_id contains the gene names;
#' \item p.values is the gene-level p.values of the DTU test;
#' \item adj.p.values is the Benjamini-Hochberg adjusted p.values (via \code{\link{p.adjust}});
#' \item p.values_inverted (only available for 2-group comparisons) is a conservative p.value, accounting for the inversion of the dominant transcript between conditions;
#' \item adj.p.values_inverted (only available for 2-group comparisons) is the Benjamini-Hochberg adjusted p.values_inverted, via \code{\link{p.adjust}};
#' \item DTU_measure (only available for 2-group comparisons) represents a measure of the intensity of changes between conditions.
#' }
#' 
#' @slot Transcript_results a \code{data.frame} containing the transcript-level results of the DTU test, structured in the following columns:
#' \itemize{
#' \item Gene_id contains the gene names;
#' \item Transcript_id contains the transcript names;
#' \item p.values is the transcript-level p.values of the DTU test;
#' \item adj.p.values is the Benjamini-Hochberg adjusted p.values (via \code{\link{p.adjust}});
#' \item Max_Gene_Tr.p.val is a conservative p.value (the maximum between the transcript p.value and corresponding gene p.value);
#' \item Max_Gene_Tr.Adj.p.val is the Benjamini-Hochberg adjusted Max_Gene_Tr.p.val (via \code{\link{p.adjust}});
#' \item Mean "group_name" indicates the posterior mean of the average relative abundance of the transcript in group "group_name".
#' }
#' 
#' @slot Convergence a \code{data.frame} containing the convercence diagnostics of the DTU test, structured in the following columns:
#' \itemize{
#' \item Gene_id contains the gene names;
#' \item converged is 1 if convergence was reached, 0 otherwise;
#' \item burn_in indicates what fraction of the chain was removed to ensure convergence 
#' (excluding the \code{burn_in} parameter specified in \code{\link{test_DTU}}.
#' }
#' 
#' @slot samples_design a \code{data.frame} containing the design of the experiment, with one row for each sample
#' and two columns with names 'sample_id' and 'group', specifying the id and group of each sample, respectively.
#' It is provided by the user to \code{\link{test_DTU}}.
#' 
#' @examples 
#' ## Preliminary information
#' 
#' # specify the directory of the internal data:
#' data_dir = system.file("extdata", package = "BANDITS")
#' data_dir
#' 
#' # load gene_to_transcript matching:
#' data("GeneTr_id", package = "BANDITS")
#' # GeneTr_id contains transcripts ids on the first column
#' # and the corresponding gene ids on the second column:
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
#' # We define the design of the study: in our case we have 2 groups, 
#' # that we call "A" and "B" of 2 samples each.
#' samples_design = data.frame(sample_id = sample_names,
#'                             group = c("A", "A", "B", "B"))
#' samples_design
#' 
#' # The groups are defined in:
#' levels(samples_design$group)
#' 
#' 
#' 
#' ## Optional (recommended): transcript pre-filtering
#' 
#' transcripts_to_keep = filter_transcripts(gene_to_transcript = GeneTr_id,
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
#' # Warning: the sample names in equiv_classes_files must have the same order
#' # as those in the design object, containted in samples_design.
#' equiv_classes_files
#' samples_design$sample_id
#' 
#' # create data and filter internally lowly abundant transcripts:
#' BANDITS_data = create_data(gene_to_transcript = GeneTr_id,
#'                            path_to_eq_classes = equiv_classes_files, eff_len = eff_len, 
#'                            n_cores = 2,
#'                            transcripts_to_keep = transcripts_to_keep)
#' 
#' # If transcripts pre-filtering is not wanted, 
#' # do not specify \code{transcripts_to_keep} parameter.
#' 
#' # Filter lowly abundant genes:
#' BANDITS_data = filter_genes(BANDITS_data, min_counts_per_gene = 20)
#' 
#' 
#' 
#' ## Optional (recommended): infer an informative prior for the precision parameter
#' 
#' # Use the same filtering criteria as in \code{\link{filter_transcripts}}; 
#' # if transcript pre-filtering is not performed, set \code{min_transcript_proportion},
#' # \code{min_transcript_counts} and \code{min_gene_counts} to 0.
#' 
#' set.seed(61217)
#' prec = prior_precision(gene_to_transcript = GeneTr_id, transcript_counts = counts,
#'                        min_transcript_proportion = 0.01, min_transcript_counts = 10,
#'                        min_gene_counts = 20, n_cores = 2)
#' 
#' # Plot the histogram of the genewise log-precision estimates.
#' # The black solid line represents the normally distributed prior distribution 
#' # for the log-precision parameter.
#' plot_precision(prec)
#' 
#' 
#' 
#' ## Test for DTU
#' set.seed(61217)
#' x = test_DTU(BANDITS_data = BANDITS_data,
#'              prior_precision = prec$prior,
#'              samples_design = samples_design,
#'              R = 10^4, burn_in = 2*10^3, n_cores = 2,
#'              gene_to_transcript = GeneTr_id)
#' x
#' 
#' # Visualize the most significant Genes, sorted by gene level significance.
#' head(top_genes(x))
#' 
#' # Alternatively, gene-level results can also be sorted according to DTU_measure, 
#' # which is a measure of the strength of the change between the 
#' # average relative abundances of the two groups.
#' head(top_genes(x, sort_by = "DTU_measure"))
#' 
#' # Visualize the most significant transcripts, sorted by transcript level significance.
#' head(top_transcripts(x, sort_by = "transcript"))
#' 
#' # Visualize the convergence output for the most significant genes, 
#' # sorted by gene level significance.
#' head(convergence(x))
#' 
#' # We can further use the \code{gene} function to gather all output for a specific gene:
#' # gene level, transcript level and convergence results.
#' top_gene = top_genes(x, n = 1)
#' gene(x, top_gene$Gene_id)
#' 
#' # Similarly we can use the \code{transcript} function to gather all output 
#' # for a specific transcript.
#' top_transcript = top_transcripts(x, n = 1)
#' transcript(x, top_transcript$Transcript_id)
#' 
#' #Finally, we can plot the estimated average transcript relative expression 
#' # in the two groups for a specific gene via \code{plot_proportions}.
#' library(ggplot2)
#' plot_proportions(x, top_gene$Gene_id)
#' 
#' @author Simone Tiberi
#'
#' @seealso \code{\link{test_DTU}}, \code{\link{create_data}}, \code{\linkS4class{BANDITS_data}}
#' 
#' @aliases BANDITS_test length top_genes top_transcripts convergence gene transcript plot_proportions
#' 
#' @export
setClass("BANDITS_test", 
         slots = representation(Gene_results = "data.frame", Transcript_results = "data.frame", 
                        Convergence = "data.frame", samples_design = "data.frame"))

#' @rdname BANDITS_test-class
#' @export
setMethod("length", "BANDITS_test", function(x){
  nrow(x@Gene_results)
})

setGeneric("top_genes", function(x, ...) 
  standardGeneric("top_genes") )

#' @rdname BANDITS_test-class
#' @export
setMethod("top_genes", "BANDITS_test", function(x, n = Inf, sort_by = "p.value"){
  if(!is.character(sort_by)){
    message("'sort_by' must be a character string, indicating either 'p.value' or 'DTU_measure'")
    return(NULL)
  }
  if(! sort_by %in% c("p.value", "DTU_measure")){
    message("'sort_by' must be either 'p.value' or 'DTU_measure'")
    return(NULL)
  }
  
  # take the min between the provided number and the size of the matrix:
  n = min(n, nrow(x@Gene_results) )
  
  if(sort_by == "p.value"){ # transcript results sorted according to gene-level p.value
    return(x@Gene_results[seq_len(n),])
  }
  
  if(is.null(x@Gene_results$DTU_measure)){
    message("Results cannot be sorted by 'DTU_measure': not available for multi-group comparisons.")
    return(NULL)
  }
  # transcript results sorted according to transcript-level p.value
  x@Gene_results[order(x@Gene_results$DTU_measure, decreasing = TRUE)[seq_len(n)],] # high DTU_measure on top.
})

setGeneric("top_transcripts", function(x, ...) 
  standardGeneric("top_transcripts") )

#' @rdname BANDITS_test-class
#' @export
setMethod("top_transcripts", "BANDITS_test", function(x, n = Inf, sort_by="gene"){
  if(!is.character(sort_by)){
    message("'sort_by' must be a character string")
    return(NULL)
  }
  if(! sort_by %in% c("gene", "transcript")){
    message("'sort_by' must be either 'gene' or 'transcript'")
    return(NULL)
  }
  
  # take the min between the provided number and the size of the matrix:
  n = min(n, nrow(x@Gene_results) )
  
  if(sort_by == "gene"){ # transcript results sorted according to gene-level p.value
    return(x@Transcript_results[seq_len(n),])
  }
  # transcript results sorted according to transcript-level p.value
  x@Transcript_results[order(x@Transcript_results$p.values)[seq_len(n)],]
})

setGeneric("convergence", function(x) 
  standardGeneric("convergence") )

#' @rdname BANDITS_test-class
#' @export
setMethod("convergence", "BANDITS_test", function(x){
  x@Convergence
})

#' @rdname BANDITS_test-class
#' @export
setMethod("show", "BANDITS_test", function(object){
  message(paste0("A 'BANDITS_test' object, with ", nrow(object@Gene_results), " genes and ", nrow(object@Transcript_results), " transcript level results."))
})

setGeneric("gene", function(x, ...)
  standardGeneric("gene") )

#' @rdname BANDITS_test-class
#' @export
setMethod("gene", "BANDITS_test", function(x, gene_id){
  if(!is.character(gene_id)){
    gene_id = as.character(gene_id)
  }
  
  if( ! all( gene_id %in% as.character( x@Gene_results$Gene_id ) ) ){
    message("One or more genes in 'gene_id' were not found in the results")
    return(NULL)
  }
  
  list( gene_results = x@Gene_results[ as.character( x@Gene_results$Gene_id ) %in% gene_id, ],
        transcript_results = x@Transcript_results[  as.character( x@Transcript_results$Gene_id ) %in% gene_id, ],
        convergence_results = x@Convergence[ as.character( x@Convergence$Gene_id ) %in% gene_id, ] )
})

setGeneric("transcript", function(x, ...) 
  standardGeneric("transcript") )

#' @rdname BANDITS_test-class
#' @export
setMethod("transcript", "BANDITS_test", function(x, transcript_id){
  if(!is.character(transcript_id)){
    transcript_id = as.character(transcript_id)
  }
  
  if( ! transcript_id %in%  as.character(x@Transcript_results$Transcript_id ) ){
    message("transcript not found in results")
    return(NULL)
  }
  
  gene_id = x@Transcript_results$Gene_id[ as.character(x@Transcript_results$Transcript_id) == transcript_id]
  list( transcript_results = x@Transcript_results[ as.character(x@Transcript_results$Transcript_id) == transcript_id, ],
        gene_results = x@Gene_results[ as.character(x@Gene_results$Gene_id) == gene_id, ],
        convergence_results = x@Convergence[ as.character(x@Convergence$Gene_id) == gene_id, ] )
})

setGeneric("plot_proportions", function(x, ...) 
  standardGeneric("plot_proportions") )

#' @rdname BANDITS_test-class
#' @export
setMethod("plot_proportions", "BANDITS_test", function(x, gene_id){
  if(!is.character(gene_id)){
    gene_id = as.character(gene_id)
  }
  
  if( ! gene_id %in% as.character( x@Transcript_results$Gene_id ) ){
    message("gene not found in results")
    return(NULL)
  }
  
  group_names = levels(x@samples_design$group)
  
  sel = x@Transcript_results$Gene_id == gene_id

  means = x@Transcript_results[sel, grep("Mean", names(x@Transcript_results))]
  n_groups = ncol(means)
  tr_names = x@Transcript_results$Transcript_id[sel]
  # impose an order to the transcripts (according to the over-all relative abudance):
  ord = order(rowSums(means), decreasing = TRUE)

  prop_samp = data.frame(feature_id = factor( rep(tr_names,n_groups), levels = tr_names[ord]), 
                         proportion = unlist(c(means)),
                         group = rep(group_names, each = nrow(means)),
                         stringsAsFactors = FALSE)

  # Plot the estimated average proportions of each groups:
  ggp <- ggplot() +
    geom_bar(data = prop_samp, aes_string(x = "feature_id", y = "proportion", 
                                          fill = "group"),
             stat = "identity", position = position_dodge(width = 0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
          axis.text=element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16), 
          legend.position = "right", 
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14)) +
    ggtitle(paste("Gene:", gene_id, sep=" ")) +
    xlab("Features") +
    ylab("Proportions")
  
  if(!is.null(prop_samp)){
    ggp <- ggp + 
      geom_point(data = prop_samp, 
                 aes_string(x = "feature_id", y = "proportion", 
                            group = "group", fill = "group"), 
                 position = position_dodge(width = 0.9), size = 3, shape = 23, 
                 alpha = 0.75)
  }
  
  ggp
})
