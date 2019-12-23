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
#' \item \code{show(object)}: prints the number of gene and transcript level results in the \code{BANDITS_test} object.
#' \item \code{top_genes(x, n = Inf,  sort_by_g = "p.value")}: returns the gene-level results of the DTU test for the top 'n' significant genes.
#' By default n = Inf and all results will be reported.
#' sort_by_g = "gene" for sorting results according to gene-level significance; sort_by_g = "DTU_measure" for sorting results according to the 'DTU_measure'.
#' \item \code{top_transcripts(x, n = Inf, sort_by_tr = "gene")}: returns the transcript-level results of the DTU test for the top 'n' significant genes.
#' By default n = Inf and all results will be reported.
#' sort_by_tr = "gene" for sorting results according to gene-level significance; sort_by_tr = "transcript" for sorting results according to transcript-level significance.
#' \item \code{convergence(x)}: returns the convergence diagnostic of the posterior MCMC chains for every gene.
#' \item \code{gene(x, gene_id)}: returns a list with all results for the gene(s) specified in 'gene_id': gene results, corresponding transcript results and convergence diagnostic.
#' \item \code{transcript(x, transcript_id)}: returns a list with all results for the trancript specified in 'transcript_id': transcript results, corresponding gene results and convergence diagnostic.
#' \item \code{plot_proportions(x, gene_id, CI = TRUE, CI_level = 0.95)}: plots the posterior means of the average transcripts
#'  relative expression (i.e., the proportions) of each condition, for the gene specified in 'gene_id'.
#'  If 'CI' is TRUE, a profile Wald type confidence interval will also be plotted for each transcript estimated proportion;
#'  the level of the confidence interval is specified by 'CI_level'.
#' }
#' 
#' @slot Gene_results a \code{data.frame} containing the gene-level results of the DTU test, structured in the following columns:
#' \itemize{
#' \item Gene_id contains the gene names;
#' \item p.values is the gene-level p.values of the DTU test;
#' \item adj.p.values is the Benjamini-Hochberg adjusted p.values (via \code{\link{p.adjust}});
#' \item p.values_inverted (only available for 2-group comparisons) is a conservative p.value, 
#' accounting for the inversion of the dominant transcript between conditions: 
#' p.values_inverted = p.values, if the dominant transcript varies between conditions,
#' and p.values_inverted = sqrt( p.values ) if both conditions have the same dominant transcript;
#' \item adj.p.values_inverted (only available for 2-group comparisons) is the Benjamini-Hochberg adjusted p.values_inverted, via \code{\link{p.adjust}};
#' \item DTU_measure (only available for 2-group comparisons) represents a measure of the intensity of changes between conditions.
#' This measure ranges between 0, when proportions are identical between groups, and 2, when an isoform is always expressed in group A and a different transcript is always chosen in group B;
#' \item Mean log-prec "group_name" indicates the posterior mean of the logarithm of the Dirichlet precision parameter in group "group_name".
#' The precision parameter models the degree of over-dispersion between samples: the higher the precision parameter (or its logarithm), the lower the sample-to-sample variability.
#' \item SD log-prec "group_name" indicates the standard deviation of the logarithm of the Dirichlet precision parameter in group "group_name".
#' }
#' 
#' @slot Transcript_results a \code{data.frame} containing the transcript-level results of the DTU test, structured in the following columns:
#' \itemize{
#' \item Gene_id contains the gene names;
#' \item Transcript_id contains the transcript names;
#' \item p.values is the transcript-level p.values of the DTU test;
#' \item adj.p.values is the Benjamini-Hochberg adjusted p.values (via \code{\link{p.adjust}});
#' \item Max_Gene_Tr.p.val is a conservative p.value and represents the maximum between the transcript p.value and corresponding gene p.value;
#' \item Max_Gene_Tr.Adj.p.val is the Benjamini-Hochberg adjusted Max_Gene_Tr.p.val (via \code{\link{p.adjust}});
#' \item Mean "group_name" indicates the posterior mean of the average relative abundance of the transcript in group "group_name"
#' (an \code{NaN} value indicates that no data was available for a group to estimate parameters);
#' \item SD "group_name" indicates the standard deviation of the average relative abundance of the transcript in group "group_name"
#' (an \code{NaN} value indicates that no data was available for a group to estimate parameters);
#' this column indicates the variability in the mean estimate and is used to plot a 
#' Wald type confidence interval for the mean relative abundance via \code{\link{plot_proportions}}.
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
#' # load the pre-computed results:
#' data("results", package = "BANDITS")
#' show(results)
#' 
#' # Visualize the most significant Genes, sorted by gene level significance.
#' head(top_genes(results))
#' 
#' # Alternatively, gene-level results can also be sorted according to DTU_measure, 
#' # which is a measure of the strength of the change between the 
#' # average relative abundances of the two groups.
#' head(top_genes(results, sort_by = "DTU_measure"))
#' 
#' # Visualize the most significant transcripts, sorted by transcript level significance.
#' head(top_transcripts(results, sort_by = "transcript"))
#' 
#' # Visualize the convergence output for the most significant genes, 
#' # sorted by gene level significance.
#' head(convergence(results))
#' 
#' # We can further use the 'gene' function to gather all output for a specific gene:
#' # gene level, transcript level and convergence results.
#' top_gene = top_genes(results, n = 1)
#' gene(results, top_gene$Gene_id)
#' 
#' # Similarly we can use the 'transcript function to gather all output 
#' # for a specific transcript.
#' top_transcript = top_transcripts(results, n = 1)
#' transcript(results, top_transcript$Transcript_id)
#' 
#' #Finally, we can plot the estimated average transcript relative expression 
#' # in the two groups for a specific gene via 'plot_proportions'.
#' plot_proportions(results, top_gene$Gene_id)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#'
#' @seealso \code{\link{test_DTU}}, \code{\link{create_data}}, \code{\linkS4class{BANDITS_data}}
#' 
#' @aliases BANDITS_test top_genes top_transcripts convergence gene transcript plot_proportions
#' 
#' @export
setClass("BANDITS_test", 
         slots = representation(Gene_results = "data.frame", Transcript_results = "data.frame", 
                                Convergence = "data.frame", samples_design = "data.frame"))

#' @rdname BANDITS_test-class
#' @param object,x a 'BANDITS_test' object.
#' @export
setMethod("show", "BANDITS_test", function(object){
  message(paste0("A 'BANDITS_test' object, with ", nrow(Gene_results(object)), " gene and ", nrow(Transcript_results(object)), " transcript level results."))
})

###############################################################################
### Set validity of the object
###############################################################################
setValidity("BANDITS_test", function(object){
  # Has to return TRUE for a valid object!
  if( !all( Transcript_results(object)$Gene_id %in% Gene_results(object)$Gene_id) ){
    return("All genes appearing in @Transcript_results slot must also appear in @Gene_results slot")
  }
  
  if( !all( Gene_results(object)$Gene_id %in% convergence(object)$Gene_id) ){
    return("All genes appearing in @Gene_results slot must also appear in @Convergence slot")
  }
  
  return(TRUE)
})


###############################################################################
### accessing methods: gene & transcript private, convergence public
###############################################################################
setGeneric("Gene_results", function(x) 
  standardGeneric("Gene_results") )
setMethod("Gene_results", "BANDITS_test", function(x) x@Gene_results)

setGeneric("Transcript_results", function(x) 
  standardGeneric("Transcript_results") )
setMethod("Transcript_results", "BANDITS_test", function(x) x@Transcript_results)

setGeneric("convergence", function(x) 
  standardGeneric("convergence") )

#' @rdname BANDITS_test-class
#' @export
setMethod("convergence", "BANDITS_test", function(x){
  x@Convergence
})

###############################################################################
### retrieve results (gene, transcript and convergence)
###############################################################################
setGeneric("top_genes", function(x, ...) 
  standardGeneric("top_genes") )

#' @rdname BANDITS_test-class
#' @param n the number of genes or transcripts to report. By default \code{n = Inf} and all results will be reported.
#' @param sort_by_g "p.value" for sorting results according to gene-level significance (i.e., p.value);
#' "DTU_measure" for sorting results according to the 'DTU_measure' (check the vignette for details).
#' @export
setMethod("top_genes", "BANDITS_test", function(x, n = Inf, sort_by_g = "p.value"){
  if(!is.character(sort_by_g)){
    message("'sort_by_g' must be a character string, indicating either 'p.value' or 'DTU_measure'")
    return(NULL)
  }
  if(! sort_by_g %in% c("p.value", "DTU_measure")){
    message("'sort_by_g' must be either 'p.value' or 'DTU_measure'")
    return(NULL)
  }
  
  # take the min between the provided number and the size of the matrix:
  n = min(n, nrow(Gene_results(x)) )
  
  if(sort_by_g == "p.value"){ # transcript results sorted according to gene-level p.value
    return(Gene_results(x)[seq_len(n),])
  }
  
  if(is.null(Gene_results(x)$DTU_measure)){
    message("Results cannot be sorted by 'DTU_measure': this feature is only available when two groups are compared.")
    return(NULL)
  }
  
  # transcript results sorted according to transcript-level p.value
  Gene_results(x)[order(Gene_results(x)$DTU_measure, decreasing = TRUE)[seq_len(n)],] # high DTU_measure on top.
})

setGeneric("top_transcripts", function(x, ...) 
  standardGeneric("top_transcripts") )

#' @rdname BANDITS_test-class
#' @param sort_by_tr "gene" for sorting results according to gene-level significance (i.e., p.value);
#' "transcript" for sorting results according to transcript-level significance (i.e., p.value).
#' @export
setMethod("top_transcripts", "BANDITS_test", function(x, n = Inf, sort_by_tr="gene"){
  if(!is.character(sort_by_tr)){
    message("'sort_by_tr' must be a character string")
    return(NULL)
  }
  if(! sort_by_tr %in% c("gene", "transcript")){
    message("'sort_by_tr' must be either 'gene' or 'transcript'")
    return(NULL)
  }
  
  # take the min between the provided number and the size of the matrix:
  n = min(n, nrow(Transcript_results(x)) )
  
  if(sort_by_tr == "gene"){ # transcript results sorted according to gene-level p.value
    return(Transcript_results(x)[seq_len(n),])
  }
  
  if(is.null(Transcript_results(x)$p.values)){
    message("Results cannot be sorted by 'transcript': this feature is only available when two or more groups are compared.")
    return(NULL)
  }
  
  # transcript results sorted according to transcript-level p.value
  Transcript_results(x)[order(Transcript_results(x)$p.values)[seq_len(n)],]
})

###############################################################################
### retrieve individual genes and transcripts
###############################################################################

setGeneric("gene", function(x, ...)
  standardGeneric("gene") )

#' @rdname BANDITS_test-class
#' @param gene_id a character string or vector indicating the gene or genes whose results should be retrieved.
#' @export
setMethod("gene", "BANDITS_test", function(x, gene_id){
  if(!is.character(gene_id)){
    gene_id = as.character(gene_id)
  }
  
  if( ! all( gene_id %in% as.character( Gene_results(x)$Gene_id ) ) ){
    message("One or more genes in 'gene_id' were not found in the results")
    return(NULL)
  }
  
  list( gene_results = Gene_results(x)[ as.character( Gene_results(x)$Gene_id ) %in% gene_id, ],
        transcript_results = Transcript_results(x)[  as.character( Transcript_results(x)$Gene_id ) %in% gene_id, ],
        convergence_results = convergence(x)[ as.character( convergence(x)$Gene_id ) %in% gene_id, ] )
})

setGeneric("transcript", function(x, ...) 
  standardGeneric("transcript") )

#' @rdname BANDITS_test-class
#' @param transcript_id a character string or vector indicating the transcript or transcripts whose results should be retrieved.
#' @export
setMethod("transcript", "BANDITS_test", function(x, transcript_id){
  if(!is.character(transcript_id)){
    transcript_id = as.character(transcript_id)
  }
  
  if( ! transcript_id %in%  as.character(Transcript_results(x)$Transcript_id ) ){
    message("transcript not found in results")
    return(NULL)
  }
  
  gene_id = Transcript_results(x)$Gene_id[ as.character(Transcript_results(x)$Transcript_id) %in% transcript_id]
  list( transcript_results = Transcript_results(x)[ as.character(Transcript_results(x)$Transcript_id) %in% transcript_id, ],
        gene_results = Gene_results(x)[ as.character(Gene_results(x)$Gene_id) %in% gene_id, ],
        convergence_results = convergence(x)[ as.character(convergence(x)$Gene_id) %in% gene_id, ] )
})

###############################################################################
### Plot results for a specific gene
###############################################################################

setGeneric("plot_proportions", function(x, ...) 
  standardGeneric("plot_proportions") )

#' @rdname BANDITS_test-class
#' @param CI a logical element indicating whether to also plot confidence boundaries (TRUE) or not (FALSE).
#' @param CI_level a number between 0 and 1, indicating the level of the confidence interval to plot.
#' @export
setMethod("plot_proportions", "BANDITS_test", function(x, gene_id,
                                                       CI = TRUE,
                                                       CI_level = 0.95){
  if(!is.logical(CI)){
    message("'CI' must be 'TRUE' or 'FALSE'")
    return(NULL)
  }

  if(!is.character(gene_id)){
    gene_id = as.character(gene_id)
  }
  
  if( length(gene_id) != 1){
    message("'gene_id' contains ", length(gene_id), " values: one gene only can be specified.")
    return(NULL)
  }
  
  if( ! gene_id %in% as.character( Transcript_results(x)$Gene_id ) ){
    message("gene not found in results")
    return(NULL)
  }
  
  group_names = levels(factor(x@samples_design$group))
  
  sel = Transcript_results(x)$Gene_id == gene_id
  
  means = as.matrix(Transcript_results(x)[sel, grep("Mean", names(Transcript_results(x)))])
  SD =  Transcript_results(x)[sel, grep("SD", names(Transcript_results(x)))]
  n_groups = ncol(means)
  tr_names = Transcript_results(x)$Transcript_id[sel]
  # impose an order to the transcripts (according to the over-all relative abudance):
  ord = order(rowSums(means), decreasing = TRUE)
  
  prop_samp = data.frame(feature_id = factor( rep(tr_names,n_groups), levels = tr_names[ord]), 
                         proportion = unlist(c(means)),
                         LB = pmax(0, unlist(c(means - qnorm(1 - (1-CI_level)/2) * SD)) ), # LB must be >= 0
                         UB = pmin(1, unlist(c(means + qnorm(1 - (1-CI_level)/2) * SD)) ), # UB must be <= 0
                         group = rep(group_names, each = nrow(means)),
                         stringsAsFactors = FALSE)
  
  # Plot the estimated average proportions of each groups:
  ggp = ggplot() +
    geom_bar(data = prop_samp, aes_string(x = "feature_id", y = "proportion", 
                                          fill = "group"),
             stat = "identity", position = position_dodge(width = 0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
          axis.text = element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16), 
          legend.position = "right", 
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14)) +
    ggtitle(paste("Gene:", gene_id, sep=" ")) +
    xlab("Features") +
    ylab("Proportions") +
    geom_point(data = prop_samp, 
               aes_string(x = "feature_id", y = "proportion", 
                          group = "group", fill = "group"), 
               position = position_dodge(width = 0.9), size = 3, shape = 1, 
               alpha = 0.75)
  
  if(CI){
    ggp = ggp +
      geom_errorbar(data = prop_samp,
                    aes_string(x = "feature_id", ymin = "LB", ymax = "UB",
                               group = "group"), 
                    position = position_dodge(width = 0.9), size = 0.5, 
                    width = 0.5,
                    alpha = 0.5)
  }
  
  ggp
})
