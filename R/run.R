#' Perform differential splicing
#'
#' \code{test_DTU} performs differential splicing, via differential transcript usage (DTU), between 2 or more groups.
#' Parameters are inferred via Markov chain Monte Carlo (MCMC) techniques and a DTU test is performed 
#' via a multivariate Wald test on the posterior densities for the average relative abundance of transcripts.
#' Warning: the samples in samples_design must have the same order
#' as those in the 'path_to_eq_classes' parameter of the \code{\link{create_data}} function.
#' 
#' 
#' @param BANDITS_data a 'BANDITS_data' object.
#' @param prior_precision a vector with the mean and standard deviation of the log-precision parameter.
#' @param samples_design a \code{data.frame} indicating the design of the experiment with one row for each sample:
#' samples_design must contain a column with the sample id and one with the group id.
#' Warning: the samples in samples_design must have the same order
#' as those in the 'path_to_eq_classes' parameter of the \code{\link{create_data}} function.
#' @param group_col_name the name of the column of 'samples_design' containing the group id.
#' By default group_col_name = "group".
#' @param R the number of iterations for the MCMC algorithm (after the burn-in).
#' Min 10^4.
#' Albeit no difference was observed in simulation studies when increasing 'R' above 10^4, 
#' we encourage users to possibly use higher values of R (e.g., double) if the computational time allows it.
#' @param burn_in the length of the burn-in to be discarded (before convergence is reached).
#' Min 2*10^3.
#' Albeit no difference was observed in simulation studies when increasing 'burn_in' above 2*10^3, 
#' we encourage users to possibly use higher values of R (e.g., double) if the computational time allows it.
#' @param n_cores the number of cores to parallelize the tasks on.
#' @param gene_to_transcript a matrix or data.frame with a list of gene-to-transcript correspondances.
#' The first column represents the gene id, while the second one contains the transcript id.
#' @param theshold_pval is a threshold between 0 and 1; when running \code{\link{test_DTU}}, if the p.value of a gene is < theshold_pval,
#' a second (independent) MCMC chain is run and the p.value is re-computed on the aggregation of the two chains.
#' By defauls theshold_pval = 0.1, while theshold_pval = 1 corresponds to running all chains twice, and theshold_pval = 0 means all chains will only run once.
#' 
#' @return A \code{\linkS4class{BANDITS_test}} object.
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
#' # Warning: the sample names in equiv_classes_files must have the same order
#' # as those in the design object, containted in samples_design.
#' equiv_classes_files
#' samples_design$sample_id
#' 
#' # create data and filter internally lowly abundant transcripts:
#' BANDITS_data = create_data(gene_to_transcript = gene_tr_id,
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
#' #set.seed(61217)
#' #prec = prior_precision(gene_to_transcript = gene_tr_id, transcript_counts = counts,
#' #                       min_transcript_proportion = 0.01, min_transcript_counts = 10,
#' #                       min_gene_counts = 20, n_cores = 2)
#' 
#' # load the pre-computes precision estimates:
#' data(prec, package = "BANDITS")
#' 
#' # Plot the histogram of the genewise log-precision estimates.
#' # The black solid line represents the normally distributed prior distribution 
#' # for the log-precision parameter.
#' plot_precision(prec)
#' 
#' 
#' 
#' ## Test for DTU
#' #set.seed(61217)
#' #x = test_DTU(BANDITS_data = BANDITS_data,
#' #             prior_precision = prec$prior,
#' #             samples_design = samples_design,
#' #             R = 10^4, burn_in = 2*10^3, n_cores = 2,
#' #             gene_to_transcript = gene_tr_id)
#' 
#' # load the pre-computed results:
#' data("DTU_results", package = "BANDITS")
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
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{create_data}}, \code{\linkS4class{BANDITS_data}}, \code{\linkS4class{BANDITS_test}}
#' 
#' @export
test_DTU = function(BANDITS_data, prior_precision = NULL, R = 10^4, burn_in = 2*10^3,
                    samples_design, group_col_name = "group", n_cores = 1, 
                    gene_to_transcript, theshold_pval = 0.1){
  #########################################################################################################
  # Give priority to "most" computationally intensive genes (total number of transcripts):
  #########################################################################################################
  if(R < 10^4){
    message("'R' must be at least 10^4")
    return(NULL)
  }
  
  if(burn_in < 2*10^3){
    message("'burn_in' must be at least 2*10^3")
    return(NULL)
  }
  
  if( !is(BANDITS_data, "BANDITS_data") ){
    message("'BANDITS_data' must be a 'BANDITS_data' object created via the create_data function")
    return(NULL)
  }
  
  if( ncol(gene_to_transcript) != 2 ){
    message("'gene_to_transcript' must be a 2 column matrix or data.frame")
    return(NULL)
  }
  
  if(is.null(prior_precision)){
    mean_log_precision = 0
    sd_log_precision = 10
  }else{
    mean_log_precision = prior_precision[1]
    sd_log_precision   = prior_precision[2]
  }
  
  if( !is.data.frame(samples_design) ){
    message("'samples_design' must be a data.frame object")
    return(NULL)
  }
  
  # select the column of samples_design which is called 'group_col_name'
  sel_col = which(group_col_name == colnames(samples_design))
  if( is.null(sel_col) ){
    message("Column ", group_col_name, " missing in 'samples_design'")
    message("'group_col_name' should specify the column name of 'samples_design' containing the group id of each sample")
    return(NULL)
  }
  if( length(sel_col) > 1.5 ){
    message( length(sel_col) , " columns from 'samples_design' are called ", group_col_name)
    message("Remove duplicated columns from 'samples_design' and provide a unique column for the group id")
    return(NULL)
  }
  
  groups = samples_design[, sel_col ]
  group_levels = levels(groups)
  N_groups = length(group_levels)
  N = vapply(group_levels,  function(g) sum(groups == g), FUN.VALUE = integer(1)) 
  # sapply( group_levels, function(g) sum(groups == g) )
  
  if( sum( N > 0.5 ) < 1.5){ # check that >=2 groups have >=1 samples
    message("At least two distinct groups must be present in 'samples_design$group'")
    return(NULL)
  }
  
  # check if data are already ordered (first A, then B, etc...)
  # if not, ORDER data (equiv classes counts) to respect the ordering in "groups"
  ord_samples = unlist(lapply( group_levels, function(g) which(groups == g) ) )
  
  #########################################################################################################
  # Give priority to "most" computationally intensive genes (total number of transcripts):
  #########################################################################################################
  # store eff_len of Unique and Together transcritps:
  eff_len_tr_Unique   = BANDITS_data@effLen[BANDITS_data@uniqueId == TRUE]
  eff_len_tr_Together = BANDITS_data@effLen[BANDITS_data@uniqueId == FALSE]
  
  K_tot_Together = vapply( eff_len_tr_Together, length, FUN.VALUE = integer(1))
  # sapply( eff_len_tr_Together, length) 
  order_Together = order(K_tot_Together, decreasing = TRUE)
  
  K_tot_Unique = vapply( eff_len_tr_Unique, length, FUN.VALUE = integer(1))
  # sapply( eff_len_tr_Unique, length)
  order_Unique = order(K_tot_Unique, decreasing = TRUE)
  
  # First analyze ALL together genes (in order), then analyze all Unique genes (in order):
  final_order = c(length(order_Unique) + order_Together, order_Unique)
  
  #########################################################################################################
  # Multi-group testing if there are 3 or more groups:
  #########################################################################################################
  if(N_groups > 2.5){ # if there are 3 or more groups, use the Multi-group test:
    return( test_DTU_multi_group(BANDITS_data = BANDITS_data,
                                 mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                 R = R, burn_in = burn_in, N = N, n_cores = n_cores,
                                 final_order = final_order, ord_samples = ord_samples, 
                                 group_levels = group_levels, samples_design = samples_design,
                                 gene_to_transcript = gene_to_transcript, theshold_pval = theshold_pval) )
  }
  
  #########################################################################################################
  # Parallelize ALL genes/groups at the same time.
  #########################################################################################################
  if(n_cores > 1){ # if n_cores > 1 (parallel computing).
    suppressWarnings({
      cl <- makeCluster(n_cores);
    })
    registerDoParallel(cl, n_cores);
  }
  
  message("Starting the MCMC")
  
  p_values_ALL = foreach(p = final_order,
                         .packages=c("BANDITS"),
                         .errorhandling = "stop") %dorng%{
                           # ORDER counts to respect the ordering in "groups" via 'ord_samples'
                           f = as.matrix(BANDITS_data@counts[[p]][,ord_samples])
                           sel_samples = colSums(f) > 0.5
                           
                           # select samples with data only:
                           f = f[, sel_samples ]
                           N1 = sum(sel_samples[seq_len(N[1])])
                           N2 = sum(sel_samples[N[1] + seq_len(N[2])])
                           
                           # only run the MCMC if there is at least 1 sample with at least 1 count in each condition:
                           cond_f = {N1 > 0.5} & {N2 > 0.5}
                           
                           # initialize results:
                           res = NULL
                           
                           # if it's NOT Unique.
                           if( cond_f & !BANDITS_data@uniqueId[[p]] ){ # for the first (Together) elements I run the together function
                             res = wald_DTU_test_Together_FULL(f = f,
                                                               l = BANDITS_data@effLen[[p]], 
                                                               exon_id = BANDITS_data@classes[[p]],
                                                               genes = BANDITS_data@genes[[p]],
                                                               transcripts = BANDITS_data@transcripts[[p]],
                                                               mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                               R = 2*R, burn_in = 2*burn_in, N_1 = N1, N_2 = N2,
                                                               theshold_pval = theshold_pval)
                             # double iterations for Together genes: they typically require more iter than Unique genes (more complex posterior space to explore).
                           }else{ # for the following elements I run the Unique function
                             if( cond_f & length(BANDITS_data@effLen[[p]]) > 1 ){ # run test only if there are at least 2 transcripts per gene.
                               res = wald_DTU_test_FULL( f = f,
                                                         l = BANDITS_data@effLen[[p]], 
                                                         exon_id = BANDITS_data@classes[[p]],
                                                         mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                         R = R, burn_in = burn_in, N_1 = N1, N_2 = N2,
                                                         theshold_pval = theshold_pval)
                               if(length(res[[1]]) > 1){
                                 res[[1]][[1]] = matrix(res[[1]][[1]], nrow = 1)
                                 rownames(res[[1]][[1]]) = BANDITS_data@genes[[p]]
                                 names(res[[1]][[2]]) = BANDITS_data@transcripts[[p]]
                               } # only if the mcmc has a return value.
                             }
                           }
                           res
                         }
  
  message("MCMC completed")
  
  stopCluster(cl) 
  stopImplicitCluster()
  
  #########################################################################################################
  # Gather together GENE level results:
  #########################################################################################################
  p_values = lapply(p_values_ALL, function(x) x[[1]][[1]])
  p_values = do.call(rbind, p_values)
  
  gene_names = rownames(p_values)
  # or simply:
  # mean( rownames(p_values) == unlist( BANDITS_data@genes[final_order] ))
  # 1, same!
  
  suppressWarnings({ p_values = apply(p_values, 2, as.numeric) })
  rownames(p_values) = gene_names
  
  # sort p.values according to their significance.
  p_values = p_values[ order(p_values[,1]), ]
  
  SEL_ALL  = (p_values[,1] != -1) | is.na(p_values[,1])  # remove -1's (genes not analyzed).
  p_values = p_values[SEL_ALL,]
  
  # COMPUTE ADJUSTED P.VALS:
  adj.p_values    = p.adjust(p_values[,1], method = "BH") # gene-test
  
  #########################################################################################################
  # Gather together TRANSCRIPT level results:
  #########################################################################################################
  p_values_tr = unlist( lapply(p_values_ALL, function(x){
    y = x[[1]]
    if(length(y) > 1){
      return(y[[2]])
    }else{
      return(NULL)
    }} ) )
  
  mode_A = unlist( lapply(p_values_ALL, function(x){
    y = x[[1]]
    if(length(y) > 1){
      return(y[[3]])
    }else{
      return(NULL)
    }} ) )
  
  mode_B = unlist( lapply(p_values_ALL, function(x){
    y = x[[1]]
    if(length(y) > 1){
      return(y[[4]])
    }else{
      return(NULL)
    }} ) )
  
  p_values_tr = cbind(p_values_tr, mode_A, mode_B)
  
  cond = !vapply(p_values_tr[,1], is.null, FUN.VALUE = logical(1))
  # sapply(p_values_tr[,1], is.null) == FALSE
  p_values_tr = p_values_tr[cond,] # filter null results
  
  cond = p_values_tr[,1] != -1
  p_values_tr = p_values_tr[cond,] # filter null results
  
  # COMPUTE ADJUSTED P.VALS:
  adj.p_values_tr = p.adjust(p_values_tr[,1],  method = "BH") # transcript-test
  
  # match GENE and TRANSCRIPT IDs.
  Gene_id     = as.character(gene_to_transcript[,1]); Tr_id = as.character(gene_to_transcript[,2]) 
  genes_in_tr = Gene_id[match(rownames(p_values_tr), Tr_id)]
  
  # Compute "conservative" transcript level test,
  # take min between gene and transcript level tests (CHECK Koen package):
  max_gene_tr_p.val     = apply( cbind(p_values_tr, p_values[match(genes_in_tr, rownames(p_values)),1]), 1, max)
  max_gene_tr_adj.p.val = apply( cbind(adj.p_values_tr, adj.p_values[match(genes_in_tr, names(adj.p_values))]), 1, max)
  
  tr_DF = data.frame(Gene_id = genes_in_tr, Transcript_id = rownames(p_values_tr), 
                     p.values = p_values_tr[,1], adj.p.values = adj.p_values_tr,
                     Max_Gene_Tr.p.val = max_gene_tr_p.val, Max_Gene_Tr.Adj.p.val = max_gene_tr_adj.p.val,
                     Mode_A = p_values_tr[,2], 
                     Mode_B = p_values_tr[,3],
                     row.names = NULL)
  
  # re-name the group names according to the groups names:
  names(tr_DF)[ seq(7,8) ] = paste("Mean", group_levels)
  
  # sort p.values according to their GENE significance (and secondly according to tr significance):
  tr_DF = tr_DF[ order(p_values[match(tr_DF$Gene_id, rownames(p_values)), 1], tr_DF$p.values) ,]
  
  #########################################################################################################
  # Return Convergence results:
  #########################################################################################################
  convergence = lapply(p_values_ALL, function(x) x[[2]])
  n_genes_convergence = vapply( BANDITS_data@genes[final_order], length, FUN.VALUE = integer(1))
  # sapply( BANDITS_data@genes[final_order], length)
  genes_convergence = unlist(BANDITS_data@genes[final_order])
  
  convergence = rep(convergence, n_genes_convergence)
  
  num = vapply(convergence, is.numeric, FUN.VALUE = logical(1))
  # sapply(convergence, class) == "numeric"
  convergence = convergence[num]
  convergence = do.call(rbind, convergence)
  rownames(convergence) = genes_convergence[num]
  
  convergence = convergence[ order(p_values[match(rownames(convergence), rownames(p_values)), 1]) ,]
  
  # remove gene ids that are not in "BANDITS_data@all_genes"; i.e., remove gene ids for (Together) genes with 1 transcript only!
  convergence = convergence[ rownames(convergence) %in% BANDITS_data@all_genes, ]
  
  #########################################################################################################
  # Return results:
  #########################################################################################################
  message("Returning results")
  
  return(new("BANDITS_test",
             Gene_results = data.frame(Gene_id = rownames(p_values), 
                                       p.values = p_values[,1], adj.p.values = adj.p_values,
                                       p.values_inverted = ifelse(p_values[,5], p_values[,1], sqrt(p_values[,1])),    # if NOT inverted, take sqrt(p.val)
                                       adj.p.values_inverted = ifelse(p_values[,5], adj.p_values, sqrt(adj.p_values)), # if NOT inverted, take sqrt(adj.p.val)
                                       DTU_measure = p_values[,7], # sum of top two differences of posterior modes between transcripts.
                                       row.names = NULL), 
             Transcript_results =  tr_DF,
             Convergence = data.frame(Gene_id = rownames(convergence),
                                      converged = ifelse(convergence[,1], TRUE, FALSE), # did the gene converge or NOT ?
                                      burn_in   = (convergence[,2]-1)/R,
                                      row.names = NULL),
             samples_design = samples_design
  )
  )
}
