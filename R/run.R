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
#' @param precision a vector with the mean and standard deviation of the log-precision parameter.
#' @param samples_design a \code{data.frame} indicating the design of the experiment with one row for each sample:
#' samples_design must contain a column with the sample id and one with the group id.
#' Warning: the samples in samples_design must have the same order
#' as those in the 'path_to_eq_classes' parameter of the \code{\link{create_data}} function.
#' @param group_col_name the name of the column of 'samples_design' containing the group id.
#' By default group_col_name = "group".
#' @param R the number of iterations for the MCMC algorithm (after the burn-in).
#' Min 10^4.
#' Albeit no difference was observed in simulation studies when increasing 'R' above 10^4, 
#' we encourage users to possibly use higher values of R (e.g., 2*10^4), if the computational time allows it,
#' particularly for comparisons between 3 or more groups.
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
#' # load gene_to_transcript matching:
#' data("gene_tr_id", package = "BANDITS")
#' 
#' # We define the design of the study
#' samples_design = data.frame(sample_id = paste0("sample", seq_len(4)),
#'                             group = c("A", "A", "B", "B"))
#' 
#' # load the pre-computed data:
#' data("input_data", package = "BANDITS")
#' input_data
#' 
#' # Filter lowly abundant genes:
#' input_data = filter_genes(input_data, min_counts_per_gene = 20)
#' 
#' # load the pre-computed precision estimates:
#' data(precision, package = "BANDITS")
#' 
#' ## Test for DTU
#' set.seed(61217)
#' results = test_DTU(BANDITS_data = input_data,
#'                    precision = precision$prior,
#'                    samples_design = samples_design,
#'                    R = 10^4, burn_in = 2*10^3, n_cores = 2,
#'                    gene_to_transcript = gene_tr_id)
#' results
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{create_data}}, \code{\linkS4class{BANDITS_data}}, \code{\linkS4class{BANDITS_test}}
#' 
#' @export
test_DTU = function(BANDITS_data, precision = NULL, R = 10^4, burn_in = 2*10^3,
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
    message("'BANDITS_data' must be a 'BANDITS_data' object created via 'create_data'")
    return(NULL)
  }
  
  # check that gene_to_transcript is a matrix or data.frame object
  if( !is.data.frame(gene_to_transcript) & !is.matrix(gene_to_transcript)  ){
    message("'gene_to_transcript' must be a matrix or data.frame")
    return(NULL)
  }
  
  if( ncol(gene_to_transcript) != 2 ){
    message("'gene_to_transcript' must be a 2 column matrix or data.frame")
    return(NULL)
  }
  
  if(is.null(precision)){
    mean_log_precision = 0
    sd_log_precision = 10
  }else{
    if( {is.list(precision)} | {length(precision) != 2} ){
      message("'precision' must be a vector of length 2")
      return(NULL)
    }
    mean_log_precision = precision[1]
    sd_log_precision   = precision[2]
  }
  
  if( !is.data.frame(samples_design) ){
    message("'samples_design' must be a data.frame object")
    return(NULL)
  }
  
  # select the column of samples_design which is called 'group_col_name'
  if( !(group_col_name %in% colnames(samples_design)) ){
    message("Column ", group_col_name, " missing in 'samples_design'")
    message("'group_col_name' should specify the column name of 'samples_design' containing the group id of each sample")
    return(NULL)
  }
  
  sel_col = which(group_col_name == colnames(samples_design))
  if( length(sel_col) > 1.5 ){
    message( length(sel_col) , " columns from 'samples_design' are called ", group_col_name)
    message("Remove duplicated columns from 'samples_design' and provide a unique column for the group id")
    return(NULL)
  }
  
  groups = samples_design[, sel_col ]
  group_levels = levels(factor(groups))
  N_groups = length(group_levels)
  N = vapply(group_levels,  function(g) sum(groups == g), FUN.VALUE = integer(1)) 
  
  # check if data are already ordered (first A, then B, etc...)
  # if not, ORDER data (equiv classes counts) to respect the ordering in "groups"
  ord_samples = unlist(lapply( group_levels, function(g) which(groups == g) ) )
  
  #########################################################################################################
  # Give priority to "most" computationally intensive genes (total number of transcripts):
  #########################################################################################################
  # store eff_len of Unique and Together transcritps:
  eff_len_tr_Unique   = effLen(BANDITS_data)[uniqueId(BANDITS_data) == TRUE]
  eff_len_tr_Together = effLen(BANDITS_data)[uniqueId(BANDITS_data) == FALSE]
  
  K_tot_Together = vapply( eff_len_tr_Together, length, FUN.VALUE = integer(1))
  # sapply( eff_len_tr_Together, length) 
  order_Together = order(K_tot_Together, decreasing = TRUE)
  
  K_tot_Unique = vapply( eff_len_tr_Unique, length, FUN.VALUE = integer(1))
  # sapply( eff_len_tr_Unique, length)
  order_Unique = order(K_tot_Unique, decreasing = TRUE)
  
  # First analyze ALL together genes (in order), then analyze all Unique genes (in order):
  final_order = c(length(order_Unique) + order_Together, order_Unique)
  
  #########################################################################################################
  # Infer model parameters from 1 group only:
  #########################################################################################################
  if(length(group_levels) == 1){
    message("One group only is present in 'samples_design$group': 
  BANDITS will infer model parameters, but will not test for DTU between groups.")
    return( infer_one_group(BANDITS_data = BANDITS_data,
                            mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                            R = R, burn_in = burn_in, n_cores = n_cores,
                            final_order = final_order, ord_samples = ord_samples,
                            group_levels = group_levels, samples_design = samples_design,
                            gene_to_transcript = gene_to_transcript) )
  }
  
  if( sum( N > 0.5 ) < 1.5){ # check that >=2 groups have >=1 samples
    message("At least two distinct groups must be present in 'samples_design$group' to perform DTU")
    return(NULL)
  }
  
  
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
  # if n_cores > 1 (parallel computing).
  suppressWarnings({
    cl = makeCluster(n_cores, setup_strategy = "sequential")
  })
  registerDoParallel(cl, n_cores);
  
  message("Starting the MCMC")
  
  p_values_ALL = foreach(p = final_order,
                         .packages=c("BANDITS"),
                         .errorhandling = "stop") %dorng%{
                           # ORDER counts to respect the ordering in "groups" via 'ord_samples'
                           f = counts(BANDITS_data)[[p]][,ord_samples]
                           # AND MIN 5 counts per group:
                           cond_min5_perGroup = {sum(f[,seq_len(N[1])]) > 4.5} & {sum(f[,N[1] + seq_len(N[2])]) > 4.5}
                           
                           # select samples with data only:
                           sel_samples = colSums(f) > 0.5
                           f = as.matrix( f[, sel_samples ] )
                           N1 = sum(sel_samples[seq_len(N[1])])
                           N2 = sum(sel_samples[N[1] + seq_len(N[2])])
                           
                           # only run the MCMC if there is at least 1 sample with at least 1 count in each condition:
                           cond_f = cond_min5_perGroup & {N1 > 0.5} & {N2 > 0.5}
                           
                           # initialize results:
                           res = NULL
                           
                           # if it's NOT Unique.
                           if( cond_f & !uniqueId(BANDITS_data)[[p]] ){ # for the first (Together) elements I run the together function
                             res = wald_DTU_test_Together_FULL(f = f,
                                                               l = effLen(BANDITS_data)[[p]], 
                                                               exon_id = classes(BANDITS_data)[[p]],
                                                               genes = genes(BANDITS_data)[[p]],
                                                               transcripts = transcripts(BANDITS_data)[[p]],
                                                               mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                               R = 2*R, burn_in = 2*burn_in, N_1 = N1, N_2 = N2,
                                                               theshold_pval = theshold_pval)
                             # double iterations for Together genes: they typically require more iter than Unique genes (more complex posterior space to explore).
                           }else{ # for the following elements I run the Unique function
                             if( cond_f & length(effLen(BANDITS_data)[[p]]) > 1 ){ # run test only if there are at least 2 transcripts per gene.
                               res = wald_DTU_test_FULL( f = f,
                                                         l = effLen(BANDITS_data)[[p]], 
                                                         exon_id = classes(BANDITS_data)[[p]],
                                                         mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                         R = R, burn_in = burn_in, N_1 = N1, N_2 = N2,
                                                         theshold_pval = theshold_pval)
                               if(length(res[[1]]) > 1){
                                 res[[1]][[1]] = matrix(res[[1]][[1]], nrow = 1)
                                 rownames(res[[1]][[1]]) = genes(BANDITS_data)[[p]]
                                 names(res[[1]][[2]]) = transcripts(BANDITS_data)[[p]]
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
  
  suppressWarnings({ p_values = apply(p_values, 2, as.numeric) })
  rownames(p_values) = gene_names
  
  SEL_ALL  = !is.na(p_values[,1])  # remove NA's (genes not converged).
  p_values = p_values[SEL_ALL,]
  
  SEL_ALL  = p_values[,1] != -1  # remove -1's (genes not analyzed).
  p_values = p_values[SEL_ALL,]
  
  # sort p.values according to their significance.
  p_values = p_values[ order(p_values[,1]), ]
  
  # COMPUTE ADJUSTED P.VALS:
  adj.p_values = p.adjust(p_values[,1], method = "BH") # gene-test
  
  gene_DF = data.frame(Gene_id = rownames(p_values), 
                       p.values = p_values[,1], 
                       adj.p.values = adj.p_values,
                       p.values_inverted = ifelse(p_values[,2], p_values[,1], sqrt(p_values[,1])),    # if NOT inverted, take sqrt(p.val)
                       adj.p.values_inverted = ifelse(p_values[,2], adj.p_values, sqrt(adj.p_values)), # if NOT inverted, take sqrt(adj.p.val)
                       DTU_measure = p_values[,3], # sum of top two differences of posterior modes between transcripts.
                       p_values[,seq.int(4,7, by = 1)],
                       row.names = NULL)
  
  names(gene_DF)[ c(7,8) ] = paste("Mean log-prec", group_levels)
  names(gene_DF)[ c(9,10) ] = paste("SD log-prec", group_levels)
  
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
  
  sd_A = unlist( lapply(p_values_ALL, function(x){
    y = x[[1]]
    if(length(y) > 1){
      return(y[[5]])
    }else{
      return(NULL)
    }} ) )
  
  sd_B = unlist( lapply(p_values_ALL, function(x){
    y = x[[1]]
    if(length(y) > 1){
      return(y[[6]])
    }else{
      return(NULL)
    }} ) )
  
  p_values_tr = cbind(p_values_tr, mode_A, mode_B, sd_A, sd_B)
  
  cond = !vapply(p_values_tr[,1], is.null, FUN.VALUE = logical(1))
  # sapply(p_values_tr[,1], is.null) == FALSE
  p_values_tr = p_values_tr[cond,] # filter null results
  
  cond = p_values_tr[,1] != -1
  p_values_tr = p_values_tr[cond,] # filter -1 results
  
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
                     p_values_tr[,seq.int(2, 5, by = 1)], 
                     row.names = NULL)
  
  # re-name the group names according to the groups names:
  names(tr_DF)[ c(7,8) ]  = paste("Mean", group_levels)
  names(tr_DF)[ c(9,10) ] = paste("SD", group_levels)
  
  # sort p.values according to their GENE significance (and secondly according to tr significance):
  tr_DF = tr_DF[ order(p_values[match(tr_DF$Gene_id, rownames(p_values)), 1], tr_DF$p.values) ,]
  
  # set rownames to 1, 2, ..., nrow(tr_DF)
  rownames(tr_DF) = seq_len(nrow(tr_DF))
  
  #########################################################################################################
  # Return Convergence results:
  #########################################################################################################
  convergence = lapply(p_values_ALL, function(x) x[[2]])
  n_genes_convergence = vapply( genes(BANDITS_data)[final_order], length, FUN.VALUE = integer(1))
  # sapply( genes(BANDITS_data)[final_order], length)
  genes_convergence = unlist(genes(BANDITS_data)[final_order])
  
  convergence = rep(convergence, n_genes_convergence)
  
  num = vapply(convergence, is.numeric, FUN.VALUE = logical(1))
  # sapply(convergence, class) == "numeric"
  convergence = convergence[num]
  convergence = do.call(rbind, convergence)
  rownames(convergence) = genes_convergence[num]
  
  convergence = convergence[ order(p_values[match(rownames(convergence), rownames(p_values)), 1]) ,]
  
  # remove gene ids that are not in "all_genes(BANDITS_data)"; i.e., remove gene ids for (Together) genes with 1 transcript only!
  convergence = convergence[ rownames(convergence) %in% all_genes(BANDITS_data), ]
  
  #########################################################################################################
  # Return results:
  #########################################################################################################
  message("Returning results")
  
  # p_values[,1] contains the p.values
  # p_values[,2] contains the inversion (0 or 1)
  # p_values[,3] contains the DTU_measure
  
  return(new("BANDITS_test",
             Gene_results = gene_DF, 
             Transcript_results =  tr_DF,
             Convergence = data.frame(Gene_id = rownames(convergence),
                                      converged = ifelse(convergence[,1], TRUE, FALSE), # did the gene converge or NOT ?
                                      burn_in   = (convergence[,2]-1)/R,
                                      row.names = NULL),
             samples_design = samples_design
  )
  )
}
