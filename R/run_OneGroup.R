infer_one_group = function(BANDITS_data, mean_log_precision, sd_log_precision,
                           R, burn_in, 
                           n_cores,
                           final_order, ord_samples,
                           group_levels, samples_design,
                           gene_to_transcript){
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
                           cond_min5_perGroup = sum(f) > 4.5

                           # select samples with data only:
                           sel_samples = colSums(f) > 0.5
                           f = as.matrix( f[, sel_samples ] )
                           N = sum(sel_samples)
                           
                           # only run the MCMC if there is at least 1 sample with at least 1 count in each condition:
                           cond_f = cond_min5_perGroup & {N > 0.5}
                           
                           # initialize results:
                           res = NULL
                           
                           # if it's NOT Unique.
                           if( cond_f & !uniqueId(BANDITS_data)[[p]] ){ # for the first (Together) elements I run the together function
                             res = infer_only_Together(f = f,
                                                       l = effLen(BANDITS_data)[[p]], 
                                                       exon_id = classes(BANDITS_data)[[p]],
                                                       genes = genes(BANDITS_data)[[p]],
                                                       transcripts = transcripts(BANDITS_data)[[p]],
                                                       mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                       R = 2*R, burn_in = 2*burn_in, N = N)
                             # double iterations for Together genes: they typically require more iter than Unique genes (more complex posterior space to explore).
                           }else{ # for the following elements I run the Unique function
                             if( cond_f & length(effLen(BANDITS_data)[[p]]) > 1 ){ # run test only if there are at least 2 transcripts per gene.
                               res = infer_only( f = f,
                                                 l = effLen(BANDITS_data)[[p]], 
                                                 exon_id = classes(BANDITS_data)[[p]],
                                                 mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                 R = R, burn_in = burn_in, N = N)
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
  p_values = p_values[ order(rownames(p_values)), ]
  
  gene_DF = data.frame(Gene_id = rownames(p_values), 
                       p_values,
                       row.names = NULL)
  
  names(gene_DF)[ 2 ] = paste("Mean log-prec", group_levels)
  names(gene_DF)[ 3 ] = paste("SD log-prec", group_levels)
  
  #########################################################################################################
  # Gather together TRANSCRIPT level results:
  #########################################################################################################
  mode_A = unlist( lapply(p_values_ALL, function(x){
    y = x[[1]]
    if(length(y) > 1){
      return(y[[2]])
    }else{
      return(NULL)
    }} ) )
  
  sd_A = unlist( lapply(p_values_ALL, function(x){
    y = x[[1]]
    if(length(y) > 1){
      return(y[[3]])
    }else{
      return(NULL)
    }} ) )
  
  p_values_tr = cbind(mode_A, sd_A)
  
  cond = !vapply(p_values_tr[,1], is.null, FUN.VALUE = logical(1))
  p_values_tr = p_values_tr[cond,] # filter null results
  
  cond = p_values_tr[,1] != -1
  p_values_tr = p_values_tr[cond,] # filter -1 results
  
  # match GENE and TRANSCRIPT IDs.
  Gene_id     = as.character(gene_to_transcript[,1]); Tr_id = as.character(gene_to_transcript[,2]) 
  genes_in_tr = Gene_id[match(rownames(p_values_tr), Tr_id)]
  
  tr_DF = data.frame(Gene_id = genes_in_tr, 
                     Transcript_id = rownames(p_values_tr), 
                     p_values_tr, 
                     row.names = NULL)
  
  # re-name the group names according to the groups names:
  names(tr_DF)[ 3 ]  = paste("Mean", group_levels)
  names(tr_DF)[ 4 ] = paste("SD", group_levels)

  # sort transcripts by GENE NAME (increasing names), and secondly according to transcript relative abundance (decreasing probs):
  tr_DF = tr_DF[ order(gene_DF$Gene_id[match(tr_DF$Gene_id, gene_DF$Gene_id)], -tr_DF[,3]), ]
  
  # set rownames to 1, 2, ..., nrow(tr_DF)
  rownames(tr_DF) = seq_len(nrow(tr_DF))
  
  #########################################################################################################
  # Return Convergence results:
  #########################################################################################################
  convergence = lapply(p_values_ALL, function(x) x[[2]])
  n_genes_convergence = vapply( genes(BANDITS_data)[final_order], length, FUN.VALUE = integer(1))
  genes_convergence = unlist(genes(BANDITS_data)[final_order])
  
  convergence = rep(convergence, n_genes_convergence)
  
  num = vapply(convergence, is.numeric, FUN.VALUE = logical(1))
  # sapply(convergence, class) == "numeric"
  convergence = convergence[num]
  convergence = do.call(rbind, convergence)
  rownames(convergence) = genes_convergence[num]
  
  # order convergence DF to keep the same order as gene names:
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
