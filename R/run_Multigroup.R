test_DTU_multi_group = function(BANDITS_data, mean_log_precision, sd_log_precision,
                                R, burn_in, N, n_cores,
                                final_order, ord_samples, 
                                group_levels, samples_design,
                                gene_to_transcript, theshold_pval){
  #########################################################################################################
  # Parallelize ALL genes/groups at the same time.
  #########################################################################################################
  # IN PARALLEL:
  if(n_cores > 1){ # if n_cores > 1 (parallel computing).
    suppressWarnings({
      cl = makeCluster(n_cores, setup_strategy = "sequential")
    })
    registerDoParallel(cl, n_cores);
  }
  
  N_tot = sum(N)
  cumulative = c(0,cumsum(N))
  splits = list()
  splits = lapply(seq_len( length(cumulative) - 1), function(i){
    {cumulative[i]+1}:cumulative[i+1]
  })
  
  message("Starting the MCMC")
  
  p_values_ALL = foreach(p = final_order,
                         .packages=c("BANDITS"),
                         .errorhandling = "stop") %dorng%{
                           # ORDER counts to respect the ordering in "groups" via 'ord_samples'
                           f = counts(BANDITS_data)[[p]][,ord_samples]
                           sel_samples = colSums(f) > 0.5
                           
                           # AND MIN 5 counts per group:
                           cond_min5_perGroup = vapply(splits, function(id) sum(f[,id]) > 4.5, FUN.VALUE = logical(1))
                           
                           # select samples with BANDITS_data only:
                           f = f[, sel_samples ]
                           N_mg = vapply(splits, function(id) sum(sel_samples[id]), FUN.VALUE = integer(1))
                           # sapply(splits, function(id) sum(sel_samples[id]))
                           sel_groups = cond_min5_perGroup & {N_mg > 0.5} # sel groups to keep (if at least 1 sample!):
                           N_mg = N_mg[ sel_groups ] # remove groups with 0 counts.
                           
                           # only run the MCMC if there are at least 2 groups with at least 1 sample with counts
                           cond_f = length(N_mg) > 1.5
                           
                           # initialize results:
                           res = NULL
                           
                           if(cond_f){ # there are at least 2 groups
                             if(length(N_mg) > 2.5){ # if there are > 2 samples, run the Multi-group comparison
                               
                               if( !uniqueId(BANDITS_data)[[p]] ){ # for the first (Together) elements I run the together function
                                 res = wald_DTU_test_MultiGroup_Together(f = f,
                                                                         l = effLen(BANDITS_data)[[p]], 
                                                                         exon_id = classes(BANDITS_data)[[p]],
                                                                         genes = genes(BANDITS_data)[[p]],
                                                                         transcripts = transcripts(BANDITS_data)[[p]],
                                                                         mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                                         R = 2*R, burn_in = 2*burn_in, N = N_mg,
                                                                         theshold_pval = theshold_pval)
                                 if( length(res[[1]]) > 1 ){
                                   sel_tmp = c(1, 1 + which(sel_groups), 1 + length(sel_groups) + which(sel_groups))
                                   res_tmp = matrix(NaN, nrow = nrow(res[[1]][[1]]), ncol = 1 + 2 * length(sel_groups)) # rows = transcripts; cols = groups
                                   res_tmp[ , sel_tmp] = res[[1]][[1]]
                                   rownames(res_tmp) = rownames(res[[1]][[1]])
                                   res[[1]][[1]] = res_tmp
                                   
                                   if(!is.null(res[[1]][[3]])){ # if trans result is not null (i.e., if the gene was analyzed):
                                     # Store the transcripts modes in a matrix, leave empty groups which were not analyzed:
                                     res[[1]][[3]] = do.call(rbind, res[[1]][[3]]) # gather together results from different genes.
                                     res[[1]][[4]] = do.call(rbind, res[[1]][[4]])
                                     
                                     # add a null column for the group which was not analyzed:
                                     trans_mode = matrix(NaN, nrow = nrow(res[[1]][[3]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                     trans_mode[, sel_groups ] = res[[1]][[3]]
                                     res[[1]][[3]] = trans_mode # rows = transcripts; cols = groups
                                     
                                     # add a null column for the group which was not analyzed:
                                     trans_sd = matrix(NaN, nrow = nrow(res[[1]][[4]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                     trans_sd[, sel_groups ] = res[[1]][[4]]
                                     res[[1]][[4]] = trans_sd # rows = transcripts; cols = groups
                                   }
                                 }                                 
                                 
                               }else{ # for the following elements I run the Unique function
                                 if( length(effLen(BANDITS_data)[[p]]) > 1 ){ # run test only if there are at least 2 transcripts per gene.
                                   res = wald_DTU_test_MultiGroup(f = as.data.frame(f),
                                                                  l = effLen(BANDITS_data)[[p]], 
                                                                  exon_id = classes(BANDITS_data)[[p]],
                                                                  mean_log_precision = mean_log_precision, 
                                                                  sd_log_precision = sd_log_precision,
                                                                  R = R, burn_in = burn_in, N = N_mg,
                                                                  theshold_pval = theshold_pval)
                                   if(length(res[[1]]) > 1){
                                     # sort p.value (pos 1) and post mean and sd of log-precision:
                                     sel_tmp = c(1, 1 + which(sel_groups), 1 + length(sel_groups) + which(sel_groups))
                                     res_tmp = matrix(NaN, nrow = 1, ncol = 1 + 2 * length(sel_groups)) # rows = transcripts; cols = groups
                                     res_tmp[ , sel_tmp] = res[[1]][[1]]
                                     res[[1]][[1]] = res_tmp
                                     rownames(res[[1]][[1]]) = genes(BANDITS_data)[[p]]
                                     
                                     names(res[[1]][[2]]) = transcripts(BANDITS_data)[[p]]

                                     # Store the transcripts modes in a matrix, leave empty groups which were not analyzed:
                                     trans_mode = matrix(NaN, nrow = length(transcripts(BANDITS_data)[[p]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                     trans_mode[ , sel_groups] = res[[1]][[3]]
                                     res[[1]][[3]] = trans_mode # rows = transcripts; cols = groups
                                     
                                     trans_sd = matrix(NaN, nrow = length(transcripts(BANDITS_data)[[p]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                     trans_sd[ , sel_groups] = res[[1]][[4]]
                                     res[[1]][[4]] = trans_sd # rows = transcripts; cols = groups
                                   } # only if the mcmc has a return value.
                                 }
                               }
                               
                             }else{  # if there are only 2 samples left (with counts), run the classical 2-group comparison:
                               f = as.matrix(f) # the 2-group comparison requires f to be a matrix
                               
                               if( !uniqueId(BANDITS_data)[[p]] ){ # for the first (Together) elements I run the together function
                                 res = wald_DTU_test_Together_FULL(f = f,
                                                                   l = effLen(BANDITS_data)[[p]], 
                                                                   exon_id = classes(BANDITS_data)[[p]],
                                                                   genes = genes(BANDITS_data)[[p]],
                                                                   transcripts = transcripts(BANDITS_data)[[p]],
                                                                   mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                                   R = 2*R, burn_in = 2*burn_in, N_1 = N_mg[1], N_2 = N_mg[2],
                                                                   theshold_pval = theshold_pval)
                                 if( length(res[[1]]) > 1 ){
                                   res[[1]][[1]] = res[[1]][[1]][,c(1, seq.int(4,7, by = 1))]
                                   # only keep the first value (the gene p.val), and 4:5 (mean log-prec) and 6:7 (sd log-prec)
                                   
                                   sel_tmp = c(1, 1 + which(sel_groups), 1 + length(sel_groups) + which(sel_groups))
                                   res_tmp = matrix(NaN, nrow = nrow(res[[1]][[1]]), ncol = 1 + 2 * length(sel_groups)) # rows = transcripts; cols = groups
                                   res_tmp[ , sel_tmp] = res[[1]][[1]]
                                   rownames(res_tmp) = rownames(res[[1]][[1]])
                                   res[[1]][[1]] = res_tmp
                                   
                                   if(!is.null(res[[1]][[3]])){ # if trans result is not null (i.e., if the gene was analyzed):
                                     # Store the transcripts modes in a matrix, leave empty groups which were not analyzed:
                                     res[[1]][[3]] = cbind( do.call(c, res[[1]][[3]]), do.call(c, res[[1]][[4]])) # gather together MEAN results from different groups
                                     res[[1]][[4]] = cbind( do.call(c, res[[1]][[5]]), do.call(c, res[[1]][[6]])) # gather together SD results from different groups
                                     
                                     # add a null column for the group which was not analyzed:
                                     trans_mode = matrix(NaN, nrow = nrow(res[[1]][[3]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                     trans_mode[, sel_groups ] = res[[1]][[3]]
                                     res[[1]][[3]] = trans_mode # rows = transcripts; cols = groups
                                     
                                     # add a null column for the group which was not analyzed:
                                     trans_sd = matrix(NaN, nrow = nrow(res[[1]][[4]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                     trans_sd[, sel_groups ] = res[[1]][[4]]
                                     res[[1]][[4]] = trans_sd # rows = transcripts; cols = groups
                                   }
                                 }
                                 # double iterations for Together genes: they typically require more iter than Unique genes (more complex posterior space to explore).
                               }else{ # for the following elements I run the Unique function
                                 if( length(effLen(BANDITS_data)[[p]]) > 1 ){ # run test only if there are at least 2 transcripts per gene.
                                   res = wald_DTU_test_FULL( f = f,
                                                             l = effLen(BANDITS_data)[[p]], 
                                                             exon_id = classes(BANDITS_data)[[p]],
                                                             mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                             R = R, burn_in = burn_in, N_1 = N_mg[1], N_2 = N_mg[2],
                                                             theshold_pval = theshold_pval)
                                   if(length(res[[1]]) > 1){
                                     res[[1]][[1]] = res[[1]][[1]][c(1, seq.int(4,7, by = 1))]
                                     # only keep the first value (the gene p.val), and 4:5 (mean log-prec) and 6:7 (sd log-prec)
                                     
                                     # sort p.value (pos 1) and post mean and sd of log-precision:
                                     sel_tmp = c(1, 1 + which(sel_groups), 1 + length(sel_groups) + which(sel_groups))
                                     res_tmp = matrix(NaN, nrow = 1, ncol = 1 + 2 * length(sel_groups)) # rows = transcripts; cols = groups
                                     res_tmp[ , sel_tmp] = res[[1]][[1]]
                                     res[[1]][[1]] = res_tmp
                                     rownames(res[[1]][[1]]) = genes(BANDITS_data)[[p]]
                                     
                                     names(res[[1]][[2]]) = transcripts(BANDITS_data)[[p]]
                                     
                                     # Store the transcripts modes in a matrix, leave empty groups which were not analyzed:
                                     trans_mode = matrix(NaN, nrow = length(transcripts(BANDITS_data)[[p]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                     trans_mode[ , which(sel_groups)[1]] = res[[1]][[3]]
                                     trans_mode[ , which(sel_groups)[2]] = res[[1]][[4]]
                                     res[[1]][[3]] = trans_mode # rows = transcripts; cols = groups
                                     
                                     trans_sd = matrix(NaN, nrow = length(transcripts(BANDITS_data)[[p]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                     trans_sd[ , which(sel_groups)[1]] = res[[1]][[5]]
                                     trans_sd[ , which(sel_groups)[2]] = res[[1]][[6]]
                                     res[[1]][[4]] = trans_sd # rows = transcripts; cols = groups
                                   } # only if the mcmc has a return value. 
                                 }
                               }
                               
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
                       p_values[,-1],
                       row.names = NULL)
  
  names(gene_DF)[ seq.int(4, 3 + length(group_levels)) ] = paste("Mean log-prec", group_levels)
  names(gene_DF)[ seq.int(4 + length(group_levels), 3 + 2*length(group_levels)) ] = paste("SD log-prec", group_levels)
  
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
  
  mode_groups = lapply(p_values_ALL, function(x){
    y = x[[1]]
    if(length(y) > 1){
      y = y[[3]]
      if(is.list(y)){ # for Together genes, I first need to wrap transcript results from a list of matrices into a unique matrix!
        y = do.call(rbind, y)
      }else{
        return(y)
      }
    }else{
      return(NULL)
    }} )
  
  sd_groups = lapply(p_values_ALL, function(x){
    y = x[[1]]
    if(length(y) > 1){
      y = y[[4]]
      if(is.list(y)){ # for Together genes, I first need to wrap transcript results from a list of matrices into a unique matrix!
        y = do.call(rbind, y)
      }else{
        return(y)
      }
    }else{
      return(NULL)
    }} )
  
  mode_groups = do.call(rbind, mode_groups)
  sd_groups = do.call(rbind, sd_groups)
  
  p_values_tr = cbind(p_values_tr, mode_groups, sd_groups)
  
  cond = !vapply(p_values_tr[,1], is.null, FUN.VALUE = logical(1))
  p_values_tr = p_values_tr[cond,] # filter null results
  
  # COMPUTE ADJUSTED P.VALS:
  adj.p_values_tr = p.adjust(p_values_tr[,1],  method = "BH") # transcript-test
  
  # match GENE and TRANSCRIPT IDs.
  Gene_id     = as.character(gene_to_transcript[,1]); Tr_id = as.character(gene_to_transcript[,2]) 
  genes_in_tr = Gene_id[match(rownames(p_values_tr), Tr_id)]
  
  # Compute "conservative" transcript level test,
  # take min between gene and transcript level tests (CHECK Koen package):
  max_gene_tr_p.val     = apply( cbind(p_values_tr, p_values[match(genes_in_tr, names(p_values))]), 1, max, na.rm = TRUE)
  max_gene_tr_adj.p.val = apply( cbind(adj.p_values_tr, adj.p_values[match(genes_in_tr, names(adj.p_values))]), 1, max, na.rm = TRUE)
  
  tr_DF = data.frame(Gene_id = genes_in_tr, Transcript_id = rownames(p_values_tr), 
                     p.values = p_values_tr[,1], adj.p.values = adj.p_values_tr,
                     Max_Gene_Tr.p.val = max_gene_tr_p.val, Max_Gene_Tr.Adj.p.val = max_gene_tr_adj.p.val,
                     mean_sd = p_values_tr[,-1],
                     row.names = NULL)
  
  # re-name the group names according to the groups names:
  names(tr_DF)[ seq(7, 6 + length(group_levels)) ] = paste("Mean", group_levels)
  names(tr_DF)[ seq(7 + length(group_levels), 6 + 2*length(group_levels)) ] = paste("SD", group_levels)
  
  # sort p.values according to their GENE significance:
  tr_DF = tr_DF[ order(p_values[match(tr_DF$Gene_id, rownames(p_values))], tr_DF$p.values) ,]
  
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
  
  convergence = convergence[order(convergence[,2], decreasing = FALSE),]
  
  # remove gene ids that are not in "all_genes(BANDITS_data)"; i.e., remove gene ids for (Together) genes with 1 transcript only!
  convergence = convergence[ rownames(convergence) %in% all_genes(BANDITS_data), ]
  
  #########################################################################################################
  # Return results:
  #########################################################################################################
  message("Returning results")
  
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
