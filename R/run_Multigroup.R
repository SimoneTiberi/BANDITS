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
      cl <- makeCluster(n_cores);
    })
    registerDoParallel(cl, n_cores);
  }
  
  N_tot = sum(N)
  cumulative = c(0,cumsum(N))
  splits = list()
  for(i in 2:length(cumulative) ){
    splits[[i-1]] = {cumulative[i-1]+1}:cumulative[i]
  }
  
  message("Starting the MCMC")
  
  p_values_ALL = foreach(p = final_order,
                         .packages=c("BANDITS"),
                         .errorhandling = "stop") %dorng%{
                           # ORDER counts to respect the ordering in "groups" via 'ord_samples'
                           f = BANDITS_data@counts[[p]][,ord_samples]
                           sel_samples = colSums(f) > 0.5
                           
                           # select samples with BANDITS_data only:
                           f = f[, sel_samples ]
                           N_mg = vapply(splits, function(id) sum(sel_samples[id]), FUN.VALUE = integer(1))
                           # sapply(splits, function(id) sum(sel_samples[id]))
                           sel_groups = N_mg > 0.5 # sel groups to keep (if at least 1 sample!):
                           N_mg = N_mg[ sel_groups ] # remove groups with 0 counts.
                           
                           # only run the MCMC if there are at least 2 groups with at least 1 sample with counts
                           cond_f = length(N_mg) > 1.5
                           
                           # initialize results:
                           res = NULL
                           
                           if(cond_f){ # there are at least 2 groups
                             if(length(N_mg) > 2.5){ # if there are > 2 samples, run the Multi-group comparison
                               
                               if( !BANDITS_data@uniqueId[[p]] ){ # for the first (Together) elements I run the together function
                                 res = wald_DTU_test_MultiGroup_Together(f = f,
                                                                         l = BANDITS_data@effLen[[p]], 
                                                                         exon_id = BANDITS_data@classes[[p]],
                                                                         genes = BANDITS_data@genes[[p]],
                                                                         transcripts = BANDITS_data@transcripts[[p]],
                                                                         mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                                         R = 2*R, burn_in = 2*burn_in, N = N_mg,
                                                                         theshold_pval = theshold_pval)
                               }else{ # for the following elements I run the Unique function
                                 if( length(BANDITS_data@effLen[[p]]) > 1 ){ # run test only if there are at least 2 transcripts per gene.
                                   res = wald_DTU_test_MultiGroup(f = as.data.frame(f),
                                                                  l = BANDITS_data@effLen[[p]], 
                                                                  exon_id = BANDITS_data@classes[[p]],
                                                                  mean_log_precision = mean_log_precision, 
                                                                  sd_log_precision = sd_log_precision,
                                                                  R = R, burn_in = burn_in, N = N_mg,
                                                                  theshold_pval = theshold_pval)
                                   if(length(res[[1]]) > 1){
                                     names(res[[1]][[2]]) = BANDITS_data@transcripts[[p]]
                                     names(res[[1]][[1]]) = BANDITS_data@genes[[p]]
                                   }
                                 }
                               }
                               
                             }else{  # if there are only 2 samples left (with counts), run the classical 2-group comparison:
                               f = as.matrix(f) # the 2-group comparison requires f to be a matrix
                               
                               if( !BANDITS_data@uniqueId[[p]] ){ # for the first (Together) elements I run the together function
                                 res = wald_DTU_test_Together_FULL(f = f,
                                                                   l = BANDITS_data@effLen[[p]], 
                                                                   exon_id = BANDITS_data@classes[[p]],
                                                                   genes = BANDITS_data@genes[[p]],
                                                                   transcripts = BANDITS_data@transcripts[[p]],
                                                                   mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                                   R = 2*R, burn_in = 2*burn_in, N_1 = N_mg[1], N_2 = N_mg[2],
                                                                   theshold_pval = theshold_pval)
                                 if(length(res[[1]]) > 1){
                                   res[[1]][[1]] = res[[1]][[1]][,1] # only keep the first value (the gene p.val)
                                   
                                   for(tr_id in seq_len(length(res[[1]][[3]]) )){ # loop over the transcripts to store results from a list into a matrix structure:
                                     if(!is.null(res[[1]][[3]][[tr_id]])){ # if trans result is not null (i.e., if the gene was analyzed):
                                       # Store the transcripts modes in a matrix, leave empty groups which were not analyzed:
                                       trans_mode = matrix(NaN, nrow = length(res[[1]][[3]][[tr_id]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                       trans_mode[ , which(sel_groups)[1]] = res[[1]][[3]][[tr_id]]
                                       trans_mode[ , which(sel_groups)[2]] = res[[1]][[4]][[tr_id]]
                                       res[[1]][[3]][[tr_id]] = trans_mod # rows = transcripts; cols = groups
                                       
                                       trans_sd = matrix(NaN, nrow = length(res[[1]][[5]][[tr_id]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                       trans_sd[ , which(sel_groups)[1]] = res[[1]][[5]][[tr_id]]
                                       trans_sd[ , which(sel_groups)[2]] = res[[1]][[6]][[tr_id]]
                                       res[[1]][[4]][[tr_id]] = trans_sd # rows = transcripts; cols = groups
                                     }
                                   }
                                 }                                 
                                 # double iterations for Together genes: they typically require more iter than Unique genes (more complex posterior space to explore).
                               }else{ # for the following elements I run the Unique function
                                 if( length(BANDITS_data@effLen[[p]]) > 1 ){ # run test only if there are at least 2 transcripts per gene.
                                   res = wald_DTU_test_FULL( f = f,
                                                             l = BANDITS_data@effLen[[p]], 
                                                             exon_id = BANDITS_data@classes[[p]],
                                                             mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                                             R = R, burn_in = burn_in, N_1 = N_mg[1], N_2 = N_mg[2],
                                                             theshold_pval = theshold_pval)
                                   if(length(res[[1]]) > 1){
                                     res[[1]][[1]] = res[[1]][[1]][1] # only keep the first value (the gene p.val)
                                     names(res[[1]][[1]]) = BANDITS_data@genes[[p]]
                                     names(res[[1]][[2]]) = BANDITS_data@transcripts[[p]]
                                     
                                     # Store the transcripts modes in a matrix, leave empty groups which were not analyzed:
                                     trans_mode = matrix(NaN, nrow = length(BANDITS_data@transcripts[[p]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
                                     trans_mode[ , which(sel_groups)[1]] = res[[1]][[3]]
                                     trans_mode[ , which(sel_groups)[2]] = res[[1]][[4]]
                                     res[[1]][[3]] = trans_mode # rows = transcripts; cols = groups
                                     
                                     trans_sd = matrix(NaN, nrow = length(BANDITS_data@transcripts[[p]]), ncol = length(sel_groups)) # rows = transcripts; cols = groups
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
  p_values = do.call(c, p_values)
  
  gene_names = names(p_values)
  
  suppressWarnings({ p_values = as.numeric(p_values) })
  names(p_values) = gene_names
  
  SEL_ALL  = !is.na(p_values)  # remove NA's (genes not converged).
  p_values = p_values[SEL_ALL]
  
  SEL_ALL  = p_values != -1  # remove -1's (genes not analyzed).
  p_values = p_values[SEL_ALL]
  
  # sort p.values according to their significance.
  p_values = p_values[ order(p_values) ]
  
  # COMPUTE ADJUSTED P.VALS:
  adj.p_values    = p.adjust(p_values, method = "BH") # gene-test
  
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
  max_gene_tr_p.val     = apply( cbind(p_values_tr, p_values[match(genes_in_tr, names(p_values))]), 1, max)
  max_gene_tr_adj.p.val = apply( cbind(adj.p_values_tr, adj.p_values[match(genes_in_tr, names(adj.p_values))]), 1, max)
  
  tr_DF = data.frame(Gene_id = genes_in_tr, Transcript_id = rownames(p_values_tr), 
                     p.values = p_values_tr[,1], adj.p.values = adj.p_values_tr,
                     Max_Gene_Tr.p.val = max_gene_tr_p.val, Max_Gene_Tr.Adj.p.val = max_gene_tr_adj.p.val,
                     mean_sd = p_values_tr[,-1],
                     row.names = NULL)
  
  # re-name the group names according to the groups names:
  names(tr_DF)[ seq(7, 6 + length(group_levels)) ] = paste("Mean", group_levels)
  names(tr_DF)[ seq(7 + length(group_levels), 6 + 2*length(group_levels)) ] = paste("sd", group_levels)
  
  # sort p.values according to their GENE significance:
  tr_DF = tr_DF[ order(p_values[match(tr_DF$Gene_id, rownames(p_values))], tr_DF$p.values) ,]
  
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
  
  convergence = convergence[order(convergence[,2], decreasing = FALSE),]
  
  # remove gene ids that are not in "BANDITS_data@all_genes"; i.e., remove gene ids for (Together) genes with 1 transcript only!
  convergence = convergence[ rownames(convergence) %in% BANDITS_data@all_genes, ]
  
  #########################################################################################################
  # Return results:
  #########################################################################################################
  message("Returning results")
  
  return(new("BANDITS_test",
             Gene_results = data.frame(Gene_id = names(p_values), 
                                       p.values = p_values, adj.p.values = adj.p_values,
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
