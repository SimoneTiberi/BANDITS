wald_DTU_test_Together_FULL = function(f, l, exon_id, N_1, N_2, R, burn_in, 
                                       mean_log_precision, sd_log_precision,
                                       genes, transcripts,
                                       theshold_pval = 0.1){
  K_tot = length(transcripts)
  gene_id = vapply(genes, function(x) names(transcripts) == x, FUN.VALUE = logical(K_tot))
  # sapply(genes, function(x) names(transcripts) == x)
  order = unlist(apply(gene_id, 2, which))
  
  # order transcripts such that the first K[1] transcripts refer to the 1st gene, and so on.
  l = l[order]
  exon_id = exon_id[order,]
  transcripts = transcripts[order]
  gene_id = gene_id[order,]
  
  # TRUE-FALSE matrix telling me what genes are associated to what transcripts!
  # gene_id[,i] refers to the i-th gene.
  
  K = vapply(genes, function(x) sum(names(transcripts) == x), FUN.VALUE = integer(1))
  # sapply(genes, function(x) sum(names(transcripts) == x))
  n_genes = length(K); # nr of genes in the group
  
  ### ### ### if K == 1, do not provide a p.value, set pi = 1 in the function above!!!
  if( all(K == 1) ){
    res_gene =  matrix(-1, ncol = 3, nrow  = length(K))
    rownames(res_gene) = genes
    res_transcript = NULL
    res = list(res_gene, res_transcript)
    return(list(p.vals = res) ) # if all genes have 1 transcript only, I cannot infer any DTU genes so I return all -1's
  }
  
  chain = MCMC_chain_Together_FULL(f = f, l = l, exon_id = exon_id, N_1 = N_1, N_2 = N_2, R = R, K = K, gene_id = gene_id, 
                                   burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  if(chain[[2]][1] == 0){ # IF the chain didn't converge (3 times), return NULL result:
    return( list(p.vals = NA, convergence = chain[[2]] ) )
  }
  
  pvals_res = pval_compute_Together_FULL( A = chain[[1]][[1]], B = chain[[1]][[2]],
                                          K = K, N = N_1+N_2, genes = genes, gene_id = gene_id, transcripts = transcripts)
  # update the p.val computation (similarly to the Unique genes!)
  #print(pvals_res)
  
  if( any(is.na(pvals_res[[1]][, 1])) == FALSE ){ # If any p.val is NA, I output a warning and redo the MCMC
    if( all( pvals_res[[1]][K>1, 1] > theshold_pval ) ){
      # if all genes tested (with at least 2 transcripts) have a p.val[55] > 0.1 I return the p.vals

      # compute the posterior mode and sd of the precision parameter
      tmp = matrix(-1, nrow = n_genes, ncol = 7)
      tmp[,c(1,2,3)] = pvals_res[[1]]
      tmp[,c(4,5)] = vapply(chain[[4]], function(x){ apply(x, 2, mean)}, FUN.VALUE = numeric(n_genes))
      tmp[,c(6,7)] = vapply(chain[[4]], function(x){ apply(x, 2, sd)}, FUN.VALUE = numeric(n_genes))
      rownames(tmp) = genes # gene id to the gene results
      
      pvals_res[[1]] = tmp
      
      return( list(p.vals = pvals_res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
    }
  }
  #  print(paste("1st:", pvals_res[,35]) )
  # If I didn't return the output yet it means either: 1) p.val is NA (never so far) 2) p.val < threshold (0.1 by default)/
  chain_2 = MCMC_chain_Together_FULL(f = f, l = l, exon_id = exon_id, N_1 = N_1, N_2 = N_2, R = R, K = K, gene_id = gene_id, 
                                     burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  # if chain_2 converged, I add it to the first one, otherwise I don't:
  if(chain_2[[2]][1] == 0){ # IF the second chain didn't converge (three times), return the result from the first one:
    # if all genes tested (with at least 2 transcripts) have a p.val[55] > 0.1 I return the p.vals

    # compute the posterior mode and sd of the precision parameter
    tmp = matrix(-1, nrow = n_genes, ncol = 7)
    tmp[,c(1,2,3)] = pvals_res[[1]]
    tmp[,c(4,5)] = vapply(chain[[4]], function(x){ apply(x, 2, mean)}, FUN.VALUE = numeric(n_genes))
    tmp[,c(6,7)] = vapply(chain[[4]], function(x){ apply(x, 2, sd)}, FUN.VALUE = numeric(n_genes))
    rownames(tmp) = genes # gene id to the gene results
    
    pvals_res[[1]] = tmp
    
    return( list(p.vals = pvals_res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
  }
  
  # merge the 2 chains together (after excluding burn-in):
  for(g in seq_along(K)){
    if(K[g] > 1){
      chain[[1]][[1]][[g]] = rbind( chain[[1]][[1]][[g]], chain_2[[1]][[1]][[g]])
      chain[[1]][[2]][[g]] = rbind( chain[[1]][[2]][[g]], chain_2[[1]][[2]][[g]])
    }
  }
  chain[[4]][[1]] = rbind(chain[[4]][[1]], chain_2[[4]][[1]])
  chain[[4]][[2]] = rbind(chain[[4]][[2]], chain_2[[4]][[2]])
  
  # I merge the two chains computed independently and return the pvals computed on the two chains merged together.
  pvals_res = pval_compute_Together_FULL( A = chain[[1]][[1]], B = chain[[1]][[2]],
                                          K = K, N = N_1+N_2, genes = genes, gene_id = gene_id, transcripts = transcripts)
  
  tmp = matrix(-1, nrow = n_genes, ncol = 7)
  tmp[,c(1,2,3)] = pvals_res[[1]]
  tmp[,c(4,5)] = vapply(chain[[4]], function(x){ apply(x, 2, mean)}, FUN.VALUE = numeric(n_genes))
  tmp[,c(6,7)] = vapply(chain[[4]], function(x){ apply(x, 2, sd)}, FUN.VALUE = numeric(n_genes))
  rownames(tmp) = genes # gene id to the gene results
  
  pvals_res[[1]] = tmp
  
  list(p.vals = pvals_res, convergence = chain[[2]]) # return the convergence result too (to check they are all converged with reasonable burn-in).
}

MCMC_chain_Together_FULL = function(f, l, exon_id, N_1, N_2, R, K, gene_id, burn_in, 
                                    mean_log_precision, sd_log_precision, FIRST_chain = 1){
  One_transcript = K ==1
  #  gene_id_Numeric = ifelse(gene_id, 1, 0) # apparently not used.
  # check what equiv classes are associated to a single transcript:
  one_transcript = colSums(exon_id) == 1
  
  n_genes = length(K); # nr of genes in the group
  pi_new_A = lapply(X = seq_len(n_genes), FUN = create_pi_new, K = K, N = N_1 ) 
  pi_new_B = lapply(X = seq_len(n_genes), FUN = create_pi_new, K = K, N = N_2 )
  
  # Matrix of posterior chains for alpha
  mcmc_alpha_A = mcmc_alpha_B = list()
  for(g in seq_len(n_genes)){
    if(One_transcript[g] == FALSE){
      mcmc_alpha_A[[g]] = matrix(NA, nrow = R+burn_in, ncol = K[g])
      mcmc_alpha_B[[g]] = matrix(NA, nrow = R+burn_in, ncol = K[g]) # hyper-parameters
    }else{
      mcmc_alpha_A[[g]] = mcmc_alpha_B[[g]]  = matrix(NA)
    }
  }
  
  alpha_new_A = alpha_new_B = alpha_new_A_original = alpha_new_B_original = list()
  chol_A = list(); chol_B = list();
  for(g in seq_len(n_genes)){
    if(One_transcript[g] == FALSE){
      chol_A[[g]] = matrix(NA, nrow = K[g], ncol = K[g])
      chol_B[[g]] = matrix(NA, nrow = K[g], ncol = K[g])
      if( mean_log_precision != 0){ # if I have a prior for the dispersion, I use it for the starting value.
        alpha_new_A[[g]] = rep( mean_log_precision - log(K[g]), K[g] ) # if not specified, I use a bigger prior for log_precision, more dispersion.
        alpha_new_B[[g]] = rep( mean_log_precision - log(K[g]), K[g] ) # if not specified, I use a bigger prior for log_precision, more dispersion.
      }else{
        alpha_new_A[[g]] = rep( log(10) - log(K[g]), K[g] ) # if not specified, I use a bigger prior for log_precision, more dispersion.
        alpha_new_B[[g]] = rep( log(10) - log(K[g]), K[g] ) # if not specified, I use a bigger prior for log_precision, more dispersion.
      }
    }else{
      chol_A[[g]] = matrix(NA)
      chol_B[[g]] = matrix(NA)
      alpha_new_A[[g]] = NA
      alpha_new_B[[g]] = NA
    }
  }
  
  res = .Call(`_BANDITS_Rcpp_Together`, R, burn_in, N_1, N_2,
              pi_new_A, pi_new_B,
              mcmc_alpha_A, mcmc_alpha_B,
              alpha_new_A, alpha_new_B,
              chol_A, chol_B,
              mean_log_precision, sd_log_precision,
              K, l, f, exon_id, 
              One_transcript, one_transcript)
  
  seq. = round( seq.int(1, R, length.out = 10^4 ) ) # thin if R > 10^4 (by construction R >= 10^4)
  convergence = my_heidel.diag(res[[3]][seq.], R = length(seq.), by. = length(seq.)/10, pvalue = 0.01)
  
  # output:
  # Stationarity test passed (1) or not (0);
  # start iteration (it'd be > burn_in);
  # p-value (for the Stationarity test).
  
  if(convergence[1] == 1){ # if it converged:
    if(convergence[2] > 1){ # remove burn-in estimated by heidel.diag (which is, AT MOST, half of the chain):
      for(g in seq_len(n_genes)){
        if(One_transcript[g] == FALSE){ # only for genes with >1 transcript
          res[[1]][[g]] = res[[1]][[g]][seq.,][-seq_len(convergence[2]-1),]
          res[[2]][[g]] = res[[2]][[g]][seq.,][-seq_len(convergence[2]-1),]
        }
      }
      res[[4]] = res[[4]][seq.,][-seq_len(convergence[2]-1),]
      res[[5]] = res[[5]][seq.,][-seq_len(convergence[2]-1),]
    }else{ # if convergence[2] == 1, seq. has already been defined above.
      if(R > 10^4){ # thin only if R > 10^4
        for(g in seq_len(n_genes)){
          if(One_transcript[g] == FALSE){ # only for genes with >1 transcript
            res[[1]][[g]] = res[[1]][[g]][seq.,]
            res[[2]][[g]] = res[[2]][[g]][seq.,]
          }
        }
        res[[4]] = res[[4]][seq.,]
        res[[5]] = res[[5]][seq.,]
      }
    }
  }else{ # IF not converged, RUN a second chain (once only):
    if(FIRST_chain < 3){ # if first or second chain re-run again:
      # message("the first chain did NOT converge, I run a second one:")
      return( MCMC_chain_Together_FULL(f = f, l = l, exon_id = exon_id, N_1 = N_1, N_2 = N_2, R = R, K = K, gene_id = gene_id, 
                                       burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision,
                                       FIRST_chain = FIRST_chain + 1) )
    }else{ # if I ran 3 chains already and none of them converged, return convergence failure message:
      return(list(NA, convergence, FIRST_chain))
    }
  }
  # thin results to return 10^4 iterations.
  # thin if R > 10^4 (to return 10^4 values).
  
  # save whether it's the first run or not (i.e. whether the convergence test failed).
  list( res[seq_len(2)], convergence, FIRST_chain, res[c(4,5)] ) 
}


pval_compute_Together_FULL = function(A, B, K, N, genes, gene_id, transcripts){
  res = matrix(-1, ncol = 3, nrow  = length(K))
  if(sum(K == 1) == length(K)){
    return( res) # if all genes have 1 transcript only, I cannot infer any DTU genes so I return all -1's
  }
  res_tr = mode_A = mode_B = sd_A = sd_B = list()
  # do usual testing procedure on each gene separately.
  
  for(g in seq_along(K)){
    if(K[g] > 1){ # if 1 transcripts per gene I keep the result equal to -1
      pval_pi_T = compute_pval_FULL( A = A[[g]], B = B[[g]], K = K[g], N = N)
      
      res[g,] = pval_pi_T[[1]]
      res_tr[[g]] = pval_pi_T[[2]] # give transcript names!!!
      mode_A[[g]] = pval_pi_T[[3]]
      mode_B[[g]] = pval_pi_T[[4]]
      sd_A[[g]] = pval_pi_T[[5]]
      sd_B[[g]] = pval_pi_T[[6]]
      names(res_tr[[g]])  = transcripts[gene_id[,g]]
      # MAKE SURE THE ORDER OF THE TRANSCRIPTS IS CORRECT!!!
    }
  }
  rownames(res) = genes # gene id to the gene results
  #  colnames(res_tr) = transcripts # transcript id to the transcript results
  list( res, res_tr, mode_A, mode_B, sd_A, sd_B)
}
