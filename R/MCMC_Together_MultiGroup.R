wald_DTU_test_MultiGroup_Together = function(f, l, exon_id, N, R, burn_in, 
                                             mean_log_precision, sd_log_precision,
                                             genes, transcripts, theshold_pval = 0.1){
  N_groups = length(N)
  n_genes = length(genes)
  
  K = vapply(genes, function(x) sum(names(transcripts) == x), FUN.VALUE = integer(1))
  # sapply(genes, function(x) sum(names(transcripts) == x))
  
  K_tot = sum(K)
  gene_id = vapply(genes, function(x) names(transcripts) == x, FUN.VALUE = logical(K_tot))
  # sapply(genes, function(x) names(transcripts) == x)
  
  order = unlist(apply(gene_id, 2, which))
  
  # order transcripts such that the first K[1] transcripts refer to the 1st gene, and so on.
  l = l[order]
  exon_id = exon_id[order,]
  transcripts = transcripts[order]
  gene_id = gene_id[order,]
  
  ### ### ### if K == 1, do not provide a p.value, set pi = 1 in the function above!!!
  if( all(K == 1) ){
    return( NULL ) # if all genes have 1 transcript only, I cannot infer any DTU genes so I return all -1's
  }
  
  chain = MCMC_chain_MultiGroup_Together(f = f, l = l, exon_id = exon_id, N = N, N_groups = N_groups,
                                         n_genes = n_genes, R = R, K = K, gene_id = gene_id,
                                         burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  if(chain[[2]][1] == 0){ # IF the chain didn't converge (3 times), return NULL result:
    return( list(p.vals = NA, convergence = chain[[2]]) )
  }
  
  pvals_res = pval_compute_Together_MultiGroup(chain[[1]], K = K, n_genes = n_genes, genes = genes, N_groups = N_groups, gene_id = gene_id, transcripts = transcripts)
  
  if( any(is.na(pvals_res[[1]])) == FALSE){
    if( all( pvals_res[[1]][K>1] > theshold_pval ) ){
      
      tmp = matrix(-1, nrow = n_genes, ncol = 1 + 2 * N_groups)
      tmp[,1] = pvals_res[[1]]
      tmp[,seq.int(2, N_groups + 1, by = 1)]                = vapply(chain[[3]], function(x){ apply(x, 2, mean)}, FUN.VALUE = numeric(N_groups))
      tmp[,seq.int(N_groups + 2, 2 * N_groups + 1, by = 1)] = vapply(chain[[3]], function(x){ apply(x, 2, sd)},   FUN.VALUE = numeric(N_groups))
      rownames(tmp) = genes # gene id to the gene results
      
      pvals_res[[1]] = tmp
      
      # if all genes tested (with at least 2 transcripts) have a p.val[55] > 0.1 I return the p.vals
      return( list(p.vals = pvals_res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
    }
  }
  #  message(paste("1st:", pvals_res[,35]) )
  # If I didn't return the output yet it means either: 1) p.val is NA (never so far) 2) p.val < threshold (0.1 by default)/
  chain_2 = MCMC_chain_MultiGroup_Together(f = f, l = l, exon_id = exon_id, N = N, N_groups = N_groups, 
                                           n_genes = n_genes, R = R, K = K, gene_id = gene_id,
                                           burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  # if chain_2 converged, I add it to the first one, otherwise I don't:
  if(chain_2[[2]][1] == 0){ # IF the second chain didn't converge (three times), return the result from the first one:
    tmp = matrix(-1, nrow = n_genes, ncol = 1 + 2 * N_groups)
    tmp[,1] = pvals_res[[1]]
    tmp[,seq.int(2, N_groups + 1, by = 1)]                = vapply(chain[[3]], function(x){ apply(x, 2, mean)}, FUN.VALUE = numeric(N_groups))
    tmp[,seq.int(N_groups + 2, 2 * N_groups + 1, by = 1)] = vapply(chain[[3]], function(x){ apply(x, 2, sd)}, FUN.VALUE = numeric(N_groups))
    rownames(tmp) = genes # gene id to the gene results
    
    pvals_res[[1]] = tmp
    
    return( list(p.vals = pvals_res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
  }
  
  for(g in seq_len(n_genes) ){
    if(K[g] > 1){
      for(n in seq_len(N_groups) ){
        chain[[1]][[g]][[n]] = rbind( chain[[1]][[g]][[n]], chain_2[[1]][[g]][[n]] )
      }
    }
  }
  
  # I merge the two chains computed independently and return the pvals computed on the two chains merged together.
  pvals_res = pval_compute_Together_MultiGroup(chain[[1]], K = K, n_genes = n_genes, genes = genes, N_groups = N_groups, gene_id = gene_id, transcripts = transcripts)
  
  tmp = matrix(-1, nrow = n_genes, ncol = 1 + 2 * N_groups)
  tmp[,1] = pvals_res[[1]]
  tmp[,seq.int(2, N_groups + 1, by = 1)]                = vapply(chain[[3]], function(x){ apply(x, 2, mean)}, FUN.VALUE = numeric(N_groups))
  tmp[,seq.int(N_groups + 2, 2 * N_groups + 1, by = 1)] = vapply(chain[[3]], function(x){ apply(x, 2, sd)}, FUN.VALUE = numeric(N_groups))
  rownames(tmp) = genes # gene id to the gene results
  
  pvals_res[[1]] = tmp
  
  list(p.vals = pvals_res, convergence = chain[[2]]) # return the convergence result too (to check they are all converged with reasonable burn-in).
}

MCMC_chain_MultiGroup_Together = function(f, l, exon_id, N, N_groups, n_genes, R, K, gene_id,burn_in, mean_log_precision, sd_log_precision,
                                          FIRST_chain = 1){
  One_transcript = K ==1
  J = ncol(exon_id);
  
  N_tot = sum(N)
  cumulative = c(0,cumsum(N))
  splits = list()
  splits = lapply(seq_len( length(cumulative) - 1), function(i){
    {cumulative[i]+1}:cumulative[i+1]
  })
  
  # define object containing the data:
  f_list = list()
  # starting values for the alpha parameters, sampled in the log-space:
  alpha_new = list() # matrix(NA, nrow = K, ncol = N_groups)
  # pi:
  pi_new = list()
  # mcmc matrices:
  mcmc_alpha = list()
  # chol matrices:
  chol_mat = list()
  # tot:
  TOT_Y_new = list()
  # log-prec:
  precision = list()
  
  for(g in seq_len(n_genes) ){
    if(!One_transcript[g]){
      precision[[g]] = matrix(NA, nrow = R, ncol = N_groups)
    }else{
      precision[[g]] = matrix(NaN, nrow = 1, ncol = N_groups)
    }
    mcmc_alpha[[g]] = list()
    for(n in seq_len(N_groups) ){
      if(!One_transcript[g]){ # only store the matrix IF more than 1 transcript in the group
        mcmc_alpha[[g]][[n]] = matrix(NA, nrow = R + burn_in, ncol = K[g]) # hyper-parameters of the DM
      }else{
        mcmc_alpha[[g]][[n]] = matrix(NaN, nrow = 1, ncol = 1)
      }
    }
  }
  
  # loop once on all objects:
  for(n in seq_len(N_groups) ){
    # define object containing the data:
    f_list[[n]] = as.matrix(f[,splits[[n]]])
    TOT_Y_new[[n]] = matrix(1, nrow = N[n], ncol = n_genes)
    alpha_new[[n]] = list()
    pi_new[[n]] = list()
    chol_mat[[n]] = list()
    
    for(g in seq_len(n_genes) ){
      if(!One_transcript[g]){
        # starting values for alpha_new (log space):
        if( mean_log_precision != 0){
          alpha_new[[n]][[g]] = rep( mean_log_precision - log(K[g]), K[g]) # delta_1, ..., delta_{K-1}, delta_{K}
        }else{
          alpha_new[[n]][[g]] = rep( log(10) - log(K[g]), K[g]) # delta_1, ..., delta_{K-1}, delta_{K}
        }
        
        # pi's:
        pi_new[[n]][[g]] = matrix( 1/K[g], nrow = N[n], ncol = K[g])
        
        # mcmc matrices:
        #mcmc_alpha[[n]][[g]] = matrix(NA, nrow = R + burn_in, ncol = K[g]) # hyper-parameters of the DM
        
        #chol matrices:
        chol_mat[[n]][[g]] = matrix(0, nrow = K[g], ncol = K[g])
      }else{
        alpha_new[[n]][[g]] = NaN
        
        # pi's:
        pi_new[[n]][[g]] = matrix( NaN )
        
        #chol matrices:
        chol_mat[[n]][[g]] = matrix(NaN)
      }
    }
  }
  
  one_transcript = colSums(exon_id) == 1
  N = as.integer(N)
  
  # Run the MCMC fully in Rcpp:
  res = .Call(`_BANDITS_Rcpp_FULL_Together_Multigroup`, R + burn_in, burn_in, N, N_groups,
              mean_log_precision, sd_log_precision, 
              pi_new, mcmc_alpha, alpha_new, chol_mat, TOT_Y_new, precision,
              K, l, f_list, exon_id, One_transcript, one_transcript)
  
  # Compute the convergence diagnostic:
  seq. = round( seq.int(1, R, length.out = 10^4 ) ) # thin if R > 10^4 (by construction R >= 10^4)
  convergence = my_heidel.diag(res[[2]][seq.], R = length(seq.), by. = length(seq.)/10, pvalue = 0.01)
  
  # output:
  # Stationarity test passed (1) or not (0);
  # start iteration (it'd be > burn_in);
  # p-value (for the Stationarity test).
  
  if(convergence[1] == 1){ # if it converged:
    if(convergence[2] > 1){ # remove burn-in estimated by heidel.diag (which is, AT MOST, half of the chain):
      for(g in seq_along(K) ){
        if(K[g] > 1){
          res[[3]][[g]] = res[[3]][[g]][seq.,][-{seq_len(convergence[2]-1)},]
          for(n in seq_len(N_groups) ){
            res[[1]][[g]][[n]] = res[[1]][[g]][[n]][seq.,][-{seq_len(convergence[2]-1)},]
          }
        }
      }
    }else{ # if convergence[2] == 1, seq. has altready been defined above.
      if(R > 10^4){ # thin if R > 10^4
        for(g in seq_along(K) ){
          if(K[g] > 1){
            res[[3]][[g]] = res[[3]][[g]][seq.,]
            for(n in seq_len(N_groups) ){
              res[[1]][[g]][[n]] = res[[1]][[g]][[n]][seq.,]
            }
          }
        }
      }
    }
  }else{ # IF not converged, RUN a second chain (once only):
    if(FIRST_chain < 3){ # if first or second chain re-run again:
      #message("the first chain did NOT converge, I run a second one:")
      return( MCMC_chain_MultiGroup_Together(f, l, exon_id, N, N_groups, n_genes, R, K, gene_id,burn_in, mean_log_precision, sd_log_precision, FIRST_chain + 1) )
    }else{ # if I ran 3 chains already and none of them converged, return convergence failure message:
      return(list(NaN, convergence, FIRST_chain))
    }
  }
  # thin results to return 10^4 iterations.
  # thin if R > 10^4 (to return 10^4 values).
  
  list( res[[1]], convergence, res[[3]] )  # I return the list of MCMC chains, excluding the burn-in, and the convergence output
}

pval_compute_Together_MultiGroup = function(mcmc, K, n_genes, genes, N_groups, gene_id, transcripts){
  res = rep(-1, n_genes)
  res_tr = mode_groups = sd_groups = list()
  
  # do usual testing procedure on each gene separately:
  for(g in seq_along(K)){
    if(K[g] > 1){ # if 1 transcripts per gene I keep the result equal to -1
      pval_pi_T = pval_compute_MultiGroup(mcmc[[g]], K = K[g], N_groups = N_groups)
      
      res[g] = pval_pi_T[[1]]
      res_tr[[g]] = pval_pi_T[[2]]
      mode_groups[[g]] = pval_pi_T[[3]]
      sd_groups[[g]] = pval_pi_T[[4]]
      names(res_tr[[g]])  = transcripts[gene_id[,g]]
    }
  }
  names(res) = genes # gene id to the gene results
  list( res, res_tr, mode_groups, sd_groups)
}
