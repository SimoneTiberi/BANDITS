infer_only_Together = function(f, l, exon_id, N, R, burn_in, 
                                             mean_log_precision, sd_log_precision,
                                             genes, transcripts){
  n_genes = length(genes)
  
  K = vapply(genes, function(x) sum(names(transcripts) == x), FUN.VALUE = integer(1))
  
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

  chain = MCMC_infer_only_Together(f = f, l = l, exon_id = exon_id, N = N,
                                   n_genes = n_genes, R = R, K = K, gene_id = gene_id,
                                   burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  if(chain[[2]][1] == 0){ # IF the chain didn't converge (3 times), return NULL result:
    return( list(res = NA, convergence = chain[[2]]) )
  }
  
  res = res_compute_infer_only_Together(chain, K = K, n_genes = n_genes, genes = genes, gene_id = gene_id, transcripts = transcripts)
  
  # if all genes tested (with at least 2 transcripts) have a p.val[55] > 0.1 I return the p.vals
  return( list(res = res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
}

MCMC_infer_only_Together = function(f, l, exon_id, N, n_genes, R, K, gene_id, burn_in, mean_log_precision, sd_log_precision,
                                    FIRST_chain = 1){
  One_transcript = K ==1
  J = ncol(exon_id);
  
  # define object containing the data:
  f_list = list()
  # starting values for the alpha parameters, sampled in the log-space:
  alpha_new = list()
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
      precision[[g]] = matrix(NA, nrow = R, ncol = 1)
    }else{
      precision[[g]] = matrix(NaN, nrow = 1, ncol = 1)
    }
    mcmc_alpha[[g]] = list()
    n = 1
    if(!One_transcript[g]){ # only store the matrix IF more than 1 transcript in the group
      mcmc_alpha[[g]][[n]] = matrix(NA, nrow = R + burn_in, ncol = K[g]) # hyper-parameters of the DM
    }else{
      mcmc_alpha[[g]][[n]] = matrix(NaN, nrow = 1, ncol = 1)
    }
  }
  
  n = 1 # group index (always 1)
  
  # define object containing the data:
  f_list[[n]] = as.matrix(f)
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
  
  one_transcript = colSums(exon_id) == 1
  N = as.integer(N)
  
  # Run the MCMC fully in Rcpp:
  res = .Call(`_BANDITS_Rcpp_FULL_Together_Multigroup`, R + burn_in, burn_in, N, 1, # 1 indicates N_groups
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
          res[[3]][[g]] = c(res[[3]][[g]])[seq.][-{seq_len(convergence[2]-1)}]
          res[[1]][[g]][[n]] = res[[1]][[g]][[n]][seq.,][-{seq_len(convergence[2]-1)},]
        }
      }
    }else{ # if convergence[2] == 1, seq. has altready been defined above.
      if(R > 10^4){ # thin if R > 10^4
        for(g in seq_along(K) ){
          if(K[g] > 1){
            res[[3]][[g]] = c(res[[3]][[g]])[seq.]
            res[[1]][[g]][[n]] = res[[1]][[g]][[n]][seq.,]
          }
        }
      }
    }
  }else{ # IF not converged, RUN a second chain (once only):
    if(FIRST_chain < 3){ # if first or second chain re-run again:
      #message("the first chain did NOT converge, I run a second one:")
      return( MCMC_infer_only_Together(f, l, exon_id, N, n_genes, R, K, gene_id, burn_in, mean_log_precision, sd_log_precision, FIRST_chain + 1) )
    }else{ # if I ran 3 chains already and none of them converged, return convergence failure message:
      return(list(NaN, convergence, FIRST_chain))
    }
  }
  # thin results to return 10^4 iterations.
  # thin if R > 10^4 (to return 10^4 values).
  
  list( res[[1]], convergence, res[[3]] )  # I return the list of MCMC chains, excluding the burn-in, and the convergence output
}

res_compute_infer_only_Together = function(mcmc, K, n_genes, genes, gene_id, transcripts){
  res = matrix(-1, nrow = n_genes, ncol = 2)

  mode_groups = sd_groups = list()
  
  # do usual testing procedure on each gene separately:
  for(g in seq_along(K)){
    if(K[g] > 1){ # if 1 transcripts per gene I keep the result equal to -1
      pval_pi_T = res_compute_infer_only( list(mcmc[[1]][[g]], NULL, mcmc[[3]][[g]]), K = K[g])
      
      res[g,] = pval_pi_T[[1]]
      mode_groups[[g]] = pval_pi_T[[2]]
      sd_groups[[g]] = pval_pi_T[[3]]
      names(mode_groups[[g]])  = transcripts[gene_id[,g]]
    }
  }
  rownames(res) = genes # gene id to the gene results
  
  list( res, mode_groups, sd_groups)
}
