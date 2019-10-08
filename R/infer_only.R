infer_only = function(f, l, exon_id, N, R, burn_in, mean_log_precision = 0, sd_log_precision = 10){
  K = nrow(exon_id) # nr of transcripts
  
  chain = MCMC_infer_only(f = f, l = l, exon_id = exon_id, N = N, R = R, K = K, 
                          burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  if(chain[[2]][1] == 0){ # IF the first chain didn't converge (3 times), return NULL result:
    return( list(NULL, convergence = chain[[2]]) )
  }
  
  res = res_compute_infer_only( chain = chain, K = K)
  
  return( list(res = res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
}

MCMC_infer_only = function(f, l, exon_id, N, R, K, burn_in, mean_log_precision, sd_log_precision,
                           FIRST_chain = 1){
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
  
  i = 1 # group index (always 1)
  
  # define object containing the data:
  f_list[[i]] = as.matrix(f)
  
  # starting values for alpha_new (log space):
  if( mean_log_precision != 0){
    alpha_new[[i]] = rep( mean_log_precision - log(K), K) # delta_1, ..., delta_{K-1}, delta_{K}
  }else{
    alpha_new[[i]] = rep( log(10) - log(K), K) # delta_1, ..., delta_{K-1}, delta_{K}
  }
  
  # pi's:
  pi_new[[i]] = matrix( 1/K, nrow = N[i], ncol = K)
  
  # mcmc matrices:
  mcmc_alpha[[i]] = matrix(NA, nrow = R + burn_in, ncol = K) # hyper-parameters of the DM
  
  #chol matrices:
  chol_mat[[i]] = matrix(0, nrow = K, ncol = K)
  
  one_transcript = colSums(exon_id) == 1
  N = as.integer(N)
  
  # Run the MCMC fully in Rcpp:
  res = .Call(`_BANDITS_Rcpp_FULL_Unique_Multigroup`, K, R + burn_in, burn_in, N, 1, # 1 indicates N_groups
              mean_log_precision, sd_log_precision, pi_new, mcmc_alpha,
              alpha_new, chol_mat, l, f_list, exon_id, one_transcript)
  
  # Compute the convergence diagnostic:
  seq. = round( seq.int(1, R, length.out = 10^4 ) ) # thin if R > 10^4 (by construction R >= 10^4)
  convergence = my_heidel.diag(res[[2]][seq.], R = length(seq.), by. = length(seq.)/10, pvalue = 0.01)
  
  # output:
  # Stationarity test passed (1) or not (0);
  # start iteration (it'd be > burn_in);
  # p-value (for the Stationarity test).
  
  n = 1 # group index (always 1)
  
  if(convergence[1] == 1){ # if it converged:
    if(convergence[2] > 1){ # remove burn-in estimated by heidel.diag (which is, AT MOST, half of the chain):
      res[[1]][[n]] = res[[1]][[n]][seq.,][-{seq_len(convergence[2]-1)},]
      res[[3]] = c(res[[3]])[seq.][-{seq_len(convergence[2]-1)}]
    }else{ # if convergence[2] == 1, seq. has altready been defined above.
      if(R > 10^4){ # thin if R > 10^4
        res[[1]][[n]] = res[[1]][[n]][seq.,]
        res[[3]] = c(res[[3]])[seq.]
      }
    }
  }else{ # IF not converged, RUN a second chain (once only):
    if(FIRST_chain < 3){ # if first or second chain re-run again:
      # message("the first chain did NOT converge, I run a second one:")
      return( MCMC_infer_only(f, l, exon_id, N, R, K, burn_in, mean_log_precision, sd_log_precision, FIRST_chain = FIRST_chain + 1) )
    }else{ # if I ran 3 chains already and none of them converged, return convergence failure message:
      return(list(NaN, convergence, FIRST_chain))
    }
  }
  # thin results to return 10^4 iterations.
  # thin if R > 10^4 (to return 10^4 values).
  
  list( res[[1]], convergence, res[[3]] )  # I return the list of MCMC chains, excluding the burn-in, and the convergence output
}

res_compute_infer_only = function(chain, K){
  mean_prec = mean(chain[[3]])
  sd_prec   = sd(chain[[3]])
  
  mode_groups = vapply(chain[[1]], function(x) colSums(x), FUN.VALUE = numeric(K) ) 
  mode_groups = apply( mode_groups, 2, function(x) x/sum(x))
  sd_groups = vapply(chain[[1]], function(x) sqrt(diag(var(x))), FUN.VALUE = numeric(K) ) 
  
  list( c(mean_prec, sd_prec), mode_groups, sd_groups )
}
