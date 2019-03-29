wald_DTU_test_MultiGroup = function(f, l, exon_id, N, R, burn_in, mean_log_precision = 0, sd_log_precision = 10, theshold_pval = 0.1){
  K = nrow(exon_id) # nr of transcripts
  N_groups = length(N)
  
  if(N_groups < 2){
    return(list(NULL, NULL))
  }
  
  chain = MCMC_chain_MultiGroup(f = f, l = l, exon_id = exon_id, N = N, N_groups = N_groups, R = R, K = K, 
                                burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  if(chain[[2]][1] == 0){ # IF the first chain didn't converge (3 times), return NULL result:
    return( list(p.vals = NA, convergence = chain[[2]]) )
  }
  
  pvals_res = pval_compute_MultiGroup( mcmc = chain[[1]], K = K, N_groups = N_groups)
  
  if(is.na(pvals_res[[1]][1]) == FALSE){
    if( pvals_res[[1]][1] > theshold_pval ){ # if p.val[55] > 0.1 I return the p.vals
      return( list(p.vals = pvals_res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
    }
  }
  
  # If I didn't return the output yet it means either: 1) p.val is NA (never so far) 2) p.val < threshold (0.1 by default)/
  chain_2 = MCMC_chain_MultiGroup(f = f, l = l, exon_id = exon_id, N = N, N_groups = N_groups, R = R, K = K, 
                                burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  if(chain_2[[2]][1] == 0){ # IF the second chain didn't converge (3 times), return the result from the first one:
    return( list(p.vals = pvals_res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
  }
  
  for(g in seq_len(N_groups) ){
    chain[[1]][[g]] = rbind( chain[[1]][[g]], chain_2[[1]][[g]])
  }
  
  # I merge the two chains computed independently and return the pvals computed on the two chains merged together.
  pvals_res = pval_compute_MultiGroup( mcmc = chain[[1]], K = K, N_groups = N_groups)
  
  return( list(p.vals = pvals_res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
}

MCMC_chain_MultiGroup = function(f, l, exon_id, N, N_groups, R, K, burn_in, mean_log_precision, sd_log_precision,
                                 FIRST_chain = 1){
  J = ncol(exon_id);

  N_tot = sum(N)
  cumulative = c(0,cumsum(N))
  splits = list()
  for(i in seq(2, length(cumulative), by = 1) ){
    splits[[i-1]] = {cumulative[i-1]+1}:cumulative[i]
  }
  
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
  
  # loop once on all objects:
  for(i in seq_len(N_groups) ){
    # define object containing the data:
    f_list[[i]] = as.matrix(f[,splits[[i]]])
    
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
  }
  
  one_transcript = colSums(exon_id) == 1
  N = as.integer(N)
  
  # Run the MCMC fully in Rcpp:
  res = .Call(`_BANDITS_Rcpp_FULL_Unique_Multigroup`, K, R + burn_in, burn_in, N, N_groups,
              mean_log_precision, sd_log_precision, pi_new, mcmc_alpha,
              alpha_new, chol_mat, l, f_list, exon_id, one_transcript)
  
  # Compute the convergence diagnostic:
  seq. = round( seq.int(1, R, length.out = 10^4 ) ) # thin if R > 10^4 (by construction R >= 10^4)
  convergence = my_heidel.diag(res[[2]][seq.], R = length(seq.), by. = length(seq.)/10, pvalue = 0.01)

  # output:
  # Stationarity test passed (1) or not (0);
  # start iteration (it'd be > burn_in);
  # p-value (for the Stationarity test).
  
  if(convergence[1] == 1){ # if it converged:
    if(convergence[2] > 1){ # remove burn-in estimated by heidel.diag (which is, AT MOST, half of the chain):
      for(n in seq_len(N_groups) ){
        res[[1]][[n]] = res[[1]][[n]][seq.,][-{seq_len(convergence[2]-1)},]
      }
    }else{ # if convergence[2] == 1, seq. has altready been defined above.
      if(R > 10^4){ # thin if R > 10^4
        for(n in seq_len(N_groups) ){
          res[[1]][[n]] = res[[1]][[n]][seq.,]
        }
      }
    }
  }else{ # IF not converged, RUN a second chain (once only):
    if(FIRST_chain < 3){ # if first or second chain re-run again:
      # message("the first chain did NOT converge, I run a second one:")
      return( MCMC_chain_MultiGroup(f, l, exon_id, N, N_groups, R, K, burn_in, mean_log_precision, sd_log_precision, FIRST_chain = FIRST_chain + 1) )
    }else{ # if I ran 3 chains already and none of them converged, return convergence failure message:
      return(list(NaN, convergence, FIRST_chain))
    }
  }
  # thin results to return 10^4 iterations.
  # thin if R > 10^4 (to return 10^4 values).
  
  list( res[[1]], convergence )  # I return the list of MCMC chains, excluding the burn-in, and the convergence output
}

pval_compute_MultiGroup = function(mcmc, K, N_groups){
  R = nrow(mcmc[[1]])
  mcmc = lapply(mcmc, function(X) X[sample.int(R, R),] )
  # Random sample to decrease the correlation between w samples!
  
  # this returns a matrix: mode_groups[,1] represents the proportions of transcript 1 in all N_groups.
  mode_groups = vapply(mcmc, function(x) apply(x, 2, sum), FUN.VALUE = numeric(K) ) 
  # sapply(mcmc, function(x) apply(x, 2, sum) ) # find.mode, adjust = 10 (mode) or sum (mean)
  mode_groups = apply( mode_groups, 2, function(x) x/sum(x))
  
  # need to remove 1 parameter to make sure I don't test it twice!
  p = (N_groups-1)*(K-1) # degrees of freedom for the Chisq.
  
  # gene level test:
  p_value = matrix(NA, nrow = N_groups, ncol = K)
  for(g in seq_len(N_groups) ){ # baseline group to compare against
    for(k in seq_len(K) ){ # transcript to be removed
      A = B = c()
      for(g_2 in {seq_len(N_groups)}[-g]){ # baseline group to compare against
        A = cbind(A, mcmc[[g_2]][,-k])
        B = cbind(B, mcmc[[g]][,-k] ) # B is repeated identically
      }
      CV     = cov(A-B)
      mode   = apply(A-B, 2, find.mode, adjust = 10)  

      # Normal (classical Wald test)
      stat = t(mode) %*% ginv(CV, tol = 0) %*% mode
      p_value[g,k] = 1-pchisq(stat, df = p)
    }
  }
  
  if(K == 2){ # if there are only 
    return(list(mean(p_value), rep(mean(p_value), K), mode_groups))
  }
  
  # transcript level test
  trancript_res = matrix(NA, nrow = N_groups, ncol = K)
  for(g in seq_len(N_groups) ){ # baseline group to compare against
    for(k in seq_len(K) ){ # transcript to be TESTED!
      A = B = c()
      for(g_2 in {seq_len(N_groups)}[-g]){ # baseline group to compare against
        A = cbind(A, mcmc[[g_2]][,k])
        B = cbind(B, mcmc[[g]][,k] ) # B is repeated identically
      }
      CV     = cov(A-B)
      mode   = apply(A-B, 2, find.mode, adjust = 10)  
      
      # Normal (classical Wald test)
      stat = t(mode) %*% ginv(CV, tol = 0) %*% mode
      trancript_res[g,k] = 1-pchisq(stat, df = N_groups-1)
    }
  }
  
  list( mean(p_value), colMeans(trancript_res), mode_groups )
}
