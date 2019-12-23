wald_DTU_test_FULL = function(f, l, exon_id, N_1, N_2, R, burn_in, mean_log_precision = 0, sd_log_precision = 10, theshold_pval = 0.1){
  K = nrow(exon_id) # nr of transcripts
  chain = MCMC_chain_FULL(f = f, l = l, exon_id = exon_id, N_1 = N_1, N_2 = N_2, R = R, K = K, 
                          burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  if(chain[[2]][1] == 0){ # IF the first chain didn't converge (3 times), return NULL result:
    return( list(p.vals = NA, convergence = chain[[2]] ) )
  }
  
  pvals_res = compute_pval_FULL( A = chain[[1]][[1]], B = chain[[1]][[2]], K = K, N = N_1 + N_2)
  
  if(is.na(pvals_res[[1]][1]) == FALSE){ # If NA, I output a warning and redo the MCMC
    if( pvals_res[[1]][1] > theshold_pval ){ # if p.val > 0.1 I return the p.vals
      mean_prec = vapply(chain[[4]], mean, FUN.VALUE = numeric(1))
      sd_prec = vapply(chain[[4]], sd, FUN.VALUE = numeric(1))
      pvals_res[[1]] = c(pvals_res[[1]], mean_prec, sd_prec)
      return( list(p.vals = pvals_res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
    }
  }
  # If I didn't return the output yet it means either: 1) p.val is NA (never so far) 2) p.val < threshold (0.1 by default)/
  chain_2 =  MCMC_chain_FULL(f = f, l = l, exon_id = exon_id, N_1 = N_1, N_2 = N_2, R = R, K = K, 
                             burn_in = burn_in, mean_log_precision = mean_log_precision, sd_log_precision = sd_log_precision)
  
  # if chain_2 converged, I add it to the first one, otherwise I don't:
  if(chain_2[[2]][1] == 0){ # IF the second chain didn't converge (3 times), return the result (already computed) from the first one:
    mean_prec = vapply(chain[[4]], mean, FUN.VALUE = numeric(1))
    sd_prec = vapply(chain[[4]], sd, FUN.VALUE = numeric(1))
    pvals_res[[1]] = c(pvals_res[[1]], mean_prec, sd_prec)
    return( list(p.vals = pvals_res, convergence = chain[[2]]) ) # return the convergence result too (to check they are all converged with reasonable burn-in).
  }
  
  # I merge the two chains computed independently and return the pvals computed on the two chains merged together.
  pvals_res = compute_pval_FULL( A = rbind( chain[[1]][[1]], chain_2[[1]][[1]]), B = rbind( chain[[1]][[2]], chain_2[[1]][[2]]), K = K, N = N_1 + N_2)
  
  chain[[4]][[1]] = c(chain[[4]][[1]], chain_2[[4]][[1]])
  chain[[4]][[2]] = c(chain[[4]][[2]], chain_2[[4]][[2]])
  
  mean_prec = vapply(chain[[4]], mean, FUN.VALUE = numeric(1))
  sd_prec = vapply(chain[[4]], sd, FUN.VALUE = numeric(1))
  pvals_res[[1]] = c(pvals_res[[1]], mean_prec, sd_prec)
  
  list(p.vals = pvals_res, convergence = chain[[2]]) # return the convergence result too (to check they are all converged with reasonable burn-in).
}

# FIRST_chain needed to keep track of the times I run the MCMC (due to convergence issues):
MCMC_chain_FULL = function(f, l, exon_id, N_1, N_2, R, K, burn_in, mean_log_precision, sd_log_precision, 
                           FIRST_chain = 1){
  one_transcript = colSums(exon_id) == 1
  
  res = .Call(`_BANDITS_Rcpp_Unique`, K, R, burn_in, N_1, N_2, mean_log_precision, sd_log_precision,
              l, f, exon_id, one_transcript)
  
  seq. = round( seq.int(1, R, length.out = 10^4 ) ) # thin if R > 10^4 (by construction R >= 10^4)
  convergence = my_heidel.diag(res[[3]][seq.], R = length(seq.), by. = length(seq.)/10, pvalue = 0.01)
  
  # output:
  # Stationarity test passed (1) or not (0);
  # start iteration (it'd be > burn_in);
  # p-value (for the Stationarity test).
  
  if(convergence[1] == 1){ # if it converged:
    if(convergence[2] > 1){ # remove burn-in estimated by heidel.diag (which is, AT MOST, half of the chain):
      res[[1]] = res[[1]][seq.,][-seq_len(convergence[2]-1),]
      res[[2]] = res[[2]][seq.,][-seq_len(convergence[2]-1),]
      res[[4]] = res[[4]][seq.][-seq_len(convergence[2]-1)]
      res[[5]] = res[[5]][seq.][-seq_len(convergence[2]-1)]
    }else{ # if convergence[2] == 1, seq. has altready been defined above.
      if(R > 10^4){ # thin if R > 10^4
        res[[1]] = res[[1]][seq.,]
        res[[2]] = res[[2]][seq.,]
        res[[4]] = res[[4]][seq.]
        res[[5]] = res[[5]][seq.]
      }
    }
  }else{ # IF not converged, RUN a second chain (once only):
    if(FIRST_chain < 3){ # if first or second chain re-run again:
      # message("the first chain did NOT converge, I run a second one:")
      return( MCMC_chain_FULL(f, l, exon_id, N_1, N_2, R, K, burn_in, mean_log_precision, sd_log_precision, FIRST_chain = FIRST_chain + 1) )
    }else{ # if I ran 3 chains already and none of them converged, return convergence failure message:
      return(list(NaN, convergence, FIRST_chain))
    }
  }
  # thin results to return 10^4 iterations.
  # thin if R > 10^4 (to return 10^4 values).
  
  # save whether it's the first run or not (i.e. whether the convergence test failed).
  list( res[seq_len(2)], convergence, FIRST_chain, res[c(4,5)] )  # I return the list of MCMC chains, excluding the burn-in
}


# sometimes R is NULL!
# check why.
compute_pval_FULL = function(A, B, K, N){
  R = nrow(A)
  A = A[sample.int(R, R),] # n indicates the nr of elements of the chain (exluded burn-in)
  
  gamma = A - B
  
  CV     = cov(gamma) # cov is 20ish times faster than posterior mode (very marginal cost).
  mode   = apply(gamma, 2, find.mode, adjust = 10)
  mode_A = colSums(A) # find.mode (mode) or sum (mean)
  mode_A = mode_A/sum(mode_A)
  sd_A = sqrt(diag(var(A)))
  mode_B = colSums(B) # find.mode (mode) or sum (mean)
  mode_B = mode_B/sum(mode_B)
  sd_B = sqrt(diag(var(B)))
  # find.mode is 20-30 % faster than posterior.moode
  
  # transcript level test:
  trancript_res = 1-pchisq(mode^2/diag(CV), df = 1)
  
  p = K-1
  
  p_value = vapply(seq_len(K), function(k){
    sel  = seq_len(K)[-k]
    # Normal (classical Wald test)
    stat = t(mode[sel]) %*% ginv(CV[sel, sel], tol = 0) %*% mode[sel]
    1-pchisq(stat, df = K-1)
  }, FUN.VALUE = numeric(1))
  
  # I return 4 versions of the p.value:
  # 1) an average of the K p.values
  # 2) the p.value obtained removing the smallest difference (min(gamma))
  # 3) the p.value obtained removing the (overall summing the two groups) most lowly expressed transcript.
  # 4) a randomly selected p_value
  #sel_1 = which.min(abs(mode)) # min diff between pi's in A and B.
  #sel_2 = which.min(mode_A + mode_B) # most lowly expressed transcript overall in A + B.
  #ran = sample.int(K, 1)
  
  # I also record if the dominant transcript is inverted between the two conditions.
  # Inverted defined w.r.t the posterior mode.
  inverted = which.max(mode_A) != which.max(mode_B)
  # In this case I can also consider less stringent constraints such as the Chi_2 maybe.
  
  # Score to highlight the impact of DS:
  # max_diff_pi_T = max(abs(mode_A - mode_B)); 
  top2_diff_pi_T = sum(sort(abs(mode_A - mode_B), decreasing = TRUE)[seq_len(2)])
  
  list( c(  mean(p_value),   # p_value[sel_1],   p_value[sel_2],   p_value[ran], 
            inverted, # max_diff_pi_T, 
            top2_diff_pi_T ),
        trancript_res,
        mode_A, mode_B, sd_A, sd_B)
}
