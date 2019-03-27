##############################################################################################################################
# General MCMC functions:
##############################################################################################################################
# initialize pi new (matrix) object
create_pi_new = function(g, K, N){
  matrix( 1/K[g], nrow = N, ncol = K[g])
}

# Fast posterior mode computation:
find.mode <- function(x, adjust, ...) {
  dx <- density(x, adjust = adjust, ...)
  dx$x[which.max(dx$y)]
}

# fast heidelberg diagnostic computation:
my_heidel.diag = function(x, R, by., pvalue = 0.01){
  start.vec <- seq(from = 1, to = R/2, by = by.)
  S0 <- my_spectrum0.ar(window(x, start = R/2), R/2+1)
  
  converged <- FALSE
  for(i in seq(along = start.vec)){
    x <- window(x, start = start.vec[i])
    n <- R + 1 - start.vec[i] # niter(x)
    B <- cumsum(x) - sum(x) * seq_len(n)/n
    Bsq <- (B * B)/(n * S0)
    I <- sum(Bsq)/n
    p = my_pcramer(I)
    if(converged <- !is.na(I) && p < 1 - pvalue) 
      break
  }
  
  if( !converged || is.na(I) ) {
    nstart <- NA
  }else {
    nstart <- start.vec[i]
  }
  return(c(converged, nstart, 1 - p))
}

my_pcramer = function(q, eps = 1e-05){
  log.eps <- log(eps)
  y = sapply(0:3, function(k){
    z <- gamma(k + 0.5) * sqrt(4 * k + 1)/(gamma(k + 1) * 
                                             pi^(3/2) * sqrt(q))
    u <- (4 * k + 1)^2/(16 * q)
    ifelse(u > -log.eps, 0, z * exp(-u) * besselK(x = u, 
                                                  nu = 1/4))
  })
  return(sum(y))
}

my_spectrum0.ar = function(x, R){
  lm.out <- lm(x ~ seq_len(R) )
  if(identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
    v0 <- 0
  }else{
    ar.out <- ar(x, aic = TRUE)
    v0 <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
  }
  #  return(list(spec = v0, order = order))
  v0
}

# char_in_char looks for string a in b.
char_in_char = function(a, b){
  len_a = nchar(a)
  len_b = nchar(b)

  # s represents the start of the sub-string of b
  for(s in seq_len(len_b - len_a + 1) ){
    if( a == substring(b, s, s+len_a-1) ){
      return(TRUE)
    }
  }
  
  FALSE
}
