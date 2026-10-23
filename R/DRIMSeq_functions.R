###############################################################################
# Internal functions adapted from the DRIMSeq package
###############################################################################
# DRIMSeq (author: Malgorzata Nowicka; License: GPL (>= 3)) is being removed
# from Bioconductor.
# BANDITS only used DRIMSeq in 'prior_precision', to obtain gene-wise
# estimates of the precision parameter of the Dirichlet-Multinomial model via:
#   d = DRIMSeq::dmDSdata(counts = counts_df, samples = samples)
#   d = DRIMSeq::dmPrecision(d, genewise_precision = TRUE,
#                            design = model.matrix(~ 1, data = samples(d)),
#                            BPPARAM = BPPARAM)
#   DRIMSeq::genewise_precision(d)
#
# The functions below reproduce exactly that computation, with all other
# 'dmPrecision' arguments at their DRIMSeq defaults.
# The code was copied from DRIMSeq (version 1.38.0) and simplified to the
# intercept-only design (i.e., one group of samples) used in BANDITS:
#   - no S4 classes: counts are stored in a named list of matrices
#     (transcripts x samples), one matrix per gene;
#   - only the "one-way" fitting with a single group is kept
#     (no regression fit, no beta-binomial model);
#   - only the "trended" precision moderation (DRIMSeq's default) is kept;
#   - progress messages and timings were removed.
# The numerical code is otherwise unchanged, hence the results (and the use of
# random numbers) are the same as those obtained with DRIMSeq.
#
# Two helper functions used by DRIMSeq come from edgeR (License: GPL (>= 2)):
# 'movingAverageByCol' (R code by Gordon Smyth, copied unchanged) and
# 'maximizeInterpolant' (C code by Aaron Lun and Lizhong Chen, as in
# edgeR >= 4.4.0, translated to R here, together with the 'fmm_spline' C
# routine that it uses, which comes from the R 'stats' package).
#
# References:
# Nowicka M, Robinson MD (2016). DRIMSeq: a Dirichlet-multinomial framework
# for multivariate count outcomes in genomics. F1000Research, 5, 1356.
# McCarthy DJ, Chen Y, Smyth GK (2012). Differential expression analysis of
# multifactor RNA-Seq experiments with respect to biological variation.
# Nucleic Acids Research, 40, 4288-4297.
###############################################################################


###############################################################################
# Data: replaces DRIMSeq::dmDSdata
###############################################################################

# Split a (transcripts x samples) count matrix into a named list of matrices,
# one per gene.
# As in DRIMSeq::dmDSdata, genes are sorted by order of first appearance
# in 'gene_id', and transcripts keep their original order within each gene.
dm_counts_by_gene <- function(counts, gene_id, feature_id){

  # same checks as in DRIMSeq::dmDSdata
  stopifnot(mode(counts) %in% "numeric")
  stopifnot(nrow(counts) == length(gene_id))
  stopifnot(nrow(counts) == length(feature_id))
  stopifnot(sum(duplicated(feature_id)) == 0)
  stopifnot(all(!is.na(gene_id)))
  stopifnot(all(!is.na(feature_id)))

  rownames(counts) <- feature_id

  gene_id <- factor(gene_id, levels = unique(gene_id))

  inds <- split(seq_len(nrow(counts)), gene_id)

  lapply(inds, function(i) counts[i, , drop = FALSE])

}


###############################################################################
# Precision estimation: replaces DRIMSeq::dmPrecision
###############################################################################

# Gene-wise precision estimates (named numeric vector, NA when the precision
# can not be estimated), as returned by DRIMSeq::genewise_precision after
# DRIMSeq::dmPrecision(..., design = model.matrix(~ 1, ...)).
# Arguments and defaults are those of DRIMSeq::dmPrecision.
# NB: as in DRIMSeq, the common precision (used to center the grid of the
# gene-wise estimates) is computed on a random subset of genes:
# call 'set.seed' beforehand to make the results reproducible.
dm_genewise_precision <- function(counts,
  prec_subset = 0.1, prec_interval = c(0, 1e+3), prec_tol = 1e+01,
  prec_grid_length = 21, prec_grid_range = c(-10, 10),
  prec_prior_df = 0, prec_span = 0.1, prop_tol = 1e-12,
  BPPARAM = BiocParallel::SerialParam()){

  # mean gene expression, used for the trended moderation:
  mean_expression <- dm_estimateMeanExpression(counts = counts)

  # common precision, estimated on a random subset of the genes:
  if(prec_subset < 1){
    genes2keep <- sample(seq_len(length(counts)),
      max(round(prec_subset * length(counts)), 1), replace = FALSE)
  }else{
    genes2keep <- seq_len(length(counts))
  }

  common_precision <- dmDS_estimateCommonPrecision(
    counts = counts[genes2keep],
    prec_interval = prec_interval, prec_tol = prec_tol,
    prop_tol = prop_tol, BPPARAM = BPPARAM)

  # gene-wise precision, using the common precision as 'prec_init':
  dmDS_estimateTagwisePrecision(counts = counts,
    mean_expression = mean_expression, prec_init = common_precision,
    prec_grid_length = prec_grid_length, prec_grid_range = prec_grid_range,
    prec_prior_df = prec_prior_df, prec_span = prec_span,
    prop_tol = prop_tol, BPPARAM = BPPARAM)

}


dm_estimateMeanExpression <- function(counts){

  # Calculate mean expression of genes
  mean_expression <- unlist(lapply(seq_along(counts), function(g){
    mean(colSums(counts[[g]]), na.rm = TRUE)
  }))

  names(mean_expression) <- names(counts)

  return(mean_expression)

}


dmDS_profileLikCommon <- function(prec, counts, prop_tol = 1e-12,
  BPPARAM = BiocParallel::SerialParam()){

  adj_lik <- dmDS_profileLik(prec = prec, counts = counts,
    prop_tol = prop_tol, BPPARAM = BPPARAM)

  adj_lik_common <- sum(adj_lik, na.rm = TRUE)

  return(adj_lik_common)

}


dmDS_estimateCommonPrecision <- function(counts,
  prec_interval = c(0, 1e+3), prec_tol = 1e+01, prop_tol = 1e-12,
  BPPARAM = BiocParallel::SerialParam()){

  optimum <- optimize(f = dmDS_profileLikCommon,
    interval = prec_interval,
    counts = counts, prop_tol = prop_tol, BPPARAM = BPPARAM,
    maximum = TRUE, tol = prec_tol)

  precision <- optimum$maximum

  return(precision)

}


dmDS_estimateTagwisePrecision <- function(counts, mean_expression,
  prec_init = 100, prec_grid_length = 21, prec_grid_range = c(-10, 10),
  prec_prior_df = 0, prec_span = 0.1, prop_tol = 1e-12,
  BPPARAM = BiocParallel::SerialParam()){

  ### Standard grid like in edgeR
  spline_pts <- seq(from = prec_grid_range[1], to = prec_grid_range[2],
    length.out = prec_grid_length)

  spline_prec <- prec_init * 2^spline_pts

  # Calculate the likelihood for each gene at the spline precision points
  loglik <- matrix(NA, nrow = length(counts), ncol = prec_grid_length)

  for(i in seq_len(prec_grid_length)){

    loglik[, i] <- dmDS_profileLik(prec = spline_prec[i], counts = counts,
      prop_tol = prop_tol, BPPARAM = BPPARAM)

  }

  not_nas <- complete.cases(loglik)
  loglik <- loglik[not_nas, , drop = FALSE]

  if(nrow(loglik) == 0){
    precision <- rep(NA, length(counts))
    names(precision) <- names(counts)
    return(precision)
  }

  # precision moderation towards the (precision versus mean expression) trend
  mean_expression <- mean_expression[not_nas]

  loglik <- dm_profileLikModeration(loglik = loglik,
    mean_expression = mean_expression,
    prec_prior_df = prec_prior_df, prec_span = prec_span)

  out <- edgeR_maximizeInterpolant(spline_pts, loglik)

  # Set NA for genes that tagwise prec could not be calculated
  precision <- rep(NA, length(counts))
  names(precision) <- names(counts)
  precision[not_nas] <- prec_init * 2^out

  return(precision)

}


###############################################################################
# Cox-Reid adjusted profile likelihood (one group of samples)
###############################################################################

dmDS_profileLik <- function(prec, counts, prop_tol = 1e-12,
  BPPARAM = BiocParallel::SerialParam()){

  fit <- dmDS_fit(counts = counts, precision = prec,
    prop_tol = prop_tol, BPPARAM = BPPARAM)

  adj <- dmDS_CRadjustment(counts = counts, prop = fit$prop,
    precision = prec, BPPARAM = BPPARAM)

  adj_lik <- fit$lik - adj

  # adj_lik has length G
  return(adj_lik)

}


dmDS_fitOneGroup_gene <- function(g, counts, prec, prop_tol){

  # equivalent to DRIMSeq's 'dm_fitManyGroups' with 1 group only
  y <- counts[[g]]
  q <- nrow(y)

  if(q < 2 || is.na(prec[g]))
    return(list(prop = rep(NA, q), lik = NA))

  f <- dm_fitOneGroup(y = y, prec = prec[g], prop_tol = prop_tol)

  if(is.na(f[["lik"]]))
    return(list(prop = rep(NA, q), lik = NA))

  return(f)

}


dmDS_fit <- function(counts, precision, prop_tol = 1e-12,
  BPPARAM = BiocParallel::SerialParam()){

  inds <- seq_along(counts)

  # Prepare precision
  if(length(precision) == 1){
    prec <- rep(precision, length(inds))
  } else {
    prec <- precision
  }

  ff <- BiocParallel::bplapply(inds, dmDS_fitOneGroup_gene,
    counts = counts, prec = prec, prop_tol = prop_tol,
    BPPARAM = BPPARAM)

  names(ff) <- names(counts)

  lik <- unlist(lapply(ff, function(f) sum(f[["lik"]])))

  prop <- lapply(ff, function(f) f[["prop"]])

  # prop is a list of G vectors of proportions (one per gene)
  # lik is a vector of length G
  return(list(prop = prop, lik = lik))

}


dmDS_CRadjustmentOneGroup_gene <- function(g, counts, prec, prop){

  # equivalent to DRIMSeq's 'dm_CRadjustmentManyGroups' with 1 group only
  if(is.na(prop[[g]][1]))
    return(NA)

  adj <- dm_CRadjustmentOneGroup(y = counts[[g]], prec = prec[g],
    prop = prop[[g]])

  if(is.na(adj))
    return(NA)

  if(abs(adj) == Inf)
    return(NA)

  return(adj)

}


dmDS_CRadjustment <- function(counts, prop, precision,
  BPPARAM = BiocParallel::SerialParam()){

  inds <- seq_along(counts)

  # Prepare precision
  if(length(precision) == 1){
    prec <- rep(precision, length(inds))
  } else {
    prec <- precision
  }

  a <- BiocParallel::bplapply(inds, dmDS_CRadjustmentOneGroup_gene,
    counts = counts, prec = prec, prop = prop, BPPARAM = BPPARAM)

  names(a) <- names(counts)
  adj <- unlist(a)

  # adj is a vector of length G
  return(adj)

}


###############################################################################
# Dirichlet-Multinomial model for one group of samples
###############################################################################

dm_fitOneGroup <- function(y, prec, prop_tol = 1e-12){
  # y matrix q x n
  # If something is wrong, return NAs

  q <- nrow(y)

  # NAs for genes with one feature
  if(q < 2 || is.na(prec))
    return(list(prop = rep(NA, q), lik = NA))

  # Check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  # There must be at least two non-zero features
  if(sum(keep_row) < 2)
    return(list(prop = rep(NA, q), lik = NA))

  # Last feature can not be zero since
  # we use the last feature as a denominator in logit
  if(keep_row[q] == 0)
    return(list(prop = rep(NA, q), lik = NA))

  y <- y[keep_row, , drop=FALSE]
  q <- nrow(y)

  # Check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]

  prop_init <- rowSums(y)/sum(y)

  # If there is only one replicate, use empirical props as output
  if(sum(keep_col) == 1){

    lik <- dm_likG(prop = prop_init[-q], prec = prec, y = y)

    keep_row[keep_row] <- prop_init
    prop <- keep_row

    return(list(prop = prop, lik = lik))
  }

  ### Must have constraint for SUM prop = 1 -->
  ### sum(prop) < 1 + eps & sum(prop) > 1 - eps
  ui <- rbind(diag(rep(1, q-1), q-1), diag(rep(-1, q-1), q-1), rep(-1, q-1))
  ci <- c(rep(0, q-1), rep(-1, q-1), -1 + .Machine$double.eps)

  # Maximization
  co <- constrOptim(theta = prop_init[-q], f = dm_likG,
    grad = dm_scoreG,
    ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol),
    method = "BFGS",
    prec = prec, y = y)

  if(co$convergence == 0){
    prop <- co$par
    prop <- c(prop, 1 - sum(prop))
    lik <- co$value
  }else{
    return(list(prop = rep(NA, length(keep_row)), lik = NA))
  }

  keep_row[keep_row] <- prop
  prop <- keep_row

  # prop numeric vector of length q
  # lik numeric of lenght 1
  return(list(prop = prop, lik = lik))

}


dm_CRadjustmentOneGroup <- function(y, prec, prop){
  # y martix q x n
  # prop vector of length q
  # If something is wrong, return NAs

  if(any(is.na(prop)))
    return(NA)

  q <- nrow(y)

  # NAs for genes with one feature
  if(q < 2 || is.na(prec))
    return(NA)

  # Check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  # There must be at least two non-zero features
  if(sum(keep_row) < 2)
    return(NA)

  # Last feature can not be zero since
  # we use the last feature as a denominator in logit
  if(keep_row[q] == 0)
    return(NA)

  y <- y[keep_row, , drop=FALSE]
  q <- nrow(y)

  # Check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]

  prop <- prop[keep_row]
  n <- ncol(y)

  H <- dm_HessianG(prop = prop[-q], prec, y)

  adj <- log(det(n * (- H) ))/2

  return(adj)

}


# DM log-likelihood (with gamma functions) for q-1 parameters
dm_likG <- function(prop, prec, y){
  # prop has length of q-1
  # prec has length 1
  # y has q rows and n columns
  # This function returns likelihhod without normalizing component,
  # but it is OK for optimization and the LR test

  n <- ncol(y)
  m <- colSums(y)

  prop <- c(prop, 1 - sum(prop))

  l <- n * lgamma(prec) - sum(lgamma(m + prec)) +
    sum( colSums( lgamma(y + prop * prec) - lgamma(prop * prec) ) )

  return(l)

}


# DM score function (with gamma functions) for q-1 parameters
dm_scoreG <- function(prop, prec, y){
  # prop has length of q-1
  # prec has length 1
  # y has q rows and n columns

  q <- nrow(y)
  n <- ncol(y)
  yq <- y[q,]
  y <- y[-q, , drop=FALSE]
  propq <- 1-sum(prop)

  S <- prec * rowSums( digamma(y + prop * prec) -
      digamma(prop * prec) - matrix(digamma(yq + prec * propq) -
          digamma(prec * propq), nrow = q-1, ncol = n, byrow = TRUE) )

  return(S)

}


# DM Hessian (with gamma functions) for q-1 parameters
dm_HessianG <- function(prop, prec, y){
  # prop has length of q-1
  # prec has length 1
  # y has q rows and n columns

  q <- nrow(y)
  n <- ncol(y)
  Djj <- rep(0, q-1)
  propq <- 1-sum(prop)
  yq <- y[q, ]
  y <- y[-q, , drop=FALSE]

  Djl <- prec^2 * sum(trigamma(yq + propq * prec) - trigamma(propq * prec))

  Djj <- prec^2 * rowSums(trigamma(y + prop * prec) - trigamma(prop * prec))

  H <- matrix(Djl, q-1, q-1)

  diag(H) <- diag(H) + Djj

  # H martix (q-1) x (q-1)
  return(H)

}


###############################################################################
# Moderation of the profile likelihood towards the trend
# (DRIMSeq's default: prec_moderation = "trended")
###############################################################################

dm_profileLikModeration <- function(loglik, mean_expression,
  prec_prior_df, prec_span){

  prec_grid_length <- ncol(loglik)

  ### Check where the grid is maximized
  grid_max <- apply(loglik, 1, which.max)

  # In the calculation of moderation, do not take into account genes
  # that have precision on the top and bottom boundry of the grid
  # (skipp 4 last grid points and 1 first grid point)
  not_boundry <- grid_max < (prec_grid_length - 3) & grid_max > 1
  boundry_last <- grid_max == prec_grid_length

  ### Calculate the span of the boundry loglikelihoods
  if(sum(boundry_last) > 1){
    loglik_span_boundry <- apply(loglik[boundry_last, , drop = FALSE], 1,
      function(x){max(x) - min(x)})
  }

  moderation <- dm_movingAverageByCol(loglik = loglik,
    mean_expression = mean_expression, not_boundry = not_boundry,
    prec_span = prec_span)

  # Estimate priorN - calculate the ratio between moderation lik span
  # and lik span of boundry genes
  if(sum(boundry_last) > 10){

    moderation_span_boundry <- apply(
      moderation[boundry_last, , drop = FALSE], 1,
      function(x){max(x) - min(x)})
    span_ratio <- moderation_span_boundry / loglik_span_boundry
    priorN <- 1/span_ratio

    ### Do loess fitting if there is enough points. Otherwise, use median
    if(length(loglik_span_boundry) > 100){

      df_priorN_loglog <- data.frame(priorN = log10(priorN),
        mean_expression = log10(mean_expression[boundry_last]))

      priorN_loess_loglog <- loess(priorN ~ mean_expression,
        df_priorN_loglog, control = loess.control(surface = "direct"))
      priorN_predict_loglog <- predict(priorN_loess_loglog,
        data.frame(mean_expression = log10(mean_expression)), se = FALSE)

      priorN <- 10 ^ priorN_predict_loglog

    }else{

      priorN <- quantile(priorN, 0.5)

    }

  }else{
    priorN <- prec_prior_df
  }

  loglik <- loglik + priorN * moderation

  return(loglik)

}


dm_movingAverageByCol <- function(loglik, mean_expression,
  not_boundry, prec_span){

  if(sum(not_boundry) == length(not_boundry)){

    o <- order(mean_expression)
    oo <- order(o)
    width <- floor(prec_span * nrow(loglik))

    moderation <- edgeR_movingAverageByCol(loglik[o,], width = width)[oo,]

  }else{

    ### Use non boundry genes for calculating the moderation
    mean_expression_not_boundry <- mean_expression[not_boundry]
    loglik_not_boundry <- loglik[not_boundry, , drop = FALSE]

    o <- order(mean_expression_not_boundry)
    oo <- order(o)

    width <- floor(prec_span * nrow(loglik_not_boundry))

    moderation_not_boundry <- edgeR_movingAverageByCol(
      loglik_not_boundry[o, , drop = FALSE], width = width)[oo, , drop = FALSE]

    ### Fill in moderation values for the boundy genes
    moderation <- matrix(NA, nrow = nrow(loglik), ncol = ncol(loglik))

    moderation[not_boundry, ] <- moderation_not_boundry

    o <- order(mean_expression)
    oo <- order(o)

    moderation <- moderation[o, , drop = FALSE]
    not_boundry <- not_boundry[o]

    ### Last value in not_boundry must be TRUE
    if(not_boundry[length(not_boundry)] == FALSE){

      last_true <- max(which(not_boundry))
      moderation[length(not_boundry), ] <- moderation[last_true, ]

      not_boundry[length(not_boundry)] <- TRUE

    }

    not_boundry_diff <- diff(not_boundry, lag = 1)

    not_boundry_cumsum <- cumsum(not_boundry)

    ### Values used for filling in the boundry NAs - swith from FALSE to TRUE
    replacement_indx <- which(not_boundry_diff == 1) + 1

    replaced_indx <- which(!not_boundry)

    replaced_freq <- as.numeric(table(not_boundry_cumsum[replaced_indx]))

    moderation_boundry  <- moderation[rep(replacement_indx,
      times = replaced_freq), , drop = FALSE]

    moderation[!not_boundry, ] <- moderation_boundry

    moderation <- moderation[oo, , drop = FALSE]

  }

  return(moderation)

}


###############################################################################
# Helper functions from edgeR
###############################################################################

# edgeR::movingAverageByCol
# (unchanged, except for the use of 'seq_len' instead of '1:...'; this gives
# identical results whenever width <= nrow(x), as is always the case here)
# Moving average smoother for columns of a matrix
# Gordon Smyth, 17 Feb 2011
edgeR_movingAverageByCol <- function(x,width=5,full.length=TRUE)
{
  x <- as.matrix(x)
  width <- as.integer(width)
  if(width<=1) return(x)
  n <- nrow(x)
  m <- ncol(x)
  if(width>n) {
    width <- n
    warning("reducing moving average width to nrow(x)")
  }
  if(full.length) {
    half1 <- ceiling(width/2)
    half2 <- floor(width/2)
    x <- rbind(matrix(0,half1,m),x,matrix(0,half2,m))
  } else {
    if(width==n) return(matrix(colMeans(x),1L,m))
    x <- rbind(matrix(0,1,m),x)
  }
  n2 <- nrow(x)
  x <- apply(x,2,cumsum)
  x <- x[(width+1):n2,,drop=FALSE]-x[seq_len(n2-width),,drop=FALSE]
  n3 <- nrow(x)
  w <- rep(width,n3)
  if(full.length) {
    if(half1>1) w[seq_len(half1-1)] <- width-rev(seq_len(half1-1))
    w[(n3-half2+1):n3] <- width-seq_len(half2)
  }
  x/w
}


# edgeR::maximizeInterpolant, translated from C to R.
# For each row of the likelihood matrix 'y' (genes x spline points), find the
# value of 'x' (ordered spline points) at which the cubic spline interpolating
# the likelihood is maximal: the spline is fitted (Forsythe, Malcolm and Moler
# end conditions, as in stats::splinefun(method = "fmm")), and its maximum is
# searched in the two segments neighbouring the grid point with the largest
# likelihood, by solving for the roots of the derivative.
# NA/NaN values are treated as in the original C code (comparisons with
# NaN are FALSE).
edgeR_maximizeInterpolant <- function(x, y){

  out <- numeric(nrow(y))

  for(g in seq_len(nrow(y))){
    out[g] <- edgeR_find_max(x = x, y = y[g, ])
  }

  return(out)

}


edgeR_find_max <- function(x, y){

  npts <- length(x)

  # Getting a good initial guess for the MLE.
  maxed <- -1
  maxed_at <- -1
  for(i in seq_len(npts)){
    if(maxed_at < 0 || isTRUE(y[i] > maxed)){
      maxed <- y[i]
      maxed_at <- i
    }
  }
  x_max <- x[maxed_at]

  coef <- fmm_spline(x = x, y = y)
  b <- coef$b
  c <- coef$c
  d <- coef$d

  # The spline coefficients are such that, in each segment, the spline is
  # equal to 'y + b*t + c*t^2 + d*t^3', where 't = x_t - x_l' and 'x_l' is the
  # left limit of the segment: the MLE is 'x_l + t', where 't' is the root of
  # the derivative 'b + 2*c*t + 3*d*t^2' corresponding to a maximum.

  # First we have a look at the segment on the left and see if it contains
  # the maximum.
  if(maxed_at > 1){
    ld <- d[maxed_at-1]
    lc <- c[maxed_at-1]
    lb <- b[maxed_at-1]

    delta <- lc*lc-3*ld*lb

    if(!isTRUE(delta < 0)){
      # Using the solution with the maximum (not minimum).
      chosen_sol <- (-lc-sqrt(delta))/(3*ld)

      if(isTRUE(chosen_sol > 0) &&
         isTRUE(chosen_sol < x[maxed_at]-x[maxed_at-1])){
        temp <- ((ld*chosen_sol+lc)*chosen_sol+lb)*chosen_sol+y[maxed_at-1]
        if(isTRUE(temp > maxed)){
          maxed <- temp
          x_max <- chosen_sol+x[maxed_at-1]
        }
      }
    }
  }

  # Repeating for the segment on the right.
  if(maxed_at < npts){
    rd <- d[maxed_at]
    rc <- c[maxed_at]
    rb <- b[maxed_at]

    delta <- rc*rc-3*rd*rb

    if(!isTRUE(delta < 0)){
      chosen_sol <- (-rc-sqrt(delta))/(3*rd)

      if(isTRUE(chosen_sol > 0) &&
         isTRUE(chosen_sol < x[maxed_at+1]-x[maxed_at])){
        temp <- ((rd*chosen_sol+rc)*chosen_sol+rb)*chosen_sol+y[maxed_at]
        if(isTRUE(temp > maxed)){
          maxed <- temp
          x_max <- chosen_sol+x[maxed_at]
        }
      }
    }
  }

  return(x_max)

}


# Spline coefficients a la Forsythe, Malcolm and Moler; translated from the
# C function 'fmm_spline' (in 'splines.c' of the R 'stats' package,
# Copyright (C) 1995-2001 Robert Gentleman, Ross Ihaka and the R Core Team;
# GPL (>= 2)), used by edgeR::maximizeInterpolant.
# The end-conditions are determined by fitting cubic polynomials to the first
# and last 4 points and matching the third derivatives of the spline at the
# end-points to the third derivatives of these cubics at the end-points.
# The spline in [x[i], x[i+1]] is y[i] + b[i]*t + c[i]*t^2 + d[i]*t^3,
# with t = x - x[i].
fmm_spline <- function(x, y){

  n <- length(x)
  b <- c <- d <- numeric(n)

  if(n < 2)
    return(list(b = b, c = c, d = d))

  if(n < 3){
    t <- (y[2] - y[1])
    b[1] <- t / (x[2]-x[1])
    b[2] <- b[1]
    return(list(b = b, c = c, d = d))
  }

  nm1 <- n - 1

  # Set up tridiagonal system
  # b = diagonal, d = offdiagonal, c = right hand side
  d[1] <- x[2] - x[1]
  c[2] <- (y[2] - y[1])/d[1]
  for(i in 2:nm1){
    d[i] <- x[i+1] - x[i]
    b[i] <- 2.0 * (d[i-1] + d[i])
    c[i+1] <- (y[i+1] - y[i])/d[i]
    c[i] <- c[i+1] - c[i]
  }

  # End conditions.
  # Third derivatives at x[0] and x[n-1] obtained from divided differences
  b[1] <- -d[1]
  b[n] <- -d[nm1]
  c[1] <- c[n] <- 0.0
  if(n > 3){
    c[1] <- c[3]/(x[4]-x[2]) - c[2]/(x[3]-x[1])
    c[n] <- c[nm1]/(x[n] - x[n-2]) - c[n-2]/(x[nm1]-x[n-3])
    c[1] <- c[1]*d[1]*d[1]/(x[4]-x[1])
    c[n] <- -c[n]*d[nm1]*d[nm1]/(x[n]-x[n-3])
  }

  # Gaussian elimination
  for(i in 2:n){
    t <- d[i-1]/b[i-1]
    b[i] <- b[i] - t*d[i-1]
    c[i] <- c[i] - t*c[i-1]
  }

  # Backward substitution
  c[n] <- c[n]/b[n]
  for(i in nm1:1)
    c[i] <- (c[i]-d[i]*c[i+1])/b[i]

  # c[i] is now the sigma[i-1] of the text
  # Compute polynomial coefficients
  b[n] <- (y[n] - y[n-1])/d[n-1] + d[n-1]*(c[n-1]+ 2.0*c[n])
  for(i in seq_len(nm1)){
    b[i] <- (y[i+1]-y[i])/d[i] - d[i]*(c[i+1]+2.0*c[i])
    d[i] <- (c[i+1]-c[i])/d[i]
    c[i] <- 3.0*c[i]
  }
  c[n] <- 3.0*c[n]
  d[n] <- d[nm1]

  return(list(b = b, c = c, d = d))

}
