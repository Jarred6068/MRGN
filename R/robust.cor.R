# This function is internal,
# A copy from library 'MRPC' (function 'RobustCor')
# Because MRPC depends on packages only available from 'bioconductor'.
# NOT to be exported: this way, we could later replace it by a dependence on
# MRPC (if it no longer depends on 'bioconductor') to attribute credit the
# authors of that package
robust.cor <- function(x,
                       Beta = 0.3, # In [0, 0.3]
                       partial = FALSE, # logical, return robust partial correlations?
                       maxiter = 500,
                       return.weights = FALSE,
                       tol = 0.005, #  Convergence control parameter
                       warn = TRUE) {
  ## Check the argument Beta
  if (!missing(Beta)) {
    Beta <- Beta[1]
    if (Beta < 0 | Beta > 0.3)
      stop("'Beta' must be in [0, 0.3].")
  }

  ## Data and sizes
  x <- as.matrix(x)
  n <- NROW(x)
  m <- NCOL(x)

  # Call cor if Beta == 0
  if (Beta == 0) {
    R <- if (partial)
      partial.cor(x, use = 'pairwise.complete.obs')
    else
      cor(x, use = "pairwise.complete.obs")
    return(list(R = R,
                mean = colMeans(x),
                sd = apply(x, MARGIN = 2, FUN = sd, na.rm = TRUE),
                weights = if(return.weights[1]) rep(1, n)))
  }

  ## Initialization
  Median <- c(apply(x, MARGIN = 2, FUN = median, na.rm = TRUE))
  Dist <- sqrt(colSums( (t(x) - Median)^2, na.rm = TRUE))
  qDist <- as.numeric(quantile(Dist, p=.5, na.rm = TRUE))
  Data0 <- sapply(1:n, FUN = function(i) {
    if (Dist[i] <= qDist) x[i,]
  })
  if (is.list(Data0))
    Data0 <- unlist(Data0, recursive = TRUE)
  if (!is.matrix(Data0))
    Data0 <- matrix(c(Data0), ncol = m, byrow = FALSE)
  else
    Data0 <- t(Data0)
  Mo <-as.numeric(colMeans(Data0, na.rm = T))
  Vo <- as.matrix(cov(Data0, use = "pairwise.complete.obs"))
  x <- replace(x, is.na(x), 0)
  DiffNorm <- 1
  Iter <- 0

  ## Iterations
  while (DiffNorm > tol && Iter <= maxiter) {
    Wx <- NULL
    for (j1 in 1:n) {
      zo <- as.numeric(x[j1,] - Mo)
      zz <- Beta*(t(zo) %*% mpinv(Vo) %*% zo)
      Wx[j1] <- exp(-zz/2)
    }
    M.new <- matrix(0, nrow = m)
    V.new <- matrix(0, nrow = m, ncol = m)
    for (j2 in 1:n) {
      M.new <- M.new + (Wx[j2] * x[j2,]) / sum(Wx)
      V.new <- V.new + (1 + Beta) * (Wx[j2] * (x[j2,] - Mo) %*% t(x[j2,] - Mo)) / sum(Wx)
    }
    DiffNorm <- sqrt(sum((M.new - Mo)^2)) / m + sqrt(sum((V.new - Vo)^2)) / m
    Mo <- M.new
    Vo <- V.new
    Iter <- Iter + 1
  }

  if (warn[1] && DiffNorm > tol) {
    warning(paste0("The current criterion value is ", round(DiffNorm, 4),
                   ", not below the tolerance value ", round(tol, 4),
                   ".\n", "     Consider increasing 'maxiter' above ", maxiter,
                   ", or decreasing 'Beta' below ", Beta, "."))
  }

  # Retain colnames from x
  dimnames(V.new) <- list(colnames(x), colnames(x))

  # Convert to partial correlations if required
  if (partial) {
    R <- - cov2cor(mpinv(V.new))
    diag(R) = 1
  }
  else
    R <- cov2cor(V.new)
  mean <- c(M.new)
  names (mean) <- colnames(x)

  return(list(R = R,
              mean = mean,
              sd = sqrt(diag(V.new)),
              weights = if(return.weights[1]) n * Wx / sum(Wx)))
}


