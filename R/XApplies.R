# Internal functions
# Apply Operations possibly using Clusters
#
# A set of apply-like functions taking an optional cluster for parallel computing (wrapper)
# The purpose is to just call this function without caring about the availability
# of a cluster or not. This function calls the appropriate non-parallel or parallel function

# Pass partout = 'matting', or 'master key', or 'skeleton key'
# matte <=> glossy ('satine' ou brillant ou luisant)

# Note that the argument 'cl' comes after '...', and must thus be explicitely named in a call.
# (same for 'chunk.size')

# cl <- parallel::makeCluster(getOption("cl.cores", 2))

# set.seed(167) ; x = rnorm(1e+8)       # user    system  elapsed
# system.time(sum(x))                   # 0.155   0.002   0.157
# system.time(sum(sapply (x, FUN=sum))) # 56.259 177.362 653.389 seconds
# system.time(sum(mrbglm:::matteSapply (x, FUN=sum))) 48.480  292.824 1397.304
# system.time(sum(mrbglm:::matteSapply (x, FUN=sum, cl = cl))) 43.485 226.251 895.329

# cl <- parallel::makeCluster(getOption("cl.cores", 8))

# set.seed(167) ; x = rnorm(1e+6)       # user    system  elapsed
# system.time(sum(x))                   # 0.001   0.000   0.002
# system.time(sum(sapply (x, FUN=sum))) # 0.299   0.010   0.308 seconds
# system.time(sum(mrbglm:::matteSapply (x, FUN=sum))) 0.288   0.005   0.293
# system.time(sum(mrbglm:::matteSapply (x, FUN=sum, cl = cl))) 0.218   0.048   0.270

matteApply <- function (X, MARGIN, FUN, ..., chunk.size = NULL, simplify = TRUE,
                        cl = parallel::getDefaultCluster()) {
  if (is.null(cl)) {
    return(
      apply(X = X, MARGIN = MARGIN, FUN = FUN, ..., simplify = simplify)
    )
  }
  else {
    if (simplify) {
      return(
        parallel::parApply(cl = cl, X = X, MARGIN = MARGIN, FUN = FUN, ...,
                           chunk.size = chunk.size)
      )
    }
    else {
      Nb.elts <- dim(X)[MARGIN]
      if (MARGIN == 2)
        X <- t(X)
      ad.args <- list(...)
      if (length(ad.args)) {
        lFUN <- function(i) {
          FUN (X[i,], ...)
        }
      }
      else {
        lFUN <- function(i) {
          FUN (X[i,])
        }
      }
      return(
        parallel::parLapply(cl = cl, X = 1:Nb.elts, fun = lFUN,
                            chunk.size = chunk.size)
      )
    }
  }
}

matteLapply <- function (X, FUN, ..., chunk.size = NULL,
                         cl = parallel::getDefaultCluster(), LB = FALSE) {
  if (is.null(cl)) {
    return(
      lapply(X = X, FUN = FUN, ...)
    )
  }
  else {
    if (LB)
      return(
        parallel::parLapplyLB(cl = cl, X = X, fun = FUN, ..., chunk.size = chunk.size))
    return(parallel::parLapply(cl = cl, X = X, fun = FUN, ..., chunk.size = chunk.size))
  }
}

matteSapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE,
          chunk.size = NULL, cl = parallel::getDefaultCluster(), LB = FALSE) {
  if (is.null(cl)) {
    return(
      sapply (X = X, FUN = FUN, ..., simplify = simplify, USE.NAMES = USE.NAMES)
    )
  }
  else {
    if (LB)
    return(
      parallel::parSapplyLB(cl = cl, X = X, FUN = FUN, ...,
                          USE.NAMES = USE.NAMES,
                          chunk.size = chunk.size))
    return(
      parallel::parSapply(cl = cl, X = X, FUN = FUN, ...,
                          USE.NAMES = USE.NAMES,
                          chunk.size = chunk.size))
  }
}

# cl <- makeCluster(getOption("cl.cores", 2))
