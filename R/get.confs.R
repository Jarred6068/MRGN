

## Replacement function for get.confs
# Extract and return the 'cors' object from get.conf.matrix
get.list.conf <- function (data = NULL,
                           p,
                           q,
                           u = NCOL(data) - p - q, # Unknown confounders pool to select from
                           measure = c("correlation", "partial_corr"), # The same for T and V nodes?
                           conditional.vars = NULL,
                           blocksize = min(p, q, 100),
                           apply.qval = TRUE,                          # The same for T and V nodes?
                           #################  get.conf.matrix does not allow 'no correction'
                           selection_fdr = 0.05,
                           adjust_by = "individual",
                           lambda = NULL,
                           pi0.method = c("smoother", "boostrap"),
                           alpha = 0.05,
                           joint = FALSE, # Jointly select V and K/U nodes? cbind V and K/U nodes as a unique pool?
                           verbose = 0L) {
  ### Check arguments
  if (is.null(data))
    return(list())
  nvars <- NCOL(data)
  stopifnot(all(c(p, q, u) <= nvars))
  m <- p + q + u
  stopifnot(m <= nvars)
  data <- data[,1:m, drop = FALSE]
  if (q == nvars)
    return(list()) # No T-node to do selection for
  if (p == nvars)
    return(list()) # No V or K/U node to select from
  stopifnot(measure %in% c("correlation", "partial_corr"))
  measure <- measure[1]
  stopifnot(pi0.method %in% c("smoother", "boostrap"))
  pi0.method <- pi0.method[1]

  if (!missing(blocksize))
    blocksize <- min(blocksize, p, q, 100)

  if (joint) {
    cov.pool <- data[, c(1:q, if (u > 0) (p + q + 1):m)]

    if (is.null(conditional.vars))
      conditional.vars <- cov.pool

    VKUconfounders <- catch.conditions({
      get.conf.matrix(data = data[,(q + 1):(p + q), drop = FALSE],
                      cov.pool = cov.pool,
                      measure = measure,
                      conditional.vars = conditional.vars,
                      blocksize = min(blocksize, p, q + u),
                      apply.qval = apply.qval,
                      selection_fdr = selection_fdr,
                      adjust_by = adjust_by,
                      lambda = lambda,
                      pi0.method = pi0.method,
                      alpha = alpha)$sig.asso.covs
    })$value
    if (any(class(VKUconfounders) %in% c("simpleError", "error", "condition"))) {
      if (identical(pi0.method, "smoother")) {
        if (verbose) {
          warning("'pi0.meth = 'smoother'' failled, trying 'pi0.meth = bootstrap''")
        }
        VKUconfounders <- catch.conditions({
          get.conf.matrix(data = data[,(q + 1):(p + q), drop = FALSE],
                          cov.pool = cov.pool,
                          measure = measure,
                          conditional.vars = conditional.vars,
                          blocksize = min(blocksize, p, q + u),
                          apply.qval = apply.qval,
                          selection_fdr = selection_fdr,
                          adjust_by = adjust_by,
                          lambda = lambda,
                          pi0.method = "boostrap",
                          alpha = alpha)$sig.asso.covs
        })$value
      }
      if (any(class(VKUconfounders) %in% c("simpleError", "error", "condition"))) {
        VKUconfounders <- list()
      }
    }

    if (!length(VKUconfounders)) {
      VKUconfounders <- list()
    }
    else {
      if (is.list(VKUconfounders)) {
        Index <- 1:(q + u)
        VKUconfounders <- lapply(VKUconfounders,
                               FUN = function(selected) {
                                 if (any(selected)) {
                                   col.indices <- Index[selected]
                                   KUnode <- col.indices > q
                                   if (any(KUnode))
                                   col.indices[KUnode] <- col.indices[KUnode] + p
                                 }
                                 else
                                   col.indices <- NULL
                                 col.indices
                               })
      }
      else {
        KUnode <- VKUconfounders > q
        VKUconfounders[KUnode] <- VKUconfounders[KUnode] + q
        VKUconfounders <- as.list(VKUconfounders)
      }
    }
  }
  else {
    ### Select V-nodes
    if (q > 0) {
      cov.pool <- data[, 1:q]

      if (is.null(conditional.vars))
        conditional.vars <- cov.pool
      Vconfounders <- catch.conditions({
        get.conf.matrix(data = data[,(q + 1):(p + q), drop = FALSE],
                        cov.pool = cov.pool,
                        measure = measure,
                        conditional.vars = conditional.vars,
                        blocksize = min(blocksize, p, q),
                        apply.qval = apply.qval,
                        selection_fdr = selection_fdr,
                        adjust_by = adjust_by,
                        lambda = lambda,
                        pi0.method = pi0.method,
                        alpha = alpha)$sig.asso.covs
      })$value
      if (any(class(Vconfounders) %in% c("simpleError", "error", "condition"))) {
        if (identical(pi0.method, "smoother")) {
          if (verbose) {
            warning("'pi0.meth = 'smoother'' failled, trying 'pi0.meth = bootstrap''")
          }
          Vconfounders <- catch.conditions({
            get.conf.matrix(data = data[,(q + 1):(p + q), drop = FALSE],
                            cov.pool = cov.pool,
                            measure = measure,
                            conditional.vars = conditional.vars,
                            blocksize = min(blocksize, p, q),
                            apply.qval = apply.qval,
                            selection_fdr = selection_fdr,
                            adjust_by = adjust_by,
                            lambda = lambda,
                            pi0.method = "boostrap",
                            alpha = alpha)$sig.asso.covs
          })$value
        }
        if (any(class(Vconfounders) %in% c("simpleError", "error", "condition"))) {
          Vconfounders <- list()
        }
      }

      if (!length(Vconfounders)) {
        Vconfounders <- list()
      }
      else {
        if (is.list(Vconfounders)) {
          Index <- 1:q
          Vconfounders <- lapply(Vconfounders,
                                 FUN = function(selected) {
                                   col.indices <- if (any(selected)) Index[selected]
                                   if (length(col.indices)) {
                                     if (is.null(names(col.indices)))
                                       names(col.indices) <- paste0("V", col.indices)
                                   }
                                   col.indices
                                 })
        }
        else {
          Vconfounders <- as.list(Vconfounders)
        }
      }
    }
    else {
      Vconfounders <- list()
    }

    ### Select K/U-nodes
    if (u > 0) {
      cov.pool <- data[, (p + q + 1):nvars]

      if (is.null(conditional.vars))
        conditional.vars <- cov.pool

      KUconfounders <- catch.conditions({
        get.conf.matrix(data = data[,(q + 1):(p + q), drop = FALSE],
                        cov.pool = cov.pool,
                        measure = measure,
                        conditional.vars = conditional.vars,
                        blocksize = min(blocksize, p, u),
                        apply.qval = apply.qval,
                        selection_fdr = selection_fdr,
                        adjust_by = adjust_by,
                        lambda = lambda,
                        pi0.method = pi0.method,
                        alpha = alpha)$sig.asso.covs
      })$value
      if (any(class(KUconfounders) %in% c("simpleError", "error", "condition"))) {
        if (identical(pi0.method, "smoother")) {
          if (verbose) {
            warning("'pi0.meth = 'smoother'' failled, trying 'pi0.meth = bootstrap''")
          }
          KUconfounders <- catch.conditions({
            get.conf.matrix(data = data[,(q + 1):(p + q), drop = FALSE],
                            cov.pool = cov.pool,
                            measure = measure,
                            conditional.vars = conditional.vars,
                            blocksize = min(blocksize, p, u),
                            apply.qval = apply.qval,
                            selection_fdr = selection_fdr,
                            adjust_by = adjust_by,
                            lambda = lambda,
                            pi0.method = "boostrap",
                            alpha = alpha)$sig.asso.covs
          })$value
        }
        if (any(class(KUconfounders) %in% c("simpleError", "error", "condition"))) {
          KUconfounders <- list()
        }
      }

      if (!length(KUconfounders)) {
        KUconfounders <- list()
      }
      else {
        if (is.list(KUconfounders)) {
          Index <- 1:u
          offset <- p + q
          KUconfounders <- lapply(KUconfounders,
                                  FUN = function(selected) {
                                    col.indices <- if (any(selected)) Index[selected] + offset
                                    col.indices
                                  })
        }
        else {
          KUconfounders <- as.list(KUconfounders)
        }
      }
    }
    else {
      KUconfounders <- list()
    }

    ### Combined all selected variables into one list
    if (!length(KUconfounders)) {
      confounders <- Vconfounders
    }
    else if (!length(Vconfounders)) {
      confounders <- KUconfounders
    }
    else {
      confounders <- mapply (function(x,y) c(x, y), Vconfounders, KUconfounders, SIMPLIFY = FALSE)
    }
  }

  if (length(confounders))
    names(confounders) <- colnames(data[,(q + 1):(p + q), drop = FALSE])

  return(confounders)
}

### Build the list of true confounding variables for T-nodes
### From an adjacency matrix
get.true.conf.set <- function(Adj, p, u, q = 0, r = 0, offset = p + q + r) { # offset = nb.VTnodes + r
  if (u == 0) {
    return(list())
  }
  if (all(dim(Adj) == c(u, p))) {
    KU_T.Adj <- Adj
  }
  else {
    KU_T.Adj <- Adj[(p + q + 1):(p + q + u), (q + 1):(p + q)]
  }
  ConfSet <- lapply(1:p, FUN = function(j) {
    TK <- KU_T.Adj[,j] * (1:u)
    if (any((Kj <- TK > 0)))
      return(TK[Kj] + offset) # True confounders for T node j
    NULL
  })
  names(ConfSet) <- colnames(KU_T.Adj)
  return(ConfSet)
}

### Build the list of parent V-nodes for T-nodes
### From an adjacency matrix
get.true.variant.set <- function(Adj, p, q) {
  if (q == 0) {
    return(list())
  }
  if (all(dim(Adj) == c(q, p))) {
    V_T.Adj <- Adj
  }
  else {
    V_T.Adj <- Adj[1:q, (q + 1):(p + q)]
  }
  VSet <- lapply(1:p, FUN = function(j) {
    VT <- V_T.Adj[,j] * (1:q)
    if (any((Vj <- VT > 0)))
      return(VT[Vj])
    NULL
  })
  names(VSet) <- colnames(V_T.Adj)
  return(VSet)
}

### Build the list of parent T-nodes for each T-node
### From an adjacency matrix
get.true.parent.genes <- function(Adj, q = 0, p = NCOL(Adj) - q) {
  if (p == 0) {
    return(list())
  }
  Adj <- Adj[(q+1):(p+q), (q+1):(p+q)]
  TSet <- lapply(1:p, FUN = function(j) {
    TT <- Adj[,j] * (1:p)
    if (any((Tj <- TT > 0)))
      return(TT[Tj])
    NULL
  })
  names(TSet) <- colnames(Adj)
  return(TSet)
}

# No need of documentation; just a proxy to be replaced by 'get.conf.matrix'
get.confs <- function (data,
                       cov.pool = NULL,
                       measure = "correlation",
                       blocksize = min(pool.size, NCOL(data), 2000),
                       apply.qval = TRUE,
                       selection_fdr = 0.05,
                       lambda = NULL,
                       pi0.method = "smoother",
                       alpha = 0.05,
                       offset = NULL) { # column number of the first confounder in 'data'minus one
  if (is.null(cov.pool))
    return(list())
  pool.size <- NCOL(cov.pool)
  data <- as.matrix(data)
  if (is.null(offset))
    offset <- NCOL(data)
  confounders <- get.conf.matrix(data = data,
                                 cov.pool = cov.pool,
                                 measure = measure,
                                 blocksize = min(blocksize, pool.size),
                                 apply.qval = apply.qval,
                                 selection_fdr = selection_fdr,
                                 lambda = lambda,
                                 pi0.method = pi0.method,
                                 alpha = alpha)$sig.asso.covs
  if (!length(confounders)) {
    return(list())
  }
  else {
    if (is.list(confounders)) {
      Index <- 1:pool.size
      confounders <- lapply(confounders,
                            FUN = function(selected) {
                              col.indices <- if (any(selected)) Index[selected] + offset
                              col.indices
                            })
    }
    else {
      confounders <- as.list(confounders + offset)
    }
    return(confounders)
  }
}
