# An internal function to generate the code to check MRGN arguments
# Also saves inputs and format some arguments for the MRGN algorithm
check.mrgn.args <- function() {
  expression({
    # ============================================================================
    # Saving inputs & Checking arguments
    # ----------------------------------------------------------------------------
    # Data
    data <- as.data.frame(data)
    n <- NROW(data)
    m <- NCOL(data)

    # Check validity and consistency of block sizes
    stopifnot(is.numeric(q), is.numeric(p), is.numeric(r), is.numeric(u))
    stopifnot(q >= 0, q <= m, p >= 0, p <= m, r >= 0, r <= m, u >= 0, u < m, p + q + u + r <= m)
    if (m > p + q + r + u)
      data <- data[, 1:(p + q + r + u)]
    m <- p + q + r + u
    nb.nodes <- p + q + r

    # Save original labels (if any) or label columns
    Labels <- colnames(data)
    if (is.null(Labels)) {
      Vlabels <- if (q > 0) paste0('V', 1:q)
      Tlabels <- if ((p + r) > 0)  paste0('T', 1:(p + r))
      Clabels <- if (u > 0) paste0('PC', 1:u)
      Labels <- c(Vlabels, Tlabels, Clabels)
      colnames(data) <- Labels
    }
    else {
      Labels <- make.unique(Labels)
      Vlabels <- if (q > 0) Labels[1:q]
      Tlabels <- if ((p + r) > 0) Labels[(q + 1):nb.nodes]
      Clabels <- if (u > 0) Labels[(nb.nodes + 1):(nb.nodes + u)]
      colnames(data) <- Labels
    }

    # Check the consistency of the argument 'confounders' with data and its blocks
    if (!missing(confounders) && !is.null(confounders)) {
      if (!is.list(confounders))
        stop ("The argument 'confounders' must be a list")
      if (!all(sapply(confounders, FUN = function(x) {is.numeric(x) || is.null(x)})))
        stop ("Each element of 'confounders' must be a numeric vector")
      n.con <- length(confounders)
      if (n.con == 0) {
        if (verbose & u > 0)
          warning(paste0("      # ", u, " confounders supplied, but none is associated with a node. \n"))
        confounders <- NULL
      }
      else if (n.con == 1) {
        stopifnot(all(!is.na(confounders[[1]])), all(floor(confounders[[1]]) == confounders[[1]]))
        # if (any(confounders[[1]] > q & confounders[[1]] <= nb.nodes) | any(confounders[[1]] > m))
        if (any(confounders[[1]] > m))
          stop ("Column index out ot bound in 'confounders'")
        confounders[1:p] <- confounders[1]
      }
      else if (n.con == p) {
        confs.range <- unlist(confounders, recursive = TRUE, use.names = FALSE)
        stopifnot(all(!is.na(confs.range)), all(floor(confs.range) == confs.range))
        confs.range <- range(confs.range)
        # if ((confs.range[1] > q & confs.range[1] <= nb.nodes) | confs.range[2] > m)
        if(confs.range[2] > m)
          stop ("Column index out of bound in 'confounders'")
        rm(confs.range)
      }
      else
        stop (paste0("The argument 'confounders' must be a list of length one or ",
                     p, ". \n"))
    }
    else {
      if (verbose & u > 0)
        warning(paste0("      # ", u, " confounders supplied, but none is associated with a node. \n"))
      confounders <- NULL
    }
    Adj <- as.matrix(Adj)

    # Check the supplied adjacency matrix
    stopifnot(is.adjacency.matrix(Adj))
    if (NCOL(Adj) %in% c(nb.nodes, m)) {
      if (NCOL(Adj) == m)
        Adj <- Adj[1:nb.nodes, 1:nb.nodes]
    }
    else
      stop("The size of argument 'Adj' is not consistent with arguments 'data', 'p', 'q', 'r', and 'u'.")

    # Orient all V -- T edges as V --> T
    Adj[(1+q):nb.nodes, 1:q] <- 0
    dimnames(Adj) <- list(Labels[1:nb.nodes], Labels[1:nb.nodes])

    # Save the input adjacency matrix
    Adj0 <- Adj

    # Check FDR control arguments if required
    stopifnot(is.character(FDRcontrol))
    FDRcontrol <- FDRcontrol[1]
    stopifnot(FDRcontrol %in% c(p.adjust.methods, "qvalue"))
    stopifnot (is.numeric(alpha))
    alpha <- alpha[1]
    stopifnot (alpha > 0, alpha < .5)
    if (identical(FDRcontrol, "qvalue")) { # qvalue method
      stopifnot (is.numeric(fdr))
      fdr <- fdr[1]
      stopifnot (fdr > 0, fdr < .5)
      stopifnot(all(pi0.meth %in% c("bootstrap", "smoother")))
      stopifnot (lambda.step > 0, lambda.step < .1)
    }
    pi0.meth <- pi0.meth[1]
  })
}
