# Wrap 'analyse.trio.i' over a set of trios
analyse.trio.set <- function(trio.set, nb.trios = NROW(trio.set),
                             alpha, FDRcontrol, fdr, lambda, lambda.step,
                             pi0.meth = "bootstrap",
                             data, confounders, p, q, use.perm, gamma, is.CNA,
                             nperms, verbose, cl = NULL, chunk.size = NULL) {
  switch(FDRcontrol,
         none = { # No correction
           trio.analysis <- matteApply(trio.set,
                                       MARGIN = 1,
                                       FUN = analyse.trio.i,
                                       data = data, confounders = confounders, p = p, q = q,
                                       use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                       alpha = alpha, nperms = nperms, half1 = FALSE, verbose = verbose > 1,
                                       cl = cl, chunk.size = chunk.size)
           trio.analysis <- do.call('rbind', trio.analysis)
           return(trio.analysis)
         },
         qvalue = { # qvalue method
           # Half trio analysis
           all.stats <- matteApply(trio.set,
                                   MARGIN = 1,
                                   FUN = analyse.trio.i,
                                   data = data, confounders = confounders, p = p, q = q,
                                   use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                   nperms = nperms, half1 = TRUE, verbose = verbose > 1,
                                   cl = cl, chunk.size = chunk.size)
           all.stats <- do.call('rbind', all.stats)
           Inferred.Model0 <- all.stats[,8]

           # FDR control: apply qvalue correction
           if (is.null(lambda)) {
             lambda <- seq(min(c(unlist(all.stats[,1:6]), lambda.step), na.rm = T),
                           max(c(unlist(all.stats[,1:6]), lambda.step), na.rm = T), lambda.step)
           }

           Qvalues <- catch.conditions (adjust.q (unlist(all.stats[,1:6]), fdr = fdr,
                                                  lambda = lambda, pi0.meth = pi0.meth))$value
           if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
             if (identical(pi0.meth, "smoother")) {
               if (verbose) {
                 cat("      #         * q value correction with 'pi0.meth = smoother' failled, trying 'pi0.meth = bootstrap'  \n")
               }
               Qvalues <- catch.conditions (adjust.q (unlist(all.stats[,1:6]), fdr = fdr,
                                                      lambda = lambda, pi0.meth = "bootstrap"))$value
             }
             if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
               if ((ll <- length(lambda)) > 1) {
                 Qvalues <- catch.conditions (adjust.q (unlist(all.stats[,1:6]),
                                                        fdr = fdr, lambda = lambda[ll],
                                                        pi0.meth = pi0.meth))$value
                 if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                   Qvalues <- catch.conditions (adjust.q (unlist(all.stats[,1:6]),
                                                          fdr = fdr, lambda = lambda[1],
                                                          pi0.meth = pi0.meth))$value
                 }
               }
               if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                 if (verbose) {
                   cat("            -- 'qvalue::qvalue' failed. \n")
                 }
                 Qvalues <- list(qvalue = c(all.stats[,1:6]),
                                 significant = c(all.stats[,1:6]) <= alpha)
               }
             }
           }

           all.stats <- cbind(matrix(Qvalues$significant, nrow = nb.trios, byrow = FALSE),
                              matrix(Qvalues$qvalue, nrow = nb.trios, byrow = FALSE),
                              Minor.freq = all.stats[,7])
           all.stats = as.data.frame(all.stats)
           colnames(all.stats) <- c("b11","b12", "b21","b22", "V1:T1", "V1:T2", "pb11",
                                    "pb12", "pb21","pb22","pV1:T1","pV1:T2", "Minor.freq")

           # Infer model structure (Second half trio analysis)
           all.stats$Inferred.Model <- c(matteApply(all.stats, MARGIN = 1, FUN = MRGN::class.vec,
                                                    cl = cl, chunk.size = chunk.size))

           # Also return the inferred model with no correction
           all.stats$Inferred.Model0 <- Inferred.Model0

           return(all.stats)
         },
         { # Bonferroni and all other Multiple Comparisons adjustment methods
           # Half trio analysis
           all.stats <- matteApply(trio.set,
                                   MARGIN = 1,
                                   FUN = analyse.trio.i,
                                   data = data, confounders = confounders, p = p, q = q,
                                   use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                   nperms = nperms, half1 = TRUE, verbose = verbose > 1,
                                   cl = cl, chunk.size = chunk.size)
           all.stats <- do.call('rbind', all.stats)
           Inferred.Model0 <- all.stats[,8]

           # Bonferroni correction
           p.adj <- stats::p.adjust(unlist(all.stats[,1:6]), method=FDRcontrol)
           p.adj <- matrix(p.adj, nrow = nb.trios, byrow = FALSE)
           all.stats <- cbind(p.adj <= alpha, p.adj,
                              Minor.freq = all.stats[,7])
           all.stats = as.data.frame(all.stats)
           colnames(all.stats) <- c("b11","b12", "b21","b22", "V1:T1", "V1:T2", "pb11",
                                    "pb12", "pb21","pb22","pV1:T1","pV1:T2", "Minor.freq")

           # Infer model structure (Second half trio analysis)
           all.stats$Inferred.Model <- c(matteApply(all.stats, MARGIN = 1, FUN = MRGN::class.vec,
                                                    cl = cl, chunk.size = chunk.size))

           # Also return the inferred model with no correction
           all.stats$Inferred.Model0 <- Inferred.Model0

           return(all.stats)
         })
}

# Call infer.trio for trio analysis
analyse.trio.i <- function (col.indices,
                            data, p, q, confounders,
                            use.perm = TRUE, gamma = 0.05, is.CNA = FALSE,
                            alpha = 0.01, nperms = 100, half1 = FALSE, verbose = FALSE) {

#  if (all(col.indices[1:3] %in% c(1, 292, 319))) {

 #   browser(text = "target trio", condition = NULL, expr = TRUE, skipCalls = 0L)

  #}

  # For each trio, take the union of the confounders for the two T nodes
  conf.set.trio <- c(if (col.indices[2] <= q + p) confounders[[col.indices[2] - q]], # Take confounder index if T-node
                     if (col.indices[3] <= q + p) confounders[[col.indices[3] - q]]) # Result is NULL if I&C node
  conf.set.trio <- unique(conf.set.trio)

  # Find and remove duplicates (V-nodes, T-nodes), if any
  VTduplicates <- conf.set.trio %in% col.indices
  if (any(VTduplicates))
    conf.set.trio <- conf.set.trio[!VTduplicates]

  # Call infer.trio
  trio.res <- infer.trio(data[, c(col.indices, conf.set.trio)],
                         use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                         alpha = alpha, nperms = nperms, verbose = verbose)
  if (half1)
    trio.res[,7:14]
  else
    trio.res
}
