# Wrap 'analyse.triplet.i' over a set of triplets
analyse.triplet.set <- function (triplet.set, nb.triplets = NROW(triplet.set),
                                 alpha, FDRcontrol, fdr, lambda, lambda.step, pi0.meth = "bootstrap",
                                 data, confounders, p, q, verbose,
                                 cl = NULL, chunk.size = NULL) {
  switch(FDRcontrol,
         qvalue = { # qvalue method
           # Half trio analysis
           all.stats <- matteApply(triplet.set,
                                   MARGIN = 1,
                                   FUN = analyse.triplet.i,
                                   data, p = p, q = q,
                                   confounders = confounders,
                                   alpha = alpha, half1 = TRUE, verbose = verbose > 1,
                                   cl = cl, chunk.size  = chunk.size)

           # FDR control: apply qvalue correction
           if (is.null(lambda)) {
             lambda <- seq(min(c(all.stats, lambda.step), na.rm = T),
                           max(c(all.stats, lambda.step), na.rm = T), lambda.step)
           }

           # There can be PROBLEM in 'pi0est': there is bug at lines 28-29 (debugging line numbers)
           Qvalues <- catch.conditions (adjust.q (c(all.stats), fdr = fdr, lambda = lambda,
                                                  pi0.meth = pi0.meth))$value

           if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
             if (identical(pi0.meth, "smoother")) {
               if (verbose) {
                 cat("      #         * q value correction with 'pi0.meth = smoother' failled, trying 'pi0.meth = bootstrap'  \n")
               }
               Qvalues <- catch.conditions (adjust.q (c(all.stats), fdr = fdr,
                                                      lambda = lambda,
                                                      pi0.meth = "bootstrap"))$value
             }
             if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
               if ((ll <- length(lambda)) > 1) {
                 Qvalues <- catch.conditions (adjust.q (c(all.stats), fdr = fdr, lambda = lambda[ll],
                                                        pi0.meth = pi0.meth))$value
                 if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                   Qvalues <- catch.conditions (adjust.q (c(all.stats), fdr = fdr, lambda = lambda[1],
                                                          pi0.meth = pi0.meth))$value
                 }
               }
               if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                 if (verbose) {
                   cat("            -- 'qvalue::qvalue' failed. \n")
                 }
                 Qvalues <- list(qvalue = c(all.stats),
                                 significant = c(all.stats) <= alpha)
               }
             }
           }

           model <- character(length = nb.triplets)
           model[!Qvalues$significant] <- "M1.1"
           if (any(Qvalues$significant)) {
             type <- triplet.set[4][Qvalues$significant]
             model[Qvalues$significant][type == 2] <- "M2.1"
             model[Qvalues$significant][type != 2] <- "Other"
           }
           all.stats <- cbind(Qvalues$significant,
                              Qvalues$qvalue,
                              model)
           all.stats = as.data.frame(all.stats)
           colnames(all.stats) <- c("rej", "pval", "Inferred.Model")

           return(all.stats)

         },
         bonferroni = { # Bonferroni method
           # Half trio analysis
           all.stats <- matteApply(triplet.set,
                                   MARGIN = 1,
                                   FUN = analyse.triplet.i,
                                   data, p = p, q = q,
                                   confounders = confounders, alpha = alpha,
                                   half1 = TRUE, verbose = verbose > 1,
                                   cl = cl, chunk.size  = chunk.size)

           # Bonferroni correction
           p.adj <- stats::p.adjust(c(all.stats), method="bonferroni")
           rej <- p.adj <= alpha
           model <- character(length = nb.triplets)
           model[!rej] <- "M1.1"
           if (any(rej)) {
             type <- triplet.set[4][rej]
             model[rej][type == 2] <- "M2.1"
             model[rej][type != 2] <- "Other"
           }
           all.stats <- cbind(rej,
                              p.adj,
                              model)

           all.stats = as.data.frame(all.stats)
           colnames(all.stats) <- c("rej", "pval", "Inferred.Model")

           return(all.stats)

         },
         {
           triplet.analysis <- matteApply(triplet.set,
                                          MARGIN = 1,
                                          FUN = analyse.triplet.i,
                                          data, p = p, q = q,
                                          confounders = confounders,
                                          alpha = alpha, half1 = FALSE,
                                          verbose = verbose > 1,
                                          cl = cl, chunk.size  = chunk.size)
           triplet.analysis <- do.call('rbind', triplet.analysis)
           return(triplet.analysis)

         })
}

# Call infer.triplet for triplet analysis
analyse.triplet.i <- function (col.indices,
                               data, p, q, confounders,
                               alpha = 0.01,
                               half1 = FALSE,
                               verbose = FALSE) {
  # For each triplet, take the union of the confounders for the three T nodes
  conf.set.triplet <- c(if (col.indices[1] <= q + p) confounders[[col.indices[1] - q]], # Take confounder index if T-node
                        if (col.indices[2] <= q + p) confounders[[col.indices[2] - q]], # Result is NULL if I&C node
                        if (col.indices[3] <= q + p) confounders[[col.indices[3] - q]])
  conf.set.triplet <- unique(conf.set.triplet)
  infer.triplet(data[, c(col.indices[-4], conf.set.triplet)],
                type = col.indices[4],
                alpha = alpha,
                half1 = half1,
                verbose = verbose)
}
