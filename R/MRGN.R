
# Keep the function simple for now
# Assume we know (as input) the set of confounding variables for each T node
# Know the type of each confounding variable (actual confounder, or one of intermediate or common child variable)
# Dealing with missing values ???

#' @name MRGN
#'
#' @title Infer a causal network using the MRGN algorithm
#'
#' @description Infer a causal network/graph with directed and undirected edges from observational data.
#'
#' @param confounders a list of length \code{p} giving the set of known confounders
#' for each phenotype, as returned by \link{get.conf} with the option
#' \code{return.for.trios = FALSE}.
#'
#'
#' @seealso \link{infer.trio} for inferring edges in small networks of one
#'  genetic variant and only two genes.

#' @export MRGN
#' @import parallel

# We have four blocks of columns in data: V, T, (I&C) and confounders

# Post filtering based on marginal tests ???


# Keep track of all added edges (directions) so that
# at each iteration, I can check for conflicts with the current adds
# but also previous adds

MRGN <- function (data, # input n-by-m data matrix: 'q' Variants, 'p' Phenotypes, r Intermediate&Common_Child and 'u' confounders
                  p, q, r, u = m - p - q - r,
                  confounders = list(), # a list of length p giving the set of known confounders for each phenotype
                  Adj, # a (p+q+r)-by-(p+q+r) adjacency matrix of an undirected acyclic graph (zeros on diagonal).
                  is.CNA = FALSE, # logical indicating if a genetic variant is a copy number alteration.
                  alpha = 0.01, # in (0, .5) type I error rate for individual tests, in the open (0, .5)
                  use.perm = TRUE, # See infer.trio
                  gamma = 0.05, #
                  nperms = 100, #
                  FDRcontrol = c("bonferroni", "qvalue", "none"),
                  fdr = 0.01, # in (0, .5) False discovery rate for all trio analyses at one iteration, in the open (0, .5)
                  lambda = NULL, # The value of the tuning parameter for q-value estimation
                  lambda.step = 0.05, # in (0, .1); Used to initialize "lambda" when lambda = NULL
                  pi0.meth = c("bootstrap", "smoother"),
                  solve.conflicts = TRUE, # Just used for testing, to be removed?
                  method = "naive",
                  parallel = FALSE,
                  cl = parallel::getDefaultCluster(),
                  chunk.size = NULL, # scalar number; number of invocations of fun or FUN in one chunk; a chunk is a unit for scheduling.
                  verbose = 0L) {
  # ============================================================================
  # Save the call to MRGN
  # ----------------------------------------------------------------------------
  mrgncl <- match.call()

  # ============================================================================
  # Parallelize computations or not?
  # ----------------------------------------------------------------------------
  stopifnot(is.logical(parallel))
  if (parallel[1]) {
    if (is.null(cl)) {
      # Use option cl.cores to choose an appropriate cluster size
      cl <- parallel::makeCluster(getOption("cl.cores", 2L))
    }
    if (!is.null(chunk.size)) {
      stopifnot(is.numeric(chunk.size))
      stopifnot(all(chunk.size > 0))
      chunk.size <- ceiling(chunk.size[1])
    }
  }
  else {
    cl <- NULL
  }

  # ============================================================================
  # Saving inputs & Checking arguments
  # ----------------------------------------------------------------------------
  eval(check.mrgn.args())

  # ============================================================================
  # Step 2: List trios involving genetic variants and direct related edges
  # ----------------------------------------------------------------------------
  # Step 2.1: Identify target trios
  # Generate trio list & Bind trios identifiers as rows
  if (verbose) {
    cat("\n      # Generating an initial list of trios involving each a genetic variant. \n")
  }
  trio.set <- enumerate.trios(i = 1:q, Adj = Adj,
                              cl = cl, chunk.size = chunk.size)
  trio.set <- do.call('rbind', trio.set)

  # Number of initial trios
  nb.trios <- NROW(trio.set)

  if (verbose) {
    cat(paste0("          * ",
               if (nb.trios == 0) "zero trio"
               else if (nb.trios == 1) "one trio"
               else paste0(nb.trios, " trios"),
               " generated. \n"))
  }

  if (nb.trios > 0) {
    # Step 2.2: Determine the directed structure for each trio
    if (verbose)
      cat("      # Inferring the causal network for each trio involving a genetic variant. \n")
    trio.analysis <- analyse.trio.set (trio.set, nb.trios = nb.trios, alpha = alpha, FDRcontrol = FDRcontrol,
                                       fdr = fdr, lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                                       data = data, confounders = confounders, p = p, q = q,
                                       use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                       nperms = nperms, verbose = verbose,
                                       cl = cl, chunk.size = chunk.size)

    # Update the adjacency matrix
    if (any(upd.trios <- trio.analysis$Inferred.Model != 'Other')) {
      if (verbose) {
        cat(if (sum(upd.trios) == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                            "direction for one trio involving a genetic variant. \n")
            else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                        sum(upd.trios), " trios involving genetic variant(s). \n"))
      }
      new.Adj <- update.adjacency.matrix (Adj = Adj,
                                          q = q,
                                          trio.set = trio.set,
                                          inferred.models = trio.analysis$Inferred.Model,
                                          solve.conflicts = solve.conflicts,
                                          method = method,
                                          cl = cl, chunk.size = chunk.size)
      added.edges <- new.Adj$added.edges
      dropped.edges <- new.Adj$dropped.edges
      if (verbose & solve.conflicts) {
        if (new.Adj$nb.conflicts)
          cat(paste0("          * ",
                     if (new.Adj$nb.conflicts == 1) "one conflict detected"
                     else paste0(new.Adj$nb.conflicts, " conflicts detected"),
                     " and handled with 'method = ", method, "'. \n"))
      }

      # Logical vector indicating trios in which new edges were added
      new.edges <- new.Adj$new.edges == 1
      nb.new_edges <- new.Adj$nb.new_edges
      new.Adj <- new.Adj$Adj

      # Check if new edges were added
      if (any(new.edges)) {
        if (verbose & nb.new_edges) {
          cat(paste0("          * ",
                     if (nb.new_edges == 1) "one new edge"
                     else paste0(nb.new_edges, " new edges"),
                     " added into the network. \n"))
        }

        # Check the presence of any undirected edge in the network
        undirected.edges <- rowSums((new.Adj + t(new.Adj)) == 2)
        if (any(undirected.edges > 1)) {
          new.trio.set <- enumerate.trios.new (Adj = new.Adj,
                                               old.Adj = Adj,
                                               new.edges = new.edges,
                                               trio.set = trio.set,
                                               p = p, q = q,
                                               cl = cl, chunk.size = chunk.size)
          if (verbose & (nb.new.trios <- NROW(new.trio.set))) {
            cat(paste0("          * ",
                       if (nb.new.trios == 1) "one new trio"
                       else paste0(nb.new.trios, " new trios"),
                       " generated. \n"))
          }
        }
        else {
          nb.new.trios <- 0
        }
      }
      else {
        nb.new.trios <- 0
      }
      Adj <- new.Adj
    }
    else {
      nb.new.trios <- 0
      added.edges <- dropped.edges <- NULL
    }

    # Run again Step 2.2 if new trios were generated
    iter <- 1
    maxiter <- choose(nb.nodes, 3)
    while (nb.new.trios > 0 & iter <= maxiter) {
      # iteration counter
      iter <- iter + 1

      # Determine the directed structure for each new trio
      if (verbose) {
        cat(paste0("      # Inferring the causal network for each new trio involving a genetic variant (iteration ",
                   iter, "). \n"))
      }
      new.trio.analysis <- analyse.trio.set (new.trio.set, nb.trios = nb.new.trios,
                                             alpha = alpha, FDRcontrol = FDRcontrol,
                                             fdr = fdr, lambda = lambda,
                                             lambda.step = lambda.step, pi0.meth = pi0.meth,
                                             data = data, confounders = confounders, p = p, q = q,
                                             use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                             nperms = nperms, verbose = verbose,
                                             cl = cl, chunk.size = chunk.size)

      # Save the additional trios in 'trio.set'
      trio.set <- rbind(trio.set, new.trio.set)

      # Save the additional trios analysis result in 'trio.analysis'
      trio.analysis <- rbind(trio.analysis, new.trio.analysis) # (not useful at this point)

      # Update the adjacency matrix
      if (any(upd.trios <- new.trio.analysis$Inferred.Model != 'Other')) {
        if (verbose) {
          cat(if (sum(upd.trios) == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                              "direction for one trio involving a genetic variant. \n")
              else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                          sum(upd.trios), " trios involving genetic variant(s). \n"))
        }
        new.Adj <- update.adjacency.matrix (Adj = Adj,
                                            q = q,
                                            trio.set = new.trio.set,
                                            inferred.models = new.trio.analysis$Inferred.Model,
                                            solve.conflicts = solve.conflicts,
                                            method = method,
                                            added.edges = added.edges,
                                            dropped.edges = dropped.edges,
                                            cl = cl, chunk.size = chunk.size)
        added.edges <- new.Adj$added.edges
        dropped.edges <- new.Adj$dropped.edges
        if (verbose & solve.conflicts) {
          if (new.Adj$nb.conflicts)
            cat(paste0("          * ",
                       if (new.Adj$nb.conflicts == 1) "one conflict detected"
                       else paste0(new.Adj$nb.conflicts, " conflicts detected"),
                       " and handled with 'method = ", method, "'. \n"))
        }

        # Logical vector indicating trios in which new edges were added
        new.edges <- new.Adj$new.edges == 1
        nb.new_edges <- new.Adj$nb.new_edges
        new.Adj <- new.Adj$Adj

        # Check if new edges were added
        if (any(new.edges)) {
          if (verbose & nb.new_edges) {
            cat(paste0("          * ",
                       if (nb.new_edges == 1) "one new edge"
                       else paste0(nb.new_edges, " new edges"),
                       " added into the network. \n"))
          }

          # Check the presence of any undirected edge in the network
          undirected.edges <- (new.Adj + t(new.Adj)) == 2
          undirected.edges <- undirected.edges[lower.tri(undirected.edges)]

          if (any(undirected.edges)) {
            new.trio.set <- enumerate.trios.new (Adj = new.Adj,
                                                 old.Adj = Adj,
                                                 new.edges = c(rep(FALSE, nb.trios),
                                                               new.edges),
                                                 trio.set = trio.set,
                                                 p = p, q = q,
                                                 cl = cl, chunk.size = chunk.size)
            nb.trios <- nb.trios + nb.new.trios
            if (verbose & (nb.new.trios <- NROW(new.trio.set))) {
              cat(paste0("          * ",
                         if (nb.new.trios == 1) "one new new trio"
                         else paste0(nb.new.trios, " new trios"),
                         " generated. \n"))
            }
          }
          else {
            nb.trios <- nb.trios + nb.new.trios
            nb.new.trios <- 0
          }
        }
        else {
          nb.trios <- nb.trios + nb.new.trios
          nb.new.trios <- 0
        }
        Adj <- new.Adj
      }
      else {
        nb.trios <- nb.trios + nb.new.trios
        nb.new.trios <- 0
      }
    }
  }
  else {
    trio.analysis <- added.edges <- dropped.edges <- NULL
  }

  # Check the presence of any undirected edge in the network
  undirected.edges <- rowSums((Adj + t(Adj)) == 2)
  if (any(undirected.edges > 0) & p + r >= 3) {

    # ======================================================================================
    # Step 3: List triplets involving strictly T-nodes, and direct related edges when possible
    # --------------------------------------------------------------------------------------
    # Step 3.1: List triplets involving strictly T-nodes and that can be updated
    if (verbose) {
      cat("      # Generating an initial list of triplets involving only T-nodes. \n")
    }
    triplet.set <- enumerate.triplets (Adj[Tlabels, Tlabels, drop = FALSE])
    # Number of initial triplets
    nb.triplets <- length(triplet.set$types)
    if (!nb.triplets) {
      if (verbose) {
        cat("      # No updateable triplet involving only T-nodes found. \n")
      }
      triplet.set <- triplet.analysis <- NULL
    }
    else {
      triplet.set <- cbind(triplet.set$triplets + q, type = triplet.set$types)

      if (verbose) {
        cat(paste0("          * ",
                   if (nb.triplets == 1) "one updateable triplet"
                   else paste0(nb.triplets, " updateable triplets"),
                   " generated. \n"))
      }

      # Step 3.2: Determine the directed structure for each triplet
      if (verbose)
        cat("      # Inferring the causal network for each triplet involving only T-nodes. \n")
      triplet.analysis <- analyse.triplet.set (triplet.set, nb.triplets = nb.triplets,
                                               alpha = alpha, FDRcontrol = FDRcontrol, fdr = fdr,
                                               lambda = lambda, lambda.step = lambda.step,
                                               pi0.meth = pi0.meth,
                                               data = data, confounders = confounders,
                                               p = p, q = q, verbose = verbose,
                                               cl = cl, chunk.size = chunk.size)

      # Update the adjacency matrix
      if (any(upd.triplets <- triplet.analysis$Inferred.Model != 'Other')) {
        if (verbose) {
          cat(if (sum(upd.triplets) == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                                 "direction for one triplet involving only T-nodes. \n")
              else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                          sum(upd.triplets), " triplets involving only T-nodes. \n"))
        }
        new.Adj <- update.adjacency.matrix (Adj = Adj,
                                            q = q,
                                            trio.set = triplet.set,
                                            inferred.models = triplet.analysis$Inferred.Model,
                                            solve.conflicts = solve.conflicts,
                                            method = method,
                                            added.edges = added.edges,
                                            dropped.edges = dropped.edges,
                                            cl = cl, chunk.size = chunk.size)

        #added.edges <- new.Adj$added.edges
        #dropped.edges <- new.Adj$dropped.edges

        # Logical vector indicating trios in which new edges were added
        new.edges <- new.Adj$new.edges == 1
        new.Adj <- new.Adj$Adj

        # RIGHT HERE Repeat triplet listing and update.adjacency.matrix

        #
        Adj <- new.Adj

      }
    }
  }
  else {
    triplet.set <- triplet.analysis <- NULL
  }

  # ======================================================================================
  list(Adj = structure(Adj, class = 'adjacency.matrix'),
       Adj0 = structure(Adj0, class = 'adjacency.matrix'),
       trio.set = trio.set,
       trio.analysis = trio.analysis,
       triplet.set = triplet.set,
       triplet.analysis = triplet.analysis,
       added.edges = added.edges,
       dropped.edges = dropped.edges,
       cl = cl, chunk.size = chunk.size,
       call = mrgncl)
}

# Alternative to check if new edges were added
#edges.new <-  (Adj.new + t(Adj.new) > 0) - (Adj + t(Adj) > 0) > 0
# any(edges.new[lower.tri(edges.new)])
