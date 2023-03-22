#'
#' @name enumerate.trios
#' @title Form trios involving a target genetic variant
#' @description Enumerate all trios involving a target genetic variant based on
#' an adjacency matrix. Not a user level function.
#'
#' @export enumerate.trios
#'
#' @usage
#' enumerate.trios (i = 1:q,
#'                  Adj,
#'                  q = length(i),
#'                  p = NROW(Adj) - q)
#'
#'
#' @param i integer vector, position of the target \code{V}-nodes (genetic
#' variants) in the adjacency matrix \code{Adj}: \code{i} must satisfy
#' \code{1 ≤ i ≤ q} (see the definition of \code{q} below). For convenience,
#' \code{i} is coerced to the integer type using the function \link{as.integer}.
#' Defaults to \code{i = 1:q} (i.e. all genetic variants).
#' @param Adj numeric, an adjacency matrix, i.e., a binary symmetric matrix
#' (graph skeleton). The first \code{q} columns of \code{Adj} represent
#' \code{V}-nodes, and the last \code{p} columns of \code{Adj} represent
#' \code{T}-nodes (expressions/phenotypes). Only the upper triangular part of
#' \code{Adj} is actually used.
#' @param p integer, number of \code{T}-nodes (non-instrumental variables).
#' @param q integer, number of \code{V}-nodes (instrumental variables).
#'
#' @details When \code{i} has more than one elements, it is first reduced to a
#' vector of unique elements using \link{unique}, and trio enumeration is performed
#' for each unique element.
#'
#' @return A list of matrices, one for each unique element of \code{i}. If the
#' input adjacency matrix \code{Adj} has non-\code{NULL} colnames (as returned
#' by \link{colnames}), then the elements of the list are named after the
#' corresponding columns (from the first \code{q} columns). Otherwise, they are
#' named as 'Vi' where 'i' are the unique elements of \code{i}.
#'
#' For each matrix, each row gives the column numbers of variables forming a
#' trio, the first giving the genetic variant (\code{i}), and the last two
#' giving phenotypes (named 'Tj', 'Tk'). If no trio can be formed for a
#' given genetic variant, the associated matrix has zero row.
#'
#' @seealso \link{enumerate.triplets}.
#
####################################################
# List all trios involving a target genetic variant
enumerate.trios <- function (i = 1:q, Adj, q = length(i), p = NROW(Adj) - q,
                             cl = NULL, chunk.size = NULL) {
  if (!missing(i)) {
    if (length(i) < 1)
      stop("'i' must be a non-NULL vector of integers.")
    i <- as.integer(unique(i))
    if (any(c(i < 1, i > q)))
      stop("The vector 'i' must satisfy '1 ≤ i ≤ q' element-wise.")
  }

  # Only use the upper triangular part of Adj
  Adj[lower.tri(Adj)] <- Adj[lower.tri(t(Adj))]

  # Labels (numbers) for all T-nodes
  Tlabels <- (1 + q):(p + q)

  # Wrap 'enumerate.trios.i' over 'i'
  res <- matteLapply(i,
                     FUN = enumerate.trios.i,
                     Adj = Adj, p = p, q = q, Tlabels = Tlabels,
                     cl = cl, chunk.size = chunk.size)

  # Name and return the obtained list
  Vnames <- colnames(Adj)
  if (!is.null(Vnames)) {
    names(res) <- Vnames[i]
  }
  else {
    names(res) <- paste0('V', i)
  }
  return(res)
}

####################################################
# 'enumerate.trios.i' is the workhorse for 'enumerate.trios'
enumerate.trios.i <- function (i, Adj, p, q, Tlabels) {
  # Binary vector indicating T-nodes associated with 'Vi'
  Tlabeli <- Adj[i, Tlabels]

  # number of T nodes that have an edge with 'Vi'
  pi <- sum(Tlabeli)

  # Terminate if 'Vi' has no edge with T-nodes
  if (pi == 0)
    return(NULL)

  # Doublets of T-nodes both associated with 'Vi'
  Apairs <- if (pi > 1) {
    t(combn(Tlabels[as.logical(Tlabeli)], 2))
  }

  # Pairs involving a T-node non-associated to Vi (with an associated one)
  Npairs <- if (pi < p) {
    form.doublets(Tlabeli = Tlabeli, Tlabels = Tlabels,
                Adj = Adj, p = p, q = q)
  }

  # Bind the two groups of doublets
  res <- rbind(Apairs, Npairs)

  # If there is in fact no doublet, return a matrix with zero row
  if (NROW(res) == 0)
    res <- matrix(NA, ncol = 3, nrow = 0)
  else { # Otherwise, format the doublets and get trios
    # Order the doublets if many
    if (NROW(res) > 1) {
      # Using order (a, b); thanks to Bandita
      res <- res[order(res[,1], res[,2]), , drop = FALSE]
    }

    # Replace column numbers (in dataset) by T numbers for T-nodes
    #res <- res - q

    # Get trios
    res <- cbind(i, res)
  }

  # Names columns and return
  colnames(res) <- c("Vi", "Tj", "Tk")
  return(res)
}

####################################################
# A routine for 'enumerate.trios.i': form all doublets involving each a T-node
# associated with Vi and a T-node non associated with 'Vi'
form.doublets <- function(Tlabeli, Tlabels, Adj, p, q) {
  # Get a list of doublets involving each T-node
  res <- lapply(Tlabels[as.logical(Tlabeli)], FUN = form.doublets.j,
                Tlabeli = Tlabeli, Tlabels=Tlabels, Adj=Adj, p=p, q=q)

  # Unlist the result
  res <- unlist(res, recursive = TRUE)

  # 'unvec' the result if any
  if (!is.null(res))
    return(unique(matrix(res, ncol = 2, byrow = TRUE)))

  # Return a matrix with zero row if no result
  return(matrix(NA, ncol = 2, nrow = 0))
}

# A child function for the T-node numbered 'j'
form.doublets.j <- function(j, Tlabeli, Tlabels, Adj, p, q) {
  # Binary vector indicating T-nodes associated with 'Tj'
  assoc <- Adj[j, (q+1):(q+p)]

  # Eliminating T-nodes already in Tlabeli (T-nodes directly associated with 'Vi')
  # to avoid duplicating T-nodes already accounted for in pairwise combinations
  assoc <- assoc * (1 - Tlabeli)

  # Terminate if no new association exists
  if (sum(assoc) == 0)
    return(NULL)

  # Form doublets (columns)
  resj <- rbind(j, Tlabels[as.logical(assoc)])

  # Sort each column
  # resj <- apply(resj, MARGIN = 2, FUN = sort) # (useful???)
  # Guess not useful, maybe counter productive

  # Return doublets involving 'Tj' (applying 'vec')
  return(c(resj))
}
