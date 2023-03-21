#' @name enumerate.triplets
#' @title List all update-able trios involving only phenotypes
#' @description Enumerate update-able triplets involving only expressions/phenotypes
#'  (no genetic variant) based on an adjacency matrix. Not a user level function.
#'
#' @export enumerate.triplets
#'
#' @usage
#' enumerate.triplets (Adj,
#'                     return.1 = FALSE)
#'
#' @param Adj numeric, a binary matrix indicating adjacency and possibly edge
#' directions (i.e. a non-directed, or a partially directed graph). All the
#' \code{p} columns of \code{Adj} represent \code{T}-nodes (phenotypes). As such,
#' \code{Adj} will in general be an extract of a larger adjacency matrix (last
#' \code{p} rows and columns in this package). See Section \code{Details} for more
#' on the interpretation of the values \code{0/1} in \code{Adj} for this function.
#'
#' @param return.1 logical, should triplets of type (1) be returned as an
#' attribute named \code{triplet.1} (of class \code{matrix}) for the output?
#' Defaults to \code{FALSE}. See Section \code{Details} for more on type (1)
#' triplets.
#'
#' @param cl a cluster object, created by one of the packages \code{parallel} and \code{snow}.
#' If \code{NULL}, the registered default cluster is used. Note that the latter can be
#' \code{NULL} too, in which case, no parallel computation is performed.
#'
#' @details The binary matrix \code{Adj} represent a non-directed, or a partially
#' directed graph. For two nodes \code{Ti} and \code{Tj}, the matrix elements
#' \code{a_ij} and \code{a_ji} define both the presence of an edge and the
#' direction of the edge if any: an edge is present between nodes \code{Ti}
#' and \code{Tj} if \code{a_ij + a_ji > 0}. If \code{a_ij + a_ji = 2}, then the
#' edge is bi-directed (or equivalently, undirected). Otherwise, if \code{a_ij = 1},
#' then the edge goes from \code{Ti} to \code{Tj}, and if \code{a_ji = 1}, the
#' edge goes from \code{Tj} to \code{Ti}.
#'
#' The function only returns trios for which graph models can be further
#' identified under the extended interpretation of the Principle of Mendelian
#' Randomization (PMR). For triplets of strictly \code{T}-nodes with at
#' least one undirected edge, two general situations are possible: (A) the trio
#' has three edges, or (B) the trio has two edges. In situation (A), no further
#' inference is possible, and such triplets are NOT returned.
#'
#' In situation (B), we have three possibilities: (1) one directed edge, and the
#' parent \code{T}-node has two edges; (2) one directed edge, and the parent
#' \code{T}-node has only one edge; and (3) two undirected edges.
#'
#' In case (1), no further inference is possible under PMR. In the absence of an
#' additional genetic variant linked to one or many of the \code{T}-nodes, we
#' have three different but Markov equivalent graph candidates, i.e., these
#' graphs share the same set of conditional and marginal independence relations
#' (see Kvamme and Fu, 2022). Such triplets are NOT returned in the outputed
#' matrix. They are nevertheless returned as an attribute named \code{triplet.1}
#' (of class \code{matrix}) for the matrix output if \code{return.1 = TRUE}
#' (note that the default is \code{return.1 = FALSE}).
#'
#' In case (2), the parent \code{T}-node can be treated as a genetic variant.
#' The second edge of the triplet can be directed using a conditional independence
#' test between the parent \code{T}-node and the non directly-related \code{T}-node.
#' Such triplets are returned in the outputed matrix.
#'
#' In case (3), one of the edges of the triplet can potentially be directed using
#' a conditional independence test between the two \code{T}-nodes not directly
#' related by an edge. Here, the ability to direct the second edge depends on
#' the test result. Such triplets are returned in the outputed matrix.
#'
#' In the outputed matrix, the columns are named 'Ti', 'Tj', and 'Tk', and each
#' row indicates three \code{T}-nodes whose edges can potentially be updated
#' (i.e. directed). The \code{T}-nodes are arranged such that 'Ti' has two
#' edges. This way, any further inference only requires regressing 'Tj' (or 'Tk')
#' on 'Ti' and 'Tk' (or 'Ti'), and confounding variables, if any. If a triplet
#' is of type \code{2}, then 'Tj' is the parent node (i.e. 'Tj' points to 'Ti').
#'
#' @return A list with two elements:
#'
#' \item{triplets}{an \code{n} by \code{3} matrix. Each row of the output gives
#' the numbers of non-instrumental variables forming a triplet with two edges.
#' The columns of the matrix are named 'Ti', 'Tj', and 'Tk'. The \code{T}-nodes
#' are arranged such that "Ti" has two edges. If no update-able triplet can be
#' formed, the returned matrix has zero row.}
#'
#' \item{types}{a \code{n}-vector of integers. The i^th element indicates the
#' type of triplet corresponding to the i^th row the matrix \code{triplets}.}
#'
#' If \code{return.1 = TRUE}, then the output matrix \code{triplets} has an
#' attribute named \code{triplet.1} giving triplets of type (1). See Section
#' \code{Details} for more on type (1) triplets. The attribute \code{triplet.1}
#' is a \code{matrix}, and it will have zero row if no triplets of type (1) is
#' found.
#'
#' @seealso \link{enumerate.trios}.
#'
#
####################################################
# List update-able triplets involving strictly T-nodes
enumerate.triplets <- function (Adj, return.1 = FALSE, cl = NULL) {

  # Set default cluster
  if (is.null(cl))
    cl <- parallel::getDefaultCluster()

  # Parallelize when possible/required
  if (!is.null(cl)) {
    enumerate.triplets.parallel (Adj = Adj, return.1 = return.1, cl = cl)
  }

  # Labels (numbers) for all T-nodes
  p <- NCOL(Adj)
  Tlabels <- 1:p

  # Adjacency matrix of the un-directed graph
  UndirAdj <- ((Adj + t(Adj)) > 0) + 0

  # An indicator for T-nodes with at least two edges
  Nedges <- rowSums(UndirAdj)
  margins <- Nedges >= 2

  # List all triplets involving a T-node with at least two edges
  triplets <- lapply (X = Tlabels[margins], FUN = find.triplets.i,
                      Adj = UndirAdj, Tlabels = Tlabels)
  triplets <- do.call('rbind', triplets)

  # Terminate if no strictly T-nodes triplet can be formed
  if (length(dim(triplets)) < 2) {
    return(triplet.null (return.1))
  }

  # Get update-able triplets (and types of triplet) if any
  triplets <- t(apply(triplets, MARGIN = 1, FUN = function(x) {
    res <- find.triplet.type (Aijk = Adj[x, x],
                              labelsijk = Tlabels[x])
    if (is.null(res$types))
      return(rep(0, 4))
    return(c(res$triplets, res$types))
  }))

  # An indicator to remove null rows
  keep <- rowSums(triplets) > 0

  # Terminate if no strictly T-nodes triplet can be formed
  if (!any(keep)) {
    return(triplet.null (return.1))
  }

  # Remove null rows and extract 'types'
  triplets <- triplets[keep, , drop = FALSE]
  colnames(triplets) <- c("Ti", "Tj", "Tk", "types")
  types <- triplets[,4]

  # Indicator of type (2) and (3) triplets
  keep <- types > 1

  # Return
  if (return.1[1]) {
    triplet.1 <- if (!all(keep))
      triplets[!keep, 1:3, drop = FALSE]
    else
      triplet.null (FALSE)$triplets
    triplets <- triplets[keep, 1:3, drop = FALSE]
    attr(triplets, 'triplet.1') <- triplet.1
  }
  else {
    triplets <- if (any(keep)) {
      triplets[keep, 1:3, drop = FALSE]
    }
    else
      triplet.null (FALSE)$triplets
  }
  return(list(triplets = triplets, types = if (any(keep)) types[keep]))
}

# Parallelized version of 'enumerate.triplets'
enumerate.triplets.parallel <- function (Adj, return.1 = FALSE, cl) {

  # Labels (numbers) for all T-nodes
  p <- NCOL(Adj)
  Tlabels <- 1:p

  # Adjacency matrix of the un-directed graph
  UndirAdj <- ((Adj + t(Adj)) > 0) + 0

  # An indicator for T-nodes with at least two edges
  Nedges <- rowSums(UndirAdj)
  margins <- Nedges >= 2

  # List all triplets involving a T-node with at least two edges
  triplets <- parallel::parLapply (cl = cl,
                                   X = Tlabels[margins], fun = find.triplets.i,
                                   Adj = UndirAdj, Tlabels = Tlabels)
  triplets <- do.call('rbind', triplets)

  # Terminate if no strictly T-nodes triplet can be formed
  if (length(dim(triplets)) < 2) {
    return(triplet.null (return.1))
  }

  # Get update-able triplets (and types of triplet) if any
  triplets <- t(parallel::parApply (cl = cl,
                                    X = triplets,
                                    MARGIN = 1, FUN = function(x) {
    res <- find.triplet.type (Aijk = Adj[x, x],
                              labelsijk = Tlabels[x])
    if (is.null(res$types))
      return(rep(0, 4))
    return(c(res$triplets, res$types))
  }))

  # An indicator to remove null rows
  keep <- rowSums(triplets) > 0

  # Terminate if no strictly T-nodes triplet can be formed
  if (!any(keep)) {
    return(triplet.null (return.1))
  }

  # Remove null rows and extract 'types'
  triplets <- triplets[keep, , drop = FALSE]
  colnames(triplets) <- c("Ti", "Tj", "Tk", "types")
  types <- triplets[,4]

  # Indicator of type (2) and (3) triplets
  keep <- types > 1

  # Return
  if (return.1[1]) {
    triplet.1 <- if (!all(keep))
      triplets[!keep, 1:3, drop = FALSE]
    else
      triplet.null (FALSE)$triplets
    triplets <- triplets[keep, 1:3, drop = FALSE]
    attr(triplets, 'triplet.1') <- triplet.1
  }
  else {
    triplets <- if (any(keep)) {
      triplets[keep, 1:3, drop = FALSE]
    }
    else
      triplet.null (FALSE)$triplets
  }
  return(list(triplets = triplets, types = if (any(keep)) types[keep]))
}

# List update-able triplets involving strictly T-nodes
# Initial version of enumerate.triplets
# Uses a greedy search
# Kept to test the effectiveness of 'enumerate.triplets'
enumerate.triplets.greedy <- function (Adj, return.1 = FALSE) {
  # Labels (numbers) for all T-nodes
  p <- NCOL(Adj)
  Tlabels <- 1:p

  # An indicator to exclude isolated T-nodes
  margins <- (rowSums(Adj) + colSums(Adj)) > 0

  # Number of non completely isolated T-nodes
  nnodes <- sum(margins)

  # Terminate if no strictly T-nodes triplet can be formed
  if (nnodes < 3) {
    return(triplet.null (return.1))
  }

  # Special case of only one possible triplet
  if (nnodes == 3) {
    # Get triplet (and type of triplet) if any
    triplets <- find.triplet.type (Aijk = Adj[margins, margins],# Extract sub-adjacency matrix
                                   labelsijk = Tlabels[margins])# Pick the 3 labels

    # Return a 0 row matrix if triplet is of type (1)
    if (triplets$types == 1) {
      # Save 'triplet.1' as an attribute if asked for
      if (return.1[1]) {
        triplet.1 <- triplets$triplets
        triplets <- triplet.null (FALSE)
        attr(triplets$triplets, 'triplet.1') <- triplet.1
        return(triplets)
      }
      return(triplet.null (FALSE))
    }

    # Give a 0 row matrix as the attribute 'triplet.1' if required
    if (return.1[1]) {
      attr(triplets$triplets, 'triplet.1') <- triplet.null (FALSE)$triplets
    }

    # Return the triplet
    return(triplets)
  }

  # List all non-isolated triplets
  triplets <- combn(Tlabels[margins], 3)

  # Get update-able triplets (and types of triplet) if any
  triplets <- t(apply(triplets, MARGIN = 2, FUN = function(x) {
    res <- find.triplet.type (Aijk = Adj[x, x],
                              labelsijk = Tlabels[x])
    if (is.null(res$types))
      return(rep(0, 4))
    return(c(res$triplets, res$types))
  }))

  # An indicator to remove null rows
  keep <- rowSums(triplets) > 0

  # Terminate if no strictly T-nodes triplet can be formed
  if (!any(keep)) {
    return(triplet.null (return.1))
  }

  # Remove null rows and extract 'types'
  triplets <- triplets[keep, , drop = FALSE]
  colnames(triplets) <- c("Ti", "Tj", "Tk", "types")
  types <- triplets[,4]

  # Indicator of type (2) and (3) triplets
  keep <- types > 1

  # Return
  if (return.1[1]) {
    triplet.1 <- if (!all(keep))
      triplets[!keep, 1:3, drop = FALSE]
    else
      triplet.null (FALSE)$triplets
    triplets <- triplets[keep, 1:3, drop = FALSE]
    attr(triplets, 'triplet.1') <- triplet.1
  }
  else {
    triplets <- if (any(keep)) {
      triplets[keep, 1:3, drop = FALSE]
    }
    else
      triplet.null (FALSE)$triplets
  }
  return(list(triplets = triplets, types = if (any(keep)) types[keep]))
}

####################################################
# A routine check if a triplet is update-able
# 'Aijk' is a 3 by 3 binary matrix
# 'labelsijk' is 3-vector of numeric labels (numbers of the 3 T-nodes)

find.triplet.type <- function (Aijk, labelsijk) {
  # Terminate if no edge is present
  if(!sum(Aijk)) {
    return(triplet.null (FALSE))
  }

  # Presence of edges per node
  n_edges <- pmax(rowSums(Aijk), colSums(Aijk))

  # Terminate if any node is isolated from the two others
  if(any(n_edges == 0)) {
    return(triplet.null (FALSE))
  }

  # Count edge indicators
  Edge_c <- Aijk[lower.tri(Aijk)] + Aijk[upper.tri(Aijk)]

  # Terminate if all present edges are already directed
  if(all(Edge_c <= 1)) {
    return(triplet.null (FALSE))
  }

  # Determine the presence of edge(s)
  #p_e <- (Edge_c > 0) + 0
  # Determine the number of edge(s)
  #n_e <- sum(p_e)
  # Terminate if n_e != 2
  if (sum(Edge_c > 0) != 2) {
    return(triplet.null (FALSE))
  }

  # Indicator for: number of edges per node = 2 or not
  n_2e <- n_edges == 2

  # Put the node with two edges in the first position if it is not yet
  if(!n_2e[1]) {
    labelsijk <- c(labelsijk[n_2e], labelsijk[!n_2e])
  }

  # Form a line matrix with node numbers
  triplets <- matrix(labelsijk, nrow = 1, ncol = 3)
  colnames(triplets) <- c("Ti", "Tj", "Tk")

  # Number of directed edges
  n_d <- sum(Edge_c == 1) # = 0 or 1 (case n_d = 2 already ruled out)

  # If n_d = 0, return triplet of type (3)
  if (n_d == 0) {
    return(list(triplets = triplets, types = 3))
  }

  # If n_d == 1, determine if we have a triplet of type (1) or (2)
  # Test if the node with two edges is the parent node
  if (sum(Aijk[n_2e,]) == 2) { # type = '1'
    return(list(triplets = triplets, types = 1))
  }

  # Otherwise, type = '2'
  if(!n_2e[1]) { # Reorganize the matrix to match 'labelsijk'
    i <- sum(n_2e * (1:3))
    Aijk <- cbind(Aijk[,n_2e], Aijk[, !n_2e])
    Aijk <- rbind(Aijk[n_2e,], Aijk[!n_2e,])
  }

  # Ensure that 'Tj' is the parent node (put it in position 2)
  if (sum(Aijk[,2]) > 0)
    triplets[,2:3] <- triplets[,3:2]
  return(list(triplets = triplets, types = 2))
}

####################################################
# A routine to return an empty matrix of triplets
# 'return.1' is a logical scalar
triplet.null <- function (return.1) {
  triplets <- matrix(NA, ncol = 3, nrow = 0)
  colnames(triplets) <- c("Ti", "Tj", "Tk")
  if (return.1[1]) {
    attr(triplets, 'triplet.1') <- triplets
  }
  return(list(triplets = triplets, types = NULL))
}
