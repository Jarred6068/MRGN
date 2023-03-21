#' From edges to adjacency matrix
#'
#' A function to convert a list of edges into an adjacency matrix.
#'
#' @param edges a \code{matrix}, \code{data.frame} or \code{list}.
#' \code{edges} must be a numeric 2-column matrix or a data.frame if
#' \code{pattern} is NULL (or missing).
#'
#' @param pattern NULL or a \code{character}, indicates the labels of nodes.
#' For instance, \code{pattern = 'G'} means that nodes in \code{edges} are
#' labelled \code{Gi} where \code{i} is a positive integer number.
#'
#' @param diag logical, should non null values be allowed on the main diagonal
#' of the adjacency matrix? Defaults to \code{FALSE}.
#'
#' @details The fitting functions in the package \code{bnlearn} return list of
#' edges. The function \code{edges2adj} is meant to derive the adjacency matrix
#' corresponding to such a list.
#'
#' @export edges2adj
#'
#' @return an object of class \code{adjacency.matrix}.
#'
#' @seealso \link{is.adjacency.matrix}.
#
edges2adj <- function(edges, pattern = NULL, diag = FALSE) {
  classedge_list <- class(edges)
  if (missing(pattern) | is.null(pattern)) {
    if (any(c('data.frame', 'matrix') %in% classedge_list)) {
    edges <- as.matrix(edges)
    if (NCOL(edges) < 2) {
      stop("argument 'edges' must have two columns or more when 'pattern = NULL'.")
    }
    stopifnot(is.numeric(edges))
    indices <- edges[,1:2]
    }
    else if (is.list(edges)) {
      if (!all(c(sapply(edges, FUN = length) == 2,
                sapply(edges, FUN = is.numeric)))) {
        stop("an argument 'edges' of class 'list' must have all 'numeric' elements, each of length 2.")
      }
      indices <- t(sapply(edges,
                          FUN = function (index) {
                            as.numeric(index[1:2])
                          }))
    }
    else
      stop("argument 'edges' must be of class 'matrix', 'data.frame' or 'list'.")
  }
  else {
    if (any(c('data.frame', 'matrix') %in% classedge_list)) {
      if (NCOL(edges) < 2) {
        stop("an argument 'edges' of class 'data.frame' or 'matrix' must have two columns or more.")
      }
      indices <- t(apply(edges,
                         MARGIN = 1,
                         FUN = function (index) {
                           c(as.numeric(stringr::str_split(index[1],
                                                           pattern = pattern)[[1]][2]),
                             as.numeric(stringr::str_split(index[2],
                                                           pattern = pattern)[[1]][2]))
                         }))
    }
    else if ('list' %in% classedge_list) {
      if (!all(c(sapply(edges, FUN = length) == 2,
                sapply(edges, FUN = is.character)))) {
        stop("an argument 'edges' of class 'list' must have all 'character' elements, each of length 2.")
      }
      indices <- t(sapply(edges,
                         FUN = function (index) {
                           c(as.numeric(stringr::str_split(index[1],
                                                           pattern = pattern)[[1]][2]),
                             as.numeric(stringr::str_split(index[2],
                                                           pattern = pattern)[[1]][2]))
                         }))
    }
    else
      stop("argument 'edges' must be of class 'matrix', 'data.frame' or 'list'.")
  }
  nb.nodes <- max(indices)
  Adj <- matrix(0, nrow = nb.nodes, ncol = nb.nodes)
  Adj[indices] <- 1
  if (!diag)
    diag(Adj) <- 0
  structure(Adj, class = 'adjacency.matrix')
}
