# A routine to check if an object is an adjacency matrix or not
#' @name is.adjacency.matrix
#' @title Test if an object is an adjacency matrix
#' @description Function to test if \code{object} is an adjacency matrix or not.
#'
#' @usage is.adjacency.matrix (object)
#'
#' @param object any \code{R} object.
#'
#' @details An adjacency matrix is restricted in this context to represent
#' acyclic graph objects (directed or not). An adjacency matrix is thus required
#' to be a square matrix with binary elements (\code{0} or \code{1}), and only
#' zeros on the main diagonal.
#'
#' @return logical indicating if \code{object} is an adjacency matrix or not.
#'
#' @examples
#'
# \dontrun{
#'
#' Adj <- matrix(c(0, 1, 1,
#'                 0, 0, 0,
#'                 1, 1, 0),
#'                 nrow = 3, byrow = TRUE)
#'
#' is.adjacency.matrix (Adj)
#'
#' is.adjacency.matrix (Adj + diag(1))
#'
# }
#'
is.adjacency.matrix <- function (object) {
  if (!is.matrix(object))
    return(FALSE)
  d <- dim(object)
  if (d[1L] != d[2L])
    return(FALSE)
  all(c(d[1] == d[2],
        all(diag(object) == 0),
        all ((1 - object) == (!object))))
}
