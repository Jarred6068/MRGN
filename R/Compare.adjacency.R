#' A generic 'compare' function
#'
#' Compare objects inheriting from a specific class.
#'
#' The default method in the package \code{MRGN}
#' compares two adjacency matrices. See \link{compare.adjacency.matrix}.
#'
#' @export compare

setGeneric("compare", function(object, Ref, ...) {
  compare.adjacency.matrix (object = object, Ref = Ref, ...)
})


#' Compare Adjacency Matrices
#'
#' Compute performance measures for an inferred adjacency matrix,
#' comparing it to a reference (TRUE) adjacency matrix.
#'
#' @export compare.adjacency.matrix
#' @exportS3Method compare adjacency.matrix
#'
#' @param object numeric matrix, the inferred adjacency matrix.
#'
#' @param Ref numeric matrix, the reference adjacency matrix.
#'
#' @param ... additional arguments passed to or from other methods.
#'
#' @param directions logical, should the comparison consider edge directions?
#' Defaults to \code{directions = FALSE} which means that the adjacency matrices
#' are considered symmetric (coerced to).
#'
#' @return A named list with elements:
#' \item{Precision}{numeric, the proportion of true edges/directions in all inferred edges/directions.}
#' \item{Recall}{numeric, the proportion of rightly inferred edges/directions in all true edges/directions.}
#' \item{TP}{numeric, number of rightly inferred edges/directions in the inferred adjacency matrix \code{object}.}
#' \item{FP}{numeric, number of wrong edges/directions in the inferred adjacency matrix \code{object}.}
#' \item{nb.true.edges}{numeric (only available when \code{directions = FALSE}), number of edges in the reference adjacency matrix \code{Ref}.}
#' \item{nb.true.directions}{numeric (only available when \code{directions = TRUE}), number of directions in the reference adjacency matrix \code{Ref}.}
#'
compare.adjacency.matrix <- function(object, Ref, directions = FALSE, ...) {
  stopifnot(is.adjacency.matrix(object),
            is.adjacency.matrix(Ref),
            is.logical(directions))
  if (!directions[1]) {
    if (!isSymmetric.matrix(object))
      object <- (object + t(object) > 0) + 0
    if (!isSymmetric.matrix(Ref))
      Ref <- (Ref + t(Ref) > 0) + 0
    scale <- 2
  }
  else {
    scale <- 1
  }
  dAdj = object - Ref

  # TRUE positives
  TP <- sum(dAdj[Ref>0] == 0)

  # FALSE positives
  FP <- sum(dAdj[Ref==0] > 0)

  # TP / inferred Positives (Precision)
  nb.inferred <- sum(object)
  Precision <- TP / nb.inferred

  # Recall
  nb.reference <- sum(Ref)
  Recall <- TP / nb.reference

  if (directions[1]) {
    list(Precision = Precision, Recall = Recall,
         TP = TP/scale, FP = FP/scale,
         nb.true.directions = nb.reference/scale)
  }
  else {
    list(Precision = Precision, Recall = Recall,
         TP = TP/scale, FP = FP/scale,
         nb.true.edges = nb.reference/scale)
  }
}


