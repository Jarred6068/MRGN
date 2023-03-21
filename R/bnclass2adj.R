#' From \code{bn} class objects to adjacency matrix
#'
#' A function to convert the list of edges from bnlearn into an adjacency matrix.
#'
#' @export bnclass2adj

bnclass2adj <- function(bnfit, pattern = "G", nodenames = NULL) {
  if (inherits(bnfit, what = "bn")) {
  Adj <- cbind(as.numeric(sapply(bnfit$arcs[,1], FUN = function(x) {
    stringi::stri_split_fixed(x, pattern = pattern)[[1]][2]})),
    as.numeric(sapply(bnfit$arcs[,2], FUN = function(x) {
      stringi::stri_split_fixed(x, pattern = pattern)[[1]][2]})))
  Adj <- edges2adj(Adj)
  if (is.null(nodenames)) {
    nodenames <- paste0(pattern, 1:NROW(Adj))
  }
  else {
    stopifnot(length(nodenames) == NROW(Adj))
  }
  dimnames(Adj) <- list(nodenames, nodenames)
  return(Adj)
  }
  else {
    stop("Only objects of class 'bn' are currently handled.")
  }
}
