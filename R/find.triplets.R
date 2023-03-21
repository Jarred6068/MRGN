# Routine for 'enumerate.triplets.seq' (which enumerates triplets sequentially)
# Find triplets involving a particular T-node
# This function is a slightly modified version of 'enumerate.trios.i'
find.triplets.i <- function(i, Adj, Tlabels) {
  # Binary vector indicating T-nodes associated with 'Ti'
  Tlabelj <- Adj[i, Tlabels]

  # number of T-nodes that have an edge with 'Ti'
  pi <- sum(Tlabelj)

  # return a matrix with zero row if 'Ti' does not have edges with at least 2 other T-nodes
  if (pi < 2) {
    res <- matrix(NA, ncol = 3, nrow = 0)
  }
  else {
    # Doublets of T-nodes both associated with 'Ti'
    res <- t(combn(Tlabels[as.logical(Tlabelj)], 2))
    # Order the doublets if many
    if (NROW(res) > 1) {
      # Using order (a, b); thanks to Bandita
      res <- res[order(res[,1], res[,2]), , drop = FALSE]
    }
    # Get triplets
    res <- cbind(i, res)
  }

  # Names columns and return
  colnames(res) <- c("Ti", "Tj", "Tk")
  return(res)
}
