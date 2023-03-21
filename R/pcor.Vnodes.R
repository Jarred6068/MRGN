
# Compute partial correlations, possibly conditioning only on V-nodes (all)
partial.cor.given.Vnodes <- function(data, p, q, onlyVnodes = TRUE,
                                     parallel = FALSE,
                                     cl =  parallel::getDefaultCluster(),
                                     chunk.size = NULL) {
  if (onlyVnodes) {
    ### Take each pair of T nodes
    resultat <- matrix(1, nrow = p, ncol = p)
    indices <- which(lower.tri(resultat, diag = FALSE), arr.ind = TRUE) + q
    if (parallel) {
      corvalues <- MRGN:::matteApply (X = indices, MARGIN = 1, FUN = function(index) {
        MRGN:::partial.cor(data[, c(index, 1:q)])[1,2]
      },
      cl = cl,
      chunk.size = chunk.size, simplify = TRUE)
    }
    else {
      corvalues <- apply(indices, MARGIN = 1, FUN = function(index) {
        MRGN:::partial.cor (data[, c(index, 1:q)])[1,2]
      })
    }
    indices <- indices - q
    resultat[indices] <- corvalues
    resultat[indices[,2:1]] <- corvalues
  }
  else {
    resultat <- MRGN:::partial.cor(data)[(q+1):(q+p), (q+1):(q+p)]
  }
  return(resultat)
}
