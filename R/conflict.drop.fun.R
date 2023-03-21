# Internal routine for 'update_adjacency.matrix_solve.conflicts'
# Deciding to drop an edge or not when there are conflicts from trio analysis results
# 'ijk' is a 3-vector of indices: c(ind_row, ind_col, ind_row <= q)
# 'Adj' is the Adjacency matrix, after adding edges based on trio analysis,
# 'method' is the handling method

conflict.drop.fun <- function(ijk, Adj, method) {
  if (ijk[3]) {
    # We have a variant, what to do?
    if (identical(method, "naive")) {
      return(FALSE) # ??????????????
    }
    else {
      return(FALSE) # ??????????????
    }
  }
  else {
    # Do we have an edge even if we drop this direction?
    if (Adj[ijk[2], ijk[1]]) {
      # Yes: leave the edge undirected (do not drop any direction)
      return(FALSE)
    }
    else {
      # No: what to do?
      if (identical(method, "naive")) {
        return(TRUE) # ??????????????
      }
      else {
        return(TRUE) # ??????????????
      }
    }
  }
}
