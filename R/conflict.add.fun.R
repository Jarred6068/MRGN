# Internal routine for 'update_adjacency.matrix_solve.conflicts'
# Deciding to add an edge or not when there are conflicts from trio analysis results
# 'ijk' is a 3-vector of indices: c(ind_row, ind_col, ind_row <= q)
# 'Adj' is the Adjacency matrix,
# 'method' is the handling method

conflict.add.fun <- function(ijk, Adj, method) {
  if (ijk[3]) {
    # We have a variant, what to do?
    if (identical(method, "naive")) {
      return(TRUE) # ??????????????
    }
    else {
      return(TRUE) # ??????????????
    }
  }
  else {
    # Do we have an edge even if we do not add this direction?
    if (Adj[ijk[2], ijk[1]]) {
      # Yes: leave the edge undirected (add the edge)
      return(TRUE)
    }
    else {
      # No: what to do?
      if (identical(method, "naive")) {
        return(FALSE) # ??????????????
      }
      else {
        return(FALSE) # ??????????????
      }
    }
  }
}
