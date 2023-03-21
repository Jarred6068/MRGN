# Internal child function of "update.adjacency.matrix
update_adjacency.matrix_triplet <- function (Adj,
                                             q,
                                             triplet.set,
                                             inferred.models,
                                             solve.conflicts = TRUE,
                                             method = "naive",
                                             added.edges = NULL,  # NOT YET USED HERE
                                             dropped.edges = NULL, # NOT YET USED HERE
                                             cl = NULL, chunk.size = NULL) {
  # Use the special function if 'solve.conflicts' is TRUE
  if (solve.conflicts[1] & FALSE) {
    # Complete (create the function 'update_adjacency.matrix_triplet_solve_conflicts')
    # or Remove
    return(
      update_adjacency.matrix_triplet_solve_conflicts (
        q = q,
        Adj = Adj,
        triplet.set = triplet.set,
        inferred.models = inferred.models,
        method = method)
    )
  }

  # Save the input Adjacency matrix
  Adj.old <- Adj

  #=============================================================================
  #------------------------------- M1.1 trios ----------------------------------
  #=============================================================================
  index <- inferred.models == 'M1.1'
  if (any(index)) {
    subtriplet.set <- triplet.set[index,, drop = FALSE]

    # Triplets of type 2
    if (any(type2 <- subtriplet.set[,4] == 2)) {
      Adj[cbind(subtriplet.set[type2,1], subtriplet.set[type2,3])] <- 1 # Yes (Ti --> Tk) # useless, Ti --> Tk is always 1 for a type 2 triplet
      Adj[cbind(subtriplet.set[type2,3], subtriplet.set[type2,1])] <- 0 # No  (Tk --> Ti)
    }
  }

  #=============================================================================
  #------------------------------- M2.1 trios ----------------------------------
  #=============================================================================
  index <- inferred.models == 'M2.1'
  if (any(index)) {
    subtriplet.set <- triplet.set[index,, drop = FALSE]
    Adj[cbind(subtriplet.set[,3], subtriplet.set[,1])] <- 1 # Yes (Ti --> Tk) # useless, Ti --> Tk is always 1 for a type 2 triplet
    Adj[cbind(subtriplet.set[,1], subtriplet.set[,3])] <- 0 # No  (Tk --> Ti)
    Adj[cbind(subtriplet.set[,2], subtriplet.set[,1])] <- 1 # Yes (Ti --> Tk) # useless, Ti --> Tk is always 1 for a type 2 triplet
    Adj[cbind(subtriplet.set[,1], subtriplet.set[,2])] <- 0 # No  (Tk --> Ti)
  }
  return(list(Adj = Adj))
}
