#' @name update.adjacency.matrix
#' @title Update an adjacency matrix
#'
#' @description This function update the adjacency matrix of a graph based on
#' trio analysis. It is used to perform Step \code{2.2} of the \code{MRGN}
#' algorithm. This is not a user level function.
#'
#' @param Adj numeric, an adjacency matrix
#'
#' @param q integer scalar, number of genetic variants in \code{Adj},
#' not required when \code{solve.conflicts = FALSE}.
#'
#' @param trio.set numeric, a matrix where each row indicates the column numbers
#' of variables forming a trio.
#'
#' A 3-column matrix \code{trio.set} corresponds to trios involving genetic
#' variant(s), i.e. each row of \code{trio.set} indicates the column number
#' of a genetic variant (first element), and the column numbers of two
#' expressions/phenotypes (second and third columns).
#'
#' A 4-column matrix \code{trio.set} corresponds to triplets involving only
#' expressions/phenotypes, no genetic variant(s), i.e. each row of \code{trio.set}
#' indicates the column numbers of three T-nodes, and the type of triplet
#' (last column). See \link{enumerate.triplets} for details on triplet types.
#'
#' @param inferred.models a character vector of length the number of rows in
#' \code{trio.set}. Each elements of \code{inferred.models} must be one of
#' \code{M0.1}, \code{M0.2}, \code{M1.1}, \code{M1.2}, \code{M2.1}, \code{M2.2},
#' \code{M3}, \code{M4}, or \code{Other}. See Badsha and Fu (2019) for the
#' definitions of these model topologies.
#'
#' @param solve.conflicts logical, should a resolution of conflicts be attempted?
#' If \code{FALSE}, edges are updated sequentially and conflicts are not noticed.
#' This makes the final inference depends on the order in which trios are analysed.
#'
#' If \code{TRUE} (the default), all edges are updated once. This put conflicts
#' into evidence and requires a method to deal with each type of errors.
#'
#' @param method character, only used if \code{solve.conflicts = TRUE}.
#' The method to be used to solve conflicts. The available \code{method} are:
#'
#' \code{naive}: conflicts about the presence of an edge are solve by including
#' an edge when at least one trio analysis inferred an edge; conflicts about an
#' edge direction are solve by letting the edge undirected.
#' \code{add_method.s}: add further method(s).
#'
#' @param added.edges character vector indicating all edges (and directions)
#' previously added into the network. Only used when \code{solve.conflicts = TRUE}.
#'
#' @param dropped.edges character vector indicating all edges (and directions)
#' previously dropped from the network (or inferred as non existent). Only used
#' when \code{solve.conflicts = TRUE}.
#'
#' @param cl a cluster object, created by one of the packages \code{parallel} and \code{snow}.
#' If \code{NULL}, the registered default cluster is used. Note that the latter can be
#' \code{NULL} too, in which case, no parallel computation is performed.
#'
#' @details \code{update.adjacency.matrix} takes an adjacency matrix, a set
#' \code{trio.set} of trios involving genetic variant(s) or not, and the underlying
#' structures under the extended Mendelian Randomization Principle, and update
#' the adjacency matrix accordingly.
#'
#' The registered default cluster is found using \code{parallel::getDefaultCluster()}.
#'
# Internal function to direct graph edges
update.adjacency.matrix <- function (Adj,
                                     q,
                                     trio.set,
                                     inferred.models,
                                     solve.conflicts = TRUE,
                                     method = "naive",
                                     added.edges = NULL,
                                     dropped.edges = NULL,
                                     cl = NULL, chunk.size = NULL) {
  # Check if 'Adj' is a valid adjacency matrix and 'trio.set' is a matrix
  stopifnot(is.adjacency.matrix(Adj), is.matrix(trio.set))
  # Un-comment the following lines if pushed to a user level function
  # stopifnot(is.character(inferred.models))
  # stopifnot(NROW(trio.set) != length(inferred.models))

  # Set default cluster
  if (is.null(cl))
    cl <- parallel::getDefaultCluster()

  # If 'trio.set' has four columns, then it contains triplets of only T-nodes
  if (NCOL(trio.set) == 4) {
    # In this case, the only allowed values in 'inferred.models'
    # are "M1.1", "M2.1", and "Other".
    # Un-comment the following line if pushed to a user level function
    # stopifnot(all(inferred.models %in% c("M1.1", "M2.1", "Other")))

    # Use the child function for triplets
    return(
      update_adjacency.matrix_triplet (Adj = Adj,
                                       q = q,
                                       triplet.set = trio.set,
                                       inferred.models = inferred.models,
                                       solve.conflicts = solve.conflicts[1],
                                       method = method[1],
                                       added.edges = added.edges,
                                       dropped.edges = dropped.edges,
                                       cl = cl, chunk.size = chunk.size)
    )
  }
  # Un-comment the following lines if pushed to a user level function
  # stopifnot(NCOL(trio.set) == 3) # Ensure 'trio.set' has only three columns
  # stopifnot(all(inferred.models %in% c('M0.1', 'M0.2', 'M1.1', 'M1.2',
  #                                   'M2.1', 'M2.2', 'M3', 'M4', 'Other')))

  # Use a special function to track conflicts if 'solve.conflicts' is TRUE
  if (solve.conflicts[1]) {
    return(
      update_adjacency.matrix_solve.conflicts.track (Adj = Adj,
                                                     q = q,
                                                     trio.set = trio.set,
                                                     inferred.models = inferred.models,
                                                     added.edges = added.edges,
                                                     dropped.edges = dropped.edges,
                                                     method = method[1],
                                                     cl = cl, chunk.size = chunk.size)
    )
  }
  else {

    ### UPDATE BELOW TO ACCOUNT FOR ARGUMENTS 'added.edges' AND 'dropped.edges'

    # Save the input Adjacency matrix
    Adj.old <- Adj

    # A vector to record the addition of new edges
    edges.new <- numeric(length(inferred.models))

    #=============================================================================
    #------------------------------- M0.1 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M0.1'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Update edge directions
      Adj[cbind(subtrio.set[,1], subtrio.set[,2])] <- 1 # Yes (V --> T1) # useless, V --> T1 is always 1 for a trio
      Adj[cbind(subtrio.set[,1], subtrio.set[,3])] <- 0 # No  (V --> T2)
      Adj[cbind(subtrio.set[,2], subtrio.set[,3])] <- 0 # No  (T1 --> T2)
      Adj[cbind(subtrio.set[,3], subtrio.set[,2])] <- 0 # No  (T2 --> T1)

      # Check for new edge(s) between nodes (V, T)
      edges.new[index] <- !Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])]
    }

    #=============================================================================
    #------------------------------- M0.2 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M0.2'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Update edge directions
      Adj[cbind(subtrio.set[,1], subtrio.set[,2])] <- 0 # No  (V --> T1)
      Adj[cbind(subtrio.set[,1], subtrio.set[,3])] <- 1 # Yes (V --> T2)
      Adj[cbind(subtrio.set[,2], subtrio.set[,3])] <- 0 # No  (T1 --> T2)
      Adj[cbind(subtrio.set[,3], subtrio.set[,2])] <- 0 # No  (T2 --> T1)

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- !Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])]
    }

    #=============================================================================
    #------------------------------- M1.1 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M1.1'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Update edge directions
      Adj[cbind(subtrio.set[,1], subtrio.set[,2])] <- 1 # Yes (V --> T1) # useless, V --> T1 is always 1 for a trio
      Adj[cbind(subtrio.set[,1], subtrio.set[,3])] <- 0 # No  (V --> T2)
      Adj[cbind(subtrio.set[,2], subtrio.set[,3])] <- 1 # Yes (T1 --> T2)
      Adj[cbind(subtrio.set[,3], subtrio.set[,2])] <- 0 # No  (T2 --> T1)

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        ((Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
            Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0)
    }

    #=============================================================================
    #------------------------------- M1.2 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M1.2'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Update edge directions
      Adj[cbind(subtrio.set[,1], subtrio.set[,2])] <- 0 # No  (V --> T1)
      Adj[cbind(subtrio.set[,1], subtrio.set[,3])] <- 1 # Yes (V --> T2)
      Adj[cbind(subtrio.set[,2], subtrio.set[,3])] <- 0 # No  (T1 --> T2)
      Adj[cbind(subtrio.set[,3], subtrio.set[,2])] <- 1 # Yes (T2 --> T1)

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
    }

    #=============================================================================
    #------------------------------- M2.1 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M2.1'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Update edge directions
      Adj[cbind(subtrio.set[,1], subtrio.set[,2])] <- 1 # Yes (V --> T1) # useless, V --> T1 is always 1 for a trio
      Adj[cbind(subtrio.set[,1], subtrio.set[,3])] <- 0 # No  (V --> T2)
      Adj[cbind(subtrio.set[,2], subtrio.set[,3])] <- 0 # No  (T1 --> T2)
      Adj[cbind(subtrio.set[,3], subtrio.set[,2])] <- 1 # Yes (T2 --> T1)

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
    }

    #=============================================================================
    #------------------------------- M2.2 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M2.2'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Update edge directions
      Adj[cbind(subtrio.set[,1], subtrio.set[,2])] <- 0 # No  (V --> T1)
      Adj[cbind(subtrio.set[,1], subtrio.set[,3])] <- 1 # Yes (V --> T2)
      Adj[cbind(subtrio.set[,2], subtrio.set[,3])] <- 1 # Yes (T1 --> T2)
      Adj[cbind(subtrio.set[,3], subtrio.set[,2])] <- 0 # No  (T2 --> T1)

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
    }

    #=============================================================================
    #-------------------------------- M3 trios -----------------------------------
    #=============================================================================
    index <- inferred.models == 'M3'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Update edge directions
      Adj[cbind(subtrio.set[,1], subtrio.set[,2])] <- 1 # Yes (V --> T1) # useless, V --> T1 is always 1 for a trio
      Adj[cbind(subtrio.set[,1], subtrio.set[,3])] <- 1 # Yes (V --> T2)
      Adj[cbind(subtrio.set[,2], subtrio.set[,3])] <- 0 # No  (T1 --> T2)
      Adj[cbind(subtrio.set[,3], subtrio.set[,2])] <- 0 # No  (T2 --> T1)

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])] == 0
    }

    #=============================================================================
    #------------------------------- M4 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M4'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Update edge directions
      Adj[cbind(subtrio.set[,1], subtrio.set[,2])] <- 1 # Yes V --> T1 # useless, V --> T1 is always 1 for a trio
      Adj[cbind(subtrio.set[,1], subtrio.set[,3])] <- 1 # No (V --> T2)
      Adj[cbind(subtrio.set[,2], subtrio.set[,3])] <- 1 # No (T1 --> T2)
      Adj[cbind(subtrio.set[,3], subtrio.set[,2])] <- 1 # No (T2 --> T1)

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0

    }
    return(list(Adj = Adj,
                nb.new_edges = 0, # TO BE DEFINED
                nb.dropped_edges = 0, # TO BE DEFINED
                new.edges = edges.new, # TO BE CORRECTED
                added.edges = added.edges,
                dropped.edges = dropped.edges,
                nb.conflicts = 0))
  }
}
