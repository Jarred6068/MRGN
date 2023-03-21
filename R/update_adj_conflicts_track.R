# Update all edges at once (not sequentially) based on trio analysis
# Internal child function of "update.adjacency.matrix
update_adjacency.matrix_solve.conflicts.track <-
  function (Adj, q, trio.set, inferred.models,
            added.edges = NULL, dropped.edges = NULL,
            method = "naive", cl = NULL, chunk.size = NULL) {
    # Save the input Adjacency matrix and added/dropped edges
    Adj.old <- Adj
    added.edges.old <- added.edges
    dropped.edges.old <- dropped.edges

    # A vector to record the addition of new edges
    edges.new <- numeric(length(inferred.models))

    # Initialize two 2-columns matrices to track edge adds (or confirmation) and edge drops
    edge_add <- edge_drop <- NULL # no need matrices; NULL is OK

    # Two columns matrices to keep track of trios associated to drops/adds in edge_add  / edge_drop
    trio_add.set <- trio_drop.set <- NULL # NOT YET USED, will need if post-filtering is included in conflict resolutions

    #=============================================================================
    #------------------------------- M0.1 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M0.1'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Set edge updating plan given trio analysis
      edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,2]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,3]),
                         cbind(subtrio.set[,2], subtrio.set[,3]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))

      # Also store the indices of the trios corresponding to each edge drop/add (may be used to confirm adds/drops)
      trio_add.set <- rbind(trio_add.set,
                            subtrio.set)
      trio_drop.set <- rbind(trio_drop.set,
                             subtrio.set,
                             subtrio.set,
                             subtrio.set)

      # Check for new edge(s) between nodes (V, T)
      edges.new[index] <- !Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])]
    }

    #=============================================================================
    #------------------------------- M0.2 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M0.2'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]
      edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,3]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,2]),
                         cbind(subtrio.set[,2], subtrio.set[,3]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))

      # Also store the indices of the trios corresponding to each edge drop/add (may be used to confirm adds/drops)
      trio_add.set <- rbind(trio_add.set,
                            subtrio.set)
      trio_drop.set <- rbind(trio_drop.set,
                             subtrio.set,
                             subtrio.set,
                             subtrio.set)

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- !Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])]
    }

    #=============================================================================
    #------------------------------- M1.1 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M1.1'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]
      edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,2]),
                        cbind(subtrio.set[,2], subtrio.set[,3]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,3]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))

      # Also store the indices of the trios corresponding to each edge drop/add (may be used to confirm adds/drops)
      trio_add.set <- rbind(trio_add.set,
                            subtrio.set,
                            subtrio.set)
      trio_drop.set <- rbind(trio_drop.set,
                             subtrio.set,
                             subtrio.set)

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
      edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,3]),
                        cbind(subtrio.set[,3], subtrio.set[,2]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,2]),
                         cbind(subtrio.set[,2], subtrio.set[,3]))

      # Also store the indices of the trios corresponding to each edge drop/add yup!
      trio_add.set <- rbind(trio_add.set,
                            subtrio.set,
                            subtrio.set)
      trio_drop.set <- rbind(trio_drop.set,
                             subtrio.set,
                             subtrio.set)

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
      edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,2]),
                        cbind(subtrio.set[,3], subtrio.set[,2]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,3]),
                         cbind(subtrio.set[,2], subtrio.set[,3]))

      # Also store the indices of the trios corresponding to each edge drop/add (may be used to confirm adds/drops)
      trio_add.set <- rbind(trio_add.set,
                            subtrio.set,
                            subtrio.set)
      trio_drop.set <- rbind(trio_drop.set,
                             subtrio.set,
                             subtrio.set)

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
      subtrio.set <- trio.set[index,, drop = FALSE]
      edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,3]),
                        cbind(subtrio.set[,2], subtrio.set[,3]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,2]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))

      # Also store the indices of the trios corresponding to each edge drop/add (may be used to confirm adds/drops)
      trio_add.set <- rbind(trio_add.set,
                            subtrio.set,
                            subtrio.set)
      trio_drop.set <- rbind(trio_drop.set,
                             subtrio.set,
                             subtrio.set)

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
      edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,2]),
                        cbind(subtrio.set[,1], subtrio.set[,3]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,2], subtrio.set[,3]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))

      # Also store the indices of the trios corresponding to each edge drop/add (may be used to confirm adds/drops)
      trio_add.set <- rbind(trio_add.set,
                            subtrio.set,
                            subtrio.set)
      trio_drop.set <- rbind(trio_drop.set,
                             subtrio.set,
                             subtrio.set)

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
      edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,2]),
                        cbind(subtrio.set[,1], subtrio.set[,3]),
                        cbind(subtrio.set[,2], subtrio.set[,3]),
                        cbind(subtrio.set[,3], subtrio.set[,2]))

      # Also store the indices of the trios corresponding to each edge drop/add (may be used to confirm adds/drops)
      trio_add.set <- rbind(trio_add.set,
                            subtrio.set,
                            subtrio.set,
                            subtrio.set,
                            subtrio.set)

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0

    }

    # Update edge directions
    ## Add edges
    if (!is.null(edge_add)) {
      # Remove duplicates
      edge_add <- unique(edge_add)

      # Check for conflicts with previous drops before adding edge directions
      # (Check that we are not adding an edge we removed earlier)
      add.edges_new <- paste0(edge_add[,1], '.', edge_add[,2])
      if (!is.null(dropped.edges)) {
        conflicts.adds <- is.element(add.edges_new, dropped.edges)
        nb.add_conflicts <- sum(conflicts.adds)

        if (nb.add_conflicts) {# If any conflict
          # Does the edge involve a genetic variant or just T-nodes?
          V.conflicts.adds <-  edge_add[conflicts.adds,1] <= q

          confirm.add <- matteApply(X = cbind(edge_add[conflicts.adds,, drop = FALSE], V.conflicts.adds),
                                    MARGIN = 1,
                                    FUN = conflict.add.fun,
                                    Adj = Adj,
                                    method = method,
                                    cl = cl, chunk.size = chunk.size)
          edge_add <- rbind(edge_add[!conflicts.adds,],
                            edge_add[conflicts.adds,,drop = FALSE][confirm.add,])

          # Update the list of newly added edges
          add.edges_new <- paste0(edge_add[,1], '.', edge_add[,2])

          # Did we re-added an edge direction we dropped?
          if (any(re.add <- is.element(dropped.edges, add.edges_new))) {
            # Update the list of dropped edge directions
            dropped.edges <- dropped.edges[!re.add]
          }
        }
      }
      else {
        nb.add_conflicts <- 0
      }

      # Update edge directions
      Adj[edge_add] <- 1

      # Update the list of added edges
      added.edges <- c(added.edges, add.edges_new)
    }
    else {
      nb.add_conflicts <- 0
    }

    ## Drop edges
    if (!is.null(edge_drop)) {
      # Remove duplicates
      edge_drop <- unique(edge_drop)

      # Check for conflicts before dropping directions
      drop.edges_new <- paste0(edge_drop[,1], '.', edge_drop[,2])
      if (!is.null(added.edges)) {
        conflicts.drops <- is.element(drop.edges_new, added.edges)
        nb.drop_conflicts <- sum(conflicts.drops)

        if (nb.drop_conflicts) {# If any conflict
          # Does the edge involve a genetic variant or just T-nodes?
          V.conflicts.drops <-  edge_drop[conflicts.drops,1] <= q
          confirm.drop <- matteApply(X = cbind(edge_drop[conflicts.drops,,drop = FALSE], V.conflicts.drops),
                                     MARGIN = 1,
                                     FUN = conflict.drop.fun,
                                     Adj = Adj,
                                     method = method,
                                     cl = cl, chunk.size = chunk.size)
          edge_drop <- rbind(edge_drop[!conflicts.drops,],
                             edge_drop[conflicts.drops,,drop = FALSE][confirm.drop,,drop = FALSE])

          # Update the list of newly dropped edges
          drop.edges_new <- paste0(edge_drop[,1], '.', edge_drop[,2])

          # Correct 'edges.new' for the changes due to conflicts
        }
      }
      else {
        nb.drop_conflicts <- 0
      }

      # Update edge directions
      Adj[edge_drop] <- 0

      # Update the list of dropped edge directions
      dropped.edges <- c(dropped.edges, drop.edges_new)

      # Update the list of added edge directions (did we drop a direction we added before)
      if (!is.null(added.edges) & sum((re.drop <- added.edges %in% drop.edges_new))) {
          added.edges <- added.edges[!re.drop]

          cat(paste0('          * Dropping previously added edge directions: ', sum(re.drop), ' drop(s).\n'))

          if (!is.null(edge_add)) { # Check this for newly added edge directions too
            if (sum((re.drop <- add.edges_new %in% drop.edges_new))) {
              edge_add <- edge_add[re.drop,, drop = FALSE]
            }
          }
      }
    }
    else {
      nb.drop_conflicts <- 0
    }

    # Count newly added edges
    if (!is.null(edge_add)) {
      ## Set of new edges (all directions)
      add.edges_new <- paste0(edge_add[,1], '.', edge_add[,2])
      add.edges_new.rv <- paste0(edge_add[,2], '.', edge_add[,1])

      ## Count new edges
      nb.new_edges <- sum((add.edges_new %in% added.edges.old) + (add.edges_new.rv %in% added.edges.old) == 0)
    }
    else
      nb.new_edges <- 0

    # Count dropped edges
    if (!is.null(edge_drop)) {
      ## Set of dropped edges (all directions)
      drop.edges_new <- paste0(edge_drop[,1], '.', edge_drop[,2])
      drop.edges_new.rv <- paste0(edge_drop[,2], '.', edge_drop[,1])

      ## Count new edges
      nb.dropped_edges <- sum((drop.edges_new %in% dropped.edges.old) + (drop.edges_new.rv %in% dropped.edges.old) == 0)
    }
    else
      nb.dropped_edges <- 0

    return(list(Adj = Adj,
                nb.new_edges = nb.new_edges,
                nb.dropped_edges = nb.dropped_edges,
                new.edges = edges.new,               # (index of trios with new edges) TO BE CORRECTED
                added.edges = added.edges,
                dropped.edges = dropped.edges,
                nb.conflicts = nb.drop_conflicts + nb.add_conflicts))
  }
