#' Generate data from a DAG
#'
#' @description Description
#'
#' @export sim.graph.data
#'
#' @import igraph
#' @import pcalg
#'
# To do: reorder variables: V, T, W, Z, K, U

sim.graph.data <- function (number.of.V,
                            number.of.T,
                            conf.num.vec, # c("K","U","W","Z")
                            graph_type = "scale-free", # not used when 'method' is supplied
                            method,
                            degree = 3,
                            connection_prob = 0.05, # not used now
                            mixed = "None", # not used now
                            theta = .5,
                            b0.1 = 0,
                            b.snp = 1,
                            b.med = 0.8 * b.snp,
                            sd.1 = 1,
                            neg.freq = 0.5,
                            conf.coef.ranges = list(K = c(0.01, 0.1),
                                                    U = c(0.15,0.5),
                                                    W = c(0.15,0.5),
                                                    Z = c(1, 1.5)),
                            sample.size) {

  # Total number of nodes
  number.nodes <-  number.of.V + number.of.T + sum(conf.num.vec)

  # Initialize the adjacency matrix
  A <- matrix(0, nrow = number.nodes, ncol = number.nodes)

  # Create all confounding variable names
  if(sum(conf.num.vec)>0){
    conf.node.names = NULL
    letter.id = c("K","U","W","Z")
    nz.letter.id = letter.id[which(conf.num.vec>0)]
    conf.num.vec2 = conf.num.vec[which(conf.num.vec>0)]
    for(i in 1:length(conf.num.vec2)){
      conf.node.names = append(conf.node.names,
                               paste0(nz.letter.id[i], c(1:conf.num.vec2[i])))
    }
  }
  else
    conf.node.names = NULL
  row.names(A) = colnames(A) = c(paste0("V", c(1:number.of.V)),
                                 paste0("T", c(1:number.of.T)),
                                 conf.node.names)

  # Create the 'method' argument for 'get.custom.graph' starting from 'graph_type'
  if (missing(method)) {
    method <- switch(graph_type,
                   `scale-free` = "power",
                   `small-world` = "watts",
                   `random graph` = "er")
  }
  stopifnot(method %in% c("regular",
                          "watts",
                          "er",
                          "power",
                          "bipartite",
                          "barabasi",
                          "geometric",
                          "interEr"))

  # Generate the graph skeleton and the effects
  A = get.custom.graph(Adj = A,
                       b.snp = b.snp,
                       b.med = b.med,
                       struct = "random",
                       number.of.V = number.of.V,
                       number.of.T = number.of.T,
                       conf.num.vec = conf.num.vec,
                       conf.coef.ranges = conf.coef.ranges,
                       neg.freq = neg.freq,
                       degree = degree,
                       method = method)

  # Save the effects
  effects.adj = A

  # Convert effects matrix to adjacency matrix
  A[A!=0] = 1

  # Create a graph object
  igraph.obj = igraph::graph_from_adjacency_matrix(A)
  graph.attr <- list(adjacency = A, effects.adj = effects.adj, igraph.obj = igraph.obj)

  # Initialize a data matrix
  X = as.data.frame(matrix(0, nrow = sample.size, ncol = dim(graph.attr$adjacency)[2]))
  colnames(X) = colnames(graph.attr$adjacency)

  # Get the topological ordering
  topo.order = colnames(graph.attr$adjacency)[as.vector(igraph::topo_sort(graph.attr$igraph.obj))]

  # Simulate data for each node
  for(i in 1:length(topo.order)){
    #generate V nodes in topo order
    location = match(topo.order[i], colnames(graph.attr$adjacency))
    parent.list = find.parents(Adjacency = graph.attr$adjacency, location = location)

    if(grepl("V", topo.order[i])){
      X[,location] = c(sample(c(0, 1, 2), size = sample.size, replace = TRUE,
                              prob = c((1 -theta)^2, 2 * theta * (1 - theta), theta^2)))
    }else{
      #catch nodes with no parents of any kind
      if(sum(unlist(lapply(parent.list, is.na)))==6){
        X[, location] = stats::rnorm(n = sample.size, mean = b0.1, sd = sd.1)
      }else{
        #simulate all other types of nodes according to parental list, topo.order, and the effects adj
        coefs = graph.attr$effects.adj[stats::na.omit(unlist(parent.list)), topo.order[i]]
        Cond.vals = as.matrix(X[, stats::na.omit(unlist(parent.list))])
        X[, location] = stats::rnorm(n = sample.size, mean = b0.1 + Cond.vals %*% coefs, sd = sd.1)
 #       X[, location] = sd.1 * X[, location] / sd(X[, location]) + b0.1
      }
    }
  }

  return(list(data = X,
              dims = list(p = number.of.T,
                         q = number.of.V,
                         r1 = conf.num.vec[3],
                         r2 = conf.num.vec[4],
                         u1 = conf.num.vec[1],
                         u2 = conf.num.vec[2]),
              sd.1 = sd.1,
              Adjacency = graph.attr$adjacency,
              Effects = graph.attr$effects.adj,
              igraph = graph.attr$igraph.obj))
}
