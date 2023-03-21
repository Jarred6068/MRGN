#' sample.graph.data
#'
#' @description sample.graph.data
#'
#' @export sample.graph.data
#'
#' @import graphsim
#' @import Matrix

# 'graph_type' specifies a limited number of graph structure
# one of "small-world", "scale-free" or "random graph"
# 'method' is passed to 'pcalg::randDAG',
# one of "regular", "watts", "er", "power", "bipartite", "barabasi", "geometric", or "interEr"

# The graphsim package is returning NaNs distances when 'geom.dist = FALSE' with large 'number.nodes'

sample.graph.data <- function (number.of.V,
                               number.of.T,
                               conf.num.vec,
                               graph_type = "scale-free", # explicit short hand for 'method', not used when 'method' is supplied
                               method,                # See '?pcalg::randDAG';
                               degree = 3,            # See '?pcalg::randDAG';
                               geom.dist = TRUE,      # Type of distance matrix to get marginal covariance matrix
                               cor.struct = "AR1",    # First Order Auto-regressive (unique available now)
                               cor.snp = sqrt(2)/2,   # Maximum absolute correlation between any Variant and any Expression
                               neg.snp.freq = 0.5,    # Frequency of negative snp effects/correlations
                               cor.med = 0.8*cor.snp, # Maximum correlation between two Expressions
                               neg.med.freq = neg.snp.freq,# Frequency of negative mediation effects/correlations
                               cor.conf.range = list(K = c(0.02, 0.10),  # Range of absolute correlations
                                                     U = c(0.15, 0.25),  # between Expressions and Confounding variables
                                                     W = c(0.15, 0.25),
                                                     Z = c(0.30, 0.50)),
                               neg.conf.freq = 0.5,   # Frequency of negative confounding effects/correlations
                               sd.expr = 0.5,         # Common marginal SD of T,K,U,W,Z nodes, Replace by a vector to sample from
                               mean.expr = 0,         # Common marginal mean of T,K,U,W,Z nodes, Replace by a vector to sample from
                               theta = 0.4,           # Minor allele frequency, Replace by a vector to sample from
                               sample.size,           # Defaults to 'number.nodes + 1'
                               verbose = 1) {         # Notice the use of nearest Positive Definite Matrix

  # Total number of nodes
  number.nodes <-  number.of.V + number.of.T + sum(conf.num.vec)

  # Default sample size
  if (missing(sample.size))
    sample.size <- number.nodes + 1

  # Initialize the adjacency matrix
  A <- matrix(0, nrow = number.nodes, ncol = number.nodes)

  # Create all confounding variable names
  if(sum(conf.num.vec)>0){
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

  # Force correlation value to be positive
  cor.snp <- abs(cor.snp)
  cor.med <- abs(cor.med)
  cor.conf.range <- lapply(cor.conf.range, FUN = abs)
  names(cor.conf.range) <- c('K', 'U', 'W', 'Z')

  # Check the support of correlation values and frequencies of negative correlations
  stopifnot(c(cor.snp, cor.med, unlist(cor.conf.range)) < 1)
  stopifnot(c(c(neg.snp.freq, neg.med.freq, neg.conf.freq) <= 1,
              c(neg.snp.freq, neg.med.freq, neg.conf.freq) >= 0))

  # Generate the graph skeleton and the effects (correlations)
  # Revise this / Remove effects sampling therein - do that herein
  A = get.custom.graph(Adj = A,
                       b.snp = cor.snp,
                       b.med = cor.med,
                       struct = "random",
                       number.of.V = number.of.V,
                       number.of.T = number.of.T,
                       conf.num.vec = conf.num.vec,
                       conf.coef.ranges = cor.conf.range,
                       neg.freq = neg.conf.freq,
                       degree = degree,
                       method = method)

  # Save the 'max' correlations
  correl.adj = A

  # Convert effects matrix to adjacency matrix
  A[A!=0] = 1

  # Allow negative 'snp' and 'med' correlations
  V.index <- 1:number.of.V
  T.index <- (number.of.V + 1):(number.of.V + number.of.T)
  if (neg.snp.freq > 0) {
    nb.effects <- sum(A[V.index, T.index])
    if (nb.effects) {
    correl.adj[V.index, T.index][c(A[V.index, T.index] == 1)] <-
      correl.adj[V.index, T.index][c(A[V.index, T.index] == 1)] *
      (1 - 2 * rbinom(nb.effects, size = 1, prob = neg.snp.freq))
    }
  }
  if (neg.med.freq > 0) {
    nb.effects <- sum(A[T.index, T.index])
    if (nb.effects) {
      correl.adj[T.index, T.index][c(A[T.index, T.index] == 1)] <-
        correl.adj[T.index, T.index][c(A[T.index, T.index] == 1)] *
        (1 - 2 * rbinom(nb.effects, size = 1, prob = neg.med.freq))
    }
  }

  # Create a graph object
  igraph.obj = igraph::graph_from_adjacency_matrix(A)
  graph.attr <- list(adjacency = A, correl.adj = correl.adj, igraph.obj = igraph.obj)

  # Obtain a distance matrix based on shortest paths
  Dist.mat <- graphsim:::make_distance_adjmat (mat = graph.attr$adjacency, directed = TRUE, absolute = !geom.dist)

  # Obtain a marginal correlation matrix
  ### initialize to the signed maximal correlation values
  Corr.mat <- (graph.attr$correl.adj + t(graph.attr$correl.adj)) / 2
  diag(Corr.mat) <- 1
  nearPD <- FALSE
  switch(cor.struct,
         AR1 = {
           Corr.mat <- Corr.mat * Dist.mat
           Val.propre <- eigen(Corr.mat)
           if (any(Val.propre$values <= 0)) { # Using Higham algorithm
             if (verbose) {
               warning("non-positive definite sampled correlation matrix; adjusting with 'nearPD'")
               nearPD < - TRUE
             }

             Corr.mat <- Matrix::nearPD(Corr.mat, corr = TRUE, base.matrix = TRUE,
                                        doSym = TRUE, ensureSymmetry = FALSE)$mat
           }
         })

  # Obtain the marginal covariance matrix (REVISE IF 'sd.expr' is allowed to be a vector)
  Cov.mat <- Corr.mat
  Cov.mat[V.index, V.index] <- 2 * theta * (1 - theta) * Corr.mat[V.index, V.index]
  Cov.mat[T.index, T.index] <- sd.expr * sd.expr * Corr.mat[T.index, T.index]
  Cov.mat[T.index, V.index] <- sqrt(2 * theta * (1 - theta)) * sd.expr * Corr.mat[T.index, V.index]
  Cov.mat[V.index, T.index] <-sqrt(2 * theta * (1 - theta)) * sd.expr * Corr.mat[V.index, T.index]

  # Build a vector of marginal expectations
  mean.vec <- c(rep(2 * theta, number.of.V),
                rep(mean.expr, number.nodes - number.of.V))

  # Initialize a data matrix
  X = as.data.frame(matrix(0, nrow = sample.size, ncol = dim(graph.attr$adjacency)[2]))
  colnames(X) = colnames(graph.attr$adjacency)

  # Initialize an effect matrix
  Effects.mat <- Cov.mat
  diag(Effects.mat) <- 0
  Effects.mat[V.index, V.index] <- 0
  Effects.mat[T.index, V.index] <- 0

  # Get the topological ordering
  topo.order = colnames(graph.attr$adjacency)[as.vector(igraph::topo_sort(graph.attr$igraph.obj))]

  # Simulate data for each node
  for(i in 1:length(topo.order)){
    #generate V nodes in topo order
    location = match(topo.order[i], colnames(graph.attr$adjacency))
    if(grepl("V", topo.order[i])){
      X[,location] = c(sample(c(0, 1, 2), size = sample.size, replace = TRUE,
                              prob = c((1 -theta)^2, 2 * theta * (1 - theta), theta^2)))
    }else{
      # Find parent nodes
      parent.list = find.parents(Adjacency = graph.attr$adjacency, location = location)
      #catch nodes with no parents of any kind
      if(sum(unlist(lapply(parent.list, is.na)))==6){
        X[, location] = stats::rnorm(n = sample.size, mean = mean.expr, sd = sd.expr)
      }else{
        #simulate all other types of nodes according to parental list, topo.order, and the effects adj
        parent.index <- stats::na.omit(unlist(parent.list))
        Xa <- t(t(X[, parent.index, drop = FALSE]) - mean.vec[parent.index])
        Sigma_a <- Cov.mat[parent.index, parent.index]
        Sigma_c <- Cov.mat[parent.index, topo.order[i]]
        Cond.var <-  c(solve(Sigma_a) %*% Sigma_c)

        # Save the actual effects
        Effects.mat[parent.index, topo.order[i]] <- Cond.var

        Cond.mean <-  mean.expr + Xa %*% as.vector(Cond.var)
        Cond.var = sd.expr^2 - sum( Sigma_c * Cond.var )
        X[, location] = stats::rnorm(n = sample.size, mean = Cond.mean, sd = sqrt(Cond.var))
      }
    }
  }

  Max.corr.mat <- graph.attr$correl.adj
  diag(Max.corr.mat) <- 1

  return(list(data = X,
              Adjacency = graph.attr$adjacency,
              Effects.mat = Effects.mat,
              Max.corr.mat = Max.corr.mat,
              Corr.mat = Corr.mat,
              Cov.mat = Cov.mat,
              igraph = graph.attr$igraph.obj,
              nearPD = nearPD))
}
