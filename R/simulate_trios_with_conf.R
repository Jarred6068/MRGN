


#' A function simulate a child node with one parent and confounding variables
#'
#' @param N Simulation sample size
#' @param P1 The data for the parent node
#' @param b0.1 the intercept
#' @param b1.1 the effect of the parent node
#' @param sd.1 the residual standard deviation (noise)
#' @param G A matrix of size N X g of confounding variables
#' @param u the indices of the confounders to be selected. If NULL confounders are omitted
#' @param gamma.u a vector of length g specifying the effects of the confounding variables
#' @return a vector of length N containing the simulated values of the child node
#' @export simData1P


#################################################
simData1P=function (N, P1, b0.1, b1.1, sd.1, G=NULL, u=NULL, gamma.u=NULL) {

  if(is.null(u)){
    t <- stats::rnorm(n = N, mean = b0.1 + b1.1 * P1, sd = sd.1)
  }else{
    U=as.matrix(G[,u])
    #gamma.u=runif(length(u), min = cs.range[1], max=cs.range[2])
    t <- stats::rnorm(n = N, mean = b0.1 + b1.1 * P1 + U%*%gamma.u, sd = sd.1)
  }

  return(t)
}

#' A function to simulate a child node with two parents and confounding variables
#'
#' @param N Simulation sample size
#' @param P1 The data for the first parent node
#' @param P2 the data for the second parent node
#' @param b0.1 the intercept
#' @param b1.1 the effect of the parent node P1
#' @param b1.2 the effect of the parent node P2
#' @param sd.1 the residual standard deviation (noise)
#' @param G A matrix of size N X g of confounding variables
#' @param u the indices of the confounders to be selected. If NULL confounders are omitted
#' @param gamma.u a vector of length g specifying the effects of the confounding variables
#' @return a vector of length N containing the simulated values of the child node
#' @export simData2P


#################################################
simData2P=function (N, P1, P2, b0.1, b1.1, b1.2, sd.1, G=NULL, u=NULL, gamma.u=NULL) {

  if(is.null(u)){
    t <- stats::rnorm(n = N, mean = b0.1 + b1.1 * P1 + b1.2 * P2, sd = sd.1)
  }else{
    U=as.matrix(G[,u])
    #gamma.u=runif(length(u), min = cs.range[1], max = cs.range[2])
    t <- stats::rnorm(n = N, mean = b0.1 + b1.1 * P1 + b1.2 * P2 + U%*%gamma.u, sd = sd.1)
  }

  return(t)
}



#################################################
# simData3P=function (N, P1, P2, P3, b0.1, b1.1, b1.2, b1.3, sd.1, cs.max) {
#
#   t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + b1.2 * P2 + b1.3 *
#                P3, sd = sd.1)
#   return(t)
# }
#################################################

#' A function to simulate a child node with no parents and confounding variables
#'
#' @param N Simulation sample size
#' @param b0.1 the intercept
#' @param sd.1 the residual standard deviation (noise)
#' @param G A matrix of size N X g of confounding variables
#' @param u the indices of the confounders to be selected. If NULL confounders are omitted
#' @param gamma.u a vector of length g specifying the effects of the confounding variables
#' @return a vector of length N containing the simulated values of the child node
#' @export simDataNP


simDataNP=function (N, b0.1, sd.1, G=NULL, u=NULL, gamma.u = NULL) {

  if(is.null(u)){
    t <- stats::rnorm(n = N, mean = b0.1, sd = sd.1)
  }else{
    U=as.matrix(G[,u])
    #gamma.u=runif(length(u), min = cs.range[1], max = cs.range[2])
    t <- stats::rnorm(n = N, mean = b0.1 + U%*%gamma.u, sd = sd.1)
  }
  return(t)
}


# simConf=function(n = NULL, means=NULL, variances=NULL){
#
#   U=rmvnorm(n = n, mean = means, sigma = diag(variances))
#   return(U)
#
# }


#' A function to simulate a trio with confounding variable effects
#'
#' This function simulates a trio from one of the model toplogies given by Badsha and Fu, 2019 with counding effects on
#' the molecular phenotypes. Note the effects of the confounding variables are selected from uniform(0.15,0.5) Yang et al. 2017
#'
#' @param theta the frequency of the alternative allele where the reference is (1-theta)
#' @param model a string specifying one of "model0","model1", "model2", "model3","model4"
#' @param b0.1 the intercept
#' @param b.snp the effect of the genetic variant
#' @param b.med the effect between the molecular phenotypes
#' @param sd.1 the residual standard deviation (noise)
#' @param G A matrix of size N X g of confounding variables
#' @param u either a numeric specifying the number of confounders to be selected at random from G or
#' the indices of the confounders to be selected.
#' @return a matrix of dimension N X (3+g) representing the simulated trio and confounding variables
#' @export simData


#################################################

simData=function (theta, model, b0.1, b.snp, b.med, sd.1, G=NULL, u=NULL){

  #generate confounders
  G.dim=dim(G)[2]
  #conf=simConf(n = N, means = conf.means, variances = conf.var)
  #uniformly sample confounders
  #upper.b=round(2*u)
  #sample.u=round(runif(1, 1, upper.b))
  #sample.b=round(runif(1, 2, upper.b))
  #sample.c=round(runif(1, 2, upper.b))
  if(length(u)==1){
    u.idx=sample(1:G.dim, u, replace = FALSE)
  }else{
    u.idx=u
  }

  #b=sample(1:conf.dim, sample.b)
  #c=sample(1:dim(conf)[2], sample.c)
  N=dim(G)[1]

  #confounder effects - allow 50% to be negative and 50% to be positive
  #different coefficients for each node U is a parent of
  ct1 = stats::rbinom(u, 1, 0.5)
  ct2 = stats::rbinom(u, 1, 0.5)
  w1 = stats::runif(u, min = 0.15, max = 0.5)
  w1[which(ct1==1)]=-1*w1[which(ct1==1)]
  w2 = stats::runif(u, min = 0.15, max = 0.5)
  w2[which(ct2==1)]=-1*w2[which(ct2==1)]

  V1 <- c(sample(c(0, 1, 2), size = N, replace = TRUE,
                 prob = c((1 -theta)^2, 2 * theta * (1 - theta), theta^2)))

  colnames(G)=paste0("U",c(1:G.dim))

  switch(model, model0 = {

    T1 <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1, G = G,
                    u = u.idx, gamma.u = w1)
    T2 <- simDataNP(N = N, b0.1 = b0.1, sd.1 = sd.1, G = G, u = u.idx,
                    gamma.u = w2)

    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U))

  }, model1 = {

    T1 <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u = w1)
    T2 <- simData1P(N = N, P1 = T1, b0.1 = b0.1, b1.1 = b.med, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u = w2)
    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U))

  }, model2 = {

    T2 <- simDataNP(N = N, b0.1 = b0.1, sd.1 = sd.1, G = G, u = u.idx,
                    gamma.u = w1)
    T1 <- simData2P(N = N, P1 = V1, P2 = T2, b0.1 = b0.1, b1.1 = b.snp, b1.2 = b.med, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u=w2)
    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U))

  }, model3 = {

    T1 <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u=w1)
    T2 <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                    G = G, u = u.idx, gamma.u=w2)

    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U))

  }, model4 = {

    T1.a <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                      G = G, u = u.idx, gamma.u=w1)
    T2.a <- simData2P(N = N, P1 = V1, P2 = T1.a, b0.1 = b0.1, b1.1 = b.snp, b1.2 = b.med, sd.1 = sd.1,
                      G = G, u = u.idx, gamma.u=w2)
    T2.b <- simData1P(N = N, P1 = V1, b0.1 = b0.1, b1.1 = b.snp, sd.1 = sd.1,
                      G = G, u = u.idx, gamma.u=w2)
    T1.b <- simData2P(N = N, P1 = V1, P2 = T2.b, b0.1 = b0.1, b1.1 = b.snp, b1.2 = b.med, sd.1 = sd.1,
                      G = G, u = u.idx, gamma.u=w1)

    coinToss <- stats::rbinom(n = N, size = 1, prob = 0.5)
    T1 <- rep(0, N)
    T1[which(coinToss == 0)] <- T1.a[which(coinToss == 0)]
    T1[which(coinToss == 1)] <- T1.b[which(coinToss == 1)]
    T2 <- rep(0, N)
    T2[which(coinToss == 0)] <- T2.a[which(coinToss == 0)]
    T2[which(coinToss == 1)] <- T2.b[which(coinToss == 1)]
    U=cbind.data.frame(G[,u.idx])
    return(data.frame(V1 = V1, T1 = T1, T2 = T2, U))

  }, stop("Model not included or missing"))
}


#################################################

#' A function to simulate effects for edges in the graph
#'
#' This function is wrapped by gen.graph.skel() and get.custom.graph()
#'
#' @param n.effects the number of effects to simulate
#' @param coef.range.list a list of length 4 containing the mininum and maximum effect size for U,K,W, and Z variables
#' (in that order)
#' @param neg.freq the frequency of negative effects for a given variable
#' @export gen.conf.coefs
#' @return a vector of simulated effects

gen.conf.coefs=function(n.effects, coef.range.list, neg.freq){

  weights = stats::runif(n.effects, coef.range.list[1], coef.range.list[2])
  cointoss = stats::rbinom(n.effects, 1, neg.freq)
  weights[which(cointoss==1)]=-1*weights[which(cointoss==1)]

  return(weights)

}
#################################################

#' A function to fill in an adjacency matrix with all connections between the confounding variable nodes
#'
#' @param A The empty adjacency matrix of the graph passed from gen.graph.skel()
#' @param b.snp either a single number or vector representing the effect(s) of the variant(s)
#' @param b.med either a single number or vector representing the effect(s) between the molecular phenotypes
#' @param struct either (1) the string "random" specifying a random DAG between the V and T nodes or (2) an adjacency matrix of dimension equal to the number of desired T and V nodes giving the connectivity between them
#' @param conf.num.vec a vector of length 4 giving the number of each type of confounding variables: specifically of K,U,W, and Z nodes
#' @param number.of.T the number of T-nodes in the graph
#' @param number.of.V the number of V-nodes in the graph
#' @param conf.coef.ranges a list of four 2-length vectors giving the minimum and maximum effect sizes of each type of confounding variable in the order of K,U,W, and Z
#' @param neg.freq the frequency for simulated effects to be negative. passed to gen.conf.coefs()
#' @param degree the expected number of neighbors (sum of in and out degree of the graph) passed to pcalg::randDAG
#' @param method the method used to generate a random DAG for the V and T nodes when struct = "random". passed to pcalg::randDAG
#' @export get.custom.graph
#' @return The adjacency matrix with all confounding variable node edges added

get.custom.graph = function(Adj, b.snp, b.med, struct, conf.num.vec, number.of.T, number.of.V,
                            conf.coef.ranges, neg.freq, degree, method = "er"){

  letter.id = c("K","U","W","Z")
  #storing info
  V.idx = 1:number.of.V
  T.idx = (1:number.of.T)+number.of.V
  #pre-structure steps
  if(struct == "random"){
    #random topology for V and T
    Adj.sub = methods::as(pcalg::randDAG(n = number.of.V+number.of.T, d = degree, method = method),"matrix")
    #replace with 1's
    Adj.sub[Adj.sub!=0] = 1
    #remove all T --> V and V --> V
    Adj.sub[,V.idx] = 0

    #handle v nodes with no edges: if no edges draw edges from bernoulli sequence
    for(i in 1:number.of.V){
      while(sum(Adj.sub[V.idx[i],])<1){
        Adj.sub[V.idx[i], T.idx] = stats::rbinom(number.of.T, 1, prob = degree/number.of.T)
      }
    }
  }else{
    #custom structure connectivity for V and T
    if(sum(colSums(Adj.sub[,V.idx]))>1){
      stop("custom graph contains edges directed towards variants!")
    }
    Adj.sub = struct
  }
  #create a second matrix for storing the simulation effects
  coefs.sub = Adj.sub
  #get the total number of V and T edges
  number.of.edges.V = sum(Adj.sub[V.idx,])
  number.of.edges.T = sum(Adj.sub[T.idx,])

  #replace edges with their simulated effects
  coefs.sub[V.idx,][coefs.sub[V.idx,]==1] = sample(b.snp, number.of.edges.V, replace = TRUE)
  coefs.sub[T.idx,][coefs.sub[T.idx,]==1] = sample(b.med, number.of.edges.T, replace = TRUE)
  #convert to igraph and get to the topological ordering
  topo.order = igraph::topo_sort(igraph::graph_from_adjacency_matrix(Adj.sub))
  Adj[1:(number.of.V+number.of.T), 1:(number.of.V+number.of.T)] = coefs.sub

  #handle confounders
  for(i in 1:2){
    loc = which(grepl(letter.id[i], colnames(Adj)))
    for(j in loc){
      weight = gen.conf.coefs(n.effects = 1, coef.range.list = conf.coef.ranges[[i]], neg.freq = neg.freq)
      Adj[j, sample(T.idx, 2, replace = F)] = weight
    }
  }

  #handle intermediate variables
  loc.int = which(grepl(letter.id[3], colnames(Adj)))
  for(i in loc.int){
    weights = gen.conf.coefs(n.effects = 2, coef.range.list = conf.coef.ranges[[3]], neg.freq = neg.freq)
    child.T = sample(topo.order[-c(1:(number.of.V+1))], 1)
    loc.in.topo.order = which(topo.order==child.T)
    poss.parents = c(topo.order[(number.of.V+1):(loc.in.topo.order-1)])
    if(length(poss.parents)>1){
      parent.T = sample(poss.parents, 1)
    }else{
      parent.T = poss.parents
    }
    # print(colnames(Adj)[topo.order])
    # print(colnames(Adj)[child.T])
    # print(colnames(Adj)[parent.T])
    Adj[i, child.T] = weights[1]
    Adj[parent.T, i] = weights[2]
  }

  loc.cc = which(grepl(letter.id[4], colnames(Adj)))
  for(i in loc.cc){
    weight = gen.conf.coefs(n.effects = 1, coef.range.list = conf.coef.ranges[[4]], neg.freq = neg.freq)
    Adj[sample(T.idx, 2, replace = F), i] = weight
  }

  return(Adj)

}

#################################################

#' A function to simulate a graph with confounder, intermediate, and common child variables
#'
#' This function simulates a graph adjacency matrix for a desired topology such as M0 - M4, custom, or random topological
#' graphs such that they include known confounder (denoted K), unkown confounder (denoted U), intermediate (denoted W), and common child
#' variables (denoted Z). This function wraps gen.conf.coefs() and get.custom.graph()
#'
#' @param model a string specifying one of "model0","model1", "model2", "model3","model4", or "custom"
#' @param b.snp a numeric or vector giving the effect(s) of the genetic variant(s).
#' @param b.med a numeric or vector giving the effect(s) between the molecular phenotypes(s).
#' @param conf.num.vec a numeric vector of length 4 containing the number of known confounders (K),
#' unknown confounders (U), intermediate (W), and common child variables (Z; in that order).
#' To exclude a variable type input a zero at the given position. Note: intermediate variables (W) cannot
#' be included for model4 and must excluded (set to zero)
#' @param number.of.T a numeric indicating the number of T variables desired for model == "custom" only
#' @param number.of.V a numeric indicating the number of V variables desired for model == "custom" only
#' @param struct For use when when model == "custom". Either (1) a sub-adjacency matrix of dimension (number.of.V + number.of.T X number.of.V + number.of.T) definining the topology of the V and T nodes or (2) the string "random" denoting a random topology
#' @param degree passed to get.custom.graph()
#' @param method passed to get.custom.graph()
#' @param plot.graph (logical) if TRUE the graph is plotted. default = TRUE
#' @param neg.freq the frequency of negative effects for simulated effects. passed to gen.conf.coefs()
#' @param conf.coef.ranges a list of length 4 representing the U,K,W, and Z (in that order) node effects where each list
#' element is a vector of length 2 giving the minimum and maximum effect sizes for the given node. default values follow
#' from the simulations of Yang et. al., 2017. passed to gen.conf.coefs()
#' @export gen.graph.skel
#' @return a list of length 3
#' \describe{
#'
#' \item{adjacency}{the indicator adjacency matrix of the graph}
#' \item{effects.Adj}{the adjacency matrix with the effects in place of the indicator values}
#' \item{igraph.obj}{an igraph object}
#' \item{true.adj}{only returned when model = "model4". It is the true adjacency matrix of the graph (the graph plotted plot.graph = TRUE) where as the
#'                 the graph in $adjacency is the graph used to emulate the model 4 network}
#'
#' }
#' @examples
#' # generate a model 1 graph with one of each type of confounding variables
#' adj = gen.graph.skel(model = "model1",
#'                      b.snp=0.8,
#'                      b.med=0.6,
#'                      conf.num.vec = c(K = 1, U = 1, W = 1, Z = 1))
#'
#' # generate a custom graph with 2 variants and 4 molecular phenotypes
#'\dontrun{
#' X=gen.graph.skel(model = "custom",
#'                  b.snp = c(0.5,0.6,0.7),
#'                  b.med = c(0.2,0.1,0.6),
#'                  conf.num.vec = c(K = 1, U = 1, W = 3, Z = 1),
#'                  number.of.T = 4,
#'                  number.of.V = 2,
#'                  struct = "random",
#'                  degree = 3,
#'                  method = "er",
#'                  plot.graph = TRUE,
#'                  neg.freq =0.5,
#'                  conf.coef.ranges=list(K=c(0.01, 0.1),
#'                                        U=c(0.15,0.5),
#'                                        W=c(0.15,0.5),
#'                                        Z=c(1, 1.5)))
#'}

#creates the graph skeleton for each model: contains u, k, w, and z variables in adj.
gen.graph.skel = function(model, b.snp, b.med, conf.num.vec, number.of.T, number.of.V, struct,
                          degree, method = "er", plot.graph = TRUE, neg.freq,
                          conf.coef.ranges=list(K=c(0.01, 0.1), U=c(0.15,0.5), W=c(0.15,0.5), Z=c(1, 1.5))){

  #preallocate names
  conf.node.names=NULL
  letter.id = c("K","U","W","Z")
  #create all confounder names
  if(sum(conf.num.vec)>0){
    nz.letter.id = letter.id[which(conf.num.vec>0)]
    conf.num.vec2 = conf.num.vec[which(conf.num.vec>0)]
    for(i in 1:length(conf.num.vec2)){
      conf.node.names = append(conf.node.names, paste0(nz.letter.id[i], c(1:conf.num.vec2[i])))
    }
  }
  #preallocate adjecency matrix:
  if(model == "custom"){
    #catch missing inputs for number of T and V nodes for custom graphs
    if(missing(number.of.T) | missing(number.of.V)){
      stop("Missing one of \"number.of.T\" or \"number.of.V\" for model = \"custom\" ")
    }
    #catch missing struct input for custom graph
    if(missing(struct)){
      stop("Missing arg: \"struct\" must enter a desired structure for V and T nodes when model = \"custom\" " )
    }
    #for custom graphs of chosen size
    B = matrix(0, nrow = sum(conf.num.vec)+number.of.T+number.of.V,
               ncol = sum(conf.num.vec)+number.of.T+number.of.V)

    row.names(B) = colnames(B) = c(paste0("V", c(1:number.of.V)),
                                   paste0("T", c(1:number.of.T)),
                                   conf.node.names)
  }else{
    #for the other M0-M4 models
    B = matrix(0, nrow = sum(conf.num.vec)+3, ncol = sum(conf.num.vec)+3)
    row.names(B) = colnames(B) = c("V1","T1","T2", conf.node.names)
  }

  #----------model-0-----------
  switch(model, model0 = {
    B[1,2] = 1*b.snp

    #----------model-1-----------
  }, model1 = {

    B[1,2] = 1*b.snp
    B[2,3] = 1*b.med

    #----------model-2-----------
  }, model2 = {

    B[1,2] = 1*b.snp
    B[3,2] = 1*b.med

    #----------model-3-----------
  }, model3 = {

    B[1,2:3] = 1*b.snp

    #----------model-4-----------
  }, model4 = {

    #catch when model 4 is entered with intermediate variables
    # if(isTRUE(conf.num.vec[3]>0)){
    #   stop("Model4 cannot include intermediate variables (W) because it produces cycles, try entering 0 for W in conf.num.vec instead")
    # }

    B[1,2:3] = 1*b.snp
    B[2,3] = 1*b.med
    #B[3,2] = 1*b.med

    #----------custom-model-----------
  }, custom = {
    #generate the custom graph::
    B = get.custom.graph(Adj = B, b.snp = b.snp, b.med = b.med, struct = struct, conf.num.vec = conf.num.vec,
                         number.of.T = number.of.T, number.of.V = number.of.V, conf.coef.ranges = conf.coef.ranges,
                         neg.freq = neg.freq, degree = degree, method = method)
    #convert effects mat to adj mat and plot
    A = B
    A[A!=0] = 1
    igraph.obj = igraph::graph_from_adjacency_matrix(A)
    if(plot.graph == TRUE){
      sub.graph.ind = 1:(number.of.V+number.of.T)
      par(mfrow = c(1,2))
      igraph::plot.igraph(igraph.obj, layout = igraph::layout_nicely, edge.arrow.size = 0.2)
      igraph::plot.igraph(igraph::graph_from_adjacency_matrix(A[sub.graph.ind,sub.graph.ind]),
                          layout = igraph::layout_nicely, edge.arrow.size = 0.2, vertex.color = "green")
    }

    return(list(adjacency = A, effects.adj = B, igraph.obj = igraph.obj))

  }, stop("Model not included or missing"))

  #handle confounders, intermediate, and common child vars for 5 basic topos
  #confounders (adding in effects to the effects matrix)

  for(i in 1:length(letter.id)){
    if(conf.num.vec[i]>0){
      if(letter.id[i] == "K" | letter.id[i] == "U"){
        #confounder
        for(j in 2:3){
          weights = gen.conf.coefs(n.effects = conf.num.vec[i], coef.range.list = conf.coef.ranges[[i]],
                                   neg.freq = neg.freq)
          B[which(grepl(letter.id[i],conf.node.names))+3, j] = 1*weights
        }
      }else if(letter.id[i] == "W"){
        #intermediate
        if(model == "model2"){
          weights1 = gen.conf.coefs(n.effects = conf.num.vec[3], coef.range.list = conf.coef.ranges[[3]],
                                    neg.freq = neg.freq)
          weights2 = gen.conf.coefs(n.effects = conf.num.vec[3], coef.range.list = conf.coef.ranges[[3]],
                                    neg.freq = neg.freq)
          B[which(grepl("W",conf.node.names))+3, 2] = 1*weights1
          B[3, which(grepl("W",conf.node.names))+3] = 1*weights2
        }else{
          weights1 = gen.conf.coefs(n.effects = conf.num.vec[3], coef.range.list = conf.coef.ranges[[3]],
                                    neg.freq = neg.freq)
          weights2 = gen.conf.coefs(n.effects = conf.num.vec[3], coef.range.list = conf.coef.ranges[[3]],
                                    neg.freq = neg.freq)
          B[which(grepl("W",conf.node.names))+3, 3] = 1*weights1
          B[2, which(grepl("W",conf.node.names))+3] = 1*weights2
        }

      }else if(letter.id[i] == "Z"){
        #common child
        for(k in 2:3){
          weights = gen.conf.coefs(n.effects = conf.num.vec[4], coef.range.list = conf.coef.ranges[[4]],
                                   neg.freq = neg.freq)
          B[k, which(grepl("Z",conf.node.names))+3] = 1*weights
        }
      }
    }
  }
  #convert effects mat to adjacency mat and plot
  A = B
  A[A!=0] = 1
  igraph.obj = igraph::graph_from_adjacency_matrix(A)
  if(plot.graph == TRUE){
    sub.graph.ind = 1:3
    par(mfrow = c(1,2))
    igraph::plot.igraph(igraph.obj, layout = igraph::layout_nicely, edge.arrow.size = 0.2)
    igraph::plot.igraph(igraph::graph_from_adjacency_matrix(A[sub.graph.ind,sub.graph.ind]),
                        layout = igraph::layout_nicely, edge.arrow.size = 0.2, vertex.color = "green")
  }

  # #edit the effects matrix if a model is model4 for later use in simData.from.graph
  # if(model == "model4"){
  #   #add in T1.a, T2.a, T2.b, T1.b into the effects matrix while retaining the other vars
  #   D = cbind.data.frame(c(B[1,1],rep(0,4),B[-c(1:3),1]),
  #                        matrix(0, nrow = dim(B)[1]+2, ncol = 4),
  #                        rbind(B[1,-c(1:3)], matrix(0, nrow = 4, ncol = dim(B)[1]-3), B[-c(1:3), -c(1:3)]))
  #   #rename cols
  #   colnames(D) = row.names(D) = c(colnames(B)[1], c("T1.a","T2.a","T1.b","T2.b"), colnames(B)[-c(1:3)])
  #   #add in effects appropriately for topo order
  #   #snp effects
  #   D[1,2:5] = b.snp
  #   #T1 - T2 effects
  #   D[2,5] = b.med
  #   D[3,4] = b.med
  #   #confounder effects
  #   D[-c(1:5), 2:5] = B[-c(1:3),2:3]
  #   #common child effects
  #   D[2:5, which(grepl("Z",colnames(D)))] = B[2:3, which(grepl("Z", colnames(B)))]
  #   #convert non-zeros in D to 1 to create adjacency matrix for use with simData.from.graph
  #   AA=as.matrix(D)
  #   AA[AA!=0]=1
  #   igraph.obj = igraph::graph_from_adjacency_matrix(AA)
  #   #return
  #   return(list(adjacency = AA, effects.adj = as.matrix(D), igraph.obj = igraph.obj, true.adj = A))
  # }
  return(list(adjacency = A, effects.adj = B, igraph.obj = igraph.obj))
}


#################################################

#' A function to find the parents of a given node in the graph
#'
#' This function uses the adjacency matrix of a graph to identify the parents of a desired node in the graph
#'
#' @param location the column index of the node of interest in the adjacency matrix
#' @param Adjacency an adjacency matrix for the graph
#' @return a list of length 6 representing the 6 nodes types, V,T,K,U,W,Z

find.parents = function(Adjacency, location){
  #define letter identifier for node types
  letter.id = c("V","T","K","U","W","Z")
  all.parents = row.names(Adjacency)[which(Adjacency[, location] == 1)]
  #allocate each type of parent indexes to list by type
  parent.list = lapply(letter.id,function(x,y,z){if(any(grepl(x,y))) match(y[which(grepl(x,y))], z) else NA},
                       y = all.parents, z = colnames(Adjacency))
  names(parent.list)=c("V","T","K","U","W","Z")

  return(parent.list)
}


#################################################

#' A function to simulate data from a desired graph topology
#'
#' This function wraps gen.graph.skel(), gen.conf.coef(), and find.parents() to simulate data from a given graph topology
#' using the topological order.
#'
#' @param model a string specifying one of "model0","model1", "model2", "model3","model4", or "custom"
#' @param theta the frequency of the minor allele of the genetic variant
#' @param b0.1 the intercept
#' @param b.snp a numeric or vector giving the effect(s) of the genetic variant(s)
#' @param b.med a numeric or vector giving the effect(s) between the molecular phenotypes(s)
#' @param sd.1 the residual standard deviation (noise)
#' @param conf.num.vec a numeric vector of length 4 containing the number of known confounders,
#' unknown confounder, intermediate, and common child variables (in that order). To exclude a variable type,
#' input a zero at the given position. Note: intermediate variables (W) cannot be included for model4 and must
#' excluded (set to 0)
#' @param number.of.T a numeric indicating the number of T variables desired for model == "custom" only
#' @param number.of.V a numeric indicating the number of V variables desired for model == "custom" only
#' @param struct For use when when model == "custom". Either (1) a sub-adjacency matrix of dimension
#' (number.of.V + number.of.T) defining the topology of the V and T nodes
#' or (2) the string "random" denoting a random topology
#' @param neg.freq the frequency of negative effects for simulated effects. passed to gen.conf.coefs()
#' @param degree passed to get.custom.graph()
#' @param method passed to get.custom.graph()
#' @param simulate.confs (logical) indicating if U and K nodes should be simulated (default = TRUE).
#' @param conf.mat For use when simulate.confs = FALSE. A named dataframe with confounders (K and U nodes) in columns
#' and observations in rows. Note that the number of observations fixes the sample size for the simulation. K nodes should
#' be first in the columns represented by column names K1, K2 ... etc followed by U nodes with column names U1, U2, ... etc
#' @param sample.size a numeric value greater than 1 specifying the sample size when simulate.confs = TRUE
#' @param plot.graph (logical) when FALSE, graph is not plotted (default is TRUE)
#' @param conf.coef.ranges a list of length 4 representing the U,K,W, and Z (in that order) node effects where each list
#' element is a vector of length 2 giving the minimum and maximum effect sizes for the given node. default values follow
#' from the simulations of Yang et. al., 2017
#' @return a list of four elements
#' \describe{
#'
#' \item{data}{a dataframe of the simulated data representing the network and its confounding variables}
#' \item{Adjacency}{the adjacency matrix of the network graph}
#' \item{Effects}{the adjacency matrix containing the effects used to simulate each node}
#' \item{igraph}{the igraph object representing the network}
#' \item{true.adj}{only returned when model = "model4". It is the true adjacency matrix of the graph (the graph plotted if plot.graph = TRUE)
#'                 where as the graph represented by \eqn{$adjacency} is the graph used to emulate the model 4 network}
#'
#' }
#' @export simData.from.graph
#' @examples
#' # simulate 1000 samples from a model 1 graph with one of each type of confounding variable and
#' # plot the graph.
#' \dontrun{
#' ## using simulated confounders:
#' X=simData.from.graph(theta=0.2,
#'                      model="model1",
#'                      b0.1=0,
#'                      b.snp=0.8,
#'                      b.med=0.6,
#'                      sd.1=0.8,
#'                      conf.num.vec = c(K = 1, U = 1, W = 1, Z = 1),
#'                      simulate.confs = TRUE,
#'                      plot.graph = TRUE,
#'                      sample.size = 1000,
#'                      conf.coef.ranges=list(K=c(0.01, 0.1), U=c(0.15,0.5),
#'                                            W=c(0.15,0.5), Z=c(1, 1.5)))
#' }
#' # simulate from a model 1 graph with one 2 "known" and 8 "unknown" confounders using a set of passed confounders and
#' # plot the graph.
#' \dontrun{
#' ## using a set of passed confounders (K and U)
#' conf.mat = WBscores[,1:10]
#' #rename the columns
#' colnames(conf.mat) = c("K1","K2", paste0("U", 1:8))
#' X=simData.from.graph(theta=0.2,
#'                      model="model1",
#'                      b0.1=0,
#'                      b.snp=0.8,
#'                      b.med=0.6,
#'                      sd.1=0.8,
#'                      conf.num.vec = c(K = 2, U = 8, W = 1, Z = 1),
#'                      simulate.confs = FALSE,
#'                      plot.graph = TRUE,
#'                      conf.mat = conf.mat,
#'                      conf.coef.ranges=list(K=c(0.01, 0.1), U=c(0.15,0.5),
#'                                            W=c(0.15,0.5), Z=c(1, 1.5)))
#' }


simData.from.graph = function(model, theta, b0.1, b.snp, b.med, sd.1, conf.num.vec, number.of.T, number.of.V,
                              struct, neg.freq = 0.5, degree = 2, method = "er", simulate.confs = TRUE,
                              conf.mat, sample.size, plot.graph = TRUE, conf.coef.ranges=list(K=c(0.01, 0.1),
                                                                                              U=c(0.15,0.5),
                                                                                              W=c(0.15,0.5),
                                                                                              Z=c(1, 1.5))){

  #generate the graph from get.graph.skel
  graph.attr = gen.graph.skel(model = model, conf.num.vec = conf.num.vec, number.of.T = number.of.T,
                              number.of.V = number.of.V, struct = struct, plot.graph = plot.graph,
                              conf.coef.ranges = conf.coef.ranges, b.med = b.med, b.snp = b.snp,
                              neg.freq = neg.freq, degree = degree, method = method)

  # #for use of real confounders or setting sample size to simulate all confounders
   if(simulate.confs==TRUE){
     if(missing(sample.size)){stop("missing sample.size: must supply sample.size when simulate.confs = TRUE")}
     N = sample.size
   }else{
     if(missing(conf.mat)){stop("missing conf.mat: must supply conf.mat when simulate.confs = FALSE")}
     N = dim(conf.mat)[1]
   }

  #preallocate data matrix
  X = as.data.frame(matrix(0, nrow = N, ncol = dim(graph.attr$adjacency)[2]))
  colnames(X) = colnames(graph.attr$adjacency)

  #get the topological ordering
  topo.order = colnames(graph.attr$adjacency)[as.vector(igraph::topo_sort(graph.attr$igraph.obj))]


  for(i in 1:length(topo.order)){

    #generate V nodes in topo order
    location = match(topo.order[i], colnames(graph.attr$adjacency))
    parent.list = find.parents(Adjacency = graph.attr$adjacency, location = location)

    if(grepl("V", topo.order[i])){
      X[,location] = c(sample(c(0, 1, 2), size = N, replace = TRUE,
                              prob = c((1 -theta)^2, 2 * theta * (1 - theta), theta^2)))
    }else{
      #catch nodes with no parents of any kind
      if(sum(unlist(lapply(parent.list, is.na)))==6){
        if( isTRUE(any(sapply(c("K", "U"),grepl, x=topo.order[i])) & !missing(conf.mat)) ){
          # add in confs from conf mat:
          X[, location] = conf.mat[,match(topo.order[i], colnames(conf.mat))]
        }else{
          #generate T nodes with no parents
          X[, location] = stats::rnorm(n = N, mean = b0.1, sd = sd.1)
        }
      }else{
        #simulate all other types of nodes according to parental list, topo.order, and the effects adj
        coefs = graph.attr$effects.adj[stats::na.omit(unlist(parent.list)), topo.order[i]]
        X[, location] = stats::rnorm(n = N, mean = b0.1 +
                                       as.matrix(X[, stats::na.omit(unlist(parent.list))])%*%coefs, sd = sd.1)
      }
    }
  }

  #shuffle the bi-directed edge and return the matrix with mixed T1 and T2
  # if(model == "model4"){
  #   coinToss <- stats::rbinom(n = N, size = 1, prob = 0.5)
  #   #shuffling
  #   T1 <- rep(0, N)
  #   T1[which(coinToss == 0)] <- X$T1.a[which(coinToss == 0)]
  #   T1[which(coinToss == 1)] <- X$T1.b[which(coinToss == 1)]
  #   T2 <- rep(0, N)
  #   T2[which(coinToss == 0)] <- X$T2.a[which(coinToss == 0)]
  #   T2[which(coinToss == 1)] <- X$T2.b[which(coinToss == 1)]
  #   #now reconfigure data
  #   X = cbind.data.frame(V1 = X[,1], T1 = T1, T2 = T2, X[,-c(1:5)])
  #   #return adjusted list adding true the adjacency with the true adjacency
  #   return(list(data = X,
  #               Adjacency = graph.attr$adjacency,
  #               Effects = graph.attr$effects.adj,
  #               igraph = graph.attr$igraph.obj,
  #               true.adj = graph.attr$true.adj))
  # }else{

    return(list(data = X,
                Adjacency = graph.attr$adjacency,
                Effects = graph.attr$effects.adj,
                igraph = graph.attr$igraph.obj))
  # }
}



