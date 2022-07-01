


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
    U=G[,u]
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
    U=G[,u]
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
    U=G[,u]
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

#' A function to simulate a graph with confounder, intermediate, and common child variables
#'
#' This function simulates a graph adjacency matrix for a desired topology such as M0 - M4, custom, or random topological
#' graphs such that they include confounder, intermediate, and common child variables.
#'
#' @param model a string specifying one of "model0","model1", "model2", "model3","model4", or "custom"
#' @param conf.num.vec a numeric vector of length 4 containing the number of confounder, known confounder, intermediate,
#' and common child variables (in that order). To exclude a variable type input a zero at the given position.
#' @param number.of.T a numeric indicating the number of T variables desired for model == "custom" only
#' @param number.of.V a numeric indicating the number of V variables desired for model == "custom" only
#' @param struct For use when when model == "custom". Either (1) a sub-adjacency matrix of dimension
#' (number.of.V + number.of.T  X  number.of.V + number.of.T) definining the topology of the V and T nodes
#' or (2) the string "random" denoting a random topology
#' @param plot.graph (logical) if TRUE the graph is plotted. default = TRUE
#' @return an adjacency matrix
#' @export gen.graph.skel
#' @examples
#' # generate a model 0 graph with one of each type of confounding variables
#' adj = gen.graph.skel(model = "model0", conf.num.vec = c(1,1,1,1))

#creates the graph skeleton for each model: contains u, k, w, and z variables in adj.
gen.graph.skel = function(model, conf.num.vec, number.of.T, number.of.V, struct, plot.graph = TRUE){

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
    A = matrix(0, nrow = sum(conf.num.vec)+number.of.T+number.of.V,
               ncol = sum(conf.num.vec)+number.of.T+number.of.V)

    row.names(A) = colnames(A) = c(paste0("V", c(1:number.of.V)),
                                   paste0("T", c(1:number.of.T)),
                                   conf.node.names)
  }else{
    #for the other M0-M4 models
    A = matrix(0, nrow = sum(conf.num.vec)+3, ncol = sum(conf.num.vec)+3)
    row.names(A) = colnames(A) = c("V1","T1","T2", conf.node.names)
  }

  #----------model-0-----------
  switch(model, model0 = {
    A[1,2] = 1

    #----------model-1-----------
  }, model1 = {

    A[1,2] = 1
    A[2,3] = 1

    #----------model-2-----------
  }, model2 = {

    A[1,2] = 1
    A[3,2] = 1

    #----------model-3-----------
  }, model3 = {

    A[1,2] = 1
    A[1,3] = 1

    #----------model-4-----------
  }, model4 = {

    A[1,2:3] = 1
    A[2,3] = 1
    A[3,2] = 1

    #----------custom-model-----------
  }, custom = {

    if(struct == "random"){
      #random topology for V and T
      A.sub = as.matrix(igraph::get.adjacency(igraph::random.graph.game(number.of.V+number.of.T, p = 0.5, directed = T)))
      #remove T --> V
      A.sub[,1:number.of.V] = 0
      #remove V --> V
      A.sub[1:number.of.V, 1:number.of.V] = 0

      #insert topology
      A[1:(number.of.V+number.of.T), 1:(number.of.V+number.of.T)] = A.sub

      #custom structure connectivity for V and T
    }else{
      A[1:(number.of.V+number.of.T), 1:(number.of.V+number.of.T)] = struct
    }

    #confounders
    for(i in 1:2){
      for(j in ((1:number.of.T)+number.of.V)){
        A[which(grepl(letter.id[i],conf.node.names))+number.of.T+number.of.V, j] = 1
      }
    }
    #intermediate
    A[which(grepl("W",conf.node.names))+number.of.T+number.of.V,
      sample(c(1:number.of.T)+number.of.V, round(runif(1,1, number.of.T)), replace = F)] = 1
    A[sample(c(1:number.of.T)+number.of.V, round(runif(1,1, number.of.T)), replace = F),
      which(grepl("W",conf.node.names))+number.of.T+number.of.V] = 1

    #common child
    for(k in ((1:number.of.T)+number.of.V)){
      A[k, which(grepl("Z",conf.node.names))+number.of.T+number.of.V] = 1
    }

    diag(A) = 0

    igraph.obj = igraph::graph_from_adjacency_matrix(A)
    if(plot.graph == TRUE){
      igraph::plot.igraph(igraph.obj, layout = igraph::layout_nicely, edge.arrow.size = 0.2)
    }

    return(list(adjacency = A, igraph.obj = igraph.obj))

  }, stop("Model not included or missing"))

  #handle confounders, intermediate, and common child vars for 5 basic topos
  #confounders

  for(i in 1:length(letter.id)){
    if(conf.num.vec[i]>0){
      if(letter.id[i] == "K" | letter.id[i] == "U"){
        for(j in 2:3){
          A[which(grepl(letter.id[i],conf.node.names))+3, j] = 1
        }
      }else if(letter.id[i] == "W"){
        #intermediate
        A[which(grepl("W",conf.node.names))+3, 3] = 1
        A[2, which(grepl("W",conf.node.names))+3] = 1
      }else if(letter.id[i] == "Z"){
        #common child
        for(k in 2:3){
          A[k, which(grepl("Z",conf.node.names))+3] = 1
        }
      }
    }
  }


  igraph.obj = igraph::graph_from_adjacency_matrix(A)
  if(plot.graph == TRUE){
    igraph::plot.igraph(igraph.obj, layout = igraph::layout_nicely, edge.arrow.size = 0.2)
  }

  return(list(adjacency = A, igraph.obj = igraph.obj))

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

#' A function to generate effects for effects in the graph
#'
#' @param parental.list a list of the parents of the given node
#' @param b.snp the desired effect for all edges with the genetic variant
#' @param b.med the desired effect for the mediation edge
#' @param num.z the number of child variables
#' @param coef.range.list a list of length 4 containing the mininum and maximum effect size for U,K,W, and Z variables
#' (in that order)
#' @param neg.freq the frequency of negative effects for a given variable
#' @return a list of length 6 containing the simulated effects for each type of variable

gen.conf.coefs=function(parental.list = NULL, b.snp=NULL, b.med=NULL, coef.range.list=NULL, neg.freq=0.5){

  #simulation of effects of variables
  #preallocate list of effects
  sim.effects = vector("list", length = 6)
  #name elements according to var type
  names(sim.effects) = c("V","T","K","U","W","Z")
  #pass the desired snp and mediation effects
  sim.effects[[1]] = b.snp
  sim.effects[[2]] = b.med
  #simulate effects of U,K,W, and Z vars from ranges in coef.range.list
  for(i in 3:6){
    if(i == 6){
      number.of.confs = length(unlist(parental.list[which(!is.na(parental.list))]))
    }else{
      number.of.confs = length(parental.list[[i]])
    }
    ct1 = rbinom(number.of.confs, 1, neg.freq)
    w1 = runif(number.of.confs, min = coef.range.list[[i-2]][1], max = coef.range.list[[i-2]][2])
    w1[which(ct1==1)]=-1*w1[which(ct1==1)]
    sim.effects[[i]]=w1
  }
  #}

  return(sim.effects)

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
#' input a zero at the given position.
#' @param number.of.T a numeric indicating the number of T variables desired for model == "custom" only
#' @param number.of.V a numeric indicating the number of V variables desired for model == "custom" only
#' @param struct For use when when model == "custom". Either (1) a sub-adjacency matrix of dimension
#' (number.of.V + number.of.T  X  number.of.V + number.of.T) definining the topology of the V and T nodes
#' or (2) the string "random" denoting a random topology
#' @param simulate.confs (logical) indicating if U and K nodes should be simulated (default = TRUE).
#' @param conf.mat For use when simulate.confs = FALSE. A dataframe with confounders (U and K nodes) in columns
#' and observations in rows. Note that the number of observations fixes the sample size for the simulation.
#' @param sample.size a numeric value greater than 1 specifying the sample size when simulate.confs = TRUE
#' @param plot.graph (logical) when FALSE, graph is not plotted (default is TRUE)
#' @param conf.coef.ranges a list of length 4 representing the U,K,W, and Z (in that order) node effects where each list
#' element is a vector of length 2 giving the minimum and maximum effect sizes for the given node. default values follow
#' from the simulations of Yang et. al., 2017
#' @return an adjacency matrix
#' @export simData.from.graph
#' @examples
#' # simulate 1000 samples from a model 1 graph with one of each confounder and
#' # plot the graph.
#' \dontrun{
#' X=simData.from.graph(theta=0.2,
#'                      model="model1",
#'                      b0.1=0,
#'                      b.snp=0.8,
#'                      b.med=0.6,
#'                      sd.1=0.8,
#'                      conf.num.vec = c(1,1,1,1),
#'                      simulate.confs = TRUE,
#'                      plot.graph = TRUE,
#'                      sample.size = 1000,
#'                      conf.coef.ranges=list(U=c(0.15,0.5), K=c(0.01, 0.1),
#'                                            W=c(0.15,0.5), Z=c(1, 1.5)))
#' }


simData.from.graph = function(model, theta, b0.1, b.snp, b.med, sd.1, conf.num.vec, number.of.T, number.of.V,
                              struct, simulate.confs = TRUE, conf.mat, sample.size, plot.graph = TRUE,
                              conf.coef.ranges=list(K=c(0.01, 0.1), U=c(0.15,0.5), W=c(0.15,0.5), Z=c(1, 1.5))){

  #generate the graph
  graph.skel = gen.graph.skel(model = model, conf.num.vec = conf.num.vec, number.of.T = number.of.T,
                              number.of.V = number.of.V, struct = struct, plot.graph = plot.graph)

  # #for use of real confounders or setting sample size to simulate all confounders
   if(simulate.confs==TRUE){
     if(missing(sample.size)){stop("missing sample.size: must supply sample.size when simulate.confs = TRUE")}
     N = sample.size
   }else{
     if(missing(conf.mat)){stop("missing conf.mat: must supply conf.mat when simulate.confs = FALSE")}
     N = dim(conf.mat)[1]
   }

  #preallocate data matrix
  X = as.data.frame(matrix(0, nrow = N, ncol = dim(graph.skel$adjacency)[2]))
  colnames(X) = colnames(graph.skel$adjacency)
  #add in the confounding variables
  if(!missing(conf.mat)){
    X[,match(colnames(conf.mat), colnames(graph.skel$adjacency))] = conf.mat
  }

  #get the topological ordering
  topo.order = colnames(graph.skel$adjacency)[as.vector(igraph::topo_sort(graph.skel$igraph.obj))]

  for(i in 1:length(topo.order)){

    #generate V nodes in topo order
    location = match(topo.order[i], colnames(graph.skel$adjacency))
    parent.list = find.parents(Adjacency = graph.skel$adjacency, location = location)

    if(grepl("V", topo.order[i])){
      X[,location] = c(sample(c(0, 1, 2), size = N, replace = TRUE,
                              prob = c((1 -theta)^2, 2 * theta * (1 - theta), theta^2)))
    }else{
      #catch nodes with no parents of any kind
      if(sum(unlist(lapply(parent.list, is.na)))==6){
        if(any(sapply(c("K", "U"),grepl, x=topo.order[i])) & missing(conf.mat)){
          #generate U,K nodes
          X[, location] = rnorm(n = N, mean = round(runif(1,0,5)), sd = 1)
        }else{
          #generate T nodes with no parents
          X[, location] = rnorm(n = N, mean = b0.1, sd = sd.1)
        }
      }else{
        #get parent idx
        parent.idx = which(unlist(lapply(parent.list, function(x) !is.na(x[1]))))
        #generate confounder effects
        coefs.list = gen.conf.coefs(parental.list = parent.list, b.snp = b.snp, b.med = b.med,
                                    coef.range.list=conf.coef.ranges)
        #find which effects to use accord to parental list
        if(any(sapply(c("W", "Z"),grepl, x=topo.order[i]))){
          which.coefs = match(substring(topo.order[i],1,1), names(coefs.list))
        }else{
          which.coefs = match(names(parent.idx), names(coefs.list))
        }
        #generate node
        X[, location] = rnorm(n = N, mean = b0.1 +
                                as.matrix(X[, unlist(parent.list[parent.idx])])%*%unlist(coefs.list[which.coefs]),
                              sd = sd.1)
      }
    }
  }

  return(X)

}



