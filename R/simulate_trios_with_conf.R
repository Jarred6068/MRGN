


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
    t <- rnorm(n = N, mean = b0.1 + b1.1 * P1, sd = sd.1)
  }else{
    U=G[,u]
    #gamma.u=runif(length(u), min = cs.range[1], max=cs.range[2])
    t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + U%*%gamma.u, sd = sd.1)
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
    t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + b1.2 * P2, sd = sd.1)
  }else{
    U=G[,u]
    #gamma.u=runif(length(u), min = cs.range[1], max = cs.range[2])
    t <- rnorm(n = N, mean = b0.1 + b1.1 * P1 + b1.2 * P2 + U%*%gamma.u, sd = sd.1)
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
    t <- rnorm(n = N, mean = b0.1, sd = sd.1)
  }else{
    U=G[,u]
    #gamma.u=runif(length(u), min = cs.range[1], max = cs.range[2])
    t <- rnorm(n = N, mean = b0.1 + U%*%gamma.u, sd = sd.1)
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
  ct1 = rbinom(u, 1, 0.5)
  ct2 = rbinom(u, 1, 0.5)
  w1 = runif(u, min = 0.15, max = 0.5)
  w1[which(ct1==1)]=-1*w1[which(ct1==1)]
  w2 = runif(u, min = 0.15, max = 0.5)
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

    coinToss <- rbinom(n = N, size = 1, prob = 0.5)
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






