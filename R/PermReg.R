

#' this function wraps PermReg.helper.fn() to preform the permuted regression. Not intended for direct use
#'
#' @param trio A dataframe with at least 3 columns and the first column containing the genetic variant
#' @param t.obs21 the t-statistic from Reg() represnting the test on beta21
#' @param t.obs22 the t-statistic from Reg() representing the test on beta22
#' @param p11 the p-value from the test on beta21
#' @param p12 the p-value from the test on beta22
#' @param m the number of permutations to perform
#' PermReg()


####################################################################
#this function uses PermReg.helper.fn() to preform the permuted regression
#analysis for rare variants - section 1.2
PermReg=function(trio=NULL, t.obs21=NULL, t.obs22=NULL, p11=NULL, p12=NULL, m=NULL){

  #preallocate a matrix of indicies ranging from 1:sample size
  #we will shuffle these numbers later to get the permutations within genotype
  #need 2: one for each regression: Section 1.2 steps 1-1.2 and eqn (3) and (4)
  mediator_perm1=matrix(c(1:dim(trio)[1]), nrow=dim(trio)[1], ncol = m)
  mediator_perm2=matrix(c(1:dim(trio)[1]), nrow=dim(trio)[1], ncol = m)

  #preallocate all permutations
  #shuffle within each genotype
  ngeno <- length (table (trio[,1]))
  geno <- unique (trio[,1])
  for (j in 1:ngeno) {
    ind <- which(trio[,1] == geno[j])
    if (length(ind) > 1) {
      mediator_perm1[ind, ] <- apply(mediator_perm1[ind, ], 2, sample)
      mediator_perm2[ind, ] <- apply(mediator_perm2[ind, ], 2, sample)
    }
  }

  #this section handles data that may not have confounding variables included
  if(dim(trio)[2]>3){
    confounders=trio[,-c(1:3)]
  }else{
    confounders=NULL
  }
  #preforms all permutations of eqn (3) in parallel
  #outputs Theta21
  Theta21=apply(mediator_perm1, 2, PermReg.help.fn,
                V=trio[,1],
                T1=trio[,2],
                T2=trio[,3],
                U=confounders,
                coln=colnames(trio),
                response="T1")
  #preforms all permutations of eqn (4) in parallel
  #outputs Theta22
  Theta22=apply(mediator_perm2, 2, PermReg.help.fn,
                V=trio[,1],
                T1=trio[,2],
                T2=trio[,3],
                U=confounders,
                coln=colnames(trio),
                response="T2")

  #Step 2.1 - calculating the nominal p-values using Theta21 and Theta22

  nominal.p21=2 * (stats::pnorm(abs((t.obs21 - mean(Theta21))/stats::sd(Theta21)), lower.tail = F))
  nominal.p22=2 * (stats::pnorm(abs((t.obs22 - mean(Theta22))/stats::sd(Theta22)), lower.tail = F))

  #concat pvalues
  pvals=c(p11, nominal.p21, p12, nominal.p22)
  #return pvalue vector
  return(pvals)
}

