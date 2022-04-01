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
  for (j in 0:2) {
    ind <- which(trio[,1] == j)
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

  nominal.p21=2 * (1 - pnorm(abs((t.obs21 - mean(Theta21))/sd(Theta21))))
  nominal.p22=2 * (1 - pnorm(abs((t.obs22 - mean(Theta22))/sd(Theta22))))

  #concat pvalues
  pvals=c(p11, nominal.p21, p12, nominal.p22)
  #return pvalue vector
  return(pvals)
}

