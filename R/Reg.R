


#' A function which preforms the standard regressions between the nodes in the trio. Not intended for direct use
#'
#' @param data A dataframe with at least 3 columns and the first column containing the genetic variant
#' @param verbose (logical) if TRUE the summary of the regressions is printed
#'
#' PermReg()

####################################################################
#A function which preforms the standard regressions from section 1.1 -step 1
Reg=function(data=NULL, verbose=FALSE){

  #data should be a n x (3+g) matrix with the variant in the first column
  #preallocate p and t vectors
  pvals=NULL
  tvals=NULL

  for(i in 2:3){

    #preform the regressions in step 1
    model=stats::lm(data[,i]~., data = data[,-i])
    if(verbose==TRUE){print(summary(model))}
    coefs=as.data.frame(summary(model)$coefficients)
    pvals=append(pvals, coefs$`Pr(>|t|)`[2:3])
    tvals=append(tvals, coefs$`t value`[2:3])

  }

  #name them according to coefficient index
  names(pvals)=c("p11","p12","p21","p22")
  names(tvals)=c("tobs11","tobs12","tobs21","tobs22")

  #output pvalues and t-stats for b11,b21,b12,b22
  return(pt.list=list(pvals=pvals,tvals=tvals))

}

