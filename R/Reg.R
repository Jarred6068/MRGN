####################################################################
#A function which preforms the standard regressions from section 1.1 -step 1
Reg=function(data=NULL, verbose=FALSE){

  #data should be a n x (3+g) matrix with the variant in the first column
  #preallocate p and t vectors
  pvals=NULL
  tvals=NULL

  for(i in 2:3){

    #preform the regressions in step 1
    model=lm(data[,i]~., data = data[,-i])
    if(verbose==TRUE){print(summary(model))}
    coefs=as.data.frame(summary(model)$coefficients)
    pvals=append(pvals, coefs$`Pr(>|t|)`[2:3])
    tvals=append(tvals, coefs$`t value`[2:3])

  }

  #name them according to coefficient index
  names(pvals)=c("p11","p21","p12","p22")
  names(tvals)=c("tobs11","tobs21","tobs12","tobs22")

  #output pvalues and t-stats for b11,b21,b12,b22
  return(pt.list=list(pvals=pvals,tvals=tvals))

}
