####################################################################
#-------------------Helper-function-for-permuted-reg---------------
PermReg.help.fn=function(perm.map=NULL, V=NULL, T1=NULL, T2=NULL, U=NULL, coln=NULL, response=NULL){

  #written for parallelization
  #permutes
  if(response=="T2"){
    if(is.null(U)){
      new.data=cbind.data.frame(V, T1[perm.map], T2)
    }else{
      new.data=cbind.data.frame(V, T1, T2[perm.map], U)
    }

  }else{
    if(is.null(U)){
      new.data=cbind.data.frame(V, T1, T2[perm.map])
    }else{
      new.data=cbind.data.frame(V, T1, T2[perm.map], U)
    }

  }

  if(isFALSE(is.null(coln))){colnames(new.data)=coln}

  #run regression
  if(response=="T2"){

    #section 1.2 step 1.1 eqn (4)
    coef.mat=as.data.frame(summary(lm(new.data[,3]~., data=new.data[,-3]))$coefficients)
    wald.stat=coef.mat$`t value`[3]
    #print(coef.mat)

  }else{

    #section 1.2 step 1.1 eqn (3)
    coef.mat=as.data.frame(summary(lm(new.data[,2]~., data=new.data[,-2]))$coefficients)
    wald.stat=coef.mat$`t value`[3]
    #print(coef.mat)

  }


  #return the wald stat: section 1.2 step 1.2
  #since it is run in parallel the output in PermReg() is a vector
  return(wald.stat)

}
