
####################################################################
#A function to classify each indicator vector returned by inf.trio()
#section 1.1 steps 4-5
class.vec=function(vec=NULL){

  #ground truth/expected test vectors
  M0=matrix(c(1,0,0,0,0,0,1,0), nrow = 2, ncol = 4, byrow = T)
  M1=matrix(c(1,1,0,1,0,1,1,1), nrow = 2, ncol = 4, byrow = T)
  M2.M4=c(1,1,1,1)
  M3=c(1,0,1,0)
  #combine gt's into matrix
  ind.mat=rbind.data.frame(M0,M1,M3,M2.M4)
  row.names(ind.mat)=c("M0.1","M0.2","M1.1","M1.2","M3","M2/M4")
  #match inferred indicator vec to gt matrix by row
  which.mod=row.names(ind.mat)[prodlim::row.match(vec[1:4], ind.mat)]
  #print(which.mod)
  #if no match is found allocate to other
  if(is.na(row.match(vec[1:4], ind.mat))){mt="Other"}

  else if(which.mod=="M2/M4"){

    #handle M2 and M4 by differentiating based on
    #marginal tests
    zp=vec[5:6]
    #if both marginals are insig or both sig we
    #allocate to M4
    if(sum(zp)==2 | sum(zp)==0){mt="M4"}

    #if only one of the marginals is sig we determine
    #which direction and allocate to M2.1 or M2.2
    else if(sum(zp==c(0,1))==2){mt="M2.1"}

    else if(sum(zp==c(1,0))==2){mt="M2.2"}

    #all other conditions allocated to class other
    else{mt="Other"}

    #if not M2/M4 and matched to gt, return model class
  }else{mt=which.mod}

  return(mt)

}