####################################################################
#a wrapper function for get.freq(), Reg(), and PermReg() to infer the trio
#combines the functions from sections 1.1-1.2
infer.trio=function(trio=NULL, gamma=0.05, alpha=0.01, nperms=10000, verbose=FALSE){

  #preallocate indicator vectors
  xp=NULL
  rp=NULL

  #preform the standard regressions and outputs t-stat and p-values
  #input is a trio with the variant in the first column followed by genes and confounders [V,T1,T2,U]
  #step 1
  pt.out=Reg(data = trio, verbose=verbose)



  #check the frequency of the minor allele using get.freq()
  #step 2
  minor=get.freq(V=trio[,1])

  #step 2.1
  if(minor<gamma){
    #preform permuted regression (section 1.2) for rare variants
    pvals=PermReg(trio = trio,
                  t.obs21 = pt.out$tvals[2],
                  t.obs22 = pt.out$tvals[4],
                  p11 = pt.out$pvals[1],
                  p12 = pt.out$pvals[3],
                  m = nperms)

  }else{
    #else return the pvals from standard reg.
    pvals=pt.out$pvals

  }

  #---steps 3-4
  #convert xp to indicator vector
  #section 1.1 step
  xp=ifelse(pvals<alpha, 1, 0)

  #preform marginal tests
  cors=c(cor.test(trio[,1], trio[,3],use="pairwise.complete.obs")$p.value,
         cor.test(trio[,1], trio[,2],use="pairwise.complete.obs")$p.value)
  #convert to indicator vector
  rp=ifelse(cors<alpha, 1, 0)
  #combine all useful stats - add indicator
  all.stats=c(append(xp, rp), append(pvals, cors), minor)

  names(all.stats)=c("b11","b21", "b12","b22", "V1:T2", "V1:T1", "pb11",
                     "pb21", "pb12","pb22","pV1:T2","pV1:T1", "Minor.freq")

  return(all.stats)
}
