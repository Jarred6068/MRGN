####################################################
#a function to recode the values for 'model' in
#SimulateData() in the MRPC package
convert.truth=function(input=NULL){

  codes=c("M0","M1","M2","M3","M4")
  unq=sort(unique(input))
  for(i in 1:length(unq)){

    input=replace(input, input==unq[i], codes[i])
  }

  return(input)

}

