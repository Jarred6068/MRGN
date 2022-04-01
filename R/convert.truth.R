

#' A function to convert the value of 'model' from MRPC::simulateData()
#'
#' This function takes in a vector of true model classes of the type 'model0', 'model1','model2'... as is input to simulateData() in the MRPC package 
#' and returns one of "M0","M1","M2"....
#' @param input a vector of true model classes
#' @export
#' convert.truth()


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

