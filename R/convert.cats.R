
#' A function to simplify inferred model classes
#'
#' This function takes in a vector of inferred classes from class.vec() and converts them to one of M0, M1, M2, M3, M4
#' @param input a vector of inferred classes as returned by class.vec()
#' @export convert.cats

####################################################
#a function to convert sub-categories to match 'model' in params
#from SimulateData() in MRPC package
convert.cats=function(input=NULL){
  input=replace(input, input=="M0.1" | input=="M0.2", "M0")
  input=replace(input, input=="M1.1" | input=="M1.2", "M1")
  input=replace(input, input=="M2.1" | input=="M2.2", "M2")

  return(input)
}
