####################################################
#a function to convert sub-categories to match 'model' in params
#from SimulateData() in MRPC package
convert.cats=function(input=NULL){
  input=replace(input, input=="M0.1" | input=="M0.2", "M0")
  input=replace(input, input=="M1.1" | input=="M1.2", "M1")
  input=replace(input, input=="M2.1" | input=="M2.2", "M2")

  return(input)
}
