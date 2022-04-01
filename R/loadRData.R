#A simple helper function for easier reading in of .Rdata files with a specific name
loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
