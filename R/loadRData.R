

#' A simple helper function to read in .Rdata files to a specific handle
#'
#' This function utilizes the load() function to read in a .Rdata or .RData file with a specific handle
#'
#' @param fileName path/to/.Rdata
#' @export loadRData

#A simple helper function for easier reading in of .Rdata files with a specific name
loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
