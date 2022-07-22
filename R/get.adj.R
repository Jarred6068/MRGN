

#' A function to get the adjacency matrix from the edge inference returned by infer.trio()
#'
#' This function takes in the output from infer.trio() and converts it to an adjacency matrix for the trio
#' @param inf a vector of length 14 as return by infer.trio()
#' @return a 3 X 3 adjacency matrix for the trio
#' @export get.adj
#' @examples
#' \dontrun{
#' inf = infer.trio(M1trio)
#' Adj = get.adj(inf)
#' }

get.adj = function(inf){
  #preallocate adjacency
  A = as.data.frame(matrix(0, nrow = 3, ncol = 3))
  colnames(A) = row.names(A) = c("V1","T1","T2")
  A[1,2:3] = inf[c(1,3)]
  A[2,3] = inf[2]
  A[3,2] = inf[4]

  return(A)
}
