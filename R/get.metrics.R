

#' A function to compute the precision and recall between two graphs
#'
#' This function takes as input true and inferred model classes as is returned by class.vec() or their
#' adjacency matrices. This function wraps MRGN::get.adj.from.class() and MRPC::RecallPrecision()
#'
#' @param Truth the true model class or adjacency matrix
#' @param Inferred the inferred model class or adjacency matrix
#' @param nV the number of genetic variants in the graph (default = 1 for trios, passed to MRPC::RecallPrecision)
#' @param includeV (logical) whether to include the edges with the genetic variant V (default = TRUE, passed to MRPC::RecallPrecision)
#' @param weight.edge.directed.present a numeric specifying the weight applied true positive edges with correct directions
#' (passed to MRPC::RecallPrecision)
#' @param weight.edge.present.only a numeric specifying the weight applied to true positive edges with incorrect directions
#' (passed to MRPC::RecallPrecision)
#' @param get.adj.inf (logical) use TRUE if inferred inputs are class labels. default = FALSE
#' @param get.adj.truth (logical) use TRUE if truth inputs are class labels. default = FALSE
#' @export get.metrics
#' @importClassesFrom graph graphNEL
#' @importFrom MRPC RecallPrecision
#' @examples
#' \dontrun{
#' #example using model class labels as input:
#' result = get.metrics(Truth = "M4", Inferred = "M3", get.adj.inf = TRUE, get.adj.truth = TRUE)
#'
#' print(result)
#'
#' #example using adjacency matrix as input
#' true = get.adj.from.class('M1')
#' pred = get.adj.from.class('M4')
#'
#' result = get.metrics(Truth = true, Inferred = pred)
#'
#' print(result)
#' }

####################################################################
# a function which returns the precision, recall, and F1-score metrics
#between the true and inferred models for testing performance in simulations
get.metrics=function(Truth, Inferred, reg.vec, nV = 1, includeV = TRUE, weight.edge.directed.present = 1,
                     weight.edge.present.only = 0.5, get.adj.inf=FALSE, get.adj.truth=FALSE){

  #get adjacency, convert to vector and remove diagonal entries
  #handle input of adjacency matrix or model class label:

  if(get.adj.inf==TRUE){
    Inf.adj=get.adj.from.class(Inferred, reg.vec = reg.vec)
  }else{
    Inf.adj = Inferred
  }

  if(get.adj.truth==TRUE){
    Truth.adj=get.adj.from.class(Truth)
  }else{
    Truth.adj = Truth
  }

  Infg = as(Inf.adj, 'graphNEL')
  Truthg = as(Truth.adj, 'graphNEL')

  RP = MRPC::RecallPrecision(g1 = Truthg,
                             g2 = Infg,
                             GV = nV,
                             includeGV = includeV,
                             edge.presence = weight.edge.directed.present,
                             edge.direction = weight.edge.present.only)

  return(c(precision=RP$Precision, recall=RP$Recall))

}
