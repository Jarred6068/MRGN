




#' A function to compute the precision and recall between two graphs
#'
#' This function takes as input true and inferred model classes as is returned by class.vec() or their
#' adjacency matrices. This function wraps MRGN::get.adj.from.class() and MRPC::RecallPrecision()
#'
#' @param Truth the true model class or adjacency matrix
#' @param Inferred the inferred model class or adjacency matrix
#' @param nV the number of genetic variants in the graph (default = 1 for trios, passed to MRPC::RecallPrecision)
#' @param includeV (logical) whether to include the edges with the genetic variant V (default = TRUE, passed to MRPC::RecallPrecision)
#' @param weight.edge.directed.present a numeric specifying the weight applied to true positive edges with correct directions
#' default value is 1 (passed to MRPC::RecallPrecision)
#' @param weight.edge.present.only a numeric specifying the weight applied to true positive edges with incorrect directions.
#' default is 0.5 (passed to MRPC::RecallPrecision)
#' @param get.adj.inf (logical) use TRUE if inferred inputs are class labels. default = FALSE
#' @param get.adj.truth (logical) use TRUE if truth inputs are class labels. default = FALSE
#' @return a named vector
#'        \describe{
#'        \item{Precision}{The edge-based precision.}
#'        \item{Recall}{The edge-based recall.}
#'        }
#' @export get.metrics
#' @importClassesFrom graph graphNEL
#' @importFrom Rdpack reprompt
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
#' @references
#' \insertRef{badsha2021mrpc}{MRGN}

####################################################################
# a function calculates the edge-based prediction performance between
# inferred graph and its corresponding true graph. Performance is given
# as edge-based precision and edge-based recall
get.metrics=function(Truth, Inferred, reg.vec, nV = 1, includeV = TRUE, weight.edge.directed.present = 1,
                     weight.edge.present.only = 0.5, get.adj.inf=FALSE, get.adj.truth=FALSE){

  #if the input for inferred is a string specifying the model:
  if(get.adj.inf==TRUE){
    #return the correct adjacency matrix for the given model
    Inf.adj=get.adj.from.class(Inferred, reg.vec = reg.vec)
  }else{
    #if input is already an adjacency
    Inf.adj = Inferred
  }

  #if the input for truth is a string specifying the model:
  if(get.adj.truth==TRUE){
    Truth.adj=get.adj.from.class(Truth)
  }else{
    #if input is already an adjacency
    Truth.adj = Truth
  }

  #convert adjacency matrices to graphNEL objects as required by MRPC::RecallPrecision
  Infg = as(Inf.adj, 'graphNEL')
  Truthg = as(Truth.adj, 'graphNEL')

  #compute edge-based score
  RP = RecallPrecision(g1 = Truthg,
                       g2 = Infg,
                       GV = nV,
                       includeGV = includeV,
                       edge.presence = weight.edge.directed.present,
                       edge.direction = weight.edge.present.only)

  return(c(precision=RP$Precision, recall=RP$Recall))

}



#'
#' This function is copied from MRPC::RecallPrecision
#'
#'


RecallPrecision <- function (g1, g2, GV, includeGV, edge.presence = 1.0, edge.direction = 0.5)
{
  if (is(g1, "pcAlgo"))
    g1 <- g1@graph
  if (is(g2, "pcAlgo"))
    g2 <- g2@graph
  if (is(g1, "graphNEL")) {
    m1 <- pcalg::wgtMatrix(g1, transpose = FALSE)  #need library(pcalg) for wgtMatrix
    m1[m1 != 0] <- 1
  }
  if (is(g2, "graphNEL")) {
    m2 <- pcalg::wgtMatrix(g2, transpose = FALSE)
    m2[m2 != 0] <- 1
  }
  #If ordering nodes
  if(any(colnames(m1)!=colnames(m2)))
  {
    Order_node <- match(colnames(m1),colnames(m2))
    m2 <- m2[Order_node,Order_node] #new
  }

  #Output matrix
  Evaluation_matrix <- matrix(0, nrow = 1, ncol = 2)
  colnames(Evaluation_matrix) <- c("TP", "FP")
  #index of evaluation matrix
  TP <- 1
  FP <- 2
  #library(combinat)
  #ind1=t(combinat::combn(ncol(m1), 2))

  if (includeGV) {     #Condition for recall and precision with GV
    #Calculate the edges from the true graph
    m11 <- m1
    for (i in 1:nrow(m11))
    {
      for (j in 1:ncol(m11))
      {
        if(m11[i,j]==m11[j,i])
        {
          m11[i,j] <- 0
        }
      }
    }
    #NTE=Number of edges in the truth graph
    NTE <- length(which(m11==1))

    #Calculate the edges from the inferred graph
    m22 <- m2
    for (i in 1:nrow(m22))
    {
      for (j in 1:ncol(m22))
      {
        if(m22[i,j]==m22[j,i])
        {
          m22[i,j] <- 0
        }
      }
    }
    #NIE=Number of edges in the inferred graph
    NIE <- length(which(m22==1))
    ind1 <- t(combn(ncol(m1), 2))
    for (i in seq_len(nrow(ind1))) {
      x <- ind1[i, 1]  #1st node
      y <- ind1[i, 2]  #2nd node
      #H0:No edge
      #Ha: edge exit
      #TP: Reject H0,when Ha is true (H0 is false)
      #Found edge in the inferred graph and edge exit in the true graph
      #FP: Reject H0,when H0 is true
      #Found edge in the inferred graph but no edge exit in the true graph

      #if x-->y  in true graph and x-->y in infrred graph then TP=1
      if ((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if y-->x  in true graph and y-->x in infrred graph then TP=1
      if ((m1[y, x]==1 & m1[x, y]!=1) & (m2[y, x]==1 & m2[x, y]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if x<-->ytrue graph and x<-->y in infrred graph then TP=1

      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if x-->y true graph and x<--y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]!=1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x<--y true graph and x-->y in infrred graph then TP=0.5
      if((m1[x, y]!=1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x<-->ytrue graph and x-->y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x<-->ytrue graph and x<--y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]!=1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x-->ytrue graph and x<-->y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x<--ytrue graph and x<-->y in infrred graph then TP=0.5
      if((m1[x, y]!=1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x y (no edge) true graph and x--y in infrred graph then FP=1
      if((m1[x, y]==0 & m1[y, x]==0) & (m2[x, y]==1 || m2[y, x]==1))
      {
        Evaluation_matrix[FP] <- Evaluation_matrix[FP] + 1.0
      }

      #Recall
      if (NTE != 0) {
        Recall <- Evaluation_matrix[TP]/NTE
      }

      else {
        Recall <-NA
      }


      if (NIE != 0 ) {
        Precision <- Evaluation_matrix[TP]/NIE
      }


      else {
        Precision <- NA
      }


    }
  }
  ####Without GV###
  else
  { #Condition for recall and precision without GV


    if(GV==0)
    {
      m11 <- m1  #since GV=0
      m22 <- m2  #since GV=0
    }
    else
    {
      m11 <- m1[-c(1:GV), -c(1:GV)]  #exclude GV
      m22 <- m2[-c(1:GV), -c(1:GV)]  #exclude GV
    }
    #Calculate the edges from the true graph
    for (i in 1:nrow(m11))
    {
      for (j in 1:ncol(m11))
      {
        if(m11[i,j]==m11[j,i])
        {
          m11[i,j] <- 0
        }
      }
    }
    #NTE=Number of edges in the truth graph
    NTE <- length(which(m11==1))

    #Calculate the edges from the inferred graph
    #m22 <- m2[-c(1:GV),-c(1:GV)]
    for (i in 1:nrow(m22))
    {
      for (j in 1:ncol(m22))
      {
        if(m22[i,j]==m22[j,i])
        {
          m22[i,j] <- 0
        }
      }
    }
    #NIE=Number of edges in the inferred graph
    NIE <- length(which(m22==1))
    if(GV==0)
    {
      m1 <- m1
      m2 <- m2
    }
    else
    {
      m1 <- m1[-c(1:GV), -c(1:GV)] #exclude GV
      m2 <- m2[-c(1:GV), -c(1:GV)] #exclude GV
    }
    #m1 <- m1[-c(1:GV),-c(1:GV)] #exclude GV
    #m2 <- m2[-c(1:GV),-c(1:GV)] #exclude GV

    ind1 <- t(combn(ncol(m1), 2))
    for (i in seq_len(nrow(ind1))) {
      x <- ind1[i, 1]  #1st node
      y <- ind1[i, 2]  #2nd node
      #H0:No edge
      #Ha: edge exit
      #TP: Reject H0,when Ha is true (H0 is false)
      #Found edge in the inferred graph and edge exit in the true graph
      #FP: Reject H0,when H0 is true
      #Found edge in the inferred graph but no edge exit in the true graph

      #if x-->y  in true graph and x-->y in infrred graph then TP=1
      if ((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if y-->x  in true graph and y-->x in infrred graph then TP=1
      if ((m1[y, x]==1 & m1[x, y]!=1) & (m2[y, x]==1 & m2[x, y]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if x<-->ytrue graph and x<-->y in infrred graph then TP=1

      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.presence
      }
      #if x-->y true graph and x<--y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]!=1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x<--y true graph and x-->y in infrred graph then TP=0.5
      if((m1[x, y]!=1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x<-->ytrue graph and x-->y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]!=1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x<-->ytrue graph and x<--y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]==1) & (m2[x, y]!=1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x-->ytrue graph and x<-->y in infrred graph then TP=0.5
      if((m1[x, y]==1 & m1[y, x]!=1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x<--ytrue graph and x<-->y in infrred graph then TP=0.5
      if((m1[x, y]!=1 & m1[y, x]==1) & (m2[x, y]==1 & m2[y, x]==1))
      {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] + edge.direction
      }

      #if x y (no edge) true graph and x--y in infrred graph then FP=1
      if((m1[x, y]==0 & m1[y, x]==0) & (m2[x, y]==1 || m2[y, x]==1))
      {
        Evaluation_matrix[FP] <- Evaluation_matrix[FP] + 1.0
      }

      #Recall

      if (NTE != 0) {
        Recall <- Evaluation_matrix[TP]/NTE
      }

      else {
        Recall <-NA
      }


      if (NIE != 0 ) {
        Precision <- Evaluation_matrix[TP]/NIE
      }


      else {
        Precision <- NA
      }


    }
  }

  #return(Evaluation_matrix)
  return(list(Matrix = Evaluation_matrix,
              TP = Evaluation_matrix[TP],
              FP = Evaluation_matrix[FP],
              NTE = NTE,
              NIE = NIE,
              Recall = Recall,
              Precision = Precision))
}
