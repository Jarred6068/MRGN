
#' A function to obtain the adjacency matrix from an inferred trio class
#'
#' This function takes in a the inferred model class from class.vec() and regression results from infer.trio() and returns
#' the correct adjacency matrix
#' @param model a character string as returned by class.vec()
#' @param reg.vec a vector of inference results as returned by infer.trio()
#' @param cnames (optional) a character vector of names of the adjacency columns and rows
#' @export
#' get.adj.from.class()


####################################################################\
# a simple function which returns the adjacency matrix based on the model classification
get.adj.from.class=function(model=NULL, reg.vec=NULL, cnames=NULL){

  #preallocate
  adj=matrix(0, nrow = 3, ncol = 3)
  if(is.null(cnames)){
    colnames(adj)=c("V1","T1","T2")
  }else{
    colnames(adj)=cnames
  }


  #directionalize adjacency
  if(sum(model=="M0.1" | model=="M0")==1){
    adj[1,2]=1
  }else if(model=="M0.2"){
    adj[1,3]=1
  }else if(sum(model=="M1.1" | model=="M1")==1){
    adj[1,2]=1
    adj[2,3]=1
  }else if(model=="M1.2"){
    adj[1,3]=1
    adj[3,2]=1
  }else if(sum(model=="M2.1" | model=="M2")==1){
    adj[1,2]=1
    adj[3,2]=1
  }else if(model=="M2.2"){
    adj[1,3]=1
    adj[2,3]=1
  }else if(model=="M3"){
    adj[1,2]=1
    adj[1,3]=1
  }else if(model=="M4"){
    adj[1,2]=1
    adj[1,3]=1
    adj[2,3]=1
    adj[3,2]=1
  }else if(model=="Other" & sum(reg.vec[1:4]-c(0,1,0,1))==0){
    adj[2,3]=1
    adj[3,2]=1
  }else if(model=="Other" & sum(reg.vec[1:4]-c(0,0,0,0))==0){
    adj=adj
  }

  return(adj)
}
