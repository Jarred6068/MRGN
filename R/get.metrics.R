####################################################################
# a function which returns the precision, recall, and F1-score metrics
#between the true and inferred models for testing performance in simulations
get.metrics=function(Truth=NULL, Inferred=NULL, reg.vec=NULL, get.adj.inf=FALSE){

  #get adjacency, convert to vector and remove diagonal entries
  Truth.adj=as.vector(get.adj.from.class(Truth))[-c(1,5,9)]
  #handle input of adjacency matrix or model class label:
  if(get.adj.inf==FALSE){
    Inf.adj=as.vector(get.adj.from.class(Inferred))[-c(1,5,9)]
  }else{
    Inf.adj=as.vector(Inferred)[-c(1,5,9)]
  }

  #get true positive edges

  if(sum(Inf.adj)==0){
    #handle the zero adj mat case of class "Other"
    return(c(precision=NA, recall=0, specificity=1, FOR=0))

  }else{
    #handle all other classes of adjacency
    #identify which entries are edges
    true.edges=which(Truth.adj==1)
    infer.edges=which(Inf.adj==1)
    #match by edge and position
    if(length(true.edges)>length(infer.edges)){
      TP=length(match(infer.edges, true.edges))
    }else{
      TP=length(match(true.edges, infer.edges))
    }

    #get false positives
    FP=sum(Inf.adj)-TP
    #get true negative edges
    true.zeros=which(Truth.adj==0)
    infer.zeros=which(Inf.adj==0)
    #match by no edge and position
    if(length(true.zeros)>length(infer.zeros)){
      TN=length(match(infer.zeros, true.zeros))
    }else{
      TN=length(match(true.zeros,infer.zeros))
    }
    #get false negatives
    FN=sum(Inf.adj==0)-TN
    specif=TN/(TN+FP)
    #calculate metrics
    prec=TP/(TP+FP)
    recall=TP/(TP+FN)
    FOR=FN/(FN+TN)

    return(c(precision=prec, recall=recall, specificity=specif, FOR=FOR))
  }

}
