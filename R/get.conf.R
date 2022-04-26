
#' A function to obtain confounders for trios using the PCs from the whole-genome expression or methylation data
#'
#' This function takes in a list of trios and the rotated variables from the PCA on whole-genome level expression
#' or methylation data and performs PC selection
#'
#' @param trios either a list of data frames containing the trios or a single data frame with samples in rows and
#' trios in the columns
#' @param PCscores the rotated variables from the whole-genome level PCA
#' @param FDR the rejection threshold for the pearson correlation test between each trio and the rotated variables
#' @param save.list (logical) if TRUE the output is saved as a .RData object
#' @param save.path string specifying the path name of the output
#' @return a list of dataframes where each element is a dataframe corresponding to a trio and its confounding PCs
#' @export get.conf
#' @examples
#' #fast example on 40 trios
#' trio.with.conf=get.conf(trios=WBtrios[1:40], PCmat=WBscores)
#'


get.conf=function(trios=NULL, PCscores=NULL, FDR=0.10, save.list=FALSE, save.path="/path/to/save/location"){

  #if data entered as list convert to dataframe
  if(typeof(trios)=="list" & is.null(dim(trios))){
    triomat=do.call("cbind.data.frame", trios)
  }else{
    triomat=trios
  }

  colnames(triomat) <- make.unique(colnames(triomat)) #duplicated colnames

  All.Pvalues <- list()
  List.significant.asso1 <- list()

  print("Getting list of PCs for each trio...")
  #find PCs that correlate with each trio
  for (m in 1:(dim(PCscores)[1]-1)) {
    corr.PCs <- psych::corr.test(PCscores[,m], triomat, use = 'pairwise.complete.obs', adjust = "none")
    # The p values
    Pvalues <- corr.PCs$p
    Pvalues.nona <- Pvalues[!is.na(Pvalues)]
    All.Pvalues [[m]] <- Pvalues.nona
    qobj <- qvalue::qvalue(Pvalues.nona, fdr.level = FDR)

    # Significant associations
    Significant.asso <- qobj$significant
    List.significant.asso1[[m]] <- which(Significant.asso,useNames = TRUE)
    print(List.significant.asso1)
  }
  print("done!")
  # columns number for trios
  start.col <- seq(1,dim(triomat)[2],3)
  end.col <- seq(3,dim(triomat)[2],3)

  trios <- dim(triomat)[2]/3
  print(trios)

  Match.PCs.col <- list()
  List.Match.significant.trios <- list()
  data.sets=list()

  print("Assembling trios with confounders...")
  for (ii in 1:trios) {
    # data
    data <- triomat[,start.col[ii]:end.col[ii]]
    # search in each PCs
    for (m1 in 1:(dim(PCscores)[1]-1)) {
      Match.PCs.col[[m1]] <- intersect(colnames(data),colnames(triomat)[List.significant.asso1[[m1]]])
    }
    List.Match.significant <- which(sapply(Match.PCs.col, function(e) length(e)!=0))
    List.Match.significant.trios[[ii]] <- List.Match.significant
    if(length(List.Match.significant)!=0)
    {
      data.withPC <- cbind(data,PCscores[,List.Match.significant])
      colnames(data.withPC)[4:ncol(data.withPC)] <- paste("PC",List.Match.significant,sep="")
    }
    else
    {
      data.withPC <- data
    }
    data.sets[[ii]]=data.withPC
  }
  print("done!")

  if(save.list==TRUE){
    print("saving...")
    save(data.sets, file = paste0(save.path,".RData"))
    print("done!")
  }

  return(data.sets)

}
