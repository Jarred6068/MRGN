
#' A function to obtain confounders for trios using the PCs from the whole-genome expression or methylation data
#'
#' This function takes in a list of trios and the scores from the PCA on whole-genome level expression
#' or methylation data and performs PC selection
#'
#' @param trios either (1) a list where each element contains a dataframe corresponding to a trio or (2) a single data
#' frame with samples in rows and trios in the columns
#' @param PCscores a dataframe of the whole-genome PCA scores with observations in rows and PCs in columns
#' @param alpha the alpha level of the pearson correlation test used by psych::corr.test() (default = 0.05)
#' @param FDR the false discovery rate (default = 0.10)
#' @param save.list (logical) if TRUE the output is saved as a .RData object (default = FALSE)
#' @param return.list (logical) if TRUE the list of the column indices of significant PCs detected for each trio
#' is returned (default = TRUE)
#' @param save.path string specifying the path name of the output
#' @return a list of the column indices of significant PCs detected for each trio
#' @export get.conf
#' @examples
#' #fast example on 40 trios
#' trio.with.conf=get.conf(trios=WBtrios[1:40], PCscores=WBscores)
#'


get.conf=function(trios=NULL, PCscores=NULL, FDR=0.10, alpha = 0.05, save.list=FALSE,
                  return.list=TRUE, save.path="/path/to/save/location"){

  #if data entered as list convert to dataframe
  if(typeof(trios)=="list" & is.null(dim(trios))){
    triomat=do.call("cbind.data.frame", trios)
  }else{
    triomat=trios
  }

  colnames(triomat) <- make.unique(colnames(triomat)) #duplicated colnames

  All.Pvalues <- list()
  List.significant.asso1 <- list()

  print("Getting list of signficant PCs for each trio...")
  #find PCs that correlate with each trio
  for (m in 1:(dim(PCscores)[1]-1)) {
    corr.PCs <- psych::corr.test(PCscores[,m], triomat, use = 'pairwise.complete.obs',
                                 adjust = "none", alpha = alpha)
    # The p values
    Pvalues <- corr.PCs$p
    Pvalues.nona <- Pvalues[!is.na(Pvalues)]
    All.Pvalues [[m]] <- Pvalues.nona
    qobj <- qvalue::qvalue(Pvalues.nona, fdr.level = FDR)

    # Significant associations
    Significant.asso <- qobj$significant
    List.significant.asso1[[m]] <- which(Significant.asso,useNames = TRUE)

  }
  print("done!")

  # columns number for trios
  start.col <- seq(1,dim(triomat)[2],3)
  end.col <- seq(3,dim(triomat)[2],3)

  number.of.trios <- dim(triomat)[2]/3
  print(paste0("number of trios = ", number.of.trios))

  Match.PCs.col <- list()
  List.Match.significant.trios <- list()
  data.sets=list()

  print("Extracting significant PCs...")
  for (ii in 1:number.of.trios) {
    # data
    data <- triomat[,start.col[ii]:end.col[ii]]
    # search in each PCs
    for (m1 in 1:(dim(PCscores)[1]-1)) {
      Match.PCs.col[[m1]] <- intersect(colnames(data),colnames(triomat)[List.significant.asso1[[m1]]])
    }
    List.Match.significant <- which(sapply(Match.PCs.col, function(e) length(e)!=0))
    List.Match.significant.trios[[ii]] <- List.Match.significant
    # if(length(List.Match.significant)!=0){
    #
    #   data.withPC <- cbind(data,PCscores[,List.Match.significant])
    #   colnames(data.withPC)[4:ncol(data.withPC)] <- paste("PC",List.Match.significant,sep="")
    #
    # }else{
    #   data.withPC <- data
    # }
    #   data.sets[[ii]]=data.withPC
  }
  print("done!")

  if(save.list==TRUE){
    print("saving...")
    save(List.Match.significant.trios, file = paste0(save.path,".RData"))
    print("done!")
  }

    return(List.Match.significant.trios)

}
