
#' A function to obtain confounders for trios using the PCs from the whole-genome expression or methylation data
#'
#' This function takes in a list of trios and the scores from the PCA on whole-genome level expression
#' or methylation data and performs PC selection
#'
#' @param trios either (1) a list where each element contains a dataframe corresponding to a trio or (2) a single data
#' frame with samples in rows and trios in the columns
#' @param PCscores a dataframe of the whole-genome PCA scores with observations in rows and PCs in columns
#' @param blocksize the number of columns to use in each block of correlation calculations passed to propagate::bigcor
#' @param apply.qval (logical) \eqn{default = TRUE} applies the qvalue adjustment to each set of correlations between a PC
#' and the columns of \eqn{trios}. If \eqn{FALSE}, the bonferroni correction is applied.
#' @param FDR the false discovery rate (default = 0.10)
#' @param lambda he value of the tuning parameter to estimate \eqn{pi_0}. Must be between 0,1. If is.null(lambda) the default
#' passed to \eqn{adjust.q()} is \eqn{seq(0.5, max(pvalues), 0.05)}
#' @param alpha the test threshold for the bonferroni correacted pvalues when \eqn{apply.qval = FALSE}
#' @param return.for.trios (logical) if \eqn{TRUE} the column indices of the PCs associated with each trio are returned. If FALSE
#' the column indices of PCs associated with each column of "trios" is returned. default=TRUE
#' @param save.list (logical) if TRUE the output is saved as a .RData object (default = FALSE)
#' @param return.list (logical) if TRUE the list of the column indices of significant PCs detected for each trio
#' is returned (default = TRUE)
#' @param save.path string specifying the path name of the output
#' @return a list of 4 elements containing:
#'   \describe{
#'   \item{sig.asso.pcs}{default: a list of length = The number of trios of the column indices of significant PCs detected for each trio
#'   Alternatively, if return.for.trios = FALSE, a list of length=ncol(trios) of the column indices of significant PCs detected
#'   for each column of trios}
#'   \item{pvalues}{a matrix of dimension \eqn{ncol(PCscores) X ncol(trios)} of the pvalues from the pearson correlation test on each set of correlations}
#'   \item{qvalues}{a matrix of dimension \eqn{ncol(PCscores) X ncol(trios)} of qvalues from each set of correlations}
#'   \item{cors}{the matrix of the calculated pairwise correlations of dimension \eqn{ncol(PCscores) X ncol(trios)}}
#'   \item{sig}{A matrix of logical values of dimension \eqn{ncol(PCscores) X ncol(trios)}}
#' }
#'  a list of length = # of trios of the column indices of significant PCs detected for each trio. Alternatively,
#' if return.for.trios = FALSE, a list of length=ncol(trios) of the column indices of significant PCs detected for each column of trios
#' @export get.conf
#' @import propagate
#' @import qvalue
#' @examples
#'
#' \dontrun{
#' #fast example on 40 trios using qvalue correction
#' trio.conf=get.conf(trios=WBtrios[1:40], PCscores=WBscores, blocksize=10)
#'
#' #fast example on 40 trios using bonferroni correction
#' trio.conf2=get.conf(trios=WBtrios[1:40], PCscores=WBscores, blocksize=10, apply.qval=FALSE)
#'}


get.conf=function(trios=NULL, PCscores=NULL, blocksize=2000, apply.qval=TRUE, FDR=0.10, lambda=NULL, alpha=0.05,
                  return.for.trios=TRUE, save.list=FALSE, return.list=TRUE, save.path="/path/to/save/location"){

  #if data entered as list convert to dataframe
  if(typeof(trios)=="list" & is.null(dim(trios))){
    triomat=do.call("cbind", trios)
  }else{
    triomat=trios
  }

  #get indicies of each trio in triomat
  trio.indices=cbind( start.col = seq(1,dim(triomat)[2],3),
                      end.col = seq(3,dim(triomat)[2],3))
  #make colnames unique to avoid error in bigcor
  colnames(triomat) <- make.unique(colnames(triomat)) #duplicated colnames
  #get the sample sizes used in each pairwise correlation calculation
  sample.sizes=apply(triomat,2, function(x) length(S4Vectors::na.omit(x)))

  #calculate the correlations of each PC with the trios
  cormat=propagate::bigcor(triomat, PCscores, verbose = T, use="pairwise.complete.obs", size = blocksize)
  #perform the pearson correlation test, apply qvalue, and return which are significant
  pr.mat=apply(as.data.frame(cormat[,]), 2, p.from.cor, n=sample.sizes)
  #extract the pvalues and correlations
  r.mat = sapply(pr.mat, function(x) x$cor)
  p.mat = sapply(pr.mat, function(x) x$pvalue)
  #perform the qvalue correction
  if(apply.qval==TRUE){
    qsig.mat = apply(p.mat, 2, adjust.q, fdr = FDR, lambda = lambda)
    sig.mat = sapply(qsig.mat, function(x) x$significant)
    q.mat = sapply(qsig.mat, function(x) x$qvalue)
    p.adj.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
  }else{
    p.adj.mat = apply(p.mat, 2, stats::p.adjust, method="bonferroni")
    sig.mat = apply(p.adj.mat, 2, function(x) x<=0.05)
    q.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
  }
  #naming
  colnames(sig.mat)=colnames(q.mat)=colnames(r.mat)=colnames(p.mat)=colnames(p.adj.mat)=paste0("PC",1:dim(PCscores)[2])
  row.names(sig.mat)=row.names(q.mat)=row.names(r.mat)=row.names(p.mat)=row.names(p.adj.mat)=make.unique(colnames(triomat))
  #find the PCs that correlated with every column of the trio matrix
  sig.asso.pcs=apply(sig.mat, 1, function(x){list(which(x))})
  #return the significant PCs for each trio or for each column in "trios"
  if(return.for.trios==TRUE){
  final.list.sig.asso.pcs=apply(trio.indices, 1,
                                function(x,y){ list(unique(unlist(y[x[1]:x[2]]))) },
                                y=sig.asso.pcs)
    if(return.list==TRUE){
      out.list = list(sig.asso.pcs=final.list.sig.asso.pcs,
                      pvalues=t(p.mat),
                      qvalues=t(q.mat),
                      cors = t(r.mat),
                      sig = t(sig.mat),
                      adj.p = t(p.adj.mat))
      return(out.list)
    }
    if(save.list==TRUE){
      save(out.list, file = paste0(save.path,".RData"))
    }


  }else{

    if(return.list==TRUE){
      out.list=list(sig.asso.pcs = sig.asso.pcs,
                    pvalues = t(p.mat),
                    qvalues = t(q.mat),
                    cors = t(r.mat),
                    sig = t(sig.mat),
                    adj.p = t(p.adj.mat))
      return(out.list)
    }
    if(save.list==TRUE){
      save(out.list, file = paste0(save.path,".RData"))
    }
  }
}

#' A function to calculate the pearson correlation p-values and apply the qvalue correction
#'
#' @param p a vector of p.values to be passed to qvalue::qvalue()
#' @param fdr the false discovery rate
#' @param lambda The value of the tuning parameter to estimate \eqn{pi_0}. Must be in \eqn{[0,1)}
#' @return an n X 3 dataframe containing the correlations, qvalues, and significance (logical)

adjust.q=function(p, fdr, lambda){
  #apply qvalue correction
  if(is.null(lambda)){
    qval.str=qvalue::qvalue(p, fdr.level = fdr, lambda=seq(0.05, max(p), 0.05))
  }else{
    qval.str=qvalue::qvalue(p, fdr.level = fdr, lambda=lambda)
  }
  #extrac significance and qvalues and return
  sig=qval.str$significant
  qval=qval.str$qvalues
  return(cbind.data.frame(significant=sig, qvalue=qval))
}

#' A function to calculate the pearson correlation p-values
#'
#' @param r either (1) a single pearson correlation coefficient or (2) a vector of correlation coefficients
#' @param n the sample size (if length(r)>1 then n must be a vector if sample sizes vary)
#' @return an n X 3 dataframe containing the correlations and pvalues

p.from.cor=function(r, n){
  #calculate t values from vector of correlations
  t=r*sqrt((n-2)/(1-r^2))
  #obtain the pvalue from two-tailed test
  p=2 * (1 -stats::pt(q = abs(t), df = n-2))
  return(cbind.data.frame(pvalue=p,cor=r))
}












