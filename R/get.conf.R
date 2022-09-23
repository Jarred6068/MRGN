
#' A function to select a set of confounding variables for an input dataset
#'
#' This function takes in a set of variables represented in the columns of \eqn{data} and a pool of potential confounding variables (e.g whole-genome expression or
#' methylation PC scores) and selects the confounding variables associated with each column \eqn{data} by testing the marginal associations of each column of \eqn{data} with
#' the covariates in \eqn{cov.pool}. This method is a simpler version of \eqn{get.conf.trios} which does not require the input data be
#' organized into trios.
#'
#' @param data a matrix or dataframe of size \eqn{n samples X p variables}.
#' @param cov.pool a matrix or dataframe of size \eqn{n samples X q covariates}whose columns represent the possible confounding variables to be selected from
#' @param blocksize the number of columns of \eqn{data} to use in each block of correlation calculations passed to propagate::bigcor
#' @param apply.qval (logical) \eqn{default = TRUE} applies the qvalue adjustment to each set of correlations between a PC
#' and the columns of \eqn{trios}. If \eqn{FALSE}, the Bonferroni correction is applied.
#' @param selection_fdr the false discovery rate (default = 0.05) for selecting confounders when apply.qval = TRUE
#' @param lambda When apply.qval = TRUE, lambda is the set of cut off points of the tuning parameter to estimate \eqn{pi_0}. Must be between 0,1. If is.null(lambda) the default
#' passed to \eqn{adjust.q()} is \eqn{seq(0.5, max(pvalues), 0.05)}
#' @param alpha the test threshold for the bonferroni correacted pvalues when \eqn{apply.qval = FALSE}
#' @param save.list (logical) if TRUE the output is saved as a .RData object (default = FALSE)
#' @param save.path string specifying the path name of the output list to be save as a .RData structure
#' @return a list of 6 elements containing:
#'   \describe{
#'   \item{sig.asso.covs}{a list of \eqn{length=ncol(data)}. Each element in the list is a vector of of indices corresponding to the columns
#'   of \eqn{cov.pool} that are significantly associated to each input variable (column) in \eqn{data}}
#'   \item{pvalues}{a matrix of dimension \eqn{ncol(cov.pool) X ncol(data)} of the pvalues from the pearson correlation test on each set of correlations}
#'   \item{qvalues}{a matrix of dimension \eqn{ncol(cov.pool) X ncol(data)} of qvalues from each set of correlations}
#'   \item{cors}{the matrix of the calculated pairwise correlations of dimension \eqn{ncol(cov.pool) X ncol(trios)}}
#'   \item{sig}{A matrix of logical values of dimension \eqn{ncol(cov.pool) X ncol(data)}}
#'   \item{adj.p}{A matrix of dimension \eqn{ncol(cov.pool) X ncol(data)} of the adjusted p-values if \eqn{apply.qval = FALSE}}
# #'   \item{filtered_covs}{A list of length = the number of trios giving the filtered covariates if \eqn{filter_int_child = TRUE}}
#' }
#' @export get.conf.matrix
#' @references
#'     \insertAllCited{}
#' @import propagate
#' @import qvalue
#' @examples
#'
#' \dontrun{
#' #fast example
#' trio.conf=get.conf.matrix(data=do.call('cbind',WBtrios[1:40]),
#'                           cov.pool=WBscores,
#'                           blocksize=10)
#'
#' # example with manual filtering
#'
#' triomat = do.call('cbind', WBtrios[1:40])
#' snp.idx = seq(1, dim(triomat)[2]-2, 3)
#' snps = triomat[,snp.idx]
#' genes = triomat[,-snp.idx]
#'
#' gene.confs=get.conf.matrix(data=genes,
#'                           cov.pool=WBscores,
#'                           blocksize=10)
#'
#' snp.confs=get.conf.matrix(data=snps,
#'                           cov.pool=WBscores,
#'                           blocksize=10)
#'
#'
#'
#'
#'
#'}


get.conf.matrix=function(data=NULL, cov.pool=NULL, blocksize=2000, apply.qval=TRUE, selection_fdr=0.05,
                         lambda=NULL, alpha=0.05, save.list=FALSE, save.path="/path/to/save/location"){

  #====================================Preprocessing====================================
  data = as.data.frame(data)
  cn.data = colnames(data)
  cov.pool = as.data.frame(cov.pool)
  cn.cov.pool = colnames(cov.pool)

  #catch columns with only 2 or less non-NA values and return their indices
  non.na.vals = apply(data, 2, function(x) length(stats::na.omit(x)))
  if(any(non.na.vals<=2)){
    stop(paste0("some columns of \"trios\" contain <= 2 non-NA values: The column(s) is/are ", paste0(which(non.na.vals<=2))))
  }

  #====================================Calculate-Correlations====================================
  #make colnames unique to avoid error in bigcor
  colnames(data) = make.unique(colnames(data)) #duplicated colnames
  #get the sample sizes used in each pairwise correlation calculation
  sample.sizes = apply(data, 2, function(x) length(na.omit(x)))
  #calculate the correlations of each PC with the trios
  cormat = propagate::bigcor(data, cov.pool, verbose = T, use = "pairwise.complete.obs", size = blocksize)

  #===================================Confounder-Selection=======================================

  message("Applying method = \"correlation\" on confounding variables...")
  #perform the pearson correlation test, apply qvalue, and return which are significant
  pr.mat=apply(as.data.frame(cormat[,]), 2, p.from.cor, n=sample.sizes)
  #extract the pvalues and correlations
  r.mat = sapply(pr.mat, function(x) x$cor)
  p.mat = sapply(pr.mat, function(x) x$pvalue)





  #-----------------apply-correction-to-pvalues-------------------
  #perform the correction
  if(apply.qval==TRUE){
    #qvalue correction
    qsig.mat = apply(p.mat, 2, adjust.q, fdr = selection_fdr, lambda = lambda)
    sig.mat = sapply(qsig.mat, function(x) x$significant)
    q.mat = sapply(qsig.mat, function(x) x$qvalue)
    p.adj.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
  }else{
    #bonferroni correction
    p.adj.mat = apply(p.mat, 2, stats::p.adjust, method="bonferroni")
    sig.mat = apply(p.adj.mat, 2, function(x) x<=alpha)
    q.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
  }



  #-------------------obtain-final-list-of-significant-covs----------------
  #find the covs that correlated with every column of the trio matrix
  sig.asso.covs=apply(sig.mat, 1, function(x){which(x)})


  #=================================Organize and return/save output====================================
  colnames(sig.mat)=colnames(q.mat)=colnames(r.mat)=colnames(p.mat)=colnames(p.adj.mat)=cn.cov.pool
  row.names(sig.mat)=row.names(q.mat)=row.names(p.mat)=row.names(p.adj.mat)=cn.data



  out.list=list(sig.asso.covs = sig.asso.covs,
                pvalues = t(p.mat),
                qvalues = t(q.mat),
                cors = t(r.mat),
                sig = t(sig.mat),
                adj.p = t(p.adj.mat))


    if(save.list==TRUE){
      save(out.list, file = paste0(save.path,".RData"))
    }

    return(out.list)

}






#' A function to select a set of confounding variables for each trio
#'
#' This function takes in a list or dataframe of trios and a matrix of potential confounding variables (e.g whole-genome expression or
#' methylation PC scores) and selects the confounding variables associated with each trio using either a correlation-based procedure of \insertCite{badsha2019learning}{MRGN}
#' or regression-based procedure of \insertCite{yang2017identifying}{MRGN}. Additionally, this function allows the user to filter common child and intermediate confounding
#' variables by removing selected confounders that are significantly associated with the genetic variant.
#'
#' @param trios either (1) a list where each element contains a dataframe corresponding to a trio or (2) a single data
#' frame with samples in rows and trios in the columns. The assumed structure of each trio is that the Genetic variant
#' is the first column followed by the two molecular phenotypes.
#' @param cov.pool a dataframe of the whole-genome PCA scores with observations in rows and PCs in columns
#' @param blocksize the number of columns to use in each block of correlation calculations passed to propagate::bigcor
#' @param selection_fdr the false discovery rate (default = 0.05) for selecting confounders
#' @param filter_int_child (logical) to determine if common child and intermediate confounding variables should be filtered out
#' (removed) from the significant confounders for each trio. Note: only used when return.for.trios==TRUE. Default = FALSE.
#' @param filter_fdr the false discovery rate for filtering common child
#' and intermediate confounding variable (only used when method = "regression").
#' @param lambda the cut off points of the tuning parameter to estimate \eqn{pi_0}. Must be between 0,1. If is.null(lambda) the default
#' passed to \eqn{adjust.q()} is \eqn{seq(0.5, max(pvalues), 0.05)}
#' @param alpha the test threshold for the bonferroni correacted pvalues when \eqn{apply.qval = FALSE}
#' @param method the method to calculate the associated pcs. A character string specifying either "correlation" or "regression".
#' Note that when \eqn{method = "regression"}, \eqn{return.for.trios} is forced to TRUE
#' @param save.list (logical) if TRUE the output is saved as a .RData object (default = FALSE)
#' @param save.path string specifying the path name of the output
#' @return a list of 6 elements containing:
#'   \describe{
#'   \item{sig.asso.pcs}{default: a list of length = The number of trios of the column indices of significant PCs detected for each trio
#'   Alternatively, if return.for.trios = FALSE, a list of length=ncol(trios) of the column indices of significant PCs detected
#'   for each column of trios}
#'   \item{pvalues}{a matrix of dimension \eqn{ncol(cov.pool) X ncol(trios)} of the pvalues from the pearson correlation test on each set of correlations}
#'   \item{qvalues}{a matrix of dimension \eqn{ncol(cov.pool) X ncol(trios)} of qvalues from each set of correlations}
#'   \item{cors}{the matrix of the calculated pairwise correlations of dimension \eqn{ncol(cov.pool) X ncol(trios)}}
#'   \item{sig}{A matrix of logical values of dimension \eqn{ncol(cov.pool) X ncol(trios)}}
#'   \item{adj.p}{A matrix of dimension \eqn{ncol(cov.pool) X ncol(trios)} of the adjusted p-values if \eqn{apply.qval = FALSE}}
#'   \item{filtered_covs}{A list of length = the number of trios giving the filtered covariates if \eqn{filter_int_child = TRUE}}
#' }
#' @export get.conf.trios
#' @references
#'     \insertAllCited{}
#' @import propagate
#' @import qvalue
#' @examples
#'
#' \dontrun{
#' #fast example on 40 trios using qvalue correction and correlation method
#' trio.conf=get.conf.trios(trios=WBtrios[1:40],
#'                          cov.pool=WBscores,
#'                          blocksize=10,
#'                          method = "correlation")
#'
#' #fast example on 40 trios using bonferroni correction and correlation method
#' trio.conf2=get.conf.trios(trios=WBtrios[1:40],
#'                           cov.pool=WBscores,
#'                           blocksize=10,
#'                           method = "regression")
#'
#' #fast example on 40 trios using the qvalue correction and the regression method
#' #Note: this method is slower than method = "correlation"
#' trio.conf3=get.conf.trio(trios=WBtrios[1:40],
#'                          cov.pool=WBscores,
#'                          blocksize=10,
#'                          filter_int_child = T,
#'                          method = "regression")
#'}


get.conf.trios=function(trios=NULL, cov.pool=NULL, blocksize=2000, selection_fdr=0.05,
                        filter_int_child = FALSE, filter_fdr = 0.1, lambda=NULL, alpha=0.05,
                        method = c("correlation, regression"), save.list=FALSE, save.path="/path/to/save/location"){

  #====================================Preprocessing====================================
  #ensure the cov pool is a dataframe
  #get the cov pool column names
  cov.pool = as.data.frame(cov.pool)
  cn.cov.pool = colnames(cov.pool)

  #this block takes the input data as list or df and converts it to one or the other.
  #the goal is to have the data in both matrix and list format for later steps
  #also calculates the total number of trios and reports it back to the user as a
  #verification
  if(typeof(trios)=="list" & is.null(dim(trios))){
    print("entered data is in list format...")
    trio.list = trios
    triomat=do.call("cbind", trios)
    num.trios = length(trios)
  }else{
    print("entered data is in df/array format...")
    triomat=trios
    num.trios = dim(trios)/3
    #get index of trios and convert to list
    trio.cols = split(cbind.data.frame(variant.idx = seq(1, dim(triomat)[2]-2, 3),
                     gene1.idx = seq(2, dim(triomat)[2]-1, 3),
                     gene2.idx = seq(3, dim(triomat)[2], 3)), seq(1, num.trios, 1))

    #obtain list of trios where each entry is a n X 3 dataframe of a trio
    trio.list = lapply(trio.cols, function(x,y) y[,unlist(x)], y = triomat)
  }
  #get a list of trio indices for later use
  trio.indices=cbind( start.col = seq(1,dim(triomat)[2],3),
                      end.col = seq(3,dim(triomat)[2],3))

  snp.idx = seq(1, dim(triomat)[2]-2, 3)

  trn = paste0("trio_", c(1:num.trios))

  print(paste0("detected ", num.trios, " trios and ", dim(cov.pool)[2], " covariates in the candidate pool..."))



  #catch columns with only 2 or less non-NA values and return their indices
  non.na.vals = apply(triomat, 2, function(x) length(stats::na.omit(x)))
  if(any(non.na.vals<=2)){
    stop(paste0("some columns of \"trios\" contain <= 2 non-NA values: The column(s) is/are ", paste0(which(non.na.vals<=2))))
  }

  #====================================Calculate/Test-Correlations====================================
  #this step is get the correlation matrix of size (3 X number of trios) x Number of covs
  #we then get the marginal test pvalues by testing the marginal relationship between each
  #variant or gene with each cov to produce a matrix of size (3 X number of trios) x Number of covs
  #containing the pvalues

  #make colnames unique to avoid error in bigcor
  colnames(triomat) = make.unique(colnames(triomat)) #duplicated colnames
  #get the sample sizes used in each pairwise correlation calculation
  sample.sizes = apply(triomat, 2, function(x) length(na.omit(x)))

  print(paste0("Calculating correlation matrix of size ", num.trios*3, " x ", dim(cov.pool)[2],
               " using ", round((num.trios*3)/blocksize), " blocks"))
  #calculate the correlations of each PC with the trios
  cormat = propagate::bigcor(triomat, cov.pool, verbose = T, use = "pairwise.complete.obs", size = blocksize)
  #perform marginal test
  pr.mat=apply(as.data.frame(cormat[,]), 2, p.from.cor, n=sample.sizes)
  #correlation matrix and pvalue matrix
  cor.r.mat= sapply(pr.mat, function(x) x$cor)
  cor.p.mat = sapply(pr.mat, function(x) x$pvalue)


#==================================filtering-and-covariate-selection========================================
  #this block of code handles both the filtering and covariate selection steps
  #the steps to the final result are different depending on method
  #the switch separates the process according to each method.

  switch(method, regression = {
    #regression method - filtering and selection covs
    result.list = regression.method(covpool = cov.pool,
                                    trio.mat = triomat,
                                    cor.pvalue.matrix = cor.p.mat,
                                    filtering = filter_int_child,
                                    filter.fdr = filter_fdr,
                                    selection.fdr = selection_fdr,
                                    lam.seq = lambda,
                                    triolist = trio.list)

    #naming
    names(result.list$final.list.covs) = names(result.list$filtered.covs) = trn
    colnames(result.list$reg.pvalues) = colnames(result.list$cor.pvalues) = colnames(result.list$q.values) = colnames(result.list$sig.mat) = cn.cov.pool
    row.names(result.list$cor.pvalues) = row.names(result.list$q.values) = row.anmes(result.list$sig.mat) = c(1:num.trios)
    row.names(result.list$reg.pvalues) = trn
    names(filt.qvalues) = trn

  }, correlation = {
    #correlation method - filtering and selection of covs
    result.list = correlation.method(pvalue.matrix = cor.p.mat,
                                     filtering = filter_int_child,
                                     selection.fdr = selection_fdr,
                                     lam.seq = lambda,
                                     trio.idx = trio.indices)
    #naming
    names(result.list$final.list.covs) = names(result.list$filtered.covs) = trn
    colnames(result.list$p.values) = colnames(result.list$q.values) = colnames(result.list$sig.mat) = cn.cov.pool
    row.names(result.list$p.values) = row.names(result.list$q.values) = row.names(result.list$sig.mat) = trn


    #stop message for missing "method"
  }, stop("Method not included or missing"))



  if(save.list==TRUE){
    save(result.list, file = paste0(save.path,".RData"))
  }


  return(result.list)

}







#' A simple function to apply the qvalue correction to a set of pvalues
#'
#' @param p a vector of p.values to be passed to qvalue::qvalue()
#' @param fdr the false discovery rate
#' @param lambda The value of the tuning parameter to estimate \eqn{pi_0}. Must be in \eqn{[0,1)}
#' @return an \eqn{n X 2} dataframe containing the qvalues, and significance (logical)
#' @export adjust.q



adjust.q=function(p, fdr, lambda){
  #apply qvalue correction
  if(is.null(lambda)){
    qval.str=qvalue::qvalue(p, fdr.level = fdr, lambda=seq(0.05, max(p, na.rm = T), 0.05))
  }else{
    qval.str=qvalue::qvalue(p, fdr.level = fdr, lambda=lambda)
  }
  #extrac significance and qvalues and return
  sig=qval.str$significant
  qval=qval.str$qvalues
  return(cbind.data.frame(significant=sig, qvalue=qval))
}



#' A simple function to apply the qvalue correction to a set of pvalues specifically from filtered regression
#'
#' @param x a two element list with the first element being a set of pvalues and the second element being
#' a set of indicies for pvalues to omit
#' @param fdr the false discovery rate
#' @param lambda The value of the tuning parameter to estimate \eqn{pi_0}. Must be in \eqn{[0,1)}
#' @return an \eqn{n X 2} dataframe containing the qvalues, and significance (logical)
#' @export adjust.q.reg



adjust.q.reg=function(x, fdr, lambda, num.of.trios){
  #preallocate significance and qvalue as vectors of NAs
  sig.vec = rep(NA, num.of.trios)
  q.vec = rep(NA, num.of.trios)
  #apply qvalue correction
  if(is.null(lambda)){
    #if cutoff values not supplied set them as:
    lambda=seq(0.05, max(x[[1]], na.rm = T), 0.05)
  }

  if(length(x[[2]]) == 0){
    #this block handles the case when there are no covariates that are significant with the snp
    #and we do not need to omit any of the pvalues
    qval.str=qvalue::qvalue(x[[1]], fdr.level = fdr, lambda = lambda)
    sig.vec = qval.str$significant
    q.vec=qval.str$qvalues
  }else{
    #this block handles the case when there are one of more covariates significant with the variant
    #and we wish to omit the pvalues corresponding them
    qval.str=qvalue::qvalue(x[[1]][-x[[2]]], fdr.level = fdr, lambda = lambda)
    #this sets omitted pvalues to NA
    sig.vec[-x[[2]]] = qval.str$significant
    q.vec[-x[[2]]]=qval.str$qvalues
  }
  #extract significance and qvalues and return

  return(cbind.data.frame(significant=sig.vec, qvalue=q.vec))
}



#' A function to calculate the pearson correlation p-values
#'
#' @param r either (1) a single pearson correlation coefficient or (2) a vector of correlation coefficients
#' @param n the sample size (if \eqn{length(r)>1} then n must be a vector of sample sizes)
#' @return an \eqn{n X 2} dataframe containing the correlations and pvalues
#' @export p.from.cor



p.from.cor=function(r, n){
  #calculate t values from vector of correlations
  t=r*sqrt((n-2)/(1-r^2))
  #obtain the pvalue from two-tailed test
  p=2 * (1 -stats::pt(q = abs(t), df = n-2))
  return(cbind.data.frame(pvalue=p,cor=r))
}

#' A function to perform the regression of a possible confounding pc aganst all pairs of cis and trans genes given
#' in \eqn{genes}
#'
#' @param pc a vector of scores from a given pc (i.e a single column of \eqn{cov.pool} in \eqn{get.conf()})
#' @param genes an \eqn{n X m} dataframe of cis and trans gene pairs of trios
#' @return an \eqn{1 X (m/2)} vector of pvalues from corresponding to the overall \eqn{F} of the regression
#' @export p.from.reg


p.from.reg=function(pc, genes){
  #extract cis and trans gene idx
  seq1=seq(1, dim(genes)[2], 2)
  seq2=seq(2, dim(genes)[2], 2)
  #bind seq into list
  seq.list=as.list(rbind.data.frame(seq1, seq2))
  #convert each pair of genes and the pc into a list of dataframes
  dfs=lapply(seq.list, function(x,y,z){cbind.data.frame(PC=z, gene1=y[,x[1]], gene2=y[,x[2]])},
             z = pc, y = genes)
  #apply the regression of the pc on each pair of genes
  fstats=sapply(dfs, function(x){summary(stats::lm(PC~., data=x))$fstatistic})
  #calculate the pvalue of the overall f statistic from each regression
  pstats=apply(fstats, 2, function(x){1-stats::pf(q = x[1], df1 = x[2], df2 = x[3])})
  #return pvalues vector
  return(pstats)
}



#' this function wraps adjust.q and adjust.q reg
#'
#' It is not recommended to use this function directly
#'
#' @param pvalues a matrix of pvalues or a list of 2-list containing the pvalues and indicies for which pvalues to exlude
#' (i.e from filtering)
#' @param fdr.level the false discovery rate passed to adjust.q or adjust.q.reg
#' @param lambda The value of the tuning parameter to estimate \eqn{pi_0}. Must be in \eqn{[0,1]} passed to adjust.q or adjust.q.reg
#' @param num.of.trios the total number of trios
#' @param input.type a string specifying one of 'array' or 'list'
#' @return of input.type = list, then the output is a list of 2-lists containing the q values and sigificance. Else if the
#' input.type == array, then the output is a 2-list where the first element is a matrix of qvalues and the second element is
#' a matrix of significance determinations
#' @export get.q.sig

get.q.sig = function(pvalues, fdr.level, lambda.seq, num.of.trios, input.type = "array"){

  if(input.type == "array"){
    #if the input is an array of pvalues
    adjust.out = apply(pvalues, 2, adjust.q, fdr = fdr.level, lambda = lambda.seq)

  }else{
    #if the input is a list --- only for filtering with regression method
    adjust.out = lapply(pvalues, adjust.q.reg, fdr = fdr.level, lambda = lambda.seq,
                        num.of.trios = num.of.trios)
  }
  #extract the qvalue and significance matrix
  cor.sig.mat = sapply(adjust.out, function(x) x$significant)
  cor.q.mat = sapply(adjust.out, function(x) x$qvalue)
  return(list(sigmat = cor.sig.mat, qmat = cor.q.mat))
}


#' this function is wrapped by get.conf.trios and executes covariate selection and filtering for the regression method
#'
#'#' It is not recommended to use this function directly
#'
#' @param covpool pool of confounding variables passed from \eqn{get.conf.trios}
#' @param trio.mat the trios in matrix form passed from \eqn{get.conf.trios}
#' @param cor.pvalue.matrix the correlation pvalue matrix passed from passed from \eqn{get.conf.trios}
#' @param filtering (logical) TRUE indicates that filtering should be done, passed from \eqn{get.conf.trios}
#' @param filter.fdr the false discovery rate for filtering common child and intermediate variables
#' @param selection.fdr the false discovery rate for selecting confounding variables
#' @param lam.seq the cutoff values for the q-value correction passed from \eqn{get.conf.trios}
#' @param triolist the trios is list format
#' @return a list
#' @export regression.method

regression.method = function(covpool, trio.mat, cor.pvalue.matrix, filtering, filter.fdr,
                             selection.fdr, lam.seq, triolist){
  #METHOD = Regression
  #pvalues from regressions on each cov
  pc.list=as.list(as.data.frame(covpool))
  #extract only genes
  genes=trio.mat[,-seq(1,dim(trio.mat)[2], 3)]
  message("Applying: method = \"regression\" on confounding variables, this step may take some time...")
  #regress each PC on each pair of genes i.e PC ~ gene1 + gene2
  p.mat=as.data.frame(sapply(pc.list, p.from.reg, genes = genes))

  #---------------filtering------------------
  if(filtering == TRUE){
    #if filtering, then filter the covariate pool for each trio
    message("For method: \"regression\" - filtering common child and intermediate variables from cov.pool")
    #get the indices for the variants
    idx = seq(1, dim(trio.mat)[2]-2, 3)
    #get filtering qvalues and significance matrices, using filtering fdr
    #only supply the values for each covariate with the genetic variants
    out = get.q.sig(pvalues = cor.pvalue.matrix[idx,], fdr.level = filter.fdr, lambda.seq = lam.seq)
    #get the indices of covs significant with the variants
    variant.sig = apply(out$sigmat, 1, which)
    #convert the regression pvalues in a list of length = ncol(p.mat)
    p.as.list = as.list(t(p.mat))

    #combine the vector of pvalues for a pc with all gene pairs with the indicators for which
    #pvalues should be omitted in the adjustment (from variant.sig) into an iterable list
    iterate.list = lapply(as.list(1:dim(covpool)[2]), function(x,y,z) list(y[[x]],z[[x]]),
                         y = p.as.list, z = variant.sig)

    #calculate the qvalues for the regression pvalues ommitted pvalues for covs significant with the snp
    #iterates over the vector of pvalues corresponding to each pc with all trios
    out.final = get.q.sig(pvalues = iterate.list, fdr.level = selection.fdr,
                          lambda.seq = lam.seq, num.of.trios = length(idx), input.type = 'list')


  }else{
    #if no filtering
    out = list(sigmat = NULL, qmat = NULL)
    out.final = get.q.sig(pvalues = p.mat, fdr.level = selection.fdr, lambda.seq = lam.seq)

  }

  #---------------------selection-of-covs---------------------
  #obtain the final list of selected covs
  final.list.sig.asso.covs = apply(out.final$sigmat, 1, function(x){which(x)})
  #obtain the final list of omitted covs
  filtered = apply(out.final$sigmat, 1, function(x){which(is.na(x))})

  return(list(final.list.covs = final.list.sig.asso.covs, filtered.covs = filtered,
              reg.pvalues = p.mat, cor.pvalues = cor.pvalue.matrix,
              q.values = out$qmat, sig.mat = out$sigmat,
              filt.qvalues = out$qmat))

}





#' this function is wrapped by get.conf.trios and executes covariate selection and filtering for the correlation method
#'
#' It is not recommended to use this function directly
#'
#' @param pvalue.matrix the correlation pvalue matrix passed from passed from \eqn{get.conf.trios}
#' @param filtering (logical) TRUE indicates that filtering should be done, passed from \eqn{get.conf.trios}
#' @param selection.fdr the false discovery rate for selecting confounding variables
#' @param lam.seq the cutoff values for the q-value correction passed from \eqn{get.conf.trios}
#' @param trio.idx an array of size \eqn{number of trios X 2} giving the column start and stop indices for each trio in the
#' trio matrix passed from \eqn{get.conf.trios}
#' @return a list
#' @export correlation.method

correlation.method = function(pvalue.matrix, filtering, selection.fdr, lam.seq, trio.idx){

  #METHOD = Correlation
  #get the qvalues and significance matrix at selection_fdr
  out = get.q.sig(pvalues = pvalue.matrix, fdr.level = selection.fdr, lambda.seq = lam.seq)
  #get the significant covariates for every column of triomat
  sig.asso.covs=apply(out$sigmat, 1, function(x){which(x)})

  #-----------------filtering-of-results------------------
  #condense to a single list of covariates for each trio
  #if filtering is specified then we omit the covs detected for the variant:
  if(filtering == TRUE){
    message("For method: \"correlation\" - filtering common child and intermediate variables from results")
    #get the covariates significant with the snp
    sig.asso.covs.trios=apply(trio.idx, 1,
                             function(x,y){ unique(unlist(y[x[1]:x[2]])) },
                             y=sig.asso.covs)

    sig.covs.snp = out$sigmat[seq(1,dim(out$sigmat)[1],3),]

    #create a list of indexes to iterate over for filtering
    iterate.idx = as.list(c(1:length(sig.asso.covs.trios)))

    which.match = lapply(iterate.idx,
                         function(x,y,z) na.omit(match(which(z[x,]), y[[x]])), y = sig.asso.covs.trios, z = sig.covs.snp )

    filtered = lapply(as.list(c(1:length(sig.asso.covs.trios))),
                      function(x,y,z) y[[x]][z[[x]]], y = sig.asso.covs.trios, z = which.match )
    #remove marginally sig covs
    final.list.sig.asso.covs = lapply(as.list(c(1:length(sig.asso.covs.trios))),
                                     function(x,y,z) y[[x]][-na.omit(z[[x]])], y = sig.asso.covs.trios, z = which.match )

    return(list(final.list.covs = final.list.sig.asso.covs, filtered.covs = filtered))
  }else{

    filtered = NULL
    #final list of pcs for trios
    final.list.sig.asso.covs=apply(trio.idx, 1,
                                  function(x,y){ unique(unlist(y[x[1]:x[2]])) },
                                  y=sig.asso.covs)

    return(list(final.list.covs = final.list.sig.asso.covs, filtered.covs = filtered,
                p.values = pvalue.matrix, q.values = out$qmat, sig.mat = out$sigmat))

  }
}


# #' This function is wrapped by get.conf() when filter_int_child = TRUE and method = "correlation". It identifies (from the selected)
# #' confounders in sig.asso.pcs which are correlated with the genetic variant and removes them from the final list of confounders.
# #' @param sig.asso.pcs a matrix of correlations between the genetic variants the covariate pool of potential confounders
# #' @param trios see get.conf() for details
# #' @param cov.pool see get.conf() for details
# #' @param filter.fdr the filtering fdr threshold
# #' @return a list of the significant associated PCs/confounders with the indices of the identified common child and intermediate
# #' confounding variables removed.
# #' @export filter_out

# filter_out = function(sig.asso.pcs, trios, cov.pool, filter.fdr){
#
#   if(is.data.frame(trios)==TRUE){
#     trio.indices=cbind( start.col = seq(1,dim(trios)[2],3),
#                         end.col = seq(3,dim(trios)[2],3))
#
#     trios = apply(trio.indices, 1, function(x,y) y[,x[1]:x[2]], y = trios)
#   }
#   #calculate the mean number of confounders selected for each trio by get.conf()
#   #if this value is less than 50 it is unlikely there will be enough pvalue to estimate pi0_est
#   #in qvalue...
#   mean.num.conf = mean(unlist(lapply(sig.asso.pcs, function(x) length(x) )))
#   #get the sample sizes for cor calc with each genetic variant (note: assumes the variant is the first column)
#   sample.sizes = unlist(lapply(trios, function(x) length(stats::na.omit(x[,1])) ))
#   # get the correlations between the selected confounders and the SNP
#   cors.list = lapply(c(1:length(trios)), function(x,y,z,w) stats::cor(y[[x]][,1], z[,w[[x]]], use = "pairwise.complete.obs"),
#                 y = trios, z = cov.pool, w = sig.asso.pcs)
#   #get the pvalues for each set of correlations
#   p.list = lapply(c(1:length(trios)), function(x,y,z) p.from.cor(r = t(unlist(y[[x]])), n = z[x])$pvalue,
#                   y = cors.list, z = sample.sizes)
#
#   #perform FDR filtering (either qvalue or BH)
#   if(mean.num.conf>=50){
#     message("in filter_out: the mean number of confounders selected for each trio is >50, applying q-value correction for filtering")
#     q.correct = lapply(p.list, qvalue::qvalue, fdr.level = filter.fdr)
#     sig.list = lapply(q.correct)$significant
#   }else{
#     message("in filter_out: the mean number of confounders selected for each trio <50: applying \" BH \" correction for filtering")
#     BH.correct = lapply(p.list, function(x) stats::p.adjust(x, method = "BH") )
#     sig.list = lapply(BH.correct, function(x) x < filter.fdr )
#   }
#
#   #remove the confoundering variables with significant correlation with genetic variant after FDR adjust
#   sig.asso.pcs.filtered = lapply(c(1:length(trios)), function(x,y,z) y[[x]][-which(z[[x]])],
#                                  y = sig.asso.pcs, z = sig.list)
#   #return filtered list
#   return(sig.asso.pcs.filtered)
#
# }
