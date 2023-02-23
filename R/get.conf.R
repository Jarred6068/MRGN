
#' A function to select a set of confounding variables for an input dataset
#'
#' This function takes in a set of variables represented in the columns of \eqn{data} and a pool of potential confounding variables (e.g whole-genome expression or
#' methylation PC scores) and selects the confounding variables associated with each column \eqn{data} by testing the marginal associations of each column of \eqn{data} with
#' the covariates in \eqn{cov.pool}. This method is a simpler version of \eqn{get.conf.trios} which does not require the input data be
#' organized into trios.
#'
#' @param data a matrix or dataframe of size \eqn{n samples X p variables}.
#' @param cov.pool a matrix or dataframe of size \eqn{n samples X q covariates}whose columns represent the possible confounding variables to be selected from
#' @param measure a string specifying one of 'correlation' or 'partial_corr'. If 'correlation', then the marginal dependence between each column in
#' data and cov.pool is tested by pearson correlation. If 'partial_corr', the conditional dependence of each column of data and each column of cov.pool is
#' tested after conditioning on all variables in 'conditional.vars' using pearson partial correlation test. Note that 'partial_corr' is not possible if the sample size
#' is less than ncol(conditional.vars)+2 .
#' @param conditional.vars a matrix or dataframe. Used when measure = 'partial_corr'. The variables for which the user wishes to use as the conditional variable(s)
#' for each pairwise partial correlation estimate.
#' @param blocksize the number of columns of \eqn{data} to use in each block of correlation calculations passed to propagate::bigcor
#' @param apply.qval (logical) \eqn{default = TRUE} applies the qvalue adjustment to each set of correlations between a PC
#' and the columns of \eqn{trios}. If \eqn{FALSE}, the Bonferroni correction is applied.
#' @param selection_fdr the false discovery rate (default = 0.05) for selecting confounders when apply.qval = TRUE
#' @param adjust_by a string specifying one of 'fwer' or 'all". If 'fwer' then the qvalue adjustment is done for the family of tests
#' corresponding to all covariates with each column of 'data'. If 'all' then the qvalue adjustment is applied to all marginal pvalues.
#' default = 'fwer'
#' @param lambda When apply.qval = TRUE, lambda is the set of cut off points of the tuning parameter to estimate \eqn{pi_0}. Must be between 0,1. If is.null(lambda) the default
#' passed to \eqn{adjust.q()} is \eqn{seq(0.5, max(pvalues), 0.05)}
#' @param pi0.method a string specifying one of 'smoother' or 'boostrap' describing the method used to estimate
#' Pi0. Passed to qvalue::qvalue(). default = 'smoother'. Note: bootstrap' can be used in place of 'smoother' if 'smoother' fails
#' @param alpha the test threshold for the bonferroni correacted pvalues when \eqn{apply.qval = FALSE}
#' @param save.list (logical) if TRUE the output is saved as a .RData object (default = FALSE)
#' @param save.path string specifying the path name of the output list to be save as a .RData structure
#' @return a list of 6 elements containing:
#'   \describe{
#'   \item{sig.asso.covs}{a list of \eqn{length=ncol(data)}. Each element in the list is a vector of of indices corresponding to the columns
#'   of \eqn{cov.pool} that are significantly associated to each input variable (column) in \eqn{data}}
#'   \item{pvalues}{a matrix of dimension \eqn{ncol(cov.pool) X ncol(data)} of the pvalues from the pearson correlation test on each set of correlations}
#'   \item{qvalues}{a matrix of dimension \eqn{ncol(cov.pool) X ncol(data)} of qvalues from each set of correlations}
#'   \item{cors}{the matrix of the calculated pairwise correlations (or partial correlations) of dimension \eqn{ncol(cov.pool) X ncol(trios)}}
#'   \item{sig}{A matrix of logical values of dimension \eqn{ncol(cov.pool) X ncol(data)}}
#'   \item{adj.p}{A matrix of dimension \eqn{ncol(cov.pool) X ncol(data)} of the adjusted p-values if \eqn{apply.qval = FALSE}}
#' }
#' @export get.conf.matrix
#' @references
#'     \insertAllCited{}
#' @import propagate
#' @import qvalue
#' @examples
#'
#' \dontrun{
#' #fast example using 'correlation' to select confounders for Whole Blood genes
#' #set the blocksize to be 1/3 the number of genes
#' result=get.conf.matrix(data=WBgenes,
#'                        cov.pool=WBscores,
#'                        measure = 'correlation'
#'                        blocksize=round((1/3)*dim(WBgenes)[2]))
#'
#'
#'
#'
#' #fast example using 'partial_corr' to select PCs for Whole Blood genes conditioning on Whole Blood SNPs
#' result=get.conf.matrix(data=WBgenes,
#'                        cov.pool=WBscores,
#'                        conditional.vars = WBsnps,
#'                        measure = 'partial_corr')
#'
#'}


get.conf.matrix=function(data=NULL, cov.pool=NULL, measure = c('correlation','partial_corr'), conditional.vars = NULL,
                         blocksize=2000, apply.qval=TRUE, selection_fdr=0.05, adjust_by = 'fwer', lambda=NULL,
                         pi0.method = 'smoother', alpha=0.05, save.list=FALSE, save.path="/path/to/save/location"){

  #====================================Preprocessing====================================
  data = as.data.frame(data)
  n1 = dim(data)[1]
  p = dim(data)[2]
  cond.p = dim(conditional.vars)[2]
  cn.data = colnames(data)
  cov.pool = as.data.frame(cov.pool)
  n2 = dim(cov.pool)[1]
  q = dim(cov.pool)[2]
  cn.cov.pool = colnames(cov.pool)
  #check to make sure both cov.pool and data are the same rowsize
  if(n1 != n2){
    stop(paste0('cov.pool and data are not the same size! data has ',n1,' samples and cov.pool has ',n2,' samples!'))
  }

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
  switch(measure, correlation = {
    message('selected measure = correlation...')
    message(paste0("Calculating correlation matrix of size ", ncol(data), " x ", ncol(cov.pool),
                   " using ", ceiling(dim(data)[2]/blocksize), " blocks"))
    #correlation matrix calc
    cormat = propagate::bigcor(data, cov.pool, verbose = T, use = "pairwise.complete.obs", size = blocksize)

    #perform the pearson correlation test, apply qvalue, and return which are significant
    pr.mat=apply(as.data.frame(cormat[,]), 2, p.from.cor, n=sample.sizes)
    #extract the pvalues and correlations
    r.mat = t(sapply(pr.mat, function(x) x$cor))
    #matrix of correlation pvalues
    p.mat = t(sapply(pr.mat, function(x) x$pvalue))

    #for partial correlations:
    }, partial_corr = {
      message('selected measure = partial correlation')
      #catch and stop if no conditional variables are specified
      if(is.null(conditional.vars)){
        stop('conditional.vars is empty!')
      }


      #catch and stop if sample size is too small to compute pcor test
      if(n1<(cond.p)){
        #check to make sure matrix is not overdetermined:
        stop(paste0('Cannot compute partial correlation test!: samples size = ',n1,'< ',((cond.p)+3),
                    ' Consider using measure = \'correlation\' instead '))
      }


      #catch and warn about NA values in the input objects
      contains.nas = c(any(is.na(data)), any(is.na(cov.pool)), any(is.na(conditional.vars)))
      if(any(contains.nas)){
        object.nm = c('data', 'cov.pool', 'conditional.vars')
        which.contain.na = object.nm[which(contains.nas)]
        #if more than one of the objects contains NAs stop
        if(length(which.contain.na)>1){
          stop(paste0(paste(which.contain.na, collapse = ' and '), ' contains NAs!...stopping'))
        }
        #if only one contains NA's just omit the NA's and proceed
        message(paste0(which.contain.na, ' contains NAs...omitting rows with missing values...'))
        #determine which object has the NA's and ID which rows to omit
        if(which(contains.nas) == 1){
          omit = which(apply(data, 2, function(x) any(is.na(x))))
        }else if(which(contains.nas) == 2){
          omit = which(apply(cov.pool, 2, function(x) any(is.na(x))))
        }else{
          omit = which(apply(conditional.vars, 2, function(x) any(is.na(x))))
        }

        #recheck the sample size condition after omitting rows
        if((n1 - length(omit)) < (cond.p+3)){
          stop('After omitting rows with missing values, sample size is too small to compute partial correlation test!...stopping')
        }
        data.omit = data[-omit,]
        cov.pool.omit = cov.pool[-omit,]
        conditional.vars.omit = conditional.vars[-omit,]
        message('Computing pairwise partial correlations...this step may take a few seconds...')
        ppout = compute.pairwise.pcors(data = data, confs = cov.pool, cond.vars = conditional.vars)
      }else{
        message('Computing pairwise partial correlations...this step may take a few seconds...')
        ppout = compute.pairwise.pcors(data = data, confs = cov.pool, cond.vars = conditional.vars)
      }
      #perform the pearson correlation test, apply qvalue, and return which are significant
      pr.mat=apply(t(ppout), 2, p.from.parcor, n=sample.sizes, S = cond.p)
      #extract the pvalues and correlations
      r.mat = t(sapply(pr.mat, function(x) x$parcor))
      #matrix of correlation pvalues
      p.mat = t(sapply(pr.mat, function(x) x$pvalue))
    }, stop('Argument \'measure\' not specified'))


  #===================================Confounder-Selection=======================================

  #-----------------apply-correction-to-pvalues-------------------
  #perform the correction
  if(apply.qval==TRUE){
    message(paste0("Applying qvalue correction to control the FDR at ", selection_fdr))
    #qvalue correction
    if(adjust_by == 'fwer'){
      qsig.mat = apply(p.mat, 2, adjust.q, fdr = selection_fdr, lambda = lambda, pi0.meth = pi0.method)
      #significance matrix (binary matrix)
      sig.mat = sapply(qsig.mat, function(x) x$significant)
      #qvalue matrix
      q.mat = sapply(qsig.mat, function(x) x$qvalue)
    }else if(adjust_by == 'all'){
      qsig.mat = adjust.q(as.vector(as.matrix(p.mat)), fdr = selection_fdr, lambda = lambda, pi0.meth = pi0.method)
      #restructure output into significance matrix of snps X covariates
      sig.mat = as.data.frame(matrix(qsig.mat$significant, nrow = nrow(p.mat),
                                    ncol = ncol(p.mat), byrow = F))
      #restructure output into qvalue matrix of snps X covariates
      q.mat = as.data.frame(matrix(qsig.mat$qval, nrow = nrow(p.mat),
                                  ncol = ncol(p.mat), byrow = F))

    }else{
      stop(paste0('Input to argument \'adjust_by \' = ', adjust_by, ' is not recognized. use one of \'fwer\' or \' all \' '))
    }
    #empty matrix
    p.adj.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
  }else{
    message(paste0("Applying bonferroni correct with threshold ", alpha))
    #bonferroni correction
    #adjustment by columns or by all pvalues
    if(adjust_by == 'fwer'){
      #adjust by columns
      #matrix of adjusted pvalues
      p.adj.mat = apply(p.mat, 2, stats::p.adjust, method="bonferroni")
      #significance matrix (binary matrix)
      sig.mat = apply(p.adj.mat, 2, function(x) x<=alpha)
      #empty matrix
      q.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))

    }else if(adjust_by == 'all'){
      #adjust by all pvalues
      #matrix of adjusted pvalues
      p.adj.mat = as.data.frame(matrix(stats::p.adjust(as.vector(as.matrix(p.mat)), method="bonferroni"),
                                          nrow = nrow(p.mat), ncol = ncol(p.mat), byrow = F))
      #significance matrix (binary matrix)
      sig.mat = apply(p.adj.mat, 2, function(x) x<=alpha)
      #empty matrix
      q.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))

    }else{
      stop(paste0('Input to argument \'adjust_by \' = ', adjust_by, ' is not recognized. use one of \'fwer\' or \' all \' '))
    }
  }



  #-------------------obtain-final-list-of-significant-covs----------------
  #find the covs that correlated with every column of the trio matrix
  message("Selecting significant covariates...")
  #extract significant covariates:
  sig.asso.covs=apply(sig.mat, 2, function(x){which(x)})


  #=================================Organize and return/save output====================================
  colnames(sig.mat)=colnames(q.mat)=colnames(r.mat)=colnames(p.mat)=colnames(p.adj.mat)=cn.data
  row.names(sig.mat)=row.names(q.mat)=row.names(p.mat)=row.names(p.adj.mat)=cn.cov.pool



  out.list=list(sig.asso.covs = sig.asso.covs,
                pvalues = p.mat,
                qvalues = q.mat,
                cors = r.mat,
                sig = sig.mat,
                adj.p = p.adj.mat)


    if(save.list==TRUE){
      message(paste0("saving output at ", save.list, ".RData"))
      save(out.list, file = paste0(save.path,".RData"))
    }

    return(out.list)

}






#' A function to select a set of confounding variables for each trio
#'
#' This function takes in a list or dataframe of trios and a matrix of potential confounding variables (e.g whole-genome expression or
#' methylation PC scores) and selects the confounding variables associated with each trio using the regression-based procedure of
#' \insertCite{yang2017identifying}{MRGN}. Additionally, this function allows the user to filter common child and intermediate confounding
#' variables by removing selected confounders that are significantly associated with the genetic variant.
#'
#' @param trios either (1) a list where each element contains a dataframe corresponding to a trio or (2) a single data
#' frame with samples in rows and trios in the columns. The assumed structure of each trio is that the Genetic variant
#' is the first column followed by the two molecular phenotypes.
#' @param cov.pool a dataframe of candidate confounding variables with observations in rows and covariates in columns
#' @param blocksize the number of variants to use in each block of correlation calculations passed to propagate::bigcor
#' @param selection_fdr the false discovery rate (default = 0.05) for selecting confounders
#' @param filter_int_child (logical) to determine if common child and intermediate confounding variables should be filtered out
#' (removed) from the significant confounders for each trio. Note: only used when return.for.trios==TRUE. Default = FALSE.
#' @param filter_fdr the false discovery rate for filtering common child
#' and intermediate confounding variable.
#' @param adjust_by a string specifying one of 'fwer' or 'all". If 'fwer' then the qvalue adjustment is done for the family of tests
#' corresponding to each pc with all genetic variants. If 'all' then the qvalue adjustment is applied to all marginal pvalues
#' @param lambda the cut off points of the tuning parameter to estimate \eqn{pi_0}. Must be between 0,1. If is.null(lambda) the default
#' passed to \eqn{adjust.q()} is \eqn{seq(0.5, max(pvalues), 0.05)}
#' @param pi0.method a string specifying one of 'smoother' or 'boostrap' describing the method used to estimate
#' Pi0. Passed to qvalue::qvalue(). default = 'smoother'. Note: bootstrap' can be used in place of 'smoother' if 'smoother' fails
#' @param save.list (logical) if TRUE the output is saved as a .RData object (default = FALSE)
#' @param save.path string specifying the path name of the output
#' @return a list of 9 elements containing:
#'   \describe{
#'   \item{sig.asso.covs}{a list of length = number of trios. each element is a set of column indices indicating selected covariates
#'   from \eqn{cov.pool}}
#'   \item{filtered.covs}{if \eqn{filter_in_child} = TRUE, then it is a list of length = number of trios giving the column indices of
#'   \eqn{cov.pool} indicating the covariates that were filtered}
#'   \item{reg.pvalues}{a matrix of dimension \eqn{number of trios X ncol(cov.pool)} of pvalues}
#'   \item{reg.qvalues}{a matrix of dimension \eqn{number of trios X ncol(cov.pool)} of qvalues}
#'   \item{reg.sig.mat}{a logical matrix of dimension \eqn{number of trios X ncol(cov.pool)} indicating the significant qvalues from \eqn{reg.q.values}}
#'   \item{filtering.pvalues}{a matrix of dimension \eqn{number of trios X ncol(cov.pool)} of filtering pvalues}
#'   \item{filtering.correlations}{a correlation matrix of dimension \eqn{number of trios X ncol(cov.pool)} representing the correlations between the genetic variants and covariates}
#'   \item{filtering.qvalues}{a matrix of dimension \eqn{number of trios X ncol(cov.pool)} of filtering qvalues}
#'   \item{filtering.sigmat}{a logical matrix of dimension \eqn{number of trios X ncol(cov.pool)} indicating the significant qvalues from \eqn{filtering.qvalues}}
#' }
#' @export get.conf.trios
#' @references
#'     \insertAllCited{}
#' @import propagate
#' @import qvalue
#' @examples
#'
#' \dontrun{
#' #fast example on 40 trios
#' trio.conf=get.conf.trios(trios=WBtrios[1:40],
#'                          cov.pool=WBscores,
#'                          blocksize=20)
#'
#' #fast example on 40 trios using common child and intermediate variable filtering
#' trio.conf2=get.conf.trios(trios=WBtrios[1:40],
#'                          cov.pool=WBscores,
#'                          blocksize=20,
#'                          filter_int_child = TRUE)
#'}


get.conf.trios=function(trios=NULL, cov.pool=NULL, blocksize=2000, selection_fdr=0.05, filter_int_child = FALSE,
                        filter_fdr = 0.1, adjust_by = 'fwer', lambda=NULL, pi0.method = 'smoother',
                        save.list=FALSE, save.path="/path/to/save/location"){

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
    num.trios = dim(trios)[2]/3
    #get index of trios and convert to list
    trio.cols = split(cbind.data.frame(variant.idx = seq(1, dim(triomat)[2]-2, 3),
                     gene1.idx = seq(2, dim(triomat)[2]-1, 3),
                     gene2.idx = seq(3, dim(triomat)[2], 3)), seq(1, num.trios, 1))

    #obtain list of trios where each entry is a n X 3 dataframe of a trio
    trio.list = lapply(trio.cols, function(x,y) y[,unlist(x)], y = triomat)
  }

  reg.qmat = reg.sigmat = matrix(rep(NA, num.trios*dim(cov.pool)[2]), nrow = num.trios, ncol = dim(cov.pool)[2])

  #make colnames unique to avoid error in bigcor
  colnames(triomat) = make.unique(colnames(triomat))
  #get the index for all genetic variants
  snp.idx = seq(1, dim(triomat)[2]-2, 3)
  #extract variants into matrix
  snp.mat = triomat[,snp.idx]
  snp.cn = colnames(triomat)[snp.idx]
  #extract genes into matrix
  gene.mat = triomat[, - snp.idx]
  gene.cn = colnames(triomat)[-snp.idx]

  #preallocate names for trios
  trn = paste0("trio_", c(1:num.trios))

  print(paste0("detected ", num.trios, " trios and ", dim(cov.pool)[2], " covariates in the candidate pool..."))



  #catch columns with only 2 or less non-NA values and return their indices
  non.na.vals = apply(triomat, 2, function(x) length(stats::na.omit(x)))
  if(any(non.na.vals<=2)){
    stop(paste0("some columns of \"trios\" contain <= 2 non-NA values: The column(s) is/are ", paste0(which(non.na.vals<=2))))
  }



  #METHOD = Regression

  #---------------filtering------------------
  if(filter_int_child == TRUE){

    #====================================Calculate/Test-Correlations====================================
    #this step is get the correlation matrix of size (3 X number of trios) x Number of covs
    #we then get the marginal test pvalues by testing the marginal relationship between each
    #variant or gene with each cov to produce a matrix of size (3 X number of trios) x Number of covs
    #containing the pvalues

    #get the sample sizes used in each pairwise correlation calculation
    sample.sizes = apply(snp.mat, 2, function(x) length(na.omit(x)))

    print(paste0("Calculating correlation matrix for variants.. size ", num.trios, " x ", dim(cov.pool)[2],
                 " using ", ceiling(num.trios/blocksize), " blocks"))
    #calculate the correlations of each PC with the trios
    cormat = propagate::bigcor(snp.mat, cov.pool, verbose = T, use = "pairwise.complete.obs", size = blocksize)
    #perform marginal test
    pr.mat=apply(as.data.frame(cormat[,]), 2, p.from.cor, n=sample.sizes)
    #correlation matrix and pvalue matrix
    cor.r.mat= sapply(pr.mat, function(x) x$cor)
    cor.p.mat = sapply(pr.mat, function(x) x$pvalue)


    #=============================================filtering===============================================
    #this block of code handles both the filtering and covariate selection steps
    #the steps to the final result are different depending on method
    #the switch separates the process according to each method.
    #if filtering, then filter the covariate pool for each trio
    message("Filtering common child and intermediate variables from cov.pool...")

    #get filtering qvalues and significance matrices, using filtering fdr
    #only supply the values for each covariate with the genetic variants
    if(adjust_by == 'fwer'){
      out = get.q.sig(pvalues = cor.p.mat, fdr.level = filter_fdr, lambda.seq = lambda,
                      pi0.method = pi0.method)
      #get the indices of covs significant with the variants
      variant.sig = apply(out$sigmat, 1, which)
      #store filtered covs and filtering information
      filtered = variant.sig
      filt.sig = out$sigmat
      filt.q = out$qmat

    }else if(adjust_by == 'all'){
      #apply the qvalue correction to all the pvalues
      out = adjust.q(p = c(t(cor.p.mat)), fdr = filter_fdr, lambda = lambda,
                     pi0.meth = pi0.method)
      #restructure output into significance matrix of snps X covariates
      sigmat = as.data.frame(matrix(out$significant, nrow = dim(cor.p.mat)[1],
                                    ncol = dim(cor.p.mat)[2], byrow = T))
      #restructure output into qvalue matrix of snps X covariates
      qmat = as.data.frame(matrix(out$qval, nrow = dim(cor.p.mat)[1],
                                  ncol = dim(cor.p.mat)[2], byrow = T))

      #get the significant covs
      variant.sig = apply(sigmat, 1, which)
      if(length(variant.sig)==0){
        stop("No covariates detected, Filtering not needed...stopping")
      }

      #store filtered covs and filtering info
      filtered = variant.sig
      filt.sig = sigmat
      filt.q = qmat
    }


    #get the indices of each gene pair
    gene.pair.idx = as.list(rbind.data.frame(seq1 = seq(1, dim(gene.mat)[2]-1, 2),
                                             seq2 = seq(2, dim(gene.mat)[2], 2)))

    #create an iteratable list to be passed to the regression pvalue calculations
    #the list combines the covariates with each gene pair and the indices of which
    #covariates should be omitted due to filtering.
    iterate.list = lapply(as.list(1:length(variant.sig)), function(x,y,z,w,u)
      if(length(z[[x]])==0) list(z[[x]], y, w[,u[[x]]]) else list(z[[x]], y[,-z[[x]]],w[,u[[x]]]),
      y = cov.pool, z = variant.sig, w = gene.mat, u = gene.pair.idx)

    #calculate the regression pvalues - NA's in the pvalue matrix represent omitted covariates for a given
    #gene pair
    message("Calculating regression p-values, this step may take some time...")
    reg.pvalues=t(as.data.frame(sapply(iterate.list, function(x) p.from.reg(inlist = x,
                                                                            cov.cols = dim(cov.pool)[2]))))

    #=========================================Covariate-selection========================================
    #calculate the qvalues for the regression pvalues for covs significant with each gene pair
    #iterates over the vector of pvalues corresponding to each pc with all trios.
    message("Selecting covariates from the filtered pool...")
    out.final = get.q.sig(reg.pvalues, fdr.level = selection_fdr, lambda.seq = lambda, pi0.method = pi0.method,
                          contains.na = any(is.na(reg.pvalues)))

    #extract the qvalues and significance matrix
    reg.sigmat = out.final$sigmat
    reg.qvalues = out.final$qmat







  }else{
    message("Calculating regression p-values, this step may take some time...")
    #if no filtering convert columns of cov pool to elements in list
    pc.list = as.list(cov.pool)
    #get the regression pvalues. simpler here as there will be no omitted/NA pvalues
    reg.pvalues=as.data.frame(sapply(pc.list, p.from.reg2, genes = gene.mat))
    message("Selecting covariates from the pool...")
    out.final = get.q.sig(pvalues = reg.pvalues, fdr.level = selection_fdr, lambda.seq = lambda,
                          pi0.method = pi0.method)
    #extract the qvalue and significance matrix
    reg.sigmat = out.final$sigmat
    reg.qvalues = out.final$qmat
    #set filtering information == NULL
    filtered = NULL
    cor.p.mat = NULL
    cor.r.mat = NULL
    filt.q = NULL
    filt.sig = NULL



  }

  #===================================-Organizing-output-=============================================
  #obtain the final list of selected covs
  final.list.sig.asso.covs = apply(reg.sigmat, 1, function(x){which(x)})

  #naming output
  trio.names = apply(cbind(seq1 = seq(1, dim(gene.mat)[2]-1, 2),
                                seq2 = seq(2, dim(gene.mat)[2], 2),
                                seq3 = seq(1, dim(snp.mat)[2], 1)),
                          1, function(x,y,z) paste(c(z[x[3]], y[x[1:2]]), collapse = ":"),
                          y = gene.cn, z = snp.cn)

  colnames(reg.sigmat) = colnames(reg.pvalues) = colnames(reg.qvalues) = colnames(cov.pool)
  row.names(reg.sigmat) = row.names(reg.pvalues) = row.names(reg.qvalues) = trn
  if(!is.null(filt.q)){
    row.names(filt.q) = row.names(filt.sig) = row.names(cor.p.mat) = row.names(cor.r.mat) = snp.cn
    colnames(filt.q) = colnames(filt.sig) = colnames(cor.p.mat) = colnames(cor.r.mat) = colnames(cov.pool)
  }

  #organize list of info
  result.list = list(sig.asso.covs = final.list.sig.asso.covs, filtered.covs = filtered,
              reg.pvalues = reg.pvalues, reg.qvalues = reg.qvalues, reg.sig.mat = reg.sigmat,
              filtering.pvalues = cor.p.mat, filtering.correlations = cor.r.mat,
              filtering.qvalues = filt.q, filtering.sigmat = filt.sig, trio.ID = trio.names)

  if(save.list==TRUE){
    #if saving true save list at path
    message(paste0("saving output at ", save.list, ".RData"))
    save(result.list, file = paste0(save.path,".RData"))
  }

  return(result.list)

}
















#' A simple function to apply the qvalue correction to a set of pvalues
#'
#' @param p a vector of p.values to be passed to qvalue::qvalue()
#' @param fdr the false discovery rate
#' @param lambda The value of the tuning parameter to estimate \eqn{pi_0}. Must be in \eqn{[0,1)}
#' @param pi0.meth a string specifying 'smoother' or 'bootstrap'. Passed to qvalue::qvalue()
#' @return an \eqn{n X 2} dataframe containing the qvalues, and significance (logical)
#' @export adjust.q



adjust.q=function(p, fdr, lambda, pi0.meth){
  #apply qvalue correction
  if(is.null(lambda)){
    #set the lambda seq
    lambda=seq(0.05, max(p, na.rm = T), 0.05)

    qval.str=qvalue::qvalue(p, fdr.level = fdr, lambda=lambda, pi0.method = pi0.meth)
  }else{
    qval.str=qvalue::qvalue(p, fdr.level = fdr, lambda=lambda, pi0.method = pi0.meth)
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
#' @param pi0.meth A string specifying one of 'smoother' or 'bootstrap'. Passed to qvalue::qvalue()
#' @return an \eqn{n X 2} dataframe containing the qvalues, and significance (logical)
#' @export adjust.q.reg

adjust.q.reg=function(x, fdr, lambda, pi0.meth, num.of.trios){
  #preallocate significance and qvalue as vectors of NAs
  which.na = which(is.na(x))
  sig.vec = rep(NA, num.of.trios)
  q.vec = rep(NA, num.of.trios)
  #apply qvalue correction
  if(is.null(lambda)){
    #if cutoff values not supplied set them as:
    lambda=seq(0.05, max(x, na.rm = T), 0.05)
  }

  if(length(which.na) == 0){
    #this block handles the case when there are no covariates that are significant with the snp
    #and we do not need to omit any of the pvalues
    qval.str=qvalue::qvalue(x, fdr.level = fdr, lambda = lambda, pi0.method = pi0.meth)
    sig.vec = qval.str$significant
    q.vec=qval.str$qvalues
  }else{
    #this block handles the case when there are one of more covariates significant with the variant
    #and we wish to omit the pvalues corresponding them
    qval.str=qvalue::qvalue(na.omit(x), fdr.level = fdr, lambda = lambda, pi0.method = pi0.meth)
    #this sets omitted pvalues to NA
    sig.vec[-which.na] = qval.str$significant
    q.vec[-which.na]=qval.str$qvalues
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


#' A function to calculate the pairwise partial correlations between variables in data and confounders in cov.pool
#' conditioned on the variables in cond.vars
#'
#' @param data a matrix or dataframe of variables for which to select confounders
#' @param confs a matrix or dataframe of covariates representing potential confounders
#' @param conf.vars a matrix or dataframe of variables on which to condition each pairwise relationship
#' @return a matrix of size ncol(confs) x ncol(data) of partial correlations
#' @export compute.pairwise.pcors

compute.pairwise.pcors = function(data, confs, cond.vars){
  #get dimensions
  p = ncol(data)
  u = ncol(confs)
  z = ncol(cond.vars)
  #create an iterable object for parallelization
  iterable = expand.grid(Uvars = c(1:u), Pvars = c(1:p))
  #compute all pairwise pcors
  pcor.mat = apply(iterable, 1, function(x,p,u,z) ppcor::pcor(cbind(u[,x[1]], p[,x[2]], z))$estimate[1,2],
                   p = data, u = confs, z = cond.vars)
  #reshape output into an array
  pcor.final = as.data.frame(matrix(pcor.mat, nrow = p, ncol = u, byrow = T))
  return(t(pcor.final))

}



#' A function to calculate the partial correlation test pvalues using fishers z-transformation
#'
#' @param r either (1) a single partial correlation coefficient or (2) a vector of partial correlation coefficients
#' @param n the sample size (if \eqn{length(r)>1} then n must be a vector of sample sizes)
#' @param S the number of variables in the conditioning set
#' @return an \eqn{n X 2} dataframe containing the partial correlations and pvalues
#' @export p.from.parcor

p.from.parcor=function(r, n, S){
  #calculate z values from vector of partial correlations
  z = (sqrt(n - S - 3)/2)*log((1 + r)/(1 - r))
  #obtain the pvalue from two-tailed test
  p = 2 * (1-pnorm(abs(z), lower.tail = T))
  return(cbind.data.frame(pvalue = p, parcor = r))
}





#' A function to perform the regression of a possible confounding pc aganst all pairs of cis and trans genes given
#' in \eqn{genes}
#'
#' @param inlist a list of 3 elements
#' @param cov.cols the column dimension of the covariate pool
#' @return a vector
#' @export p.from.reg

p.from.reg=function(inlist, cov.cols){

  which.removed = inlist[[1]]
  pc.list = as.list(inlist[[2]])
  genes = inlist[[3]]
  pstats.final = rep(NA, cov.cols)
  #convert each pair of genes and the pc into a list of dataframes
  dfs=lapply(pc.list, function(x,y) {cbind.data.frame(PC=x, gene1 = y[,1], gene2 = y[,2])}, y = genes)
  #apply the regression of the pc on each pair of genes
  fstats=sapply(dfs, function(x){summary(stats::lm(PC~., data=x))$fstatistic})
  #calculate the pvalue of the overall f statistic from each regression
  pstats = apply(fstats, 2, function(x){stats::pf(q = x[1], df1 = x[2], df2 = x[3], lower.tail = F)})
  #return pvalues vector
  if(length(which.removed) == 0){
    pstats.final = pstats
  }else{
    pstats.final[-which.removed] = pstats
  }
  return(pvalue = pstats.final)
}










#' A function to perform the regression of a possible confounding pc aganst all pairs of cis and trans genes given
#' in \eqn{genes} - when no filtering is required
#'
#' @param pc a vector of scores from a given pc (i.e a single column of \eqn{cov.pool} in \eqn{get.conf()})
#' @param genes an \eqn{n X m} dataframe of cis and trans gene pairs of trios
#' @return an \eqn{1 X (m/2)} vector of pvalues from corresponding to the overall \eqn{F} of the regression
#' @export p.from.reg
#'
p.from.reg2=function(pc, genes){
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
  pstats=apply(fstats, 2, function(x){stats::pf(q = x[1], df1 = x[2], df2 = x[3], lower.tail = F)})
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
#' @param pi0.meth A string specifying 'smoother' or 'bootstrap'. Passed to qvalue::qvalue()
#' @return of input.type = list, then the output is a list of 2-lists containing the q values and sigificance. Else if the
#' input.type == array, then the output is a 2-list where the first element is a matrix of qvalues and the second element is
#' a matrix of significance determinations
#' @export get.q.sig

get.q.sig = function(pvalues, fdr.level, lambda.seq, pi0.method, contains.na = FALSE){

  if(contains.na == FALSE){
    adjust.out = apply(pvalues, 2, adjust.q, fdr = fdr.level, lambda = lambda.seq, pi0.meth = pi0.method)
  }else{
    adjust.out = apply(pvalues, 2, adjust.q.reg, fdr = fdr.level, lambda = lambda.seq, pi0.meth = pi0.method,
                       num.of.trios = dim(pvalues)[1])
  }

  #extract the qvalue and significance matrix
  cor.sig.mat = sapply(adjust.out, function(x) x$significant)
  cor.q.mat = sapply(adjust.out, function(x) x$qvalue)
  return(list(sigmat = cor.sig.mat, qmat = cor.q.mat))
}
















