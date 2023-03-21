#' @name conf.control
#' @title Auxiliary for controlling the selection of confounding variables for trios
#' @description Auxiliary function for (generating and) checking control parameters
#' for \link{get.conf}. Typically only used internally by \link{MRGN}, but may be
#' used to construct control arguments for \link{get.conf}.
#'
#' @usage
#' conf.control (blocksize = 2000,
#'               apply.qval = TRUE,
#'               selection_fdr = 0.05,
#'               filter_int_child = FALSE,
#'               filter_fdr = 0.1,
#'               lambda = NULL,
#'               alpha = 0.05,
#'               method = c('correlation', 'regression'),
#'               return.for.trios = TRUE,
#'               save.list = FALSE,
#'               save.path = '/path/to/save/location')
#'
#' @param blocksize a positive integer, the number of columns to use in each block
#' of correlation calculations. It is passed to \link{propagate::bigcor}.
#' Defaults to \code{2000}.
#'
#' @param apply.qval logical, should the \code{qvalue} adjustment be applied to
#' each set of correlations between a confounder and the variables forming each
#' trio? Defaults to \code{TRUE}. If \code{FALSE}, Bonferroni correction is used.
#'
#' @param selection_fdr real in the open (0, 1), the false discovery rate for
#' selecting confounders. Defaults to \code{0.05}
#'
#' @param filter_int_child logical, should common child and intermediate
#' confounders be filtered out (removed) from the significant confounders for
#' each trio? Only used when \code{return.for.trios = TRUE}. Defaults to
#' \code{FALSE}.
#'
#' @param filter_fdr a real in the open (0, 1), the false discovery rate for
#' filtering common child and intermediate confounding variables. Only used
#' when \code{method = 'regression'}. Defaults to \code{0.1}.
#'
#' @param lambda \code{NULL} or a real in the open (0, 1), the cut off points
#' of the tuning parameter to estimate \eqn{pi_0}. The default value \code{NULL}
#' means that \eqn{seq(0.5, max(pvalues), 0.05)} is passed to \link{adjust.q}.
#'
#' @param alpha a real in the open (0, 1), the type I error rate for the
#' Bonferroni corrected pvalues. Only used when \code{apply.qval = FALSE}.
#' Defaults to \code{0.05}.
#'
#' @param method a character string indicating the method to evaluate the
#' association between potential confounders and each trio. Either "correlation"
#' or "regression". Note that when \code{method = 'regression'}, the argument
#' \code{return.for.trios} is forced to \code{TRUE}.
#'
#' @param return.for.trios logical, should the column indices of the confounders
#' associated with each trio be returned? If \code{FALSE}, the column indices of
#' confounders associated with each variable in each \code{trios} is
#' returned. Defaults to \code{TRUE}.
#'
#' @param save.list logical, should the output be saved as a \code{.RData}
#' object? Defaults to \code{FALSE}.
#'
#' @param save.path a character string indicating the location (and the path to
#' this location) to save the output. Ignored when \code{save.list = FALSE}.
#'
#' @return A list with components named as the arguments.
#'
#' @export conf.control
#'
#' @seealso \link{get.conf}.
#

conf.control <- function(blocksize = 2000,
                         apply.qval = TRUE,
                         selection_fdr = 0.05,
                         filter_int_child = FALSE,
                         filter_fdr = 0.1,
                         lambda = NULL,
                         alpha = 0.05,
                         method = c('correlation', 'regression'),
                         return.for.trios = TRUE,
                         save.list = FALSE,
                         save.path = '/path/to/save/location') {
  if (any(c(!is.numeric(blocksize), blocksize <= 0)))
    stop("Argument 'blocksize' must be a positive integer.")
  if ((!is.numeric(apply.qval) && !is.logical(apply.qval)))
    stop("The value of 'apply.qval' must numeric or logical.")
  apply.qval <- apply.qval[1] > 1
  if (any(c(!is.numeric(selection_fdr), selection_fdr <= 0, selection_fdr >= 1)))
    stop("The value of 'selection_fdr' must be in the open (0, 1).")
  if ((!is.numeric(filter_int_child) && !is.logical(filter_int_child)))
    stop("The value of 'filter_int_child' must numeric or logical.")
  filter_int_child <- filter_int_child[1] > 1
  if (any(c(!is.numeric(filter_fdr), filter_fdr <= 0, filter_fdr >= 1)))
    stop("The value of 'filter_fdr' must be in the open (0, 1).")
  if (!is.null(lambda)) {
    if (any(c(!is.numeric(lambda), lambda <= 0, lambda >= 1)))
      stop("The value of 'lambda' must be in the open (0, 1).")
  }
  if (any(c(!is.numeric(alpha), alpha <= 0, alpha >= 1)))
    stop("The value of 'alpha' must be in the open (0, 1).")
  if (!all(c(is.character(method), method %in% c('correlation', 'regression'))))
    stop("The value of 'method' must be a character, one of: 'correlation' or 'regression'.")
  if ((!is.numeric(return.for.trios) && !is.logical(return.for.trios)))
    stop("The value of 'return.for.trios' must numeric or logical.")
  return.for.trios <- return.for.trios[1] > 1
  if ((!is.numeric(save.list) && !is.logical(save.list)))
    stop("The value of 'save.list' must numeric or logical.")
  save.list <- save.list[1] > 1
  if (save.list) {
    if (missing(save.path))
      save.list <- FALSE
    if (!is.character(save.path))
      stop("Argument 'save.path' must be a character in the for '/path/to/save/location'.")
    save.path <- save.path[1]
  }
  list(blocksize = blocksize[1],
       apply.qval = apply.qval,
       selection_fdr = selection_fdr[1],
       filter_int_child = filter_int_child,
       filter_fdr = filter_fdr[1],
       lambda = lambda[1],
       alpha = alpha[1],
       method = method[1],
       return.for.trios = return.for.trios,
       save.list = save.list,
       save.path = save.path)
}
