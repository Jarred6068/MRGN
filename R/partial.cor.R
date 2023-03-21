#' @name partial.cor
#' @aliases marginal2partial.cor
#' @title Calculate partial correlations (matrix)
#'
#' @description \code{partial.cor} calculates the partial correlation matrix
#' (i.e. the correlation between each pair of variables given all other
#' variables) starting from a covariance matrix. \code{marginal2partial.cor}
#' converts a marginal correlation matrix (with the associated vector of
#' standard deviations) or the corresponding covariance matrix into a partial
#' correlation matrix.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param R numeric \code{p}-by-\code{p} symmetric matrix (marginal correlation matrix).
#' @param sd numeric vector of positive scalars, standard deviations of the
#' variables for which \code{R} is the correlation matrix.
#' @param cov numeric matrix (marginal covariance matrix). When supplied and not
#' \code{NULL}, \code{cov} is used whereas \code{R} and \code{sd} are ignored.
#' @param ... additional arguments passed to \link{cov}. See Section \code{Details}.
#'
#' @export partial.cor
#' @export marginal2partial.cor
#'
#' @usage
#' partial.cor (x,
#'              ...)
#'
#' marginal2partial.cor (R,
#'                       sd = rep(1, p),
#'                       cov = NULL)
#'
#' @details The arguments \code{...} include 'use' (which can be one of the
#' strings "everything", "all.obs", "complete.obs", "na.or.complete", or
#' "pairwise.complete.obs", or an abbreviation of these strings) and 'method'
#' (which can here be one of "pearson", or "spearman"). Note that these
#' arguments are not checked before they are passed to \link{cov}.
#'
#' @return
#' Either \code{partial.cor} or \code{marginal2partial.cor} returns a square
#' matrix of partial correlations.
#'
#' @seealso \link{pairwise.cor.test}.
#'

# A routine to calculate partial correlations
# ... are passed to cov
# they can be: use = "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs",
#              method = "pearson", "kendall" or "spearman"
partial.cor <- function(x, ...) {
  H <- - cov2cor(MRGN:::mpinv(cov(x, ...)))
  diag(H) <- 1
  dimnames(H) <- list(colnames(x), colnames(x))
  H
}

marginal2partial.cor <- function (R, sd = rep(1, p), cov = NULL) {
  if (missing(cov) || is.null(cov)) {
    p <- (d <- dim(R))[1L]
    if (!is.numeric(R) || length(d) != 2L || p != d[2L])
      stop("'R' is not a square numeric matrix")
    if (!is.numeric(sd) || length(sd) != p || any(sd <= 0))
      stop("'sd' is not a numeric vector of positive reals, or of a size consistent with 'R'")
    if (p <= 2)
      return(R)
    # Convert R into a covariance matrix
    cov <- sd * R * rep(sd, each = p)
  }
  else {
    p <- (d <- dim(cov))[1L]
    if (!is.numeric(cov) || length(d) != 2L || p != d[2L])
      stop("'cov' is not a square numeric matrix")
    if (p <= 2)
      return(cov)
  }

  # Get the partial correlations
  H <- - cov2cor(mpinv(cov))
  diag(H) <- 1
  # dimnames(H) <- dimnames(cov)
  H
}

# A routine to calculate Moore-Penrose Generalized inverse
# Not exported
mpinv <- function(X, eps = NULL) {
  # Perform singular value decomposition
  s <- svd(X)
  d <- s$d

  # Number of singular values (min(dim(X)))
  m <- length(d)

  # remove singular values ~ zero
  if (missing(eps) | is.null(eps))
    eps <- .Machine$double.eps * max(NROW(X), NCOL(X)) * max(d)
  d <- d[d > eps]
  n <- length(d)

  # Inverse of positive singular values
  inv <- if (n > 1) diag(1/d) else 1/d

  # Add rows, columns and rows of zeros if X has null singular values
  if (n != m) {
    inv <- cbind(inv, matrix(0, nrow=n, ncol=(m - n)))
    inv <- rbind(inv, matrix(0, nrow=(m-n), ncol=m))
  }

  # compute the Moore-Penrose inverse
  inv <- s$v %*% inv %*% t(s$u)

  # set very small values to zero
  inv[abs(inv) < eps] <- 0
  return(inv)
}
