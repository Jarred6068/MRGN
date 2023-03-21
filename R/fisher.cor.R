#' @name fisher.cor
#' @title Test for (partial) correlation between variables
#'
#' @description A function to (optionally) calculate the (partial) pearson
#' correlations coefficients and the associated p-values using Fisher's
#' z-transformation.
#'
#' If a type error level \code{alpha} is specified, the test
#' conclusions (accept/reject = 0/1) are returned. In this case, the overall
#' type I error rate can (optionally) be controlled, instead of controlling
#' just individual type I error rates.
#'
#' @param cor.mat numeric (optional), \code{m}-by-\code{m} correlation matrix.
#' @param n integer (optional), size of the samples used to calculate
#' \code{cor.mat}. The argument \code{n} must be specified for \code{cor.mat}
#' to be used (specifying \code{cor.mat} but not \code{n} and \code{data} wiil
#' result in an error).
#' @param data a numeric (optional), \code{n}-by-\code{m} matrix or data.frame
#' (\code{m} > 1). The argument \code{data} must be supplied when either of
#' \code{cor.mat} or \code{n} is missing. When \code{cor.mat} and \code{n} are
#' given, \code{data} is ignored.
#' @param cols (optional) a numeric or logical vector indicating the target
#' columns in \code{data}. Defaults to all columns.
#' @param partial logical, should partial correlations be considered? This is
#' used only when \code{cor.mat} or \code{n} is missing. To test partial
#' correlations with a marginal correlation matrix at hand, first convert the
#' marginal correlation matrix into a partial correlation matrix.
#' @param alpha (optional) a scalar in the open (0, 1), significance level. When
#' set to \code{NULL} or not supplied, individual p-values are returned for all
#' pairwise correlation coefficients.
#' @param overall logical, should the overall type I error rate be controlled?
#' Only used when the type I error rate \code{alpha} is supplied. Defaults to
#' \code{TRUE}. See Section \code{Details} for more on the control process.
#'
#' @usage
#' fisher.cor (data,
#'             cols,
#'             cor.mat = NULL,
#'             n = NULL,
#'             partial = FALSE,
#'             alpha = NULL,
#'             overall = TRUE)
#'
#' @export fisher.cor
#'
#' @details
#' To conclude the results of the tests for significance of the correlations
#' coefficients, the argument \code{alpha} must be supplied. When \code{alpha}
#' is given and \code{overall = FALSE}, each test is decided at the nominal
#' level \code{alpha}, and a matrix of the test results is also returned. If
#' \code{overall = TRUE} with \code{alpha} supplied, the the LOND (Level On the
#' Number of Discoveries) method (Javanmard and Montanari, 2015), which is a
#' sequential hypothesis testing procedure and sets value of alpha for each test
#' based on the number of rejections to control the overall type I error rate.
#'
#' @return A \code{list} with at least two elements:
#' \item{pvalue}{a numeric vector, containing the lower triangular part of the
#' symmetric matrix of pairwise p-values.}
#' \item{cor}{a symmetric matrix with pairwise partial correlation coefficients.}
#'
#' If the level of significance \code{alpha} is specified, the returned object
#' additionally has element: \item{Adj}{adjacency matrix of an undirected graph
#' where the target variables are nodes.}
#'
#' @references Javanmard A and Montanari A (2015). On Online Control of False Discovery Rate. arXiv:150206197 [statME].
#'

fisher.cor <- function (data, cols,
                        cor.mat = NULL, n = NULL,
                        partial = FALSE,
                        alpha = NULL, overall = TRUE) {
  if (any(missing(cor.mat), is.null(cor.mat), is.null(n), !is.numeric(n))) {
    if (missing(data))
      stop ("One of arguments 'data' and 'cor.mat' must be specified.")

    # Data
    data <- as.data.frame(data)
    n <- NROW(data)
    if (!missing(cols))
      data <- data[cols, cols]
    p <- NCOL(data)
    if (p == 1)
      stop("At least two target columns are required in 'data'.")
    if (n <= p + 1)
      stop(paste0("For ", p, " target variables (columns), at least ",
                  p + 2, " observations (rows) are required in 'data'. /n"))

    # Calculate correlation Matrix
    if (partial)
      cor.mat <- partial.cor (data, use = 'pairwise.complete.obs')
    else
      cor.mat <- cor (data, use = 'pairwise.complete.obs')
  }
  else {
    stopifnot(isSymmetric(cor.mat))
    p <- NCOL(cor.mat)
    if (p == 1)
      stop("At least two target columns are required in 'data'.")
    n <- n[1]
    if (n <= p + 1)
      stop(paste0("For ", p, " target variables (columns), at least ",
                  p + 2, " observations (rows) are required in 'data'. /n"))
    cor.mat <- cov2cor(cor.mat)
  }

  # Calculate z-values from correlation values
  targ <- lower.tri(cor.mat, diag = FALSE)
  df <- if (partial) p - 2 else 0
  zp <- .5 * (log1p(cor.mat[targ]) - log1p(-cor.mat[targ])) * sqrt(n - df - 3)

  # Calculate each p-value for a two-tailed test from Z ~ N(0,1)
  zp <- 2 * stats::pnorm (q = abs(zp), lower.tail = FALSE)

  # Return the results if 'alpha = NULL'
  if (is.null(alpha)) {
    res <- list(pvalue = zp, cor = cor.mat)
    if (missing(alpha))
      return(res)
  }

  # Check the argument 'alpha'
  stopifnot(is.numeric(alpha))
  alpha <- alpha[1]
  if (alpha <= 0 | alpha >= 1)
    stop("Invalid argument 'alpha': only values in the open (0, 1) are allowed.")

  # Initialize the adjacency matrix
  A <- matrix(0, nrow = p, ncol = p)
  dimnames(A) <- dimnames(cor.mat)

  if (!overall) { # Running individual tests are at level a fixed level
    A[targ] <- zp <= alpha
  }
  else { # LOND method to control the overall type I error rate
    C <- alpha * 6 / (pi^2)

    # Initialize a decision vector
    Decisions <- zp

    # A counter for rejections
    R <- 0

    # Setting an arbitrary testing order based on distance from average p-value
    ordp <- order(abs(zp - mean(zp)))

    # Decide each test
    for (i in ordp) {
      alpha_i <- (1 + sum(R)) * C / (i^2)
      Decisions[i] <- zp[i] <= alpha_i
      R <- R + 1 - Decisions[i]
    }
    A[targ] <- Decisions
  }
  return(list(Adj = structure(A + t(A), class = 'adjacency.matrix'),
              pvalue = zp, cor = cor.mat))
}

