#' @name pairwise.cor.test
#' @title Pairwise tests for significant (partial) correlations
#'
#' @description The function calculates the (partial) pearson correlations
#' coefficients (optionally) and the associated p-values. If a type error level
#' \code{alpha} is specified, the test conclusions (accept = 0, or reject = 1)
#' are returned. In this case, the overall type I error rate can (optionally)
#' be controlled with the argument \code{overall}, instead of controlling just
#' individual type I error rates.
#'
#' @param data a numeric (optional), \code{n}-by-\code{m} matrix or data.frame
#' (\code{m} > 1). The argument \code{data} must be supplied when either of
#' \code{cor.mat} or \code{n} is missing. When \code{cor.mat} and \code{n} are
#' given, \code{data} is ignored.
#' @param cols (optional) a numeric or logical vector indicating the target
#' columns in \code{data}. Defaults to all columns.
#' @param cor.mat numeric (optional), \code{m}-by-\code{m} correlation matrix.
#' @param n integer (optional), size of the samples used to calculate
#' \code{cor.mat}. The argument \code{n} must be specified for \code{cor.mat}
#' to be used (specifying \code{cor.mat} but not \code{n} and \code{data} wiil
#' result in an error).
#' @param partial logical, should partial correlations be considered? This is
#' considered when \code{data} is used (and \code{cor.mat} or \code{n} is
#' missing). To test partial correlations with a marginal correlation matrix at
#' hand (and not the original data), first convert the marginal correlation
#' matrix into a partial correlation matrix (see \link{marginal2partial.cor}).
#' @param alternative logical, indicates the alternative hypothesis and must be
#' one of "two.sided", "greater" or "less". You can specify just the initial
#' letter. "greater" corresponds to positive association, "less" to negative
#' association.
#' @param alpha (optional) a scalar in the open (0, 1), significance level. When
#' set to \code{NULL} or not supplied, individual p-values are returned for all
#' pairwise correlation coefficients.
#' @param overall logical, should the overall type I error rate be controlled?
#' Only used when the type I error rate \code{alpha} is supplied. Defaults to
#' \code{TRUE}. See Section \code{Details} for more on the control process.
#' @param method a character string indicating the test to perform. One of
#' 'z.test' (using the asymptotic normal distribution to compute p-values,
#' after Fisher's z-transformation) or 't.test' (using Student's t distribution
#' to compute p-values).
#' @param cor.method a character string indicating which correlation coefficient
#' is to be computed for the tests. One of "pearson", or "spearman". Only used
#' along with \code{data}.
#'
#' @usage
#' pairwise.cor.test (data,
#'                    cols,
#'                    cor.mat = NULL,
#'                    n = NULL,
#'                    partial = TRUE,
#'                    alternative = c('two.sided', 'less', 'greater'),
#'                    method = c('z.test', 't.test'),
#'                    cor.method = c("pearson", "spearman"),
#'                    alpha = NULL,
#'                    overall = TRUE)
#'
#' @export pairwise.cor.test
#'
#' @details The method \code{z.test} is for very large samples and \code{t.test}
#' is for small samples.
#'
#' To conclude the tests for significance of the correlations coefficients, the
#' argument \code{alpha} must be supplied. When \code{alpha} is given and
#' \code{overall = FALSE}, each test is decided at the nominal level \code{alpha},
#' and a matrix of the test results is returned. If \code{overall = TRUE} with
#' \code{alpha} supplied, the LOND (Level On the Number of Discoveries) method
#' is used. LOND is a sequential hypothesis testing procedure and sets value of
#' alpha for each test based on the number of rejections to control the overall
#' type I error rate (Javanmard and Montanari, 2015).
#'
#' @return A \code{list} with at least two elements:
#' \item{pvalue}{a numeric vector, containing the lower triangular part (as
#' returned by \link{lower.tri}) of the symmetric matrix of pairwise p-values.}
#' \item{cor}{a symmetric matrix with pairwise partial correlation coefficients.}
#'
#' If the level of significance \code{alpha} is specified, the returned object
#' additionally has element: \item{Adj}{adjacency matrix of an undirected graph
#' where the target variables are nodes.}
#'
#' @references Javanmard A and Montanari A (2015). On Online Control of False
#' Discovery Rate. arXiv:150206197 [statME].
#'
#' @seealso \link{partial.cor}.
#'

pairwise.cor.test <- function (data, cols,
                               cor.mat = NULL, n = NULL,
                               partial = TRUE,
                               alternative = c('two.sided', 'less', 'greater'),
                               method = c('z.test', 't.test'),
                               cor.method = c('pearson', 'spearman'),
                               alpha = NULL, overall = TRUE) {
  if (any(missing(cor.mat), is.null(cor.mat), is.null(n), !is.numeric(n))) {
    if (missing(data))
      stop ("One of the argument sets ('data') and ('cor.mat' and 'n') must be specified.")

    # Data
    data <- as.data.frame(data)
    n <- NROW(data)
    if (!missing(cols))
      data <- data[, cols, drop = FALSE]
    p <- NCOL(data)
    if (p == 1)
      stop("At least two target columns are required in 'data'.")
    if (n <= p + 1)
      stop(paste0("For ", p, " target variables (columns), at least ",
                  p + 2, " observations (rows) are required in 'data'. /n"))

    # Check the argument 'cor.method'
    if (!missing(cor.method)) {
      stopifnot(is.character(cor.method), cor.method %in% c('pearson', 'spearman'))
    }

    # Calculate correlation Matrix
    if (partial)
      cor.mat <- partial.cor (data, use = 'pairwise.complete.obs',
                              method = cor.method[1])
    else
      cor.mat <- cor (data, use = 'pairwise.complete.obs',
                      method = cor.method[1])
  }
  else {
    stopifnot(isSymmetric(cor.mat, check.attributes = FALSE))
    p <- NCOL(cor.mat)
    if (p == 1)
      stop("At least two target columns are required in 'data'.")
    n <- n[1]
    if (n <= p + 1)
      stop(paste0("For ", p, " target variables (columns), at least ",
                  p + 2, " observations (rows) are required in 'data'. /n"))
    cor.mat <- cov2cor(cor.mat)
  }

  # Check the arguments 'method' and 'alternative'
  if (!missing(method)) {
    stopifnot(is.character(method), method %in% c('z.test', 't.test'))
  }
  if (!missing(alternative)) {
    stopifnot(is.character(alternative), alternative %in%
                c('two.sided', 't', 'less', 'l', 'greater', 'g'))
  }

  # An indicator for the lower triangular part of 'cor.mat'
  targ <- lower.tri(cor.mat, diag = FALSE)

  if (identical(method[1], 'z.test')) {
    # Degrees of freedom for the tests
    df <- if (partial) p - 2 else 0

    # Calculate z-values from correlation values
    zp <- .5 * sqrt(n - df - 3) * (log1p(cor.mat[targ]) - log1p(-cor.mat[targ]))

    # Calculate each p-value for a two-tailed test from Z ~ N(0,1)
    zp <- switch(alternative[1],
                 two.sided = 2 * stats::pnorm (q = abs(zp), lower.tail = FALSE),
                 t         = 2 * stats::pnorm (q = abs(zp), lower.tail = FALSE),
                 less      = stats::pnorm (q = zp, lower.tail = TRUE),
                 l         = stats::pnorm (q = zp, lower.tail = TRUE),
                 greater   = stats::pnorm (q = zp, lower.tail = FALSE),
                 g         = stats::pnorm (q = zp, lower.tail = FALSE))
  }
  else {
    # Degrees of freedom for the tests
    df <- if (partial) p - 2 else 0

    # Calculate t-values from correlation values
    zp <- sqrt(n - df - 2) * cor.mat[targ] / sqrt(1 -cor.mat[targ]^2)

    # Calculate each p-value for a two-tailed test from Z ~ N(0,1)
    zp <- switch(alternative[1],
                 two.sided = 2 * stats::pt (q = abs(zp), df = n - df - 2, lower.tail = FALSE),
                 t         = 2 * stats::pt (q = abs(zp), df = n - df - 2, lower.tail = FALSE),
                 less      = stats::pt (q = zp, df = n - df - 2, lower.tail = TRUE),
                 l         = stats::pt (q = zp, df = n - df - 2, lower.tail = TRUE),
                 greater   = stats::pt (q = zp, df = n - df - 2, lower.tail = FALSE),
                 g         = stats::pt (q = zp, df = n - df - 2, lower.tail = FALSE))
  }

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

  if (!overall) {
    # Running individual tests each at the fixed level alpha
    A[targ] <- zp <= alpha
  }
  else {
    # LOND method to control the overall type I error rate
    # The constant C to control the error rate
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
