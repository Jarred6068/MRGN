% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.conf.R
\name{p.from.reg2}
\alias{p.from.reg2}
\title{A function to perform the regression of a possible confounding pc aganst all pairs of cis and trans genes given
in \eqn{genes} - when no filtering is required}
\usage{
p.from.reg2(pc, genes)
}
\arguments{
\item{pc}{a vector of scores from a given pc (i.e a single column of \eqn{cov.pool} in \eqn{get.conf()})}

\item{genes}{an \eqn{n X m} dataframe of cis and trans gene pairs of trios}
}
\value{
an \eqn{1 X (m/2)} vector of pvalues from corresponding to the overall \eqn{F} of the regression
}
\description{
A function to perform the regression of a possible confounding pc aganst all pairs of cis and trans genes given
in \eqn{genes} - when no filtering is required
}
