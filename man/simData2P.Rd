% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_trios_with_conf.R
\name{simData2P}
\alias{simData2P}
\title{A function to simulate a child node with two parents and confounding variables}
\usage{
simData2P(
  N,
  P1,
  P2,
  b0.1,
  b1.1,
  b1.2,
  sd.1,
  G = NULL,
  u = NULL,
  gamma.u = NULL
)
}
\arguments{
\item{N}{Simulation sample size}

\item{P1}{The data for the first parent node}

\item{P2}{the data for the second parent node}

\item{b0.1}{the intercept}

\item{b1.1}{the effect of the parent node P1}

\item{b1.2}{the effect of the parent node P2}

\item{sd.1}{the residual standard deviation (noise)}

\item{G}{A matrix of size N X g of confounding variables}

\item{u}{the indices of the confounders to be selected. If NULL confounders are omitted}

\item{gamma.u}{a vector of length g specifying the effects of the confounding variables}
}
\value{
a vector of length N containing the simulated values of the child node
}
\description{
A function to simulate a child node with two parents and confounding variables
}
