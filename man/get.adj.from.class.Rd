% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.adj.from.class.R
\name{get.adj.from.class}
\alias{get.adj.from.class}
\title{A function to obtain the adjacency matrix from an inferred trio class}
\usage{
get.adj.from.class(model = NULL, reg.vec = NULL, cnames = NULL)
}
\arguments{
\item{model}{a character string as returned by class.vec()}

\item{reg.vec}{a vector of inference results as returned by infer.trio()}

\item{cnames}{(optional) a character vector of names of the adjacency columns and rows}
}
\description{
This function takes in a the inferred model class from class.vec() and regression results from infer.trio() and returns
the correct adjacency matrix
}
