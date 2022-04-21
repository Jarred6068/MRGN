#'
#' A function to perform PCA on the whole-genome expression or methylation data
#'
#' This function takes in a matrix representing the whole-genome expression or methylation data and returns the
#' rotated variables from the Principle Component Analysis
#'
#' @param x a matrix of df with samples in rows and genes in columns
#' @param center (logical) default = TRUE
#' @param scale (logical) default = FALSE
#' @return returns a matrix/df containing the rotated variables of the PCA
#' @export get.PCs
