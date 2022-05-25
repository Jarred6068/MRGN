#' A function to calculate the minor allele frequency
#'
#' This function takes in a vector V representing the counts of alternative alleles at an eQTL.
#' where 0 = homozygous for the reference allele, 1 = heterozygous for the reference allele,
#' and 2 = homozygous for the alternative allele. It returns the min(alternative, reference)
#' @param V a vector of alternative allele counts at a given eQTL
#' @return min(alternative,reference)
#' @examples
#' get.freq(M1trio$V1)
#' @export get.freq


####################################################################
#function To Get Minor Variant Frequency - section 1.1 -step 2
get.freq=function(V=NULL){
  #used in step 2
  #remove missing values
  V=S4Vectors::na.omit(V)
  #count the alternative alleles
  alternative = sum(V)/(2*length(V))
  reference=1-alternative
  #return the minor allele freq.
  return(min(alternative, reference))
}
