#' @title 8,823 trios formed from the GTEx tissue WholeBlood
#' @references The GTEx Consortium. GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2). https:// gtexportal.org/home/datasets, 2021. [Online; accessed 8-June-2020]
#' @description A list of dataframes where each dataframe represents a single trio from one of the 5 model topologies
#'
#' @format a list of length = 8,823 where each element is an 670 X 3 dataframe where the first
#'  column is an eQTL representing V1 and columns 2:3 are the cis and trans targets representing
#'  T1 and T2 (respectively)
"WBtrios"

#' @title A simulated M1 trio with 1000 observations
#'
#' @description A 1000 X 3 dataframe where V1 is a simulated eQTL and T1 and T2 are simulated molecular phenotypes
#'
#' @format dataframe 650 X 3
"M1trio"

#' @title The rotated variables from the whole-genome level expression PCA for the GTEx tissue Whole Blood
#'
#' @description A 670 X 670 dataframe of PCA scores where the rows are observations and the columns are the PCs
#'
#' @format dataframe 670 X 670
"WBscores"
