
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRGNtrio: Mendelian Randomization Genomic Network for Trios

<!-- ![](./static/MRGN-logos.jpeg width="800" height="600" style="display: block; margin: 0 auto") -->

<p align="center">
    <img src="static/MRGN-logos.jpeg" width="60%" height="60%" />
</p>
<!--  -->
<!-- badges: start --> <!-- badges: end -->

## Overview

MRGN (Mendelian Randomization for Genomic Networks) is a novel software tool designed to infer potential causal relationships between two biological variables, such as the expression of two genes or between molecular phenotypes like gene expression and DNA methylation. The software is built on the principle of Mendelian randomization, which is a robust approach to assess causality in genomic networks. MRGNtrio utilizes a genetic variant as an instrumental variable under a regression framework to infer 5 mutually exclusive causal models for a genomic trio. 

\[ T_1 = \alpha_1 +\beta_{11}V+\beta_{12}T_2+{{\bf \Gamma}_1 {\bf U}}+ \epsilon_1;\]

One of the key challenges in causal inference is accounting for confounding variables $\bf U$, which can significantly impact the results. MRGN addresses this challenge by integrating a regression-based method that can handle a large number of confounding variables effectively.

## Installation

You can install MRGNtrio from github release with

``` r
# install.packages("devtools")
devtools::install_github("Jarred6068/MRGN")
```

## Example

``` r
library(MRGN)
#inference on a single eQTL trio
result=infer.trio(M1trio)
print(result)

## Not run: 
#fast example on 10 eQTL trios from the built in dataset WBtrios
results = sapply(WBtrios[1:10], function(x) infer.trio(x))
print(results)
#return just the inferred model topology
models = sapply(WBtrios[1:10], function(x) infer.trio(x)$Inferred.Model)
print(models)
#fast example on 10 CNA trios from the built in dataset CNAtrios using permutation
models = sapply(CNAtrios[1:10], function(x) infer.trio(x, is.CNA = T)$Inferred.Model)
print(models)

## End(Not run)
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
