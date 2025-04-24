
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRGNtrio: Mendelian Randomization Genomic Network for Trios

![](./static/MRGN-logos.jpeg width="800" height="600" style="display: block; margin: 0 auto")
<!-- <img src="static/MRGN-logos.jpeg" width="50%" height="50%" /> -->
<!-- badges: start --> <!-- badges: end -->


MRGN (Mendelian Randomization for Genomic Networks) is a novel software tool designed to infer potential causal relationships between two biological variables, such as the expression of two genes or between molecular phenotypes like gene expression and DNA methylation. By utilizing a genetic variant as an instrumental variable, MRGN forms a trio that facilitates the causal inference process.

The software is built on the principle of Mendelian randomization (MR), which is a robust approach to assess causality in genomic networks. One of the key challenges in causal inference is accounting for confounding variables, which can significantly impact the results. MRGN addresses this challenge by integrating a regression-based method that can handle a large number of confounding variables effectively.

MRGNtrio, the core component of the software, allows for the detection of diverse causal models for genomic trios using individual-level data. It offers powerful inference capabilities while maintaining computational efficiency, making it suitable for for high throughput analyses. 

## Installation

current version `v1.3.7.3`

You can install MRGNtrio from with

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
