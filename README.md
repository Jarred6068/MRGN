
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRGNtrio: Mendelian Randomization Genomic Network for Trios

<!-- ![](./static/MRGN-logos.jpeg width="800" height="600" style="display: block; margin: 0 auto") -->

<p align="center">
    <img src="static/MRGN-logos.jpeg" width="50%" height="50%" />
</p>
<!--  -->
<!-- badges: start --> <!-- badges: end -->

## Overview

MRGNtrio is a novel software tool designed to infer potential causal relationships between two biological variables, such as the expression of two genes or between molecular phenotypes like gene expression and DNA methylation. The software is built on the principle of Mendelian randomization, which is a robust approach to assess causality in genomic networks. MRGNtrio utilizes a genetic variant as an instrumental variable to infer 5 mutually exclusive causal models for a genomic trio. 

<p align="center">
    <img src="static/5 causal models.drawio.png" width="60%" height="60%" />
</p>

One of the key challenges in causal inference is accounting for confounding variables $\bf U$, which can significantly impact the results. MRGNtrio addresses this challenge by integrating a regression-based framework for conditional dependence testing that can handle a large number of confounding variables effectively.

$$ T_1 = \alpha_1 +\beta_{11}V+\beta_{12}T_2+{{\bf \Gamma}_1 {\bf U}}+ \epsilon_1; \ \ (1)$$

$$ T_2 = \alpha_2 +\beta_{21}V+\beta_{22}T_1+{{\bf \Gamma}_2} {\bf U}+\epsilon_2  \ \ (2)$$

**Table: Conditional and marginal dependence tests utilized by MRGNtrio**  
*A reference table of the conditional and marginal dependence tests used to determine which model is supported by the data. The first three models each have two possible configurations. Coefficients β₁₁ and β₁₂ are from regression (Eq. 1), and coefficients β₂₁ and β₂₂ from regression (Eq. 2). These coefficients test conditional independence. Correlations $r_{V, T_1}$ and $r_{V, T_2}$ test marginal independence.*

| Model | Configuration                                | β₁₁    | β₁₂    | β₂₁    | β₂₂    | $r_{V, T_1}$ | $r_{V,T_2}$ |
|-------|----------------------------------------------|--------|--------|--------|--------|-----------|-----------|
|       | **Dependence Relation**                      | T₁ ⫫ V | T₁ ⫫ T₂ | T₂ ⫫ V | T₂ ⫫ T₁ | V ⫫ T₁    | V ⫫ T₂    |
|       | **Conditioned on**                           | T₂, **U** | V, **U** | T₁, **U** | V, **U** | --        | --        |
| M₀    | V → T₁; T₂ is singleton                      | ≠ 0   | = 0    | = 0    | = 0    | --        | --        |
|       | V → T₂; T₁ is singleton                      | = 0    | = 0    | ≠ 0   | = 0    | --        | --        |
| M₁    | V → T₁ → T₂                                  | ≠ 0   | ≠ 0   | = 0    | ≠ 0   | --        | --        |
|       | V → T₂ → T₁                                  | = 0    | ≠ 0   | ≠ 0   | ≠ 0   | --        | --        |
| M₂    | V → T₁ ← T₂                                  | ≠ 0   | ≠ 0   | ≠ 0   | ≠ 0   | ≠ 0       | = 0       |
|       | V → T₂ ← T₁                                  | ≠ 0   | ≠ 0   | ≠ 0   | ≠ 0   | = 0       | ≠ 0       |
| M₃    | T₁ ← V → T₂                                  | ≠ 0   | = 0    | ≠ 0   | = 0    | --        | --        |
| M₄    | T₁ ← V → T₂; T₁ ↔ T₂                         | ≠ 0   | ≠ 0   | ≠ 0   | ≠ 0   | ≠ 0       | ≠ 0       |




## Installation

You can install MRGNtrio from GitHub release with

``` r
install.packages("devtools")
devtools::install_github("Jarred6068/MRGN")
```

## Example

``` r
library(MRGN)
#inference on a single eQTL trio
result=infer.trio(M1trio)
print(result)
 
#fast example on 10 eQTL trios from the built in dataset WBtrios
results = sapply(WBtrios[1:10], function(x) infer.trio(x))
print(results)
#return just the inferred model topology
models = sapply(WBtrios[1:10], function(x) infer.trio(x)$Inferred.Model)
print(models)
#fast example on 10 CNA trios from the built in dataset CNAtrios using permutation
models = sapply(CNAtrios[1:10], function(x) infer.trio(x, is.CNA = T)$Inferred.Model)
print(models)

```

