# highcor Package

## Introduction

`highcor` is an R package for direct estimation and inference of higher-level correlations from lower-level variables. Examples of higher-level variables include protein abundances and gene-pathway profiles; the corresponding lower-level variables are peptides and individual gene expression levels, respectively. 

## Installation

To install the latest version of `highcor` from GitHub, use:

```R
# install.packages("devtools")  # Uncomment this line if you haven't installed 'devtools' package
devtools::install_github("taryue/highcor")
```
## Usage 

Below are some basic examples of how to use the `highcor` package:

Example 1: Direct Estimation
This example shows how to perform a direct estimation using `direct_est()` function.

```R
library(highcor)
data(binding_mat)
data(gene_expression_stage1)
fit.est <- direct_est(Z=gene_expression_stage1, A=binding_mat)
```

Example 2: Optimal Shrinkage Estimation
This example shows how to perform an optimal shrinkage estimation using `opt_thres()` function based on the output of `direct_est()`.

```R
data(toy_lower_dat)
data(toy_binding)
fit.cv <- cv_opt_thres(A=toy_binding, Z=toy_lower_dat, n.split=20)
Sigma_hat <- direct_est(Z=toy_lower_dat, A=toy_binding)$cov.est
fit <- opt_thres(Sigma=Sigma_hat, A=toy_binding, Z=toy_lower_dat, kappa = fit.cv$kappa.min)
```

Example 3: Direct Inference

This example shows how to perform direct inference of the higher-level correlations using `direct_inf()` function based on the output of `direct_est()`.

```R
data(toy_lower_dat)
data(toy_binding)
Sigma_hat <- direct_est(Z=toy_lower_dat, A=toy_binding)$cov.est
inf <- direct_inf(Sigma=Sigma_hat, A=toy_binding, Z=toy_lower_dat, xi=c(0,0.1), parallel=FALSE)
```
