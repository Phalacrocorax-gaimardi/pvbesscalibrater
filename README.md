
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pvbesscalibrater

<!-- badges: start -->
<!-- badges: end -->

The goal of pvbesscalibrater is to create a set of weights and empirical
partial utilities used by pvbessmicrosimr. This is the micro-calibration
step of the ABM. Some of the tasks that pvbesscalibrater handles include

- mapping survey adoption likert scores to partial utilities using a
  parameter (epsilon) that describes hypothetical bias
- choice of financial variable(s) in ABM (usually highest bill q14,
  annual bill q_ab or highest and lowest bills q14 and q15)
- finds a set of regularised weights for financial, social and barrier
  terms

## Installation

You can install the development version of pvbesscalibrater like so:

``` r
remotes::install_github(")
```

## Example

This explains the workflow.

``` r
library(pvbesscalibrater)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
