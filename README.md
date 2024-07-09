
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geocomplexity <img src="man/figures/logo.png" align="right" height="120"/>

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/geocomplexity)](https://CRAN.R-project.org/package=geocomplexity)
[![r-universe](https://spatlyu.r-universe.dev/badges/geocomplexity)](https://spatlyu.r-universe.dev/geocomplexity)
<!-- badges: end -->

The goal of **geocomplexity** is to *support geocomplexity computation
and applications*.

## Overview

Full document of the most recent release of **geocomplexity** is online:
<https://spatlyu.github.io/geocomplexity/>

## Installation

Please install **Rcpp** and the corresponding **C++** compilation
environment in advance.

You can install the development version of **geocomplexity** from
[github](https://github.com/SpatLyu/geocomplexity) with:

``` r
# install.packages("devtools")
devtools::install_github("SpatLyu/geocomplexity",
                         dep = TRUE)
```

You can also install the binary version of **geocomplexity** from
[r-universe](https://spatlyu.r-universe.dev/geocomplexity):

``` r
install.packages("geocomplexity", 
                 repos = c("https://spatlyu.r-universe.dev",
                           "https://cran.rstudio.com/"),
                 dep = TRUE)
```