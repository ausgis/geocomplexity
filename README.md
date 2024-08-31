
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geocomplexity <img src="man/figures/logo.png" align="right" height="120"/>

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/geocomplexity)](https://CRAN.R-project.org/package=geocomplexity)
[![r-universe](https://ausgis.r-universe.dev/badges/geocomplexity)](https://ausgis.r-universe.dev/geocomplexity)
<!-- badges: end -->

The goal of **geocomplexity** is to *support geocomplexity computation
and applications*.

## Overview

Full document of the most recent release of **geocomplexity** is online:
<https://ausgis.github.io/geocomplexity/>

## Installation

- Install development binary version from
  [r-universe](https://ausgis.r-universe.dev/geocomplexity) with:

``` r
install.packages('geocomplexity',
                 repos = c("https://ausgis.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install development source version from
  [GitHub](https://github.com/ausgis/geocomplexity) with:

``` r
# install.packages("devtools")
devtools::install_github("ausgis/geocomplexity",
                         build_vignettes = TRUE,
                         dep = TRUE)
```

Please ensure that appropriate **C++** compilation environment is
configured and **Rcpp** & **RcppArmadillo** is properly installed in
advance if you want to install **geocomplexity** from source. See the
next topic on C++ Settings for more information.

## Set up to use **C++** compilation environment

**Windows**:

- Download the Rtools installer that matches your version of R from
  <https://cran.r-project.org/bin/windows/Rtools/>
- Run the installer, `Rtools.exe`, keeping the default settings.

**Linux**

Debian/Ubuntu:

``` sh
apt-get update
apt-get install r-base-dev
```

Fedora/RedHat: should be set up already.

**MacOS**

Option 1

- [Register as an Apple developer (for
  free)](https://apps.apple.com/cn/app/apple-developer/id640199958/)
  - Then, in the terminal:

    ``` sh
    xcode-select --install
    ```

Option 2

- Install the current release of full [Xcode from the Mac App
  Store](https://apps.apple.com/ca/app/xcode/id497799835?mt=12)
- Within XCode go to Preferences -\> Downloads and install the Command
  Line Tools
- More convenient but installs a lot you don’t need
