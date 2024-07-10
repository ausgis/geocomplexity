
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

#### From source:

You can install the development version of **geocomplexity** from
[github](https://github.com/SpatLyu/geocomplexity) with:

``` r
# install.packages("devtools")
devtools::install_github("SpatLyu/geocomplexity",
                         dep = TRUE)
```

Please ensure that **Rcpp** is properly installed and the appropriate
**C++** compilation environment is configured in advance if you want to
install **geocomplexity** from
[github](https://github.com/SpatLyu/geocomplexity). See the next topic
on C++ Settings for more information.

#### From binary

You can also install the binary version of **geocomplexity** from
[r-universe](https://spatlyu.r-universe.dev/geocomplexity):

``` r
install.packages("geocomplexity", 
                 repos = c("https://spatlyu.r-universe.dev",
                           "https://cran.rstudio.com/"),
                 dep = TRUE)
```

## Set up to use **C++** compilation environment

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
  free)](https://developer.apple.com/programs/register/)
  - Then, in the terminal:

    ``` sh
    xcode-select --install
    ```

Option 2

- Install the current release of full [Xcode from the Mac App
  Store](https://itunes.apple.com/ca/app/xcode/id497799835?mt=12)
- Within XCode go to Preferences -\> Downloads and install the Command
  Line Tools
- More convenient but installs a lot you don’t need

**Windows**:

- Download the Rtools installer that matches your version of R from
  <https://cran.r-project.org/bin/windows/Rtools/>
- Run the installer, `Rtools.exe`, keeping the default settings.
