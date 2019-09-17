
<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://tidyverse.org/lifecycle/#maturing)
[![Travis-CI Build
Status](https://travis-ci.org/thierrygosselin/radiator.svg?branch=master)](https://travis-ci.org/thierrygosselin/radiator)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/radiator.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/radiator)
[![packageversion](https://img.shields.io/badge/Package%20version-1.1.2-orange.svg)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2019--09--17-brightgreen.svg)](/commits/master)
<!-- badges: end -->

# radiator: an R package for RADseq Data Exploration, Manipulation and Visualization

Most genomic analysis look for patterns and trends with various
statistics. Bias, noise and outliers can have bounded influence on
estimators and interfere with polymorphism discovery. Avoid bad data
exploration and control the impact of filters on your downstream genetic
analysis. Use radiator to: import, explore, manipulate, visualize,
filter, impute and export your GBS/RADseq data.

**radiator** is designed and optimized for fast computations using
Genomic Data Structure [GDS](https://github.com/zhengxwen/gdsfmt) file
format and data science packages in
[tidyverse](https://www.tidyverse.org). **radiator** handles VCF files
with millions of SNPs and files of several GB.

## Installation

To try out the dev version of **radiator**, copy/paste the code below:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("thierrygosselin/radiator")
library(radiator)
```

Web site with additional info:
<https://thierrygosselin.github.io/radiator/>

  - [Computer setup and
    troubleshooting](https://thierrygosselin.github.io/radiator/articles/rad_genomics_computer_setup.html)
  - [Learning
    radiator](https://thierrygosselin.github.io/radiator/articles/get_started.html)
    and [Overview of
    features](https://thierrygosselin.github.io/radiator/articles/get_started.html#overview)
  - [Vignettes](https://thierrygosselin.github.io/radiator/articles/index.html)

## [Life cycle](https://thierrygosselin.github.io/radiator/articles/life_cycle.html)

radiator is maturing, but in order to make the package better, changes
are inevitable. Experimental functions will change, argument names will
change. Your codes and workflows might break from time to time until
radiator is stable. Consequently, depending on your tolerance to change,
radiator might not be for you.

  - Philosophy, major changes and deprecated functions/arguments are
    documented in life cycle section of functions.
  - The latest changes are documented
    ([here](https://thierrygosselin.github.io/radiator/articles/life_cycle.html))
    and in [changelog, versions, new features and bug
    history](https://thierrygosselin.github.io/radiator/news/index.html)
  - [issues](https://github.com/thierrygosselin/radiator/issues/new/choose)
    and
    [contributions](https://github.com/thierrygosselin/radiator/issues/new/choose)
