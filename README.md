<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://tidyverse.org/lifecycle/#maturing)
[![Travis-CI Build
Status](https://travis-ci.org/thierrygosselin/radiator.svg?branch=master)](https://travis-ci.org/thierrygosselin/radiator)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/radiator?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/radiator)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/radiator)](http://cran.r-project.org/package=radiator)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/radiator.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/radiator)
[![packageversion](https://img.shields.io/badge/Package%20version-1.0.0-orange.svg)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2019--04--28-brightgreen.svg)](/commits/master)

------------------------------------------------------------------------

radiator: an R package for RADseq Data Exploration, Manipulation and Visualization
==================================================================================

This is the development page of the **radiator**, if you want to help,
see [contributions
section](https://github.com/thierrygosselin/radiator#contributions)

Most genomic analysis look for patterns and trends with various
statistics. Bias, noise and outliers can have bounded influence on
estimators and interfere with polymorphism discovery. Avoid bad data
exploration and control the impact of filters on your downstream genetic
analysis. Use radiator to: import, explore, manipulate, visualize,
filter, impute and export your GBS/RADseq data.

**radiator** is designed and optimized for fast computations using
Genomic Data Structure [GDS](http://zhengxwen.github.io/gdsfmt) file
format and data science packages in
[tiverse](https://www.tidyverse.org). **radiator** handles VCF files
with millions of SNPs and files of several GB.

Installation
------------

To try out the dev version of **radiator**, copy/paste the code below:

``` r
if (!require("pak")) install.packages("pak")
pak::pkg_install("thierrygosselin/radiator")
library(radiator)
```

To minimize dependencies, just the basic required packages are installed
with the command above. If you want the full suits of functions and
don’t want to be preoccupied, run:

``` r
radiator::radiator_pkg_install() # that's it. It will update, when necessary, radiator.
```

[Computer setup and troubleshooting
vignette](http://thierrygosselin.github.io/assigner/articles/rad_genomics_computer_setup.html)

[Learning radiator](https://thierrygosselin.github.io/radiator/articles/getting_started.html)
---------------------------------------------------------------------------------------------

[Overview of features](https://thierrygosselin.github.io/radiator/articles/getting_started.html#overview)
---------------------------------------------------------------------------------------------------------

[New features and changelog](https://thierrygosselin.github.io/radiator/news/index.html)
----------------------------------------------------------------------------------------

Change log, version, new features and bug history lives in the [NEWS.md
file](https://thierrygosselin.github.io/radiator/news/index.html)

Life cycle
----------

Philosophy, major changes and deprecated functions/arguments.

[Vignettes](https://thierrygosselin.github.io/radiator/articles/index.html)
---------------------------------------------------------------------------

Citation:
---------

To get the citation, inside R:

``` r
citation("radiator")
```
