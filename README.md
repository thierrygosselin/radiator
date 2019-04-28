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

[Overview](https://thierrygosselin.github.io/radiator/articles/getting_started.html)
------------------------------------------------------------------------------------

New features
------------

Change log, version, new features and bug history lives in the [NEWS.md
file](https://thierrygosselin.github.io/radiator/news/index.html)

Life cycle
----------

**DArT users**:

-   `filter_dart`: is now deprecated. Please use `filter_rad`.
-   `tidy_dart` and `tidy_silico_dart`: are now deprecated. Please use
    `read_dart` for all the 4 DArT files recognized by radiator.

**Missing data: visualization and imputations**

Visualizing missing data and it’s imputations requires special attention
that fall outside the scope of **radiator**. Inside my package called
[grur](https://github.com/thierrygosselin/grur), users can **visualize
patterns of missingness** associated with different variables (lanes,
chips, sequencers, populations, sample sites, reads/samples,
homozygosity, etc). Several **Map-independent imputations** of missing
genotypes are available: **Random Forests** (on-the-fly-imputations or
predictive modeling), **Extreme Gradient Tree Boosting**, Strawman
imputations (\~ max/mean/mode: the most frequently observed, non-missing
genotypes is used). Imputations can be conducted **overall samples** or
**by populations/strata/grouping**. `radiator::genomic_converter` is
integrated with the imputation function of **grur**.

Prerequisite - Suggestions - Troubleshooting
--------------------------------------------

-   **Parallel computing**: follow the steps in this [notebook
    vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0)
    to install the packages with OpenMP-enabled compiler and conduct
    imputations in parallel.
-   [Installation
    problems.](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0)
-   **Windows users**: Install
    [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
-   The R GUI is unstable with functions using parallel ([more
    info](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html)),
    so I recommend using
    [RStudio](https://www.rstudio.com/products/rstudio/download/) for a
    better experience.
-   Using my R Notebook: use the option to [run chunks of codes in
    console, not
    inline](https://bookdown.org/yihui/rmarkdown/notebook.html#fig:notebook-console).

Vignettes, R Notebooks and examples
-----------------------------------

**Vignettes (in development, check periodically for updates):**

-   Vignettes with real data for example in the form of R Notebooks take
    too much space to be included in package, without CRAN complaining.
    Consequently, vignettes are gradually being excluded from the
    package and distributed separately, follow the links below.
-   **installation problems** [notebook
    vignette](https://www.dropbox.com/s/1kz59xpolb5y52m/rad_genomics_computer_setup.nb.html?dl=0)
-   **parallel computing during imputations** [notebook
    vignette](https://www.dropbox.com/s/1kz59xpolb5y52m/rad_genomics_computer_setup.nb.html?dl=0)
-   **vcf2dadi**
    [Rmd](https://www.dropbox.com/s/bl0mv6kavz97ibz/vignette_vcf2dadi.Rmd?dl=0)
    or
    [html](https://www.dropbox.com/s/qo0ujxmye7g7ora/vignette_vcf2dadi.html?dl=0)

**R Notebooks:**

-   Missing data visualization and analysis [(html
    vignette)](https://www.dropbox.com/s/btw1jos6yfck407/vignette_missing_data_analysis.nb.html?dl=0)
    and
    [(Rmd)](https://www.dropbox.com/s/tjjld6jczefyrj2/vignette_missing_data_analysis.Rmd?dl=0)

Citation:
---------

To get the citation, inside R:

``` r
citation("radiator")
```

Roadmap of future developments:
-------------------------------

-   Updated filters: more efficient, interactive and visualization
    included: *in progress*.
-   Workflow tutorial that links functions and points to specific
    vignettes to further explore some problems: *in progress*
-   Use Shiny and ggvis (when subplots and/or facets becomes available
    for ggvis).
-   Until publication **radiator** will change rapidly, stay updated
    with releases and contribute with bug reports.
-   Suggestions ?

Contributions:
--------------

This package has been developed in the open, and it wouldn’t be nearly
as good without your contributions. There are a number of ways you can
help me make this package even better:

-   If you don’t understand something, please let me know.
-   Your feedback on what is confusing or hard to understand is
    valuable.
-   If you spot a typo, feel free to edit the underlying page and send a
    pull request.

New to pull request on github ? The process is very easy:

-   Click the edit this page on the sidebar.
-   Make the changes using github’s in-page editor and save.
-   Submit a pull request and include a brief description of your
    changes.
-   “Fixing typos” is perfectly adequate.
