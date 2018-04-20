---
title: "Problems"
output: github_document
---

## Problems during installation
  
Sometimes you'll get warnings while installing dependencies required for **stackr** or other R packages.
```r
Warning: cannot remove prior installation of package ‘stringi’
```

To solve this problem, delete manually the problematic package in the installation folder (on mac: `/Library/Frameworks/R.framework/Resources/library`) or in the `Terminal`:
```r
sudo rm -R /Library/Frameworks/R.framework/Resources/library/package_name
# Changing 'package_name' to the problematic package.
# Reinstall the package.
```

## Problem running functions

My packages are still in developpement, so it's very likely you'll find bugs,
before contacting me, try these tricks:

* try to reproduce the problem after turning off/on your computer
* a clean R/RStudio session works best
* re-install the problematic package
* update [R CRAN](https://cran.r-project.org)
* update [RStudio](https://www.rstudio.com/products/rstudio/download/)
* raise an issue on the package github page
    * give the output of `devtools::session_info()`
    * my help will be limited if I can't reproduce your problem
    * output the full command and highlight problem with a small dataset
    * if you prefer sharing your data directly with me,
    know that it will remain anonymous
