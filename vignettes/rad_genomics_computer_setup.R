## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----eval=FALSE---------------------------------------------------------------
#  if (!require("devtools")) install.packages("devtools") # to install
#  install.packages("tidyverse")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("gsl")

## ----eval=FALSE---------------------------------------------------------------
#  devtools::install_github("thierrygosselin/grur")
#  devtools::install_github("thierrygosselin/assigner")

## ----eval=FALSE---------------------------------------------------------------
#  usethis::edit_r_makevars()

## ----eval=FALSE---------------------------------------------------------------
#  #Warning: cannot remove prior installation of package ‘stringi’

## ----eval=FALSE---------------------------------------------------------------
#  file.exists("~/.Renviron")

## ---- eval=FALSE--------------------------------------------------------------
#  Sys.getenv("LD_LIBRARY_PATH")
#  # [1]""

## ---- eval=FALSE--------------------------------------------------------------
#  # in R:
#  Sys.setenv(LD_LIBRARY_PATH="/usr/local/lib/")
#  # For Linux you could use: /usr/local/lib/:/usr/lib64

## ---- eval=FALSE--------------------------------------------------------------
#  Sys.setenv(LD_LIBRARY_PATH="/usr/local/lib64/R/lib:/lib:/usr/local/lib64:/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.222.b10-0.amzn2.0.1.x86_64/jre/lib/amd64/server:/usr/local/lib/:/usr/lib64")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("rgeos", repos="http://R-Forge.R-project.org", type="source")
#  install.packages("rgdal", repos="http://R-Forge.R-project.org", type="source")
#  devtools::install_github("r-spatial/sf", configure.args = "--with-proj-lib=/usr/local")
#  install.packages("adegenet")

## ---- eval=FALSE--------------------------------------------------------------
#  Error: C stack usage  7971092 is too close to the limit

## ----eval=FALSE---------------------------------------------------------------
#  usethis::edit_r_makevars()

## ----eval=FALSE---------------------------------------------------------------
#  #fstcore fst data.table
#  CC=/usr/local/bin/gcc -fopenmp
#  CXX=/usr/local/bin/g++ -fopenmp
#  CXX11=/usr/local/bin/g++ -fopenmp
#  CXX14=/usr/local/bin/g++ -fopenmp
#  CXX17=/usr/local/bin/g++ -fopenmp
#  CXX1X=/usr/local/bin/g++ -fopenmp
#  CXX98=/usr/local/bin/g++ -fopenmp

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("fstcore", type = "source")
#  install.packages("fst", type = "source")
#  install.packages("data.table", type = "source")

## ----eval=FALSE---------------------------------------------------------------
#  usethis::edit_r_makevars()

## ----eval=FALSE---------------------------------------------------------------
#  usethis::edit_r_makevars()

## ----eval=FALSE---------------------------------------------------------------
#  require(xgboost)
#  x <-  matrix(rnorm(100*10000), 10000, 100)
#  y <-  x %*% rnorm(100) + rnorm(1000)
#  
#  system.time({bst = xgboost(data = x, label = y, nthread = 1, nround = 100, verbose = FALSE)})
#  system.time({bst = xgboost(data = x, label = y, nthread = 4, nround = 100, verbose = FALSE)})

## ----eval=FALSE---------------------------------------------------------------
#  usethis::edit_r_makevars()

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("ranger")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("missRanger")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("miceRanger")

## ----eval=FALSE---------------------------------------------------------------
#  BiocManager::install("pcaMethods")

## ----terminal_folder_shortcut.png, echo=FALSE---------------------------------
knitr::include_graphics("terminal_folder_shortcut.png")

## ----copy_path.png, echo=FALSE------------------------------------------------
knitr::include_graphics("copy_path.png")

