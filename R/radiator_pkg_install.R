# radiator dependencies install ------------------------------------------------
#' @name radiator_pkg_install
#' @title radiator packages install helper
#' @description This function help to install packages for the full suits of
#' radiator functions or the minimal install to lower packages dependencies installation.
#' It also re-install the latest version of radiator when available.
#' @rdname radiator_pkg_install
#' @param check (logical, optional) Check the packages installed and system info.
#' Default: \code{check = TRUE}
#' @param minimal.install (logical, optional) The minimal install does nothing,
#' and the function just check what's installed.
#' Default: \code{check = FALSE}, does the full installation of packages you might
#' need using radiator.
#' @export
radiator_pkg_install <- function(
  check = TRUE,
  minimal.install = FALSE
) {
  required.pkg = c(
    "devtools", "tidyverse",
    "adegenet", "BiocManager", "UpSetR",
    "ggtern", "HardyWeinberg",
    "knitr", "rmarkdown", "stringdist", "quantreg",
    "gdsfmt", "SeqArray", "SeqVarTools", "SNPRelate",
    "radiator"
  )

  # List the packages currently installed --------------------------------------
  list.pkg <- utils::installed.packages()[,"Package"]

  # Required package for the function to work out-of-the-box -------------------
  if (Sys.info()[['sysname']] == "Windows") {
    message("You'll need Rtools ... this is not an R package!")
  }
  if (!"devtools" %in% list.pkg) utils::install.packages("devtools")
  if (!"tidyverse" %in% list.pkg) utils::install.packages("tidyverse")


  # System info ----------------------------------------------------------------
  info <- devtools::session_info()
  message("R Version: ", info$platform$version)
  message("Computer OS: ", info$platform$os)
  message("Number of CPU available: ", parallel::detectCores())
  # message("Number of RAM available: ", )

  # Function -------------------------------------------------------------------
  check_installation <- function(required.pkg, list.pkg) {
    i <- required.pkg %in% list.pkg
    res <- tibble::tibble(PACKAGES = required.pkg, INSTALLED = i)
    return(res)
  }

  # Tibble of packages installation status--------------------------------------
  tpis <- purrr::map_dfr(.x = required.pkg, .f = check_installation, list.pkg = list.pkg)
  pi <- paste0(tpis$PACKAGES[tpis$INSTALLED], collapse = "\n")
  pm <- tpis$PACKAGES[!tpis$INSTALLED]
  names(pm) <- pm
  message("\n\nPackages installed:\n", pi)
  message("\n\nPackages missing:\n", paste0(pm, collapse = "\n"))

  # If missing and wants to install --------------------------------------------
  if (!minimal.install) {
    message("\n\nInstalling missing packages...")
    # CRAN packages---------------------------------------------------------------
    cran.p <- c("BiocManager", "UpSetR", "adegenet", "ggtern",
                "HardyWeinberg", "knitr", "rmarkdown", "stringdist", "quantreg")
    cran.p <- purrr::keep(.x = cran.p, cran.p %in% pm)
    purrr::walk(.x = cran.p, .f = utils::install.packages)

    # BiocManager packages -----------------------------------------------------
    if ("BiocManager" %in% pm) utils::install.packages("BiocManager")
    if (any(c("gdsfmt", "SeqArray", "SeqVarTools") %in% pm)) {
      BiocManager::install("SeqVarTools")
    }
    if ("SNPRelate" %in% pm) BiocManager::install("SNPRelate")

    # Github packages ----------------------------------------------------------
    devtools::install_github("thierrygosselin/radiator")

    # check and report a second time
    list.pkg <- utils::installed.packages()[,"Package"]
    tpis <- purrr::map_dfr(.x = required.pkg, .f = check_installation, list.pkg = list.pkg)
    pi <- paste0(tpis$PACKAGES[tpis$INSTALLED], collapse = "\n")
    pm <- tpis$PACKAGES[!tpis$INSTALLED]
    names(pm) <- pm

    if (length(pm) > 0) {
      message("\n\nPackages missing:\n", paste0(pm, collapse = "\n"))
      message("\nRe-run the function please")
    } else {
      message("\n\nYou're all good!")
    }
  }
  return(list(pkg.info = tpis, system.info = info))
}#End radiator_pkg_install
