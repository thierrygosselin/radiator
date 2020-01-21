# write a LDna object from a tidy data frame

#' @name write_ldna
#' @title Write a LDna object from a tidy data frame
#' @description Write a \href{https://github.com/petrikemppainen/LDna}{LDna}
#' object from a biallelic tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.
#' \strong{The genotypes are biallelic.}


#' @param filename (optional) The file name of the LDna lower matrix file.
#' Radiator will append \code{.ldna.rds} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_datetime.ldna.rds} is chosen.
#' With default, \code{filename = NULL}, no file is generated, only an object in
#' the Global Environment.
#' To read the data back into R, use readRDS("filename.ldna.rds").

#' @inheritParams tidy_genomic_data

#' @param ... (optional) To pass further argument for fine-tuning the
#' function (see details).

#' @export
#' @rdname write_ldna
#' @references Kemppainen P, Knight CG, Sarma DK et al. (2015)
#' Linkage disequilibrium network analysis (LDna) gives a global view of
#' chromosomal inversions, local adaptation and geographic structure.
#' Molecular Ecology Resources, 15, 1031-1045.

#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' (2012) A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606

#' @details The function requires \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}
#' to prepare the data for LDna.
#'
#' To install SNPRelate:
#' install.packages("BiocManager")
#' BiocManager::install("SNPRelate")
#'
#' To install LDna:
#' devtools::install_github("petrikemppainen/LDna")
#'

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_ldna <- function(data,
                       filename = NULL,
                       parallel.core = parallel::detectCores() - 1,
                       ...
) {


  # testing
  # data <- unfiltered.data
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1

  # Check that snprelate is installed
  radiator_packages_dep(package = "SNPRelate", cran = FALSE, bioc = TRUE)

  # dotslist -------------------------------------------------------------------
  radiator.dots <- list(...)


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    write.ldna <- FALSE
    filename <- stringi::stri_join("radiator_", file.date, ".ldna")
  } else {
    write.ldna <- TRUE
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".ldna")
    } else {
      filename <- stringi::stri_join(filename, ".ldna")
    }
  }

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  markers <- unique(data$MARKERS)

  # Check if data is biallelic -------------------------------------------------
  biallelic <- radiator::detect_biallelic_markers(data = data)
  if (!biallelic) rlang::abort("LDna requires biallelic genotypes")

  # Generating SNPRelate data --------------------------------------------------
  message("Generating SNPRelate data")
  data <- radiator::write_snprelate(
    data = data,
    biallelic = TRUE,
    filename = filename,
    verbose = FALSE)
  message("SNPRelate GDS file generated: ", filename, ".gds")
  message("To close the connection use SNPRelate::snpgdsClose(filename)")


  # Compute LD -----------------------------------------------------------------
  message("Computing LD matrix...")
  LD <- NULL
  long.distance.ld <- SNPRelate::snpgdsLDMat(
    gdsobj = data,
    snp.id = NULL,
    sample.id = NULL,
    slide = -1,
    mat.trim = FALSE,
    method = "r", #composite and corr option are the same with 0, 1, 2 gt coding
    num.thread = parallel.core,
    verbose = TRUE) %$%
    LD

  # anyNA(long.distance.ld)
  # long.distance.ld.bk <- long.distance.ld
  # long.distance.ld <- long.distance.ld.bk

  # fill diagonal and upper matrix with NA
  message("Preparing the data for LDna")
  long.distance.ld[upper.tri(long.distance.ld, diag = TRUE)] <- rlang::na_dbl
  # long.distance.ld[is.nan(long.distance.ld)] <- NA # removes NaN
  long.distance.ld <- long.distance.ld^2 # R-square required by LDna

  # dim(long.distance.ld)
  colnames(long.distance.ld) <- rownames(long.distance.ld) <- markers
  markers <- NULL

  # write object ----------------------------------------------------------------
  if (write.ldna) {
    filename <- stringi::stri_join(filename, ".rds")
    message("Writing LDna file: ", filename)
    saveRDS(object = long.distance.ld, file = filename)
  }
  return(long.distance.ld)
} # End write_ldna
