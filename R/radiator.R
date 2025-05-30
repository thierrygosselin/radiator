# radiator package startup message ---------------------------------------------
# .onAttach <- function(libname, pkgname) {
#   radiator.version <- utils::packageDescription("radiator", fields = "Version")
#   radiator.build <- utils::packageDescription("radiator", fields = "Built")
#   startup.message <- stringi::stri_join(
#     "******************************* IMPORTANT NOTICE *******************************\n",
#     "radiator v.", radiator.version, " was modified heavily.\n",
#     "Read functions documentation and available vignettes.\n\n",
#     "For reproducibility:\n",
#     "    radiator version: ", radiator.version,"\n",
#     "    radiator build date: ", radiator.build,"\n",
#     "    Keep zenodo DOI.\n",
#     "********************************************************************************",
#     sep = "")
#   packageStartupMessage(startup.message)
# }

# radiator common arguments ----------------------------------------------------
#' @name radiator_common_arguments
#' @title radiator common arguments
#' @description radiator common arguments
#' @rdname radiator_common_arguments
#' @export
#' @keywords internal

#' @param interactive.filter (optional, logical) Do you want the filtering session to
#' be interactive. Figures of distribution are shown before asking for filtering
#' thresholds.
#' Default: \code{interactive.filter = TRUE}.

#' @param gds (2 options) A Genomic Data Structure (GDS) file or object
#'
#' \emph{How to get GDS ?}
#' Look into \code{\link{tidy_genomic_data}},
#' \code{\link{read_vcf}},
#' \code{\link{tidy_vcf}} or \code{\link{write_gds}}.


#' @param data (4 options) A file or object generated by radiator:
#' \itemize{
#' \item tidy data
#' \item Genomic Data Structure (GDS)
#' }
#'
#' \emph{How to get GDS and tidy data ?}
#' Look into \code{\link{tidy_genomic_data}},
#' \code{\link{read_vcf}} or
#' \code{\link{tidy_vcf}}.


#' @param verbose (optional, logical) When \code{verbose = TRUE}
#' the function is a little more chatty during execution.
#' Default: \code{verbose = TRUE}.

#' @param parallel.core (optional) The number of core used for parallel
#' execution during import.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @param random.seed (integer, optional) For reproducibility, set an integer
#' that will be used inside the function that requires randomness. With default,
#' a random number is generated and printed in the appropriate output.
#' Default: \code{random.seed = NULL}.



#' @param ... (optional) Advance mode that allows to pass further arguments
#' for fine-tuning the function. Also used for legacy arguments (see details or
#' special section)


# @inheritParams radiator_common_arguments


radiator_common_arguments <- function(
  interactive.filter = TRUE,
  gds,
  data,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  random.seed = NULL,
  ...) {
  data = NULL
  parallel.core <- NULL
  verbose <- NULL
  random.seed <- NULL
}#End radiator_common_arguments


# radiator_question ---------------------------------------------------------
#' @title radiator_question
#' @description Ask to enter a word or number
#' @rdname radiator_question
#' @keywords internal
#' @export
radiator_question <- function(x, answer.opt = NULL, minmax = NULL) {
  # Note to myself: tryCatch might be simpler here... investigate
  message(x)
  question <- function(x, answer.opt) {

    if (!is.null(answer.opt)) {
      answer.type <- (unique(class(answer.opt)))
      if (answer.type == "character") {
        x <-   match.arg(arg = readLines(n = 1), choices = answer.opt)
      }
      if (answer.type == "numeric") {
        x <-   match.arg(arg = as.numeric(readLines(n = 1)), choices = answer.opt)
      }
      if (answer.type == "integer") {
        x <-   match.arg(arg = as.integer(readLines(n = 1)), choices = answer.opt)
      }
    }
    answer.type <- NULL
    return(x)
  }
  safe_question <- purrr::safely(.f = question, otherwise = FALSE)
  answers.ok <- FALSE
  while (!answers.ok) {

    if (!is.null(minmax)) {
      # answer <- as.numeric(readLines(n = 1))
      answer <- readLines(n = 1)

      check.answer <- stringi::stri_detect_regex(str = answer, pattern = "[0-9]")
      if (check.answer) {
        answer <- as.numeric(answer)
        good.value <- (answer >= minmax[1] & answer <= minmax[2])
        if (!good.value) {
          answers.ok <- FALSE
          message("Please try again: ")
        } else {
          answers.ok <- TRUE
        }
        good.value <- NULL
      } else {
        answers.ok <- FALSE
        message("Please try again: ")
      }
    } else {
      answer <- safe_question(x, answer.opt = answer.opt)
      if (is.null(answer$error)) {
        answer <- answer$result
        answers.ok <- TRUE
      } else {
        answers.ok <- FALSE
        message("Please try again, options are: ",
                stringi::stri_join(answer.opt, collapse = " or "))
      }
    }
  }
  return(answer)
}#End radiator_question


# radiator_parameters-------------------------------------------------------------
#' @title radiator_parameters
#' @description Generate or update a filters parameters file and object.
#' Used internally in radiator, not usefull outside the package.
#' @rdname radiator_parameters
#' @export
#' @keywords internal

# Note to myself: might be able to increase timing here by reading a whitelist
# instead of markers.meta for gds file...
# Then figure out what to do with individuals and strata...

radiator_parameters <- function(
  generate = FALSE,
  initiate = FALSE,
  update = TRUE,
  parameter.obj = NULL,
  data = NULL,
  filter.name = "",
  param.name = "",
  values = paste(NULL, NULL, sep = " / "),
  units = "individuals / strata / chrom / locus / markers",
  comments = "",
  path.folder = NULL,
  file.date = NULL,
  internal = FALSE,
  verbose = TRUE
) {
  if (internal && !verbose) return(NULL)
  res <- list()# initiate the list to store the results
  if (is.null(file.date)) file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # check for existing file
  if (is.null(path.folder)) path.folder <- getwd()
  if (!is.null(parameter.obj) && generate && !initiate) {
    generate <- initiate <- update <- FALSE
    res <- parameter.obj
  }
  if (!is.null(parameter.obj) && generate && initiate) generate <- FALSE
  if (is.null(parameter.obj) && update) rlang::abort("parameter.obj = NULL not accepted")
  if (internal) verbose <- FALSE

  # GENERATE filters parameters file
  if (generate) {
    filters.parameters.name <- stringi::stri_join("filters_parameters_", file.date, ".tsv")
    parameter.obj$filters.parameters.path <- res$filters.parameters.path <- file.path(path.folder, filters.parameters.name)
    res$filters.parameters <- tibble::tibble(
      FILTERS = as.character(),
      PARAMETERS = as.character(),
      VALUES = as.character(),
      BEFORE = as.character(),
      AFTER = as.character(),
      BLACKLIST = as.integer(),
      UNITS = as.character(),
      COMMENTS = as.character())

    write_rad(
      data = res$filters.parameters,
      filename = res$filters.parameters.path,
      tsv = TRUE,
      internal = internal,
      append = FALSE,
      write.message = NULL,
      verbose = verbose)
    if (verbose) message("Filters parameters file generated: ", filters.parameters.name)
  }#End generate

  # INITIATE filters parameters file
  if (initiate) {
    if (is.null(data)) rlang::abort("GDS or tidy data object required")
    res$info <- parameter.obj$info <- data_info(data)
    res$filters.parameters.path <- parameter.obj$filters.parameters.path
  }#End initiate

  # UPDATE filters parameters file
  if (update) {
    if (is.null(data)) rlang::abort("GDS or tidy data object required")
    info <- parameter.obj$info
    info.new <- data_info(data) # updating parameters

    res$filters.parameters <- tibble::tibble(
      FILTERS = filter.name,
      PARAMETERS = param.name,
      VALUES = if (!is.null(values)) {
        values
      } else {
        "not filtering"
      },
      BEFORE = paste(info$n.ind, info$n.pop, info$n.chrom, info$n.locus, info$n.snp, sep = " / "),
      AFTER = paste(info.new$n.ind, info.new$n.pop, info.new$n.chrom, info.new$n.locus, info.new$n.snp, sep = " / "),
      BLACKLIST = paste(info$n.ind - info.new$n.ind, info$n.pop - info.new$n.pop, info$n.chrom - info.new$n.chrom, info$n.locus - info.new$n.locus, info$n.snp - info.new$n.snp, sep = " / "),
      UNITS = units,
      COMMENTS = comments
    )

    write_rad(
      data = res$filters.parameters,
      filename = parameter.obj$filters.parameters.path,
      tsv = TRUE,
      internal = internal,
      append = TRUE, col.names = FALSE,
      write.message = NULL,
      verbose = verbose)
    # update info
    res$info <- info.new
    res$filters.parameters.path <- parameter.obj$filters.parameters.path
  }#End update

  # messages
  # if (initiate && update) {
  #   if (verbose) message("Filters parameters file: initiated and updated")
  # }

  # if (initiate && !update) {
  #   if (verbose) message("Filters parameters file: initiated")
  # }

  # if (!initiate && update) {
  #   if (verbose) message("Filters parameters file: updated")
  # }

  return(res)
}#End radiator_parameters

# data.info -------------------------------------------------------------
#' @title data_info
#' @description function generate tidy data main info
#' @rdname data_info
#' @keywords internal
#' @export
data_info <- function(x, verbose = FALSE) {
  res <- list()

  data.type <- class(x)[1]

  if (data.type == "tbl_df") {
    if (rlang::has_name(x, "POP_ID") || rlang::has_name(x, "STRATA")) {

      if (rlang::has_name(x, "POP_ID")) {
        res$n.pop <- length(unique(x$POP_ID))
      } else {
        res$n.pop <- length(unique(x$STRATA))
      }
    } else {
      res$n.pop <- NA_integer_
    }

    if (rlang::has_name(x, "INDIVIDUALS")) {
      res$n.ind <- length(unique(x$INDIVIDUALS))
    } else {
      res$n.ind <- NA_integer_
    }


    if (rlang::has_name(x, "MARKERS")) {
      res$n.snp <- length(unique(x$MARKERS))
    } else {
      res$n.snp <- NA_integer_
    }

    if (rlang::has_name(x, "LOCUS")) {
      res$n.locus <- length(unique(x$LOCUS))
    } else {
      res$n.locus <- NA_integer_
    }

    if (rlang::has_name(x, "CHROM")) {
      res$n.chrom <- length(unique(x$CHROM))
    } else {
      res$n.chrom <- NA_integer_
    }
  } else {
    i <- extract_markers_metadata(
      gds = x,
      markers.meta.select = c("CHROM", "LOCUS", "MARKERS"),
      whitelist = TRUE)
    res$n.chrom <- length(unique(i$CHROM))
    res$n.locus <- length(unique(i$LOCUS))
    res$n.snp <- length(unique(i$MARKERS))
    i <- extract_individuals_metadata(
      gds = x,
      ind.field.select = c("STRATA", "INDIVIDUALS"),
      whitelist = TRUE)
    res$n.pop <- length(unique(i$STRATA))
    res$n.ind <- length(unique(i$INDIVIDUALS))
    # res$n.chrom <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/markers.meta/CHROM", silent = TRUE))))
    # res$n.locus <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/markers.meta/LOCUS", silent = TRUE))))
    # res$n.snp <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/markers.meta/MARKERS", silent = TRUE))))
    # res$n.pop <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/individuals/STRATA", silent = TRUE))))
    # res$n.ind <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/individuals/INDIVIDUALS", silent = TRUE))))
    res[is.null(res)] <- NA_integer_
  }

  if (verbose) {
    message("Number of chrom: ", res$n.chrom)
    message("Number of locus: ", res$n.locus)
    message("Number of SNPs: ", res$n.snp)
    message("Number of strata: ", res$n.pop)
    message("Number of individuals: ", res$n.ind)
  }
  return(res)
}

# radiator_results_message ------------------------------------------------------------
#' @title radiator_results_message
#' @description Message printed at the end of most radiator functions
#' @keywords internal
#' @export
radiator_results_message <- function(
  rad.message = NULL,
  filters.parameters,
  internal = FALSE,
  verbose = TRUE
) {
  if (!internal) {
    if (verbose) cat("################################### RESULTS ####################################\n")
    if (!is.null(rad.message)) message(rad.message)
    message("Number of individuals / strata / chrom / locus / SNP:")
    if (verbose) message("    Before: ", filters.parameters$filters.parameters$BEFORE)
    message("    Blacklisted: ", filters.parameters$filters.parameters$BLACKLIST)
    if (verbose) message("    After: ", filters.parameters$filters.parameters$AFTER)
  }
}#End radiator_results_message

# radiator_folder--------------------------------------------------------------------
#' @title radiator_folder
#' @description Generate the rad folders
#' @param path.folder path of the folder
#' @param prefix_int Use an integer prefix padded left with 0.
#' Default: \code{prefix_int = TRUE}.
#' @inheritParams folder_short
#' @keywords internal
#' @export
#' @rdname radiator_folder
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

radiator_folder <- function(f, path.folder = NULL, prefix_int = TRUE) {
  if (is.null(path.folder)) path.folder <- getwd()
  if (prefix_int) {
    n.dir <- list.dirs(path = path.folder, full.names = FALSE)[-1]
    n.dir <- length(n.dir) - sum(stringi::stri_count_fixed(
      str = n.dir, pattern = "/")
    ) + 1L

    f <- stringi::stri_join(
      stringi::stri_pad_left(str = n.dir, width = 2, pad = 0),
      "_",
      f
    )
    n.dir <- NULL

    # old version for bk
    # f <- stringi::stri_join(stringi::stri_pad_left(
    #   str = length(list.dirs(path = path.folder, full.names = FALSE)[-1]) + 1L,
    #   width = 2,
    #   pad = 0
    # ), "_", f)
  }
  folder.prefix <- file.path(path.folder, f)
  return(folder.prefix)
}#End radiator_folder


# generate_squeleton_folders----------------------------------------------------
#' @title generate_squeleton_folders
#' @description Generate squeleton folders
#' @keywords internal
#' @export
generate_squeleton_folders <- function(
  fp = 0L,
  path.folder = NULL,
  interactive.filter = TRUE,
  ...
) {

  # test
  # fp = 0L
  # file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  # interactive.filter = TRUE

  if (is.null(path.folder)) path.folder <- getwd()
  folders.labels <- c(
    "filter_dart_reproducibility",
    "filter_individuals", "filter_individuals", "filter_individuals",
    "filter_common_markers",
    "filter_ma",
    "filter_coverage",
    "filter_genotyping",
    "filter_snp_position_read",
    "filter_snp_number",
    "filter_ld", "filter_ld",
    "detect_mixed_genomes",
    "detect_duplicate_genomes",
    "filter_hwe")

  if (!interactive.filter) {
    get.filters <- ls(envir = as.environment(1))
    need <- c(
      "filter.reproducibility",
      "filter.individuals.missing",
      "filter.individuals.heterozygosity",
      "filter.individuals.coverage.total",
      "filter.common.markers",
      "filter.ma",
      "filter.coverage",
      "filter.genotyping",
      "filter.snp.position.read",
      "filter.snp.number",
      "filter.short.ld",
      "filter.long.ld",
      "detect.mixed.genomes",
      "detect.duplicate.genomes",
      "filter.hwe")
    folders <- purrr::keep(.x = get.filters, .p = get.filters %in% need)
    wanted_filters <- function(x) {
      !is.null(rlang::eval_tidy(rlang::parse_expr(x)))
    }
    folders <- purrr::keep(.x = folders, .p = wanted_filters)
    folders <- factor(
      x = folders,
      levels = need,
      labels = folders.labels,
      ordered = TRUE
    ) %>%
      droplevels(.) %>%
      unique %>%
      sort %>%
      as.character
  } else {
    folders <- unique(folders.labels)
  }

  folders <- c("radiator", folders)

  res <- list()
  fp.loop <- fp
  temp <- NULL
  for (f in folders) {
    # message("Processing: ", f)
    temp <- folder_prefix(
      prefix_int = fp.loop,
      prefix.name = f,
      path.folder = path.folder)
    res[[f]] <- temp$folder.prefix
    fp.loop <- temp$prefix_int
  }
  return(res)
}#End generate_squeleton_folders
# generate_filename-------------------------------------------------------------
#' @title Filename radiator
#' @description Generate a filename object
#' @name generate_filename
#' @rdname generate_filename
#' @keywords internal
#' @export
generate_filename <- function(
  name.shortcut = NULL,
  path.folder = NULL,
  date = TRUE,
  extension = c(
    "tsv", "gds.rad", "rad", "gds", "gen", "dat",
    "genind", "genlight", "gtypes", "vcf", "colony",
    "bayescan", "gsisim", "hierfstat", "hzar", "ldna",
    "pcadapt", "related", "stockr", "structure", "arlequin",
    "arrow.parquet"
  )
) {

  if (is.null(path.folder)) path.folder <- getwd()

  # date and time-
  if (date) {
    file.date <- stringi::stri_join("_", format(Sys.time(), "%Y%m%d@%H%M"))
  } else {
    file.date <- ""
  }

  # path.folder
  if (!dir.exists(path.folder)) dir.create(path.folder)

  # Extension
  want <- c("tsv", "gds.rad", "rad", "gds", "gen", "dat", "genind", "genlight", "gtypes",
            "vcf", "colony", "bayescan", "gsisim", "hierfstat", "hzar", "ldna",
            "pcadapt", "plink", "related", "stockr", "structure", "arlequin",
            "arrow.parquet")
  extension <- match.arg(extension, want)

  # note to myself: currently excluded output :
  # "fineradstructure", "maverick", "plink", "betadiv"


  # with same extension
  # extension <- "tsv"
  if (extension %in% c("tsv", "gds.rad", "rad", "gds", "vcf", "colony", "ldna",
                       "arrow.parquet")) {
    extension <- stringi::stri_join(file.date, ".", extension)
  }


  # Radiator saveRDS
  # extension <- "genind"
  if (extension %in% c("genind", "genlight", "gtypes", "stockr")) {
    extension <- stringi::stri_join("_", extension, file.date, ".RData")
  }

  # Radiator tsv
  if (extension %in% c("tsv")) {
    extension <- stringi::stri_join("_", extension, file.date, ".tsv")
  }

  # Radiator txt
  if (extension %in% c("bayescan", "pcadapt", "related")) {
    extension <- stringi::stri_join("_", extension, file.date, ".txt")
  }

  # Radiator csv
  if (extension %in% c("hzar", "arlequin")) {
    extension <- stringi::stri_join("_", extension, file.date, ".csv")
  }

  # custom
  if (extension == "gen") extension <- stringi::stri_join("_genepop", file.date, ".gen")
  if (extension == "dat") extension <- stringi::stri_join("_fstat", file.date, ".dat")
  if (extension == "hierfstat") extension <- stringi::stri_join("_hierfstat", file.date, ".dat")
  if (extension == "structure") extension <- stringi::stri_join("_structure", file.date, ".str")

  # Filename
  if (is.null(name.shortcut)) {
    filename <- stringi::stri_join("radiator", extension)
  } else {
    filename.problem <- file.exists(stringi::stri_join(name.shortcut, extension))
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_radiator", extension)
    } else {
      filename <- stringi::stri_join(name.shortcut, extension)
    }
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join("duplicated_", filename)
    }
  }


  # Include path.folder in returned object
  return(res = list(filename.short = filename, filename = file.path(path.folder, filename)))
}#End generate_filename


# generate_folder---------------------------------------------------------------

#' @title generate_folder
#' @description Generate a folder based on ...
#' @name generate_folder
#' @param rad.folder Name of the rad folder
#' @param internal (optional, logical) Is the function internal or not
#' @param append.date Include the date and time with the folder.
#' Default: \code{append.date = TRUE}.
#' @param file.date The file date included in as argument/value or with
#' default \code{file.date = NULL}, generated by the fucntion.
#' @inheritParams radiator_folder
#' @inheritParams radiator_common_arguments
#' @inheritParams folder_short
#' @keywords internal
#' @export
#' @rdname generate_folder
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

generate_folder <- function(
  f,
  rad.folder = NULL,
  internal = FALSE,
  append.date = TRUE,
  file.date = NULL,
  prefix_int = TRUE,
  verbose = FALSE
) {

  if (internal) {
    rad.folder <- NULL
  }
  if (!is.null(rad.folder)) f <- radiator_folder(rad.folder, f, prefix_int = prefix_int)


  f.temp <- f
  if (is.null(file.date)) {
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")# Date and time
  }

  if (is.null(f)) {
    f <- getwd()
  } else {
    #working directory in the path?
    wd.present <- TRUE %in% unique(stringi::stri_detect_fixed(str = f, pattern = c(getwd(), paste0(getwd(), "/"))))
    date.present <- TRUE %in% unique(stringi::stri_detect_fixed(str = f, pattern = "@"))
    if (!date.present && append.date) f <- stringi::stri_join(f, file.date, sep = "_")
    if (!wd.present) f <- file.path(getwd(), f)
    if (verbose && !identical(f.temp, f)) message("Folder created: ", folder_short(f))
  }
  if (!dir.exists(f)) dir.create(f)
  return(f)
}#End generate_folder

# folder_short------------------------------------------------------------------
#' @title folder_short
#' @description Extract the short name of the folder generated by radiator
#' @name folder_short
#' @param f Folder name
#' @rdname folder_short
#' @keywords internal
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
folder_short <- function(f) {
  # remove / if found last
  if (stringi::stri_sub(str = f, from = -1, length = 1) == "/") {
    f <- stringi::stri_replace_last_regex(
      str = f,
      pattern = "[/$]",
      replacement = "")
  }

  # detect the presence of /
  if (stringi::stri_detect_fixed(str = f, pattern = "/")) {
    f %<>% stringi::stri_sub(
      str = .,
      from = stringi::stri_locate_last_fixed(
        str = .,
        pattern = "/")[2] + 1,
      to = stringi::stri_length(str = .)
    )
  }
  return(f)
}#End folder_short

# folder_prefix-----------------------------------------------------------------
#' @title folder_prefix
#' @description Generate a seq and folder prefix
#' @name folder_prefix
#' @rdname folder_prefix
#' @keywords internal
#' @export
folder_prefix <- function(
  prefix_int = NULL,
  prefix.name = NULL,
  path.folder = NULL
) {
  if (is.null(path.folder)) {
    path.folder <- getwd()
  } else {
    if (stringi::stri_sub(str = path.folder, from = -1, length = 1) == "/") {
      path.folder <- stringi::stri_replace_last_regex(
        str = path.folder,
        pattern = "[/$]",
        replacement = "")
    }
  }

  if (is.null(prefix_int)) {
    prefix_int <- 0L
  } else {
    if (is.list(prefix_int)) {
      prefix_int <- as.integer(prefix_int$prefix_int) + 1L
    } else {
      prefix_int <- as.integer(prefix_int) + 1L
    }
  }

  if (is.null(prefix.name)) {
    folder.prefix <- stringi::stri_join(
      stringi::stri_pad_left(
        str = prefix_int, width = 2, pad = 0
      ), "_"
    )
  } else {
    folder.prefix <- stringi::stri_join(
      stringi::stri_pad_left(
        str = prefix_int, width = 2, pad = 0
      ),
      prefix.name,
      sep = "_"
    )
  }
  folder.prefix <- file.path(path.folder, folder.prefix)
  res = list(prefix_int = prefix_int, folder.prefix = folder.prefix)
}#End folder_prefix


# radiator_snakecase------------------------------------------------------------
#' @title radiator_snakecase
#' @description Transform CamelCase to snake_cases
#' @name radiator_snakecase
#' @rdname radiator_snakecase
#' @keywords internal
#' @export
radiator_snakecase <- function(x) {
  x <- base::gsub(pattern = "([A-Za-z])([A-Z])([a-z])", replacement = "\\1_\\2\\3", x = x) %>%
    base::gsub(pattern = ".", replacement = "_", x =  ., fixed = TRUE) %>%
    base::gsub(pattern = "([a-z])([A-Z])", replacement = "\\1_\\2", x = .) %>%
    stringi::stri_trans_toupper(str = .)
  return(x)
}#End radiator_snakecase


# radiator_packages_dep---------------------------------------------------------
#' @title radiator_packages_dep
#' @description Verify required packages
#' @rdname radiator_packages_dep
#' @keywords internal
#' @export
radiator_packages_dep <- function(package, cran = TRUE, bioc = FALSE) {
  if (cran) bioc <- FALSE
  if (bioc) cran <- FALSE
  installer <- "devtools::install_github"
  if (cran) installer <- "install.packages"
  if (bioc) installer <- "BiocManager::install"
  how.to <- stringi::stri_join(installer, "('", package, "')")
  if (suppressPackageStartupMessages(!requireNamespace(package, quietly = TRUE))) {
    rlang::abort(
      paste0(paste0("Package required: ", package),
             paste0("\n       Install with: ", how.to))
    )
  }

}#End radiator_packages_dep

# radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)
# requireNamespace
# installed.packages

# radiator_function_header -----------------------------------------------------
#' @title radiator_function_header
#' @description Generate function header
#' @rdname radiator_function_header
#' @keywords internal
#' @export
radiator_function_header <- function(f.name = NULL, start = TRUE, verbose = TRUE) {
  if (is.null(f.name)) invisible(NULL)
  if (start) {
    if (verbose) {
      # cat("################################################################################\n")
      cat(paste0(stringi::stri_pad_both(str = "", width = 80L, pad = "#"), "\n"))
      cat(paste0(stringi::stri_pad_both(str = paste0(" radiator::", f.name, " "), width = 80L, pad = "#"), "\n"))
      cat(paste0(stringi::stri_pad_both(str = "", width = 80L, pad = "#"), "\n"))
      # cat("################################################################################\n")
    }
  } else {
    if (verbose) {
      cat(paste0(stringi::stri_pad_both(str = paste0(" completed ", f.name, " "), width = 80L, pad = "#"), "\n"))
    }
  }
}# End radiator_function_header

# radiator_clock ---------------------------------------------------------------
#' @title radiator_tic
#' @description radiator tictoc function
#' @rdname radiator_tic
#' @keywords internal
#' @export
radiator_tic <- function(timing = proc.time()) {
    invisible(timing)
}# End radiator_tic

#' @title radiator_toc
#' @description radiator tictoc function
#' @rdname radiator_toc
#' @keywords internal
#' @export
radiator_toc <- function(
  timing,
  end.message = "Computation time, overall:",
  verbose = TRUE
) {
  if (verbose) message("\n", end.message, " ", round((proc.time() - timing)[[3]]), " sec")
}# End radiator_toc


