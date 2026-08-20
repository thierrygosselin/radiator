# radiator_dots ---------------------------------------------------------------
#' @title radiator_dots
#' @description
#' Parse and manage `...` for radiator functions:
#'   - assign supported arguments ("keepers") into a target environment
#'   - fill in internal defaults for missing keepers
#'   - detect deprecated and unknown arguments
#'   - return a compact summary tibble for logging/debugging.
#'
#' @name radiator_dots
#' @rdname radiator_dots
#'
#' @param func.name Name of the calling function
#'   (used only in messages).
#'   Default: \code{as.list(sys.call())[[1]]}.
#'
#' @param fd (optional) The formal argument names of the calling function.
#'   Default: \code{rlang::fn_fmls_names()}.
#'
#' @param args.list (optional) Named list of current arguments from the calling
#'   environment. Default: \code{as.list(environment())}.
#'
#' @param dotslist Captured `...` arguments, usually from
#'   \code{rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)}.
#'
#' @param keepers Character vector of argument names that are valid inside
#'   \code{...} and should be assigned into the calling environment.
#'
#' @param deprecated Character vector of deprecated argument names.
#'
#' @param env Environment where keeper arguments and defaults should be
#'   assigned. Default: \code{parent.frame()}.
#'
#' @param assign (logical) If \code{TRUE}, assign keeper and default values
#'   into \code{env}. If \code{FALSE}, no assignment is performed and the
#'   function only returns the summary tibble.
#'   Default: \code{assign = TRUE}.
#'
#' @param verbose (logical) When \code{TRUE}, messages about function call
#'   arguments, `...` resolution, deprecated and unknown arguments are printed.
#'   Default: \code{verbose = TRUE}.
#'
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{ARGUMENTS}: argument name
#'     \item \code{VALUES}: character representation of the value or class
#'     \item \code{GROUPS}: one of \code{"fct.call.args"}, \code{"fct.call..."},
#'       \code{"default..."}, \code{"deprecated..."}, \code{"unknowned..."}.
#'   }
#'
#' The side-effect (when \code{assign = TRUE}) is to assign resolved keeper
#' arguments and their defaults into \code{env}.
#'
#' @keywords internal
#' @export
radiator_dots <- function(
    func.name  = as.list(sys.call())[[1]],
    fd         = NULL,
    args.list  = NULL,
    dotslist   = NULL,
    keepers    = c(
      "subsample.markers.stats", "force.stats", "id.stats", "subsample",
      "filter.reproducibility",
      "filter.individuals.missing",
      "filter.individuals.heterozygosity",
      "filter.individuals.coverage.total",
      "filter.individuals.coverage.median",
      "filter.individuals.coverage.iqr",
      "filter.common.markers", "filter.monomorphic",
      "filter.ma",
      "ma.stats",
      "filter.coverage", "dp",
      "filter.genotyping",
      "filter.snp.position.read",
      "filter.snp.number",
      "filter.short.ld", "filter.long.ld", "long.ld.missing", "ld.method", "ld.figures",
      "detect.mixed.genomes", "ind.heterozygosity.threshold",
      "detect.duplicate.genomes",
      "filter.hwe",
      "filter.strands",
      "random.seed",
      "path.folder", "filename",
      "parameters",
      "blacklist.genotypes", "erase.genotypes",
      "gt", "gt.bin", "gt.vcf", "gt.vcf.nuc",
      "pop.levels", "pop.labels", "pop.select", "blacklist.id",
      "markers.info", "keep.allele.names", "keep.gds", "calibrate.alleles",
      "vcf.metadata", "vcf.stats", "wide",
      "whitelist.markers",
      "write.tidy",
      "missing.memory",
      "dart.sequence",
      "internal",
      "heatmap.fst",
      "tidy.check", "tidy.vcf", "tidy.dart",
      "species",
      "population",
      "tau",
      "threshold.y.markers",
      "threshold.y.silico.markers",
      "sex.id.input",
      "threshold.x.markers.qr",
      "threshold.x.markers.RD",
      "threshold.x.markers.RD.silico",
      "mis.threshold.data",
      "mis.threshold.silicodata",
      "zoom.data",
      "zoom.silicodata",
      "sex.id.input",
      "het.qr.input"
    ),
    deprecated = c(
      "maf.thresholds",
      "common.markers",
      "max.marker",
      "monomorphic.out",
      "snp.ld",
      "filter.call.rate",
      "filter.markers.coverage",
      "filter.markers.missing",
      "number.snp.reads",
      "mixed.genomes.analysis",
      "duplicate.genomes.analysis",
      "ref.calibration"
    ),
    env     = parent.frame(),
    assign  = TRUE,
    verbose = TRUE
) {

  opt.change <- getOption("width")
  options(width = 70)
  on.exit(options(width = opt.change), add = TRUE)

  if (is.null(fd)) {
    fd <- rlang::fn_fmls_names()
  }
  if (is.null(args.list)) {
    args.list <- as.list(environment())
  }
  if (is.null(dotslist)) {
    dotslist <- list()
  }

  # ---------------------------------------------------------------------------
  # 1. Function call arguments
  # ---------------------------------------------------------------------------
  args.list <- purrr::map(.x = args.list, .f = check_args_class)

  func.call <- tibble::tibble(
    ARGUMENTS = names(args.list),
    VALUES    = args.list
  ) %>%
    dplyr::filter(ARGUMENTS %in% fd) %>%
    dplyr::mutate(GROUPS = "fct.call.args")

  if (verbose) {
    message("\n", func.name, " function call arguments:")
    purrr::walk2(
      .x = func.call$ARGUMENTS,
      .y = func.call$VALUES,
      .f = message_func_call,
      verbose = verbose
    )
  }

  res <- dplyr::mutate(func.call, VALUES = as.character(VALUES))

  # ---------------------------------------------------------------------------
  # 2. Dots processing
  # ---------------------------------------------------------------------------
  deprecated <- sort(unique(deprecated))
  keepers    <- sort(unique(keepers))

  want <- c(keepers, deprecated)

  unknowned_param <- setdiff(names(dotslist), want) |> sort()
  unk <- length(unknowned_param) > 0L

  dots.keepers <- dotslist[names(dotslist) %in% keepers]
  dots.keepers <- dots.keepers[sort(names(dots.keepers))]
  rdk <- length(dots.keepers) > 0L

  dots.deprecated <- dotslist[names(dotslist) %in% deprecated]
  rdd <- length(dots.deprecated) > 0L

  dots.defaults <- purrr::keep(
    .x = keepers,
    .p = ~ !.x %in% unique(c(deprecated, names(dotslist)))
  ) |> sort()
  rdf <- length(dots.defaults) > 0L

  if (unk || rdk || rdd) {
    if (verbose) message("\ndots-dots-dots ... arguments")
  }

  # 2a. Keepers actually supplied in ...
  if (rdk) {
    if (verbose) {
      message("\nArguments inside \"...\" assigned in ", func.name, ":")
    }

    res.df <- purrr::map2_df(
      .x = names(dots.keepers),
      .y = dots.keepers,
      .f = extract_dots,
      env.arg = env,
      assign  = assign,
      verbose = verbose
    )
    res <- dplyr::bind_rows(res, res.df)
  }

  # 2b. Defaults for keepers not present in ... nor deprecated
  if (rdf) {
    if (verbose) {
      message("\nDefault \"...\" arguments assigned in ", func.name, ":")
    }

    res.df <- purrr::map_df(
      .x = dots.defaults,
      .f = assign_defaults,
      env.arg = env,
      assign  = assign,
      verbose = verbose
    )
    res <- dplyr::bind_rows(res, res.df)
  }

  # 2c. Deprecated arguments
  if (rdd) {
    message("\nDeprecated arguments identified inside \"...\": ")
    message("    ", stringi::stri_join(sort(names(dots.deprecated)),
                                       collapse = "\n    "))
    res <- dplyr::bind_rows(
      res,
      tibble::tibble(ARGUMENTS = names(dots.deprecated)) |>
        dplyr::mutate(
          VALUES = "NA",
          GROUPS = "deprecated..."
        )
    )

    if (verbose) {
      check.strata <- c("pop.levels", "pop.labels", "pop.select", "blacklist.id")
      if (any(check.strata %in% names(dots.deprecated))) {
        message("\nNote: manipulating strata related arguments\n",
                "is best done inside radiator::read_strata\n")
      }
    }
  }

  # 2d. Unknown arguments
  if (unk) {
    message("\nUnknowned arguments identified inside \"...\": ")
    message("    ", stringi::stri_join(unknowned_param, collapse = "\n    "))
    res <- dplyr::bind_rows(
      res,
      tibble::tibble(ARGUMENTS = unknowned_param) |>
        dplyr::mutate(
          VALUES = "NA",
          GROUPS = "unknowned..."
        )
    )
  }

  if (rdd || unk) {
    message("\nRead documentation for latest changes, and modify your codes!\n")
  }

  if (verbose) message("\n")

  res
} # End radiator_dots

# Internal nested functions ----------------------------------------------------

#' @title message_func_call
#' @description Message the function call arguments
#' @name message_func_call
#' @keywords internal
#' @export
message_func_call <- function(n, v, verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  v.txt <- check_args_class(v)
  message("    ", n, " = ", v.txt)
}#END message_func_call


#' @title extract_dots
#' @description Extract and assign `...` arguments (keepers)
#' @name extract_dots
#' @keywords internal
#' @export
extract_dots <- function(n, v, env.arg, assign = TRUE, verbose = TRUE) {

  # mutate value for nicer printing in some special cases ----------------------
  if (identical(n, "path.folder") && !is.null(v)) {
    v.print <- basename(v)
  } else if (identical(n, "subsample")) {
    v.print <- length(v)
  } else if (identical(n, "pop.levels")) {
    v.print <- length(v)
  } else if (identical(n, "pop.labels")) {
    v.print <- length(v)
  } else if (identical(n, "quantiles.ci")) {
    v.print <- paste(v, collapse = "-")
  } else {
    v.print <- v
  }

  if (assign) {
    assign(x = n, value = v, envir = env.arg)
  }

  v.txt <- check_args_class(v.print)
  if (verbose) message("    ", n, " = ", v.txt)

  tibble::tibble(
    ARGUMENTS = n,
    VALUES    = as.character(v.txt),
    GROUPS    = "fct.call..."
  )
}#END extract_dots

#' @title assign_defaults
#' @description Assign the default values for `...` arguments
#' @name assign_defaults
#' @keywords internal
#' @export
assign_defaults <- function(n, env.arg, assign = TRUE, verbose = TRUE) {

  v <- NULL  # default: NULL

  # Arguments with default TRUE
  dots.true <- c(
    "keep.gds", "vcf.stats", "vcf.metadata",
    "filter.common.markers", "filter.monomorphic",
    "ld.figures", "dart.sequence", "force.stats"
  )

  # Arguments with default FALSE
  dots.false <- c(
    "keep.allele.names", "calibrate.alleles", "long.ld.missing",
    "detect.mixed.genomes", "detect.duplicate.genomes",
    "dp", "internal", "heatmap.fst", "wide", "filter.hwe",
    "gt", "gt.bin", "gt.vcf", "gt.vcf.nuc",
    "filter.haplotype.format"
  )

  if (n %in% dots.true)  v <- TRUE
  if (n %in% dots.false) v <- FALSE

  # Specific defaults
  if (identical(n, "filter.strands"))         v <- "blacklist"
  if (identical(n, "ld.method"))             v <- "r2"
  if (identical(n, "iteration.subsample"))   v <- 1L

  if (assign) {
    assign(x = n, value = v, envir = env.arg)
  }

  v.txt <- check_args_class(v)
  if (verbose) message("    ", n, " = ", v.txt)

  tibble::tibble(
    ARGUMENTS = n,
    VALUES    = as.character(v.txt),
    GROUPS    = "default..."
  )
}#End assign_defaults

#' @title check_args_class
#' @description Check the class of the argument/parameter value and return a
#'   compact printable representation (value or class name).
#' @name check_args_class
#' @keywords internal
#' @export
check_args_class <- function(x) {
  cls <- class(x)[1]
  if (!cls %in% c("logical", "character", "numeric", "double", "integer")) {
    # Non-simple objects are printed by class name
    out <- cls
  } else {
    out <- x
  }
  if (length(out) > 1) {
    out <- paste(out, collapse = ", ")
  }
  out
}# End check_args_class


# dots_keepers_core ------------------------------------------------------------
# R/radiator_dots_defaults.R
# Core global keepers/deprecated used across radiator functions


#' @title Global dot-argument "keepers" for radiator
#' @description
#' Internal helper returning the core set of supported `...` arguments
#' used across most radiator functions.
#' @keywords internal
#' @export
dots_keepers_core <- function() {
  c(
    "subsample.markers.stats", "force.stats", "id.stats", "subsample",
    "filter.reproducibility",
    "filter.individuals.missing",
    "filter.individuals.heterozygosity",
    "filter.individuals.coverage.total",
    "filter.individuals.coverage.median",
    "filter.individuals.coverage.iqr",
    "filter.common.markers", "filter.monomorphic",
    "filter.ma",
    "ma.stats",
    "filter.coverage", "dp",
    "filter.genotyping",
    "filter.snp.position.read",
    "filter.snp.number",
    "filter.short.ld", "filter.long.ld", "long.ld.missing",
    "ld.method", "ld.figures",
    "detect.mixed.genomes", "ind.heterozygosity.threshold",
    "detect.duplicate.genomes",
    "filter.hwe",
    "filter.strands",
    "random.seed",
    "path.folder", "filename",
    "parameters",
    "blacklist.genotypes", "erase.genotypes",
    "gt", "gt.bin", "gt.vcf", "gt.vcf.nuc",
    "pop.levels", "pop.labels", "pop.select", "blacklist.id",
    "markers.info", "keep.allele.names", "keep.gds", "calibrate.alleles",
    "vcf.metadata", "vcf.stats", "wide",
    "whitelist.markers",
    "write.tidy",
    "missing.memory",
    "dart.sequence",
    "internal",
    "heatmap.fst",
    "tidy.check", "tidy.vcf", "tidy.dart",
    "species",
    "population",
    "tau",
    "threshold.y.markers",
    "threshold.y.silico.markers",
    "sex.id.input",
    "threshold.x.markers.qr",
    "threshold.x.markers.RD",
    "threshold.x.markers.RD.silico",
    "mis.threshold.data",
    "mis.threshold.silicodata",
    "zoom.data",
    "zoom.silicodata",
    "het.qr.input"
  ) |> unique()
}#END dots_keepers_core


# dots_deprecated_core --------------------------------------------------------
#' @title Global deprecated dot-arguments for radiator
#' @description
#' Internal helper returning deprecated argument names for `...`.
#' @keywords internal
#' @export
dots_deprecated_core <- function() {
  c(
    "maf.thresholds",
    "common.markers",
    "max.marker",
    "monomorphic.out",
    "snp.ld",
    "filter.call.rate",
    "filter.markers.coverage",
    "filter.markers.missing",
    "number.snp.reads",
    "mixed.genomes.analysis",
    "duplicate.genomes.analysis",
    "ref.calibration"
  ) |> unique()
}#END dots_deprecated_core
