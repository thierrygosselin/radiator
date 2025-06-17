# radiator_dots-----------------------------------------------------------------
#' @title radiator_dots
#' @description Extract and assign the dots-dots-dots
#' @name radiator_dots
#' @rdname radiator_dots
#' @param func.name Default: \code{as.list(sys.call())[[1]]}.
#' @param fd (optional) Default: \code{rlang::fn_fmls_names()}.
#' @param args.list (optional) Default:\code{args.list = as.list(environment())}.
#' @param dotslist The dots dots dots arguments captured.
#' Default: \code{
#' dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)}.
#' @param keepers (optional) The arguments that are used inside the function.
#' Default: \code{keepers = c(
#' "subsample.markers.stats", "subsample",
#' "filter.reproducibility", "filter.individuals.missing",
#' "filter.individuals.heterozygosity", "filter.individuals.coverage.total",
#' "filter.individuals.coverage.median", "filter.individuals.coverage.iqr",
#' "filter.common.markers", "filter.monomorphic", "filter.ma", "ma.stats"
#' "filter.coverage", "filter.genotyping", "filter.snp.position.read",
#' "filter.snp.number", "filter.short.ld", "filter.long.ld", "long.ld.missing",
#' "ld.method", "detect.mixed.genomes", "ind.heterozygosity.threshold",
#' "detect.duplicate.genomes",
#' "filter.hwe", "filter.strands", "random.seed", "path.folder", "filename",
#' "blacklist.genotypes", "erase.genotypes",
#' "gt", "gt.bin", "gt.vcf", "gt.vcf.nuc",
#' "pop.levels", "pop.labels", "pop.select", "blacklist.id",
#' "markers.info", "keep.allele.names", "keep.gds", "calibrate.alleles",
#' "vcf.metadata", "vcf.stats", "id.stats", "dp",
#' "whitelist.markers",
#' "write.tidy",
#' "dart.sequence",
#' "missing.memory",
#' "internal", "heatmap.fst", "tidy.check", "tidy.vcf", "tidy.dart",
#' "species", "population", "tau", "threshold.y.markers", "threshold.y.silico.markers",
#' "sex.id.input", "threshold.x.markers.qr", "threshold.x.markers.RD", "threshold.x.markers.RD.silico",
#' "mis.threshold.data", "mis.threshold.silicodata", "zoom.data", "zoom.silicodata",
#' "sex.id.input", "het.qr.input")}.
#' @param deprecated (optional) radiator's deprecated arguments.
#' Default: \code{deprecated = c("maf.thresholds", "common.markers",
#' "max.marker","monomorphic.out", "snp.ld", "filter.call.rate",
#' "filter.markers.coverage", "filter.mac","filter.markers.missing", "number.snp.reads",
#' "mixed.genomes.analysis", "duplicate.genomes.analysis", "maf.data", "ref.calibration")}.
#' @keywords internal
#' @inheritParams radiator_common_arguments
#'
#' @export

radiator_dots <- function(
  func.name = as.list(sys.call())[[1]],
  fd = NULL,
  args.list = NULL,
  dotslist = NULL,
  keepers = c(
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
  verbose = TRUE
) {
  opt.change <- getOption("width")
  options(width = 70)
  env.arg <- parent.frame()
  res <- tibble::tibble(
    ARGUMENTS = character(0),
    VALUES = character(0),
    GROUPS = character(0)
    )

  # function call --------------------------------------------------------------
  args.list <- purrr::map(.x = args.list, .f = check_args_class)

  func.call <- tibble::tibble(
    ARGUMENTS = names(args.list),
    VALUES = args.list
  ) %>%
    dplyr::filter(ARGUMENTS %in% fd) %>%
    dplyr::mutate(GROUPS = "fct.call.args")#,VALUES = VALUES)
  # print(func.call)

  if (verbose) if (verbose) message("\n", func.name, " function call arguments:")
  purrr::walk2(
    .x = func.call$ARGUMENTS,
    .y = func.call$VALUES,
    .f = message_func_call,
    verbose = verbose
  )
  res %<>% dplyr::bind_rows(dplyr::mutate(func.call, VALUES = as.character(VALUES)))

  # Dots dots dots -------------------------------------------------------------
  deprecated <- sort(deprecated)
  keepers <- sort(keepers)

  want <- c(keepers, deprecated)
  unknowned_param <- setdiff(names(dotslist), want)
  unknowned_param <- sort(unknowned_param)

  unk <- length(unknowned_param) > 0

  dots.keepers <- dotslist[names(dotslist) %in% keepers]
  dots.keepers <- dots.keepers[sort(names(dots.keepers))]
  rdk <- length(dots.keepers) > 0

  dots.deprecated <- dotslist[names(dotslist) %in% deprecated]
  rdd <- length(dots.deprecated) > 0

  dots.defaults <- purrr::keep(
    .x = keepers,
    .p = !keepers %in% unique(c(deprecated, names(dotslist)))
  ) %>% sort
  rdf <- length(dots.defaults) > 0


  if (unk || rdk || rdd)
    if (verbose) message("\ndots-dots-dots ... arguments")

  # The args present:
  if (rdk) {
    if (verbose) message("\nArguments inside \"...\" assigned in ", func.name, ":")
    res.df <- purrr::map2_df(
      .x = names(dots.keepers),
      .y = dots.keepers,
      .f = extract_dots,
      env.arg = env.arg,
      verbose = verbose
      )
    res %<>% dplyr::bind_rows(res.df)
  }

  # defaults
  if (rdf) {
    if (verbose) message("\nDefault \"...\" arguments assigned in ", func.name, ":")
    res.df <- purrr::map_df(
      .x = dots.defaults,
      .f = assign_defaults,
      env.arg = env.arg,
      verbose = verbose
      )
    res %<>% dplyr::bind_rows(res.df)
  }


  # The deprecated args
  if (rdd) {
    message("\nDeprecated arguments identified inside \"...\": ")
    message("    ", stringi::stri_join(sort(names(dots.deprecated)),
                                                    collapse = "\n    "))
    res %<>% dplyr::bind_rows(
      tibble::tibble(ARGUMENTS = names(dots.deprecated)) %>%
        dplyr::mutate(
          VALUES = "NA", GROUPS = "deprecated..."
        )
    )

    if (verbose) {
      check.strata <- c("pop.levels", "pop.labels", "pop.select", "blacklist.id")
      if (TRUE %in% (check.strata %in% names(dots.deprecated))) {
        message("\nNote: manipulating strata related arguments\nis best done inside the function radiator::read_strata\n")
      }
    }
  }
  if (unk) {
    message("\nUnknowned arguments identified inside \"...\": ")
    message("    ", stringi::stri_join(unknowned_param, collapse = "\n    "))
    res %<>% dplyr::bind_rows(
      tibble::tibble(ARGUMENTS = unknowned_param) %>%
        dplyr::mutate(
          VALUES = "NA", GROUPS = "unknowned..."
        )
    )
  }

  if (rdd || unk) {
    message("\nRead documentation, for latest changes, and modify your codes!\n")
  }
  options(width = opt.change)
  if (verbose) message("\n")
  return(res)
}#End radiator_dots


# Internal nested functions ----------------------------------------------------

#' @title message_func_call
#' @description Message the function call
#' @name message_func_call
#' @keywords internal
#' @export
message_func_call <- function(n, v, verbose = TRUE) {
  if (verbose) {
    message(
      "    ",
      stringi::stri_join(n, " = ", paste(rlang::quo_name(v), collapse = "," ))
    )
  }
}# End message_func_call


#' @title extract_dots
#' @description Extract and Assign ...
#' @name extract_dots
#' @keywords internal
#' @export

extract_dots <- function(n, v, env.arg, verbose = TRUE) {
  assign(x = n, value = v, pos = env.arg, envir = env.arg)
  if (n == "path.folder" && !is.null(v)) v <- basename(v)
  if (n == "subsample") v <- length(n)
  if (n == "pop.levels") v <- length(n)
  if (n == "pop.labels") v <- length(n)
  if (n == "quantiles.ci") v <- paste(n, collapse = "-")

  v <- check_args_class(x = v)
  if (verbose) message("    ", n, " = ", v)
  res <- tibble::tibble(
    ARGUMENTS = n,
    VALUES = rlang::quo_name(v),
    GROUPS = "fct.call...")
  return(res)
}#End extract_dots

#' @title assign_defaults
#' @description Assign the default values for...
#' @name assign_defaults
#' @keywords internal
#' @export
assign_defaults <- function(n, env.arg, verbose = TRUE) {
  v <- NULL # by defaults all NULL

  # Specifics...
  # Arguments that default value is TRUE
  dots.true <- c(
    "keep.gds", "vcf.stats", "vcf.metadata",
    "filter.common.markers", "filter.monomorphic",
    "ld.figures", "dart.sequence", "force.stats"
  )
  # Arguments that default value is FALSE
  dots.false <- c(
    "keep.allele.names", "calibrate.alleles", "long.ld.missing",
    "detect.mixed.genomes", "detect.duplicate.genomes",
    "dp", "internal", "heatmap.fst", "wide", "filter.hwe",
    "gt", "gt.bin", "gt.vcf", "gt.vcf.nuc"
  )
  if (n %in% dots.true) v <- TRUE
  if (n %in% dots.false) v <- FALSE

  # Specific values...
  if (n == "filter.strands") v <- "blacklist"
  if (n == "ld.method") v <- "r2"
  if (n == "iteration.subsample") v <- 1L

  # assignment
  assign(rlang::quo_name(n), v, pos = env.arg, envir = env.arg)
  if (verbose) message("    ", n, " = ", rlang::quo_name(v))
  v <- check_args_class(x = v)
  res <- tibble::tibble(
    ARGUMENTS = n,
    VALUES = rlang::quo_name(v),
    GROUPS = "default...")
  return(res)
}#End assign_defaults

#' @title check_args_class
#' @description Check the class of the argument/parameter value
#' @name check_args_class
#' @keywords internal
#' @export
check_args_class <- function(x) {
  y <- class(x)[1]
  if (!y %in% c("logical", "character", "numeric", "double", "integer")) {
    x <- y
  } else {
    x
  }
  if (length(x) > 1) x <- paste(x, collapse = ", ")
  return(x)
}# End check_args_class

