# radiator package startup message ---------------------------------------------
#' @importFrom utils packageDescription
#' @importFrom stringi stri_join

.onAttach <- function(libname, pkgname) {
  radiator.version <- utils::packageDescription("radiator", fields = "Version")
  radiator.build <- utils::packageDescription("radiator", fields = "Built")
  startup.message <- stringi::stri_join(
    "******************************* IMPORTANT NOTICE *******************************\n",
    "radiator v.", radiator.version, " was modified heavily.\n",
    "Read functions documentation and available vignettes.\n\n",
    "For reproducibility:\n",
    "    radiator version: ", radiator.version,"\n",
    "    radiator build date: ", radiator.build,"\n",
    "    Keep zenodo DOI.\n",
    "********************************************************************************",
    sep = "")
  packageStartupMessage(startup.message)
}

# magrittr ---------------------------------------------------------------------
#' @title Forward-pipe operator
#' @description magrittr forward-pipe operator
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

# Exposition pipe-operator
#' @title Exposition pipe-operator
#' @description magrittr Exposition pipe-operator
#' @name %$%
#' @rdname Exposition_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %$%
#' @usage lhs \%$\% rhs
NULL

# compound assignment pipe operator
#' @title compound assignment pipe operator
#' @description magrittr compound assignment pipe operator
#' @name %<>%
#' @rdname compound_assignment_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL

# dplyr n ----------------------------------------------------------------------
# The number of observations in the current group.
#' @title The number of observations in the current group.
#' @description Check dplyr
#' @name n
#' @rdname n
#' @keywords internal
#' @export
#' @importFrom dplyr n
#' @usage n()
NULL
# split_vec_row ----------------------------------------------------------------
#' @title split_vec_row
#' @description Split input into chunk for parallel processing
#' @rdname split_vec_row
#' @keywords internal
#' @export
split_vec_row <- function(x, cpu.rounds, parallel.core = parallel::detectCores() - 1) {
  if (!is.integer(x)) {
    n.row <- nrow(x)
  } else {
    n.row <- x
  }
  split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.row - 1) / n.row) + 1))
  return(split.vec)
}#End split_vec_row

# replace_by_na ----------------------------------------------------------------
#' @title replace_by_na
#' @description Fast removal of NA
#' @rdname replace_by_na
#' @keywords internal
#' @export
replace_by_na <- function(data, what = ".") {
  replace(data, which(data == what), NA)
}#End replace_by_na

# separate_gt ------------------------------------------------------------------
#' @title separate_gt
#' @description Separate genotype field
#' @rdname separate_gt
#' @keywords internal
#' @export
separate_gt <- function(
  x,
  sep = "/",
  gt = "GT_VCF_NUC",
  gather = TRUE,
  exclude = c("LOCUS", "INDIVIDUALS", "POP_ID"),
  cpu.rounds = 10,
  parallel.core = parallel::detectCores() - 1
) {
  # sep <-  "/"
  # x <- input2
  n.row <- nrow(x)
  split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.row - 1) / n.row) + 1))

  separate_genotype <- function(x, sep, gt, exclude){
    res <- tidyr::separate(
      data = x,
      col = gt, into = c("ALLELE1", "ALLELE2"),
      sep = sep,
      extra = "drop", remove = TRUE
    )

    if (gather) {
      res <- tidyr::gather(data = res, key = ALLELE_GROUP,
                           value = HAPLOTYPES, -dplyr::one_of(exclude))
    }
    return(res)
  }

  res <- split(x = x, f = split.vec) %>%
    .radiator_parallel_mc(
      X = .,
      FUN = separate_genotype,
      mc.cores = parallel.core,
      sep = sep, gt = gt, exclude = exclude) %>%
    dplyr::bind_rows(.)
  return(res)
}#End separate_gt

# melt the dist matrice into a tibble df --------------------------------------
#' @title distance2tibble
#' @description melt the distance matrix into a tibble
#' @rdname distance2tibble
#' @export
#' @keywords internal
distance2tibble <- function(
  x,
  remove.diag = TRUE,
  na.diag = FALSE,
  remove.lower = TRUE,
  relative = TRUE,
  pop.levels = NULL,
  distance.class.double = TRUE
) {
  # x <- dist.computation
  x <- as.matrix(x)
  if (remove.diag || na.diag) diag(x) <- NA
  if (remove.lower) x[lower.tri(x)] <- NA
  x <- dplyr::bind_cols(tibble::tibble(ID1 = rownames(x)),
                        tibble::as_tibble(x)) %>%
    data.table::as.data.table(.) %>%
    data.table::melt.data.table(
      data = ., id.vars = "ID1", variable.name = "ID2", value.name = "DISTANCE",
      variable.factor = FALSE) %>%
    tibble::as_tibble(.)

  if (!na.diag) x  %<>% dplyr::filter(!is.na(DISTANCE))

  if (distance.class.double) {
    x %<>% dplyr::mutate(DISTANCE = as.double(as.character(DISTANCE)))
  }

  x %<>% dplyr::arrange(DISTANCE)

  if (relative) {
    x  %<>% dplyr::mutate(DISTANCE_RELATIVE = DISTANCE/max(DISTANCE))
  }
  if (!is.null(pop.levels)) {
    x  %<>% dplyr::mutate(
      ID1 = factor(x = ID1, levels = pop.levels, ordered = TRUE),
      ID2 = factor(x = ID2, levels = pop.levels, ordered = TRUE)
    )
  }

  return(x)
}#End distance2tibble


# count ind_total_reads -------------------------------------------------------
#' @title ind_total_reads
#' @description Counts the total number of reads per samples
#' @rdname ind_total_reads
#' @keywords internal
#' @export
ind_total_reads <- function(x, path.folder = NULL) {
  # x <- unfiltered
  # path.folder <- path.folder.coverage
  if (is.null(path.folder)) path.folder <- getwd()
  # generate the stats
  read.info <- dplyr::group_by(x, INDIVIDUALS, POP_ID) %>%
    dplyr::summarise(TOTAL_READ_COUNTS = sum(READ_DEPTH, na.rm = TRUE))

  if (is.factor(read.info$POP_ID)) {
    read.pop.levels <- levels(read.info$POP_ID)
  } else {
    read.pop.levels <- unique(read.info$POP_ID)
  }
  read.pop.levels <- c(read.pop.levels, "OVERALL")

  overall.read <- dplyr::mutate(read.info, POP_ID = "OVERALL")


  read.info <- suppressWarnings(dplyr::bind_rows(read.info, overall.read)) %>%
    dplyr::mutate(POP_ID = factor(POP_ID, levels = read.pop.levels))

  read.counts.stats <- read.info %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      MEAN = mean(TOTAL_READ_COUNTS, na.rm = TRUE),
      MEDIAN = stats::median(TOTAL_READ_COUNTS, na.rm = TRUE),
      SD = stats::sd(TOTAL_READ_COUNTS, na.rm = TRUE),
      Q25 = stats::quantile(TOTAL_READ_COUNTS, 0.25, na.rm = TRUE),
      Q75 = stats::quantile(TOTAL_READ_COUNTS, 0.75, na.rm = TRUE),
      IQR = stats::IQR(TOTAL_READ_COUNTS, na.rm = TRUE),
      MIN = min(TOTAL_READ_COUNTS, na.rm = TRUE),
      MAX = max(TOTAL_READ_COUNTS, na.rm = TRUE),
      OUTLIERS_LOW = Q25 - (1.5 * IQR),
      OUTLIERS_HIGH = Q75 + (1.5 * IQR),
      OUTLIERS_LOW_N = length(TOTAL_READ_COUNTS[TOTAL_READ_COUNTS < OUTLIERS_LOW]),
      OUTLIERS_HIGH_N = length(TOTAL_READ_COUNTS[TOTAL_READ_COUNTS > OUTLIERS_HIGH]),
      OUTLIERS_TOTAL = OUTLIERS_HIGH_N + OUTLIERS_LOW_N,
      OUTLIERS_PROP = round(OUTLIERS_TOTAL / length(x), 3)
    )

  # plots
  element.text <- ggplot2::element_text(size = 10,
                                        family = "Helvetica", face = "bold")
  n.pop <- dplyr::n_distinct(x$POP_ID)
  ind.plot <- suppressWarnings(
    ggplot2::ggplot(
      read.info, ggplot2::aes(x = POP_ID, y = TOTAL_READ_COUNTS, na.rm = TRUE)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(x = "Populations",
                    y = "Individuals total read counts",
                    title = "Individuals total read counts") +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold", hjust = 0.5),
        # plot.subtitle = ggplot2::element_text(size = 10, family = "Helvetica", hjust = 0.5),
        axis.title.y = element.text,
        axis.title.x = element.text,
        axis.text.x = element.text))

  suppressWarnings(
    ggplot2::ggsave(
      plot = ind.plot,
      filename = file.path(path.folder, "plot.ind.total.reads.pdf"),
      width = n.pop * 2, height = 10, dpi = 300, units = "cm",
      useDingbats = FALSE)
  )

}#End ind_coverage


# interactive_question ---------------------------------------------------------
#' @title interactive_question
#' @description Ask to enter a word or number
#' @rdname interactive_question
#' @keywords internal
#' @export
interactive_question <- function(x, answer.opt = NULL, minmax = NULL) {
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
}#End interactive_question

# markers_genotyped_helper------------------------------------------------------
#' @title markers_genotyped_helper
#' @description Help individual's genotyped threshold
#' @rdname markers_genotyped_helper
#' @export
#' @keywords internal
markers_genotyped_helper <- function(x, y, overall.only = FALSE) {
  # x <- res$missing.genotypes.markers.pop
  # Set the breaks for the figure
  max.markers <- dplyr::n_distinct(y$MARKERS)

  threshold.helper.overall <- y %>%
    dplyr::ungroup(.) %>%
    dplyr::summarise(
      `0` = length(PERCENT[PERCENT == 0]),
      `10` = length(PERCENT[PERCENT <= 10]),
      `20` = length(PERCENT[PERCENT <= 20]),
      `30` = length(PERCENT[PERCENT <= 30]),
      `40` = length(PERCENT[PERCENT <= 40]),
      `50` = length(PERCENT[PERCENT <= 50]),
      `60` = length(PERCENT[PERCENT <= 60]),
      `70` = length(PERCENT[PERCENT <= 70]),
      `80` = length(PERCENT[PERCENT <= 80]),
      `90` = length(PERCENT[PERCENT <= 90]),
      `100` = length(PERCENT[PERCENT <= 100])
    ) %>%
    tidyr::gather(data = ., key = GENOTYPED_THRESHOLD, value = NUMBER_MARKERS) %>%
    dplyr::mutate(POP_ID = rep("OVERALL", n()))

  if (!overall.only){
    threshold.helper.pop <- x %>%
      dplyr::group_by(POP_ID) %>%
      dplyr::summarise(
        `0` = length(PERCENT[PERCENT == 0]),
        `10` = length(PERCENT[PERCENT <= 10]),
        `20` = length(PERCENT[PERCENT <= 20]),
        `30` = length(PERCENT[PERCENT <= 30]),
        `40` = length(PERCENT[PERCENT <= 40]),
        `50` = length(PERCENT[PERCENT <= 50]),
        `60` = length(PERCENT[PERCENT <= 60]),
        `70` = length(PERCENT[PERCENT <= 70]),
        `80` = length(PERCENT[PERCENT <= 80]),
        `90` = length(PERCENT[PERCENT <= 90]),
        `100` = length(PERCENT[PERCENT <= 100])
      ) %>%
      tidyr::gather(data = ., key = GENOTYPED_THRESHOLD, value = NUMBER_MARKERS, -POP_ID)

    mean.pop <- threshold.helper.pop %>%
      dplyr::group_by(GENOTYPED_THRESHOLD) %>%
      dplyr::summarise(
        NUMBER_MARKERS = round(mean(NUMBER_MARKERS), 0)
      ) %>%
      dplyr::mutate(POP_ID = rep("MEAN_POP", n()))

    # Check if x$POP_ID is a factor
    if (is.factor(x$POP_ID)) {
      x.pop.levels <- levels(x$POP_ID)
    } else {
      x.pop.levels <- unique(x$POP_ID)
    }

    threshold.helper <- suppressWarnings(
      dplyr::bind_rows(threshold.helper.pop, mean.pop, threshold.helper.overall) %>%
        dplyr::mutate(
          GENOTYPED_THRESHOLD = as.numeric(GENOTYPED_THRESHOLD),
          POP_ID = factor(POP_ID, levels = c(x.pop.levels, "MEAN_POP", "OVERALL"), ordered = TRUE)
        ))
    threshold.helper.pop <- mean.pop <- threshold.helper.overall <- x <- y <- NULL
  } else {
    threshold.helper <- threshold.helper.overall %>%
      dplyr::mutate(
        GENOTYPED_THRESHOLD = as.numeric(GENOTYPED_THRESHOLD)
      )
    threshold.helper.overall <- x <- y <- NULL
  }
  #Function to replace package plyr round_any
  rounder <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }

  if (max.markers >= 1000) {
    y.breaks.by <- rounder(max.markers / 10, 100, ceiling)
    y.breaks.max <- rounder(max.markers, 1000, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  } else {
    y.breaks.by <- rounder(max.markers / 10, 10, ceiling)
    y.breaks.max <- rounder(max.markers, 100, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  }


  plot.markers.geno.threshold <- ggplot2::ggplot(
    threshold.helper,
    ggplot2::aes(x = GENOTYPED_THRESHOLD, y = NUMBER_MARKERS)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
    ggplot2::scale_x_continuous(name = "Marker's missing genotype threshold (percent)", breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
    ggplot2::scale_y_continuous(name = "Markers\n(whitelisted number)", breaks = y.breaks, limits = c(0, y.breaks.max)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica"),#, angle = 90, hjust = 1, vjust = 0.5),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    ) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(~POP_ID)
  # plot.markers.geno.threshold
  return(plot.markers.geno.threshold)
}#End markers_genotyped_helper

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
      VALUES = values,
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
    res$filters.parameters.path <-parameter.obj$filters.parameters.path
  }#End update

  # messages
  if (initiate && update) {
    # if (verbose && verbose.on) message("Filters parameters file: initiated and updated")
    if (verbose) message("Filters parameters file: initiated and updated")
  }

  if (initiate && !update) {
    # if (verbose && verbose.on) message("Filters parameters file: initiated")
    if (verbose) message("Filters parameters file: initiated")
  }

  if (!initiate && update) {
    # if (verbose && verbose.on) message("Filters parameters file: updated")
    if (verbose) message("Filters parameters file: updated")
  }

  return(res)
}#End radiator_parameters
# data.info -------------------------------------------------------------
#' @title data_info
#' @description function generate tidy data main info
#' @rdname data_info
#' @keywords internal
#' @export
data_info <- function(x, print.info = FALSE) {
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
    res$n.chrom <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/markers.meta/CHROM", silent = TRUE))))
    res$n.locus <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/markers.meta/LOCUS", silent = TRUE))))
    res$n.snp <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/markers.meta/MARKERS", silent = TRUE))))
    res$n.pop <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/individuals/STRATA", silent = TRUE))))
    res$n.ind <- length(unique(gdsfmt::read.gdsn(gdsfmt::index.gdsn(node = x, path = "radiator/individuals/INDIVIDUALS", silent = TRUE))))
    res[is.null(res)] <- NA_integer_
  }

  if (print.info) {
    message("Number of chrom: ", res$n.chrom)
    message("Number of locus: ", res$n.locus)
    message("Number of SNPs: ", res$n.snp)
    message("Number of populations: ", res$n.pop)
    message("Number of individuals: ", res$n.ind)
  }
  return(res)
}

# tibble_stats-----------------------------------------------------------------
#' @title tibble_stats
#' @description Generate a tibble of statistics
#' @rdname tibble_stats
#' @keywords internal
#' @export
tibble_stats <- function(x, group, subsample = NULL) {
  if (is.null(subsample)) subsample <- 1L
  Q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  res <- tibble::tibble(
    GROUP = group,
    MIN = min(x, na.rm = TRUE),
    # Q25 = stats::quantile(x, 0.25, na.rm = TRUE),
    Q25 = Q[1],
    MEDIAN = stats::median(x, na.rm = TRUE),
    Q75 = Q[2],
    MAX = max(x, na.rm = TRUE),
    IQR = abs(diff(Q)),
    # IQR = stats::IQR(depth, na.rm = TRUE),
    OUTLIERS_LOW = Q25 - (1.5 * IQR),
    OUTLIERS_HIGH =  Q75 + (1.5 * IQR)) %>%
    dplyr::mutate_if(.tbl = ., .predicate = is.integer, .funs = as.numeric) %>%
    dplyr::mutate(
      OUTLIERS_LOW = dplyr::if_else(OUTLIERS_LOW < MIN, MIN, OUTLIERS_LOW),
      OUTLIERS_HIGH = dplyr::if_else(OUTLIERS_HIGH > MAX, MAX, OUTLIERS_HIGH),
      SUBSAMPLE = subsample
    )
  Q <- NULL
  return(res)
}#End tibble_stats

# boxplot_stats-----------------------------------------------------------------
# boxplot of stats
#' @title boxplot_stats
#' @description Generate a boxplot
#' @rdname boxplot_stats
#' @keywords internal
#' @export
boxplot_stats <- function(
  data,
  title,
  subtitle = NULL,
  x.axis.title = NULL,
  y.axis.title,
  facet.columns = FALSE,
  facet.rows = FALSE,
  bp.filename = NULL,
  path.folder = NULL
) {
  # data <- test
  # x.axis.title = NULL
  # x.axis.title <- "SNP position on the read groupings"
  # title <- "Individual's QC stats"
  # subtitle = "testing"
  # y.axis.title <- "Statistics"
  # y.axis.title <- "SNP position (base pair)"
  # bp.filename <- "vcf.snp.position.read.pdf"
  # bp.filename <- "test.pdf"
  # facet.columns = TRUE
  # facet.rows = FALSE
  # path.folder = NULL

  if (is.null(path.folder)) path.folder <- getwd()

  n.group <- dplyr::n_distinct(data$GROUP)
  element.text <- ggplot2::element_text(size = 10,
                                        family = "Helvetica", face = "bold")

  if (facet.columns) {
    data <- dplyr::mutate(data, X = "1")
    fig.boxplot <- ggplot2::ggplot(data = data, ggplot2::aes(X))
  } else {
    fig.boxplot <- ggplot2::ggplot(data = data, ggplot2::aes(GROUP))
  }


  fig.boxplot <- fig.boxplot +
    ggplot2::geom_boxplot(
      ggplot2::aes(ymin = OUTLIERS_LOW, lower = Q25, middle = MEDIAN, upper = Q75,
                   ymax = OUTLIERS_HIGH), stat = "identity") +
    ggplot2::labs(y = y.axis.title, title = title)

  if (!is.null(subtitle)) fig.boxplot <- fig.boxplot + ggplot2::labs(subtitle = subtitle)

  # Draw upper outliers
  if (facet.columns) {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = "1", xend = "1",
          y = OUTLIERS_HIGH, yend = MAX),
        linetype = "dashed")
  } else {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = GROUP, xend = GROUP,
          y = OUTLIERS_HIGH, yend = MAX),
        linetype = "dashed")
  }

  # Draw lower outliers
  if (facet.columns) {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = "1", xend = "1",
          y = OUTLIERS_LOW, yend = MIN),
        linetype = "dashed")
  } else {
    fig.boxplot <- fig.boxplot +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = GROUP, xend = GROUP,
          y = OUTLIERS_LOW, yend = MIN),
        linetype = "dashed")
  }

  fig.boxplot <- fig.boxplot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, family = "Helvetica",
                                         face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.title.y = element.text,
      axis.text.y = element.text
    ) +
    ggplot2::theme_bw()

  if (is.null(x.axis.title)) {
    fig.boxplot <- fig.boxplot +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  } else {
    fig.boxplot <- fig.boxplot +
      ggplot2::xlab(x.axis.title) +
      ggplot2::theme(
        axis.title.x = element.text,
        # axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  if (!is.null(subtitle)) {
    fig.boxplot <- fig.boxplot +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_text(size = 10, family = "Helvetica"))
  }

  if (facet.columns) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(GROUP ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 10
    height <- 10 + (4 * n.group)
  }

  # else {
  # width <-  10 + (5 * n.group) + 1
  # height <-  10
  # }

  if (facet.rows) {
    fig.boxplot <- fig.boxplot + ggplot2::facet_grid(FACET_ROWS ~ ., scales = "free")
    n.facet <- n.group * 2
    width <- 10
    height <- 10 + (4 * n.group)
  }

  if (!facet.rows && !facet.columns) {
    width <-  10 + (5 * n.group) + 1
    height <-  10
  }

  print(fig.boxplot)
  if (!is.null(bp.filename)) {
    suppressMessages(ggplot2::ggsave(
      filename = file.path(path.folder, bp.filename),
      plot = fig.boxplot,
      width = width,
      height = height,
      dpi = 300, units = "cm", useDingbats = FALSE))
  }
  return(fig.boxplot)
}#Endboxplot_stats


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
    "filter_mac",
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
      "filter.mac",
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
  path.folder = getwd(),
  date = TRUE,
  extension = c(
    "tsv", "gds.rad", "rad", "gds", "gen", "dat",
    "genind", "genlight", "gtypes", "vcf", "colony",
    "bayescan", "gsisim", "hierfstat", "hzar", "ldna",
    "pcadapt", "related", "stockr", "structure", "arlequin"
  )
) {

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
            "pcadapt", "plink", "related", "stockr", "structure", "arlequin")
  extension <- match.arg(extension, want)

  # note to myself: currently excluded output : "fineradstructure", "maverick", "plink", "betadiv"


  # with same extension
  # extension <- "tsv"
  if (extension %in% c("tsv", "gds.rad", "rad", "gds", "vcf", "colony", "ldna")) {
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

# folder_short------------------------------------------------------------------
#' @title folder_short
#' @description Extract the short name of the folder generated by radiator
#' @name folder_short
#' @param f Folder name
#' @rdname folder_short
# @keywords internal
#' @export
#' @rdname folder_short
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

# generate_folder---------------------------------------------------------------

#' @title generate_folder
#' @description Generate a folder based on ...
#' @name generate_folder
#' @rdname generate_folder
#' @param rad.folder Name of the rad folder
#' @param internal (optional, logical) Is the function internal or not
#' @param file.date The file date included
#' @inheritParams radiator_folder
#' @inheritParams radiator_common_arguments
#' @inheritParams folder_short
#' @export
#' @rdname generate_folder
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

generate_folder <- function(
  f,
  rad.folder = NULL,
  internal = FALSE,
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
    if (!date.present) f <- stringi::stri_join(f, file.date, sep = "_")
    if (!wd.present) f <- file.path(getwd(), f)
    if (verbose && !identical(f.temp, f)) message("Folder created: ", folder_short(f))
  }
  if (!dir.exists(f)) dir.create(f)
  return(f)
}#End generate_folder

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
# radiator_folder--------------------------------------------------------------------
#' @title radiator_folder
#' @description Generate the rad folders
#' @param path.folder path of the folder
#' @param prefix_int Use an integer prefix padded left with 0.
#' Default: \code{prefix_int = TRUE}.
#' @inheritParams folder_short
# @keywords internal
#' @export
#' @rdname radiator_folder
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

radiator_folder <- function(f, path.folder = NULL, prefix_int = TRUE) {
  if (is.null(path.folder)) path.folder <- getwd()
  if (prefix_int) {
    f <- stringi::stri_join(stringi::stri_pad_left(
      str = length(list.dirs(path = path.folder, full.names = FALSE)[-1]) + 1L,
      width = 2,
      pad = 0
    ), "_", f)
  }
  folder.prefix <- file.path(path.folder, f)
  return(folder.prefix)
}#End radiator_folder

#' @title rad_write
#' @description Generate the rad folders
#' @keywords internal
#' @export



# tidy2wide --------------------------------------------------------------------
#' @title tidy2wide
#' @description tidy2wide
#' @keywords internal
#' @export
tidy2wide <- function(x = NULL, gds = NULL, individuals = NULL, markers = NULL, tidy = TRUE, wide = TRUE, wide.markers = TRUE) {
  res <- list()
  if (is.null(markers) && !is.null(gds)) {
    markers <- extract_markers_metadata(gds = gds, markers.meta.select = "MARKERS") %$% MARKERS
  }
  if (is.null(individuals) && !is.null(gds)) {
    individuals <- extract_individuals(gds = gds, ind.field.select = "INDIVIDUALS") %$% INDIVIDUALS
  }

  n.markers <- length(markers)
  n.ind <- length(individuals)

  res$data.tidy <- suppressWarnings(
    tibble::as_tibble(
      matrix(data = NA, nrow = n.markers, ncol = n.ind)
    ) %>%
      magrittr::set_colnames(x = ., value = individuals) %>%
      magrittr::set_rownames(x = ., value = markers) %>%
      data.table::as.data.table(x = ., keep.rownames = "MARKERS") %>%
      data.table::melt.data.table(
        data = .,
        id.vars = "MARKERS",
        variable.name = "INDIVIDUALS",
        value.name = "GT",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::select(-GT) %>%
      dplyr::mutate(
        MARKERS = factor(x = MARKERS,
                         levels = markers, ordered = TRUE),
        INDIVIDUALS = factor(x = INDIVIDUALS,
                             levels = individuals,
                             ordered = TRUE)) %>%
      dplyr::arrange(MARKERS, INDIVIDUALS) %>%
      dplyr::bind_cols(x)
  )

  if (wide) {
    if (wide.markers) {
      res$data.wide <- data.table::as.data.table(res$data.tidy) %>%
        data.table::dcast.data.table(
          data = .,
          formula = INDIVIDUALS ~ MARKERS,
          value.var = names(x)
        ) %>%
        tibble::as_data_frame(.)
    } else {
      res$data.wide <- data.table::as.data.table(res$data.tidy) %>%
        data.table::dcast.data.table(
          data = .,
          formula = MARKERS ~ INDIVIDUALS,
          value.var = names(x)
        ) %>%
        tibble::as_data_frame(.)
    }
  }
  if (!tidy) res$data.tidy <- NULL
  return(res)
}# End tidy2wide
