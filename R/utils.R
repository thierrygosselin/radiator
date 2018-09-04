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


# .onUnload <- function(libpath) {
#   library.dynam.unload("radiator", libpath)
# }


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


#' @title replace_by_na
#' @description Fast removal of NA
#' @rdname replace_by_na
#' @keywords internal
#' @export
replace_by_na <- function(data, what = ".") {
  replace(data, which(data == what), NA)
}#End replace_by_na

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


#' @title compute_maf
#' @description Compute MAF
#' @rdname compute_maf
#' @keywords internal
#' @export
compute_maf <- function(x, biallelic) {
  if (tibble::has_name(x, "GT_BIN") && biallelic) {
    x <- x %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        NN = as.numeric(2 * n()),
        PP = as.numeric(2 * length(GT_BIN[GT_BIN == 0])),
        PQ = as.numeric(length(GT_BIN[GT_BIN == 1])),
        QQ = as.numeric(2 * length(GT_BIN[GT_BIN == 2]))
      ) %>%
      # need this step because seen cases where the minor allele is not minor
      dplyr::mutate(
        PP = PP + PQ,
        QQ = QQ + PQ,
        PQ = NULL) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        NN_G = sum(NN),
        PP_G = sum(PP),
        QQ_G = sum(QQ)) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::mutate(
        ALT = dplyr::if_else(PP_G < QQ_G, PP, QQ),
        MAF_LOCAL = (ALT / NN),
        PP = NULL,
        # QQ = NULL,
        NN = NULL) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        ALT = dplyr::if_else(PP_G < QQ_G, PP_G, QQ_G),
        MAF_GLOBAL = (ALT / NN_G),
        ALT = NULL,
        PP_G = NULL,
        # QQ_G = NULL,
        NN_G = NULL) %>%
      dplyr::ungroup(.) %>%
      dplyr::rename(ALT_LOCAL = QQ, ALT_GLOBAL = QQ_G)
  } else {
    if (!tibble::has_name(x, "GT_VCF_NUC")) {
      x <- x %>%
        dplyr::select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(GT, 1, 3),
          A2 = stringi::stri_sub(GT, 4,6)
        ) %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID))

      maf.local <- x %>%
        dplyr::group_by(MARKERS, POP_ID, GT) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_LOCAL = n / n.al.tot, ALT_LOCAL = n) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, ALT_LOCAL)

      x <- x %>%
        dplyr::group_by(MARKERS, GT) %>%
        dplyr::tally(.) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_GLOBAL = n / n.al.tot, ALT_GLOBAL = n) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(MARKERS, MAF_GLOBAL, ALT_GLOBAL) %>%
        dplyr::left_join(maf.local, by = c("MARKERS")) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, ALT_LOCAL, MAF_GLOBAL, ALT_GLOBAL)
      maf.local <- NULL
    } else {
      x <- x %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_VCF_NUC) %>%
        tidyr::separate(
          data = .,
          col = GT_VCF_NUC, into = c("A1", "A2"),
          sep = "/",
          extra = "drop", remove = TRUE
        ) %>%
        tidyr::gather(
          data = ., key = ALLELE_GROUP, value = HAPLOTYPES,
          -dplyr::one_of(c("MARKERS", "INDIVIDUALS", "POP_ID"))) %>%
        dplyr::select(-ALLELE_GROUP) %>%
        dplyr::group_by(MARKERS, HAPLOTYPES, POP_ID) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, HAPLOTYPES), fill = list(n = 0)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(N_LOCAL = sum(n)) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(N_GLOBAL = sum(n)) %>%
        dplyr::arrange(MARKERS, POP_ID) %>%
        dplyr::group_by(MARKERS, POP_ID, HAPLOTYPES) %>%
        dplyr::mutate(MAF_LOCAL = n / N_LOCAL) %>%
        dplyr::group_by(MARKERS, HAPLOTYPES) %>%
        dplyr::mutate(
          ALT_GLOBAL = sum(n),
          MAF_GLOBAL = ALT_GLOBAL / N_GLOBAL,
          N_LOCAL = NULL,
          N_GLOBAL = NULL
        ) %>%
        dplyr::rename(ALT_LOCAL = n) %>%
        dplyr::ungroup(.)

      ref.info <- dplyr::distinct(x, MARKERS, HAPLOTYPES, MAF_GLOBAL) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::filter(MAF_GLOBAL == max(MAF_GLOBAL)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(REF = rep("REF", n()), MAF_GLOBAL = NULL) %>%
        dplyr::bind_rows(
          dplyr::distinct(x, MARKERS, HAPLOTYPES, MAF_GLOBAL) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::filter(MAF_GLOBAL == min(MAF_GLOBAL)) %>%
            dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
            dplyr::ungroup(.) %>%
            dplyr::mutate(REF = rep("MAF", n()), MAF_GLOBAL = NULL)
        )

      x <- dplyr::left_join(x, ref.info, by = c("MARKERS", "HAPLOTYPES")) %>%
        dplyr::mutate(REF = stringi::stri_replace_na(REF, replacement = "ALT"))
      ref.info <- NULL
    }
  }
  return(x)
}#End compute_maf

# update data.info
#' @title data_info
#' @description function generate tidy data main info
#' @rdname data_info
#' @keywords internal
#' @export
data_info <- function(x, print.info = FALSE) {

  if (tibble::has_name(x, "POP_ID")) {
    x.pop.ind <- dplyr::distinct(x, POP_ID, INDIVIDUALS)
    n.pop <- dplyr::n_distinct(x.pop.ind$POP_ID)
    n.ind <- dplyr::n_distinct(x.pop.ind$INDIVIDUALS)
  } else {
    n.pop <- NA_integer_
    n.ind <- NA_integer_
  }
  if (tibble::has_name(x, "MARKERS")) {
    n.snp <- dplyr::n_distinct(x$MARKERS)
  } else {
    n.snp <- NA_integer_
  }

  if (tibble::has_name(x, "LOCUS")) {
    n.locus <- dplyr::n_distinct(x$LOCUS)
  } else {
    n.locus <- NA_integer_
  }

  if (tibble::has_name(x, "CHROM")) {
    n.chrom <- dplyr::n_distinct(x$CHROM)
  } else {
    n.chrom <- NA_integer_
  }

  res <- list(
    n.pop = n.pop,
    n.ind = n.ind,
    n.chrom = n.chrom,
    n.locus = n.locus,
    n.snp = n.snp
  )
  if (print.info) {
    message("Number of chrom: ", res$n.chrom)
    message("Number of locus: ", res$n.locus)
    message("Number of SNPs: ", res$n.snp)
    message("Number of populations: ", res$n.pop)
    message("Number of individuals: ", res$n.ind)
  }
  return(res)
}


# update ind_total_reads
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


# interactive_question
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
