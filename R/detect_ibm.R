#' @rdname detect_ibm
#' @title Detect identity-by-missingness structure
#' @description Computes a binary identity-by-missingness matrix directly from a
#' GDS object and returns a heatmap together with marker- and individual-level
#' missingness summaries.
#'
#' @param gds A \code{SeqVarGDSClass} object or a path to a GDS file.
#'
#' @param strata (optional) Path to a strata file with a minimum of 2 columns:
#' \code{INDIVIDUALS} and the grouping column selected in \code{strata.select}.
#' Default: \code{strata = NULL}.
#'
#' @param strata.select (character) Column used to facet individuals in the
#' heatmap.
#' Default: \code{strata.select = "POP_ID"}.
#'
#' @param sort.individuals (character) Sorting method for individuals.
#' Choices are: \code{"input"}, \code{"missingness"}, \code{"strata"}.
#' Default: \code{sort.individuals = "input"}.
#'
#' @param sort.markers (character) Sorting method for markers.
#' Choices are: \code{"input"}, \code{"missingness"}, \code{"position"}.
#' Default: \code{sort.markers = "input"}.
#'
#' @param sample.max (optional, integer) Maximum number of individuals plotted.
#' Default: \code{sample.max = NULL}.
#'
#' @param marker.max (optional, integer) Maximum number of markers plotted.
#' Default: \code{marker.max = 50000}.
#'
#' @param facet (logical) Should the heatmap be faceted by \code{strata.select}?
#' Default: \code{facet = TRUE}.
#'
#' @param parallel.core (optional, integer) Number of cores.
#' Default: \code{parallel::detectCores() - 1}.
#'
#' @param ... (optional) Further arguments passed to advanced mode.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{heatmap}: a \code{ggplot2} object
#'   \item \code{ibm.matrix}: integer matrix, 0 = missing, 1 = genotyped
#'   \item \code{individuals.missingness}: tibble with per-individual missingness
#'   \item \code{markers.missingness}: tibble with per-marker missingness
#'   \item \code{plot.data}: long tibble used to generate the heatmap
#' }
#'
#' @export
detect_ibm <- function(
    gds,
    strata = NULL,
    strata.select = "POP_ID",
    sort.individuals = "input",
    sort.markers = "input",
    sample.max = NULL,
    marker.max = 50000,
    facet = TRUE,
    parallel.core = parallel::detectCores() - 1,
    ...
) {

  if (missing(gds)) rlang::abort("Input gds missing")

  dots <- radiator::radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = NULL,
    deprecated = NULL,
    verbose = FALSE
  )

  valid.sort.ind <- c("input", "missingness", "strata")
  valid.sort.mark <- c("input", "missingness", "position")

  if (!sort.individuals %in% valid.sort.ind) {
    rlang::abort(
      message = stringi::stri_join(
        "sort.individuals must be one of: ",
        stringi::stri_join(valid.sort.ind, collapse = ", ")
      )
    )
  }

  if (!sort.markers %in% valid.sort.mark) {
    rlang::abort(
      message = stringi::stri_join(
        "sort.markers must be one of: ",
        stringi::stri_join(valid.sort.mark, collapse = ", ")
      )
    )
  }

  gds.opened.here <- FALSE

  if (inherits(gds, "character")) {
    gds <- SeqArray::seqOpen(gds)
    gds.opened.here <- TRUE
  }

  if (!inherits(gds, "SeqVarGDSClass")) {
    rlang::abort("gds must be a SeqVarGDSClass object or a path to a GDS file")
  }

  if (gds.opened.here) {
    on.exit(SeqArray::seqClose(gds), add = TRUE)
  }

  on.exit(SeqArray::seqResetFilter(gds, verbose = FALSE), add = TRUE)

  sample.ids <- SeqArray::seqGetData(gds, "sample.id")
  marker.ids <- SeqArray::seqGetData(gds, "variant.id")

  chrom <- tryCatch(
    SeqArray::seqGetData(gds, "chromosome"),
    error = function(e) NULL
  )

  pos <- tryCatch(
    SeqArray::seqGetData(gds, "position"),
    error = function(e) NULL
  )

  n.samples <- length(sample.ids)
  n.markers <- length(marker.ids)

  individuals.missing.count <- integer(n.samples)

  markers.missing.count <- SeqArray::seqApply(
    gds,
    "genotype",
    FUN = function(x) {
      sum(colSums(is.na(x)) == nrow(x))
    },
    margin = "by.variant",
    as.is = "integer"
  )

  invisible(
    SeqArray::seqApply(
      gds,
      "genotype",
      FUN = function(x) {
        individuals.missing.count <<- individuals.missing.count +
          as.integer(colSums(is.na(x)) == nrow(x))
        NULL
      },
      margin = "by.variant",
      as.is = "none"
    )
  )

  individuals.missingness <- tibble::tibble(
    INDIVIDUALS = sample.ids,
    MISSING_GENOTYPE = individuals.missing.count,
    MARKER_NUMBER = n.markers,
    MISSING_GENOTYPE_PROP = MISSING_GENOTYPE / MARKER_NUMBER,
    PERCENT = round(MISSING_GENOTYPE_PROP * 100, 2)
  )

  markers.missingness <- tibble::tibble(
    MARKERS = marker.ids,
    CHROM = chrom,
    POS = pos,
    MISSING_GENOTYPE = markers.missing.count,
    INDIVIDUALS_NUMBER = n.samples,
    MISSING_GENOTYPE_PROP = MISSING_GENOTYPE / INDIVIDUALS_NUMBER,
    PERCENT = round(MISSING_GENOTYPE_PROP * 100, 2)
  )

  if (is.null(strata)) {
    strata.df <- tibble::tibble(
      INDIVIDUALS = sample.ids,
      POP_ID = "all"
    )
    if (strata.select != "POP_ID") strata.df[[strata.select]] <- "all"
  } else {
    strata.df <- radiator::read_strata(
      strata = strata,
      tidy = FALSE,
      verbose = FALSE
    )

    if (!"INDIVIDUALS" %in% names(strata.df)) {
      rlang::abort("strata must contain an INDIVIDUALS column")
    }

    if (!strata.select %in% names(strata.df)) {
      rlang::abort(
        message = stringi::stri_join(
          "Column '", strata.select, "' not found in strata"
        )
      )
    }

    strata.df <- dplyr::semi_join(
      strata.df,
      tibble::tibble(INDIVIDUALS = sample.ids),
      by = "INDIVIDUALS"
    )
  }

  sample.order <- sample.ids
  marker.order <- marker.ids

  if (sort.individuals == "missingness") {
    sample.order <- individuals.missingness %>%
      dplyr::arrange(MISSING_GENOTYPE_PROP, INDIVIDUALS) %>%
      dplyr::pull(INDIVIDUALS)
  }

  if (sort.individuals == "strata") {
    sample.order <- strata.df %>%
      dplyr::select(INDIVIDUALS, dplyr::all_of(strata.select)) %>%
      dplyr::distinct() %>%
      dplyr::arrange(.data[[strata.select]], INDIVIDUALS) %>%
      dplyr::pull(INDIVIDUALS)
  }

  if (sort.markers == "missingness") {
    marker.order <- markers.missingness %>%
      dplyr::arrange(MISSING_GENOTYPE_PROP, MARKERS) %>%
      dplyr::pull(MARKERS)
  }

  if (sort.markers == "position" && !is.null(chrom) && !is.null(pos)) {
    marker.order <- markers.missingness %>%
      dplyr::arrange(CHROM, POS, MARKERS) %>%
      dplyr::pull(MARKERS)
  }

  if (!is.null(sample.max)) {
    sample.order <- sample.order[seq_len(min(sample.max, length(sample.order)))]
  }

  if (!is.null(marker.max)) {
    marker.order <- marker.order[seq_len(min(marker.max, length(marker.order)))]
  }

  SeqArray::seqSetFilter(
    gds,
    sample.id = sample.order,
    variant.id = marker.order,
    verbose = FALSE
  )

  ibm.blocks <- SeqArray::seqApply(
    gds,
    "genotype",
    FUN = function(x) {
      as.integer(!(colSums(is.na(x)) == nrow(x)))
    },
    margin = "by.variant",
    as.is = "list"
  )

  ibm.matrix <- do.call(cbind, ibm.blocks)
  rownames(ibm.matrix) <- SeqArray::seqGetData(gds, "sample.id")
  colnames(ibm.matrix) <- SeqArray::seqGetData(gds, "variant.id")

  plot.data <- as.data.frame(ibm.matrix) %>%
    tibble::rownames_to_column(var = "INDIVIDUALS") %>%
    tidyr::pivot_longer(
      cols = -INDIVIDUALS,
      names_to = "MARKERS",
      values_to = "GT_MISSING_BINARY"
    ) %>%
    dplyr::left_join(
      dplyr::select(strata.df, INDIVIDUALS, dplyr::all_of(strata.select)),
      by = "INDIVIDUALS"
    ) %>%
    dplyr::mutate(
      Missingness = dplyr::if_else(
        GT_MISSING_BINARY == 0L,
        "missing",
        "genotyped"
      ),
      Missingness = factor(Missingness, levels = c("genotyped", "missing")),
      INDIVIDUALS = factor(INDIVIDUALS, levels = sample.order),
      MARKERS = factor(MARKERS, levels = marker.order)
    )

  heatmap <- ggplot2::ggplot(
    plot.data,
    ggplot2::aes(
      x = INDIVIDUALS,
      y = MARKERS,
      fill = Missingness
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c("grey", "black")) +
    ggplot2::labs(
      x = "Individuals",
      y = "Markers",
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        size = 10, family = "Helvetica", face = "bold"
      ),
      axis.title.y = ggplot2::element_text(
        size = 10, family = "Helvetica", face = "bold"
      ),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(
        size = 10, family = "Helvetica", face = "bold"
      )
    )

  if (facet) {
    heatmap <- heatmap +
      ggplot2::facet_grid(
        stats::as.formula(stringi::stri_join("~", strata.select)),
        scales = "free_x",
        space = "free_x"
      )
  }

  return(list(
    heatmap = heatmap,
    ibm.matrix = ibm.matrix,
    individuals.missingness = individuals.missingness,
    markers.missingness = markers.missingness,
    plot.data = plot.data
  ))
}#END detect_ibm
