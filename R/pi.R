## Nucleotide diversity

#' @title Nucleotide diversity
#' @description Calculates the nucleotide diversity (Nei & Li, 1979).
#'
#' To get an estimate with the consensus reads, use the
#' function \emph{summary_haplotypes} found in the package
#' [stackr](https://github.com/thierrygosselin/stackr). The estimate in
#' \emph{summary_haplotypes} integrates the consensus markers found in
#' [STACKS](http://catchenlab.life.illinois.edu/stacks/)
#' populations.haplotypes.tsv file.
#' Both radiator and stackr functions requires \code{stringdist} package.
#'
#' The \code{read.length} argument below is used directly in the calculations.
#' To be correctly estimated, the reads obviously need to be of identical size...


#' @inheritParams tidy_genomic_data

#' @param monomorphic.out (optional) Set by default to remove
#' monomorphic markers that might have avoided filters.
#' Default: \code{monomorphic.out = TRUE}.

#' @param common.markers (optional) Logical. The argument for common markers
#' between populations is set by default to maximize genome coverage of
#' individuals and populations.
#' Default: \code{common.markers = FALSE}

#' @param read.length (number) The length in nucleotide of your reads
#' (e.g. \code{read.length = 100}).

# @importFrom stringdist stringdist
#' @importFrom utils combn count.fields
#' @importFrom stats lm na.omit
#' @importFrom stringi stri_replace_all_fixed stri_replace_na stri_join stri_count_fixed
#' @importFrom tibble as_data_frame data_frame add_column add_row
#' @importFrom dplyr select rename n_distinct distinct mutate summarise group_by ungroup arrange left_join full_join semi_join anti_join bind_rows bind_cols if_else
#' @importFrom readr write_tsv read_tsv
#' @importFrom tidyr separate gather
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid stat_smooth
#' @importFrom purrr flatten_chr map_df flatten_dbl

#' @return The function returns a list with the function call and:
#' \enumerate{
#' \item $pi.individuals: the pi estimated for each individual
#' \item $pi.populations: the pi statistics estimated per populations and overall.
#' \item $boxplot.pi: showing the boxplot of Pi for each populations and overall.
#' }
#' use $ to access each #' objects in the list.

#' @examples
#' \dontrun{
#' require(stringdist)
#' # The simplest way to run the function:
#' sum <- radiator::pi(
#' data = "batch_1.vcf",
#' strata = "strata_brook_charr.tsv",
#' read.length = 90)
#' }


#' @rdname pi
#' @export

#' @references Nei M, Li WH (1979)
#' Mathematical model for studying genetic variation in terms of
#' restriction endonucleases.
#' Proceedings of the National Academy of Sciences of
#' the United States of America, 76, 5269–5273.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

pi <- function(
  data,
  strata,
  read.length,
  monomorphic.out = TRUE,
  common.markers = FALSE,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  max.marker = NULL,
  snp.ld = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    stop("stringdist needed for this function to work
         Install with install.packages('stringdist')", call. = FALSE)
  }

  cat("#######################################################################\n")
  cat("############################## radiator::pi #############################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  res <- list() # to store results

  # manage missing arguments -----------------------------------------------------
  if (missing(data)) stop("Input file missing")
  if (missing(read.length)) stop("read.length argument is required")
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(
      pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

  # store function call
  res$function.call <- match.call()


  # import data ----------------------------------------------------------------
  message("Importing and tidying the data...")
  input <- suppressMessages(radiator::tidy_genomic_data(
    data = data,
    vcf.metadata = FALSE,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    max.marker = max.marker,
    snp.ld = snp.ld,
    common.markers = common.markers,
    strata = strata,
    pop.select = pop.select,
    pop.levels = pop.labels,
    pop.labels = pop.labels,
    filename = NULL,
    verbose = FALSE))

  if (!"MARKERS" %in% colnames(input) & "LOCUS" %in% colnames(input)) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }

  # Nei & Li 1979 Nucleotide Diversity -----------------------------------------
  message("Nucleotide diversity (Pi):")
  message("    Read length used: ", read.length)

  # Pi: by individuals----------------------------------------------------------
  message("    Pi calculations: individuals...")
  separate_gt <- function(x) {
    res <- x  %>%
      tidyr::separate(
        col = GT, into = c("ALLELE1", "ALLELE2"),
        sep = 3, extra = "drop", remove = TRUE
      ) %>%
      dplyr::mutate(ALLELE2 = ifelse(is.na(ALLELE2), ALLELE1, ALLELE2))
    return(res)
  }#End separate_haplo

  # keep genotyped info
  input <- dplyr::filter(input, GT != "000000")

  # split genotypes in 2
  n.row <- nrow(input)
  split.vec <- as.integer(floor((parallel.core * 20 * (1:n.row - 1) / n.row) + 1))
  n.row <- NULL

  input <- split(x = input, f = split.vec) %>%
    .radiator_parallel_mc(
      X = ., FUN = separate_gt, mc.cores = parallel.core) %>%
    dplyr::bind_rows(.)

  split.vec <- NULL

  res$pi.individuals <- input %>%
      dplyr::mutate(
        PI = (stringdist::stringdist(a = ALLELE1, b = ALLELE2, method = "hamming"))/read.length
      ) %>%
      dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      dplyr::summarise(PI = mean(PI))

  # Pi function ----------------------------------------------------------------
  pi <- function(data, read.length) {
    y <- dplyr::select(data, ALLELES) %>% purrr::flatten_chr(.)

    if (length(unique(y)) <= 1) {
      pi <- tibble::data_frame(PI = as.numeric(0))
    } else {

      #1 Get all pairwise comparison
      allele.pairwise <- utils::combn(unique(y), 2)

      #2 Calculate pairwise nucleotide mismatches
      pairwise.mismatches <- apply(allele.pairwise, 2, function(z) {
        stringdist::stringdist(a = z[1], b = z[2], method = "hamming")
      })

      #3 Calculate allele frequency
      allele.freq <- table(y)/length(y)

      #4 Calculate nucleotide diversity from pairwise mismatches and allele frequency
      pi <- apply(allele.pairwise, 2, function(y) allele.freq[y[1]] * allele.freq[y[2]])
      pi <- tibble::data_frame(PI = sum(pi * pairwise.mismatches) / read.length)
    }
    return(pi)
  }#End pi

  pi_pop <- function(data, read.length, parallel.core) {
    pop <- unique(data$POP_ID)
    message("    Pi calculations for pop: ", pop)
    # data <- df.split.pop[["DD"]]
    # data <- df.split.pop[["SKY"]]

    pi.pop <- data %>%
      split(x = ., f = .$MARKERS) %>%
      .radiator_parallel(
        X = .,
        FUN = pi,
        mc.cores = parallel.core,
        read.length = read.length
      ) %>%
      dplyr::bind_rows(.) %>%
      dplyr::summarise(PI_NEI = mean(PI)) %>%
      tibble::add_column(.data = ., POP_ID = pop, .before = "PI_NEI")

    return(pi.pop)
  }#End pi_pop

  # Pi: by pop------------------------------------------------------------------
  message("    Pi calculations: populations...")
  input <- dplyr::select(input, POP_ID, INDIVIDUALS, MARKERS, ALLELE1, ALLELE2) %>%
    tidyr::gather(ALLELE_GROUP, ALLELES, -c(POP_ID, INDIVIDUALS, MARKERS))

  res$pi.populations <- input %>%
    split(x = ., f = .$POP_ID) %>%
    purrr::map_df(
      .x = ., .f = pi_pop,
      read.length = read.length, parallel.core = parallel.core
    )

  # Pi: overall  ---------------------------------------------------------------
  message("    Pi calculations: overall")
  res$pi.populations <- tibble::add_row(
    .data = res$pi.populations,
    POP_ID = "OVERALL",
    PI_NEI = input %>%
      split(x = ., f = .$MARKERS) %>%
      .radiator_parallel(
        X = .,
        FUN = pi,
        mc.cores = parallel.core,
        read.length = read.length
      ) %>%
      dplyr::bind_rows(.) %>%
      dplyr::summarise(PI_NEI = mean(PI)) %>% purrr::flatten_dbl(.)
  )

  # figure ---------------------------------------------------------------------
  res$boxplot.pi <- dplyr::filter(res$pi.individuals, POP_ID != "OVERALL") %>%
    ggplot2::ggplot(data = ., ggplot2::aes(x = factor(POP_ID), y = PI, na.rm = TRUE)) +
    ggplot2::geom_violin(trim = F) +
    ggplot2::geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
    ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
    ggplot2::labs(x = "Sampling sites") +
    ggplot2::labs(y = "Individual nucleotide diversity (Pi)") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )

  # results --------------------------------------------------------------------
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)
}
