#' @name ibdg_fh
#' @title FH measure of IBDg
#' @description FH is a proxy mesure of IBDg based on the excess in the observed
#' number of homozygous genotypes within an individual,
#' relative to the mean number of homozygous genotypes expected under random mating
#' (Keller et al., 2011; Kardos et al., 2015; Hedrick & Garcia-Dorado, 2016).
#'
#' \strong{IBDg} is the realized proportion of the individual genome
#' that is identical by descent by reference to the current population
#' under hypothetical random mating
#' (Keller et al., 2011; Kardos et al., 2015; Hedrick & Garcia-Dorado, 2016).
#'
#' This function is using a modified version of the FH measure
#' (constructed using \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK}
#' \code{-het} option) described in (Keller et al., 2011; Kardos et al., 2015).
#'
#' The novelties are:
#'
#' \itemize{
#' \item \strong{population-wise:} the individual's observed homozygosity is
#' contrasted against the expected homozygosity.
#' Two estimates of the expected homozygosity are provided based
#' on the population and/or the overall expected homozygosity
#' averaged across markers.
#' \item \strong{tailored for RADseq:} instead of using the overall number
#' of markers, the population and the overall expected homozygosity
#' are averaged with the same markers the individual's are genotyped for.
#' This reduces the bias potentially introduced by comparing the individual's
#' observed homozygosity (computed from non-missing genotypes) with
#' an estimate computed with more markers found at the population or at the
#' overall level.
#' }
#'
#' The FH measure is also computed in
#' [stackr](https://github.com/thierrygosselin/stackr) \emph{summary_haplotypes}
#' function and [grur](https://github.com/thierrygosselin/grur)
#' \emph{missing_visualization} functions.
#' See \strong{note} below for the equations.

#' @inheritParams tidy_genomic_data

#' @param monomorphic.out (optional) Set by default to remove
#' monomorphic markers that might have avoided filters.
#' Default: \code{monomorphic.out = TRUE}.

#' @param common.markers (optional) Logical. The argument for common markers
#' between populations is set by default to maximize genome coverage of
#' individuals and populations.
#' Default: \code{common.markers = FALSE}

#' @param filename (optional) Name of the tidy data set,
#' written to the working directory.

#' @return A list is created with 6 objects (function call, tables, manhattan,
#' boxplot and distribution plot).
#' FH measure is on average negative when the parents are less related than
#' expected by random mating. The distribution \code{fh.distribution.plot}
#' should be centered around 0 in samples of non-inbred individuals.
#' The first table, \code{$fh}, gives the individual's value
#' while the second table, \code{$fh.stats}, show the population and overall averaged.

#' @note
#'
#' \strong{Modified FH:}
#' \deqn{F_{h_i} = \frac{\overline{Het}_{obs_{ij}} - \overline{Het}_{exp_j}}{\sum_{i}snp_{ij} - \overline{Het}_{exp_j}}}
#'
#' \strong{Individual Observed Heterozygosity averaged across markers:}
#' \deqn{\overline{Het}_{obs_i} = \frac{\sum_iHet_{obs_i}}{\sum_i{snp_i}}}
#'
#' \strong{Population expected Heterozygosity (under Hardy-Weinberg) and
#' tailored by averaging for each individual using his genotyped markers:}
#' #\deqn{\overline{Het}_{exp_j} = \frac{\sum_jHet_{exp_j}}{\sum_j{snp_j}}}


#' @examples
#' \dontrun{
#' # Using a  VCF file, the simplest for of the function:
#' fh <- ibdg_fh(
#' data = "batch_1.vcf",
#' strata = "strata.panda.tsv"
#' )
#' # To see what's inside the list
#' names(fh)
#' # To view the manhattan plot:
#' fh$fh.manhattan.plot
#' # To view the boxplot:
#' fh$fh.boxplot
#' # To view the distribution of FH values:
#' fh$fh.distribution.plot
#' }

#' @references Keller MC, Visscher PM, Goddard ME (2011)
#' Quantification of inbreeding due to distant ancestors and its detection
#'  using dense single nucleotide polymorphism data. Genetics, 189, 237–249.
#' @references Kardos M, Luikart G, Allendorf FW (2015)
#' Measuring individual inbreeding in the age of genomics: marker-based
#' measures are better than pedigrees. Heredity, 115, 63–72.
#' @references Hedrick PW, Garcia-Dorado A. (2016)
#' Understanding Inbreeding Depression, Purging, and Genetic Rescue.
#' Trends in Ecology and Evolution. 2016; 31: 940-952.

#' @export
#' @rdname ibdg_fh

#' @importFrom dplyr distinct rename arrange mutate select summarise group_by ungroup filter inner_join left_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_replace_all_regex
#' @importFrom utils count.fields
#' @importFrom readr read_tsv write_tsv
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom ape pcoa
#' @importFrom stats dist
#' @importFrom tibble data_frame
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid geom_histogram


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

ibdg_fh <- function(
  data,
  strata = NULL,
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
  filename = NULL,
  verbose = TRUE
) {
  if (verbose) {
    cat("#######################################################################\n")
    cat("########################### radiator::ibdg_fh ###########################\n")
    cat("#######################################################################\n")
    timing <- proc.time()
  }
  # manage missing arguments -----------------------------------------------------
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")

  # store function call
  function.call <- match.call()

  # import data ----------------------------------------------------------------
  input <- radiator::tidy_genomic_data(
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
    filename = filename,
    verbose = FALSE
  )

  if (!"MARKERS" %in% colnames(input) & "LOCUS" %in% colnames(input)) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }

  # population names if pop.levels/pop.labels were request
  # input <- radiator::change_pop_names(data = input, pop.levels = pop.labels, pop.labels = pop.labels)

  # Detect if biallelic --------------------------------------------------------
  biallelic <- radiator::detect_biallelic_markers(input)

  # IBDg computations ----------------------------------------------------------
  message("Genome-Wide Identity-By-Descent calculations using FH...")
  if (tibble::has_name(input, "GT_VCF") & biallelic) {
    # freq.full <- input %>%
    #   dplyr::filter(GT_VCF != "./.") %>%
    #   dplyr::group_by(MARKERS, POP_ID) %>%
    #   dplyr::summarise(
    #     N = n(),
    #     HOM_REF = length(GT_VCF[GT_VCF == "0/0"]),
    #     HOM_ALT = length(GT_VCF[GT_VCF == "1/1"]),
    #     HOM = HOM_REF + HOM_ALT,
    #     HET = length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])    ) %>%
    #   dplyr::mutate(
    #     FREQ_ALT = ((HOM_ALT * 2) + HET) / (2 * N),
    #     FREQ_REF = 1 - FREQ_ALT,
    #     HET_O = HET / N,
    #     HOM_O = HOM / N,
    #     HOM_REF_O = HOM_REF / N,
    #     HOM_ALT_O = HOM_ALT / N,
    #     HOM_E = (FREQ_REF^2) + (FREQ_ALT^2),
    #     # HET_E2 = 1 - HOM_E2,
    #     HET_E = 2 * FREQ_REF * FREQ_ALT
    #   )

    # Remove missing
    input <- dplyr::filter(.data = input, GT_VCF != "./.")

    freq <- input %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        N = n(),
        HOM_ALT = length(GT_VCF[GT_VCF == "1/1"]),
        HET = length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])
      ) %>%
      dplyr::mutate(
        FREQ_ALT = ((HOM_ALT * 2) + HET) / (2 * N),
        FREQ_REF = 1 - FREQ_ALT,
        HOM_E = (FREQ_REF^2) + (FREQ_ALT^2)
      ) #%>% dplyr::group_by(POP_ID) %>% dplyr::summarise(HOM_E = mean(HOM_E, na.rm = TRUE))

    hom.e <- dplyr::full_join(
      input,
      dplyr::select(.data = freq, MARKERS, POP_ID, HOM_E)
      , by = c("MARKERS", "POP_ID")
    ) %>%
      dplyr::select(MARKERS, POP_ID, INDIVIDUALS, HOM_E) %>%
      dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      dplyr::summarise(HOM_E = mean(HOM_E, na.rm = TRUE)) #%>% dplyr::group_by(POP_ID) %>% dplyr::summarise(HOM_E = mean(HOM_E, na.rm = TRUE))

    fh <- input %>%
      dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      dplyr::summarise(
        N = n(),
        HOM_REF = length(GT_VCF[GT_VCF == "0/0"]),
        HOM_ALT = length(GT_VCF[GT_VCF == "1/1"]),
        HOM = HOM_REF + HOM_ALT
        # HET = length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])
      ) %>%
      dplyr::mutate(
        # FREQ_ALT = ((HOM_ALT * 2) + HET) / (2 * N),
        # FREQ_REF = 1 - FREQ_ALT,
        # HET_O = HET / N,
        HOM_O = HOM / N #, HOM_REF_O = HOM_REF / N, HOM_ALT_O = HOM_ALT / N
      ) %>%
      dplyr::full_join(dplyr::select(.data = hom.e, INDIVIDUALS, POP_ID, HOM_E), by = c("POP_ID", "INDIVIDUALS")) %>%
      dplyr::mutate(FH = ((HOM_O - HOM_E)/(N - HOM_E))) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)

    ind.levels <- fh$INDIVIDUALS
    fh <- dplyr::mutate(.data = fh, INDIVIDUALS = factor(INDIVIDUALS, levels = ind.levels, ordered = TRUE))
  } else {
    # not biallelic
    input.alleles <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
      dplyr::filter(GT != "000000") %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(GT, 1, 3),
        A2 = stringi::stri_sub(GT, 4,6)
      ) %>%
      dplyr::select(-GT)

    freq <- input.alleles %>%
      tidyr::gather(data = ., key = ALLELE_GROUP, value = ALLELES, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
      dplyr::group_by(MARKERS, ALLELES, POP_ID) %>%
      dplyr::tally(.) %>%
      dplyr::ungroup(.) %>%
      tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, ALLELES), fill = list(n = 0)) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::mutate(
        FREQ = n/sum(n),
        HOM_E = FREQ^2
      ) %>%
      dplyr::summarise(HOM_E = sum(HOM_E)) #%>% dplyr::mutate(HET_E = 1 - HOM_E)

    hom.e <- dplyr::full_join(
      dplyr::filter(.data = input, GT != "000000"),
      dplyr::select(.data = freq, MARKERS, POP_ID, HOM_E)
      , by = c("MARKERS", "POP_ID")
    ) %>%
      dplyr::select(MARKERS, POP_ID, INDIVIDUALS, HOM_E) %>%
      dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      dplyr::summarise(HOM_E = mean(HOM_E, na.rm = TRUE)) #%>% dplyr::group_by(POP_ID) %>% dplyr::summarise(HOM_E = mean(HOM_E, na.rm = TRUE))

    fh <- input.alleles %>%
      dplyr::group_by(POP_ID, INDIVIDUALS) %>%
      dplyr::summarise(
        N = n(),
        HOM = length(INDIVIDUALS[A1 == A2])
      ) %>%
      dplyr::mutate(HOM_O = HOM / N) %>%
      dplyr::full_join(dplyr::select(.data = hom.e, INDIVIDUALS, POP_ID, HOM_E), by = c("POP_ID", "INDIVIDUALS")) %>%
      dplyr::mutate(FH = ((HOM_O - HOM_E)/(N - HOM_E))) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)

    ind.levels <- fh$INDIVIDUALS
    fh <- dplyr::mutate(.data = fh, INDIVIDUALS = factor(INDIVIDUALS, levels = ind.levels, ordered = TRUE))
  }

  # FH statistics per pop
  fh.stats <- fh %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(FH = mean(FH))

  # per pop and overall combined
  fh.stats <- tibble::add_row(
    .data = fh.stats,
    POP_ID = "OVERALL",
    FH = unlist(dplyr::summarise(.data = fh.stats, FH = mean(FH)))
    )

  # plots ----------------------------------------------------------------------
  message("Generating plots")
  # manhattan
  fh.manhattan.plot <- ggplot2::ggplot(data = fh, ggplot2::aes(x = INDIVIDUALS, y = FH, colour = POP_ID)) +
    ggplot2::geom_jitter() +
    ggplot2::labs(y = "Individual IBDg (FH)") +
    ggplot2::labs(x = "Individuals") +
    ggplot2::labs(colour = "Populations") +
    # theme_minimal() +
    ggplot2::theme_classic() +
    # theme_dark() +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
    )
  # fh.manhattan.plot

  fh.boxplot <- ggplot2::ggplot(data = fh, ggplot2::aes(x = POP_ID, y = FH)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(y = "Individual IBDg (FH)") +
    ggplot2::labs(x = "Populations") +
    # theme_bw() +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
    )
  # fh.boxplot


  # Histogram
  fh.distribution.plot <- ggplot2::ggplot(data = fh, ggplot2::aes(x = FH)) +
    ggplot2::geom_histogram() +
    ggplot2::labs(x = "Individual IBDg (FH)") +
    ggplot2::labs(y = "Markers (number)") +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  # fh.distribution.plot


  # Results --------------------------------------------------------------------
  if (verbose) {
    timing <- proc.time() - timing
    message(stringi::stri_join("Computation time: ", round(timing[[3]]), " sec"))
    cat("############################## completed ##############################\n")
  }
  res = list(call = function.call, fh = fh, fh.stats = fh.stats, fh.manhattan.plot = fh.manhattan.plot, fh.boxplot = fh.boxplot, fh.distribution.plot = fh.distribution.plot)
  return(res)
}
