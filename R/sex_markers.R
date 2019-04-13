
#' @name sexy_markers
#' @title sexy_markers finds sex-linked markers and re-assigns sex
#'
#' @description This function identifies sex-linked markers putatively located on
#' homogametic or heterogametic chromosomes and re-assign the sex in a dataset
#' based on findings. The function work best in: DArT silico (counts) >
#' DArT counts or RADseq with allele read depth > DArT silico (genotypes) >
#' RADseq (genotypes) and DArT (1-row, 2-rows genotypes).

#' @param data (object or file) DArT file \code{.csv or .tsv}, VCF file \code{.vcf},
#' GDS file or object (\code{.gds}).
#' See Input data section data for more details.

#' @param silicodata (optional, file) A silico DArT file (.csv or .tsv). This can be
#' count or genotyped data. Note that both \code{data} and \code{silicodata} can
#' be used at the same time.
#' Default: \code{silicodata = NULL}.

#' @param boost.analysis (logical, optional) This method uses machine learning
#' approaches to find sex markers and re-assign samples in sex group.
#' The approach is currently been tested and will be available for uses soon.#'
#' Default: \code{boost.analysis = FALSE}.


#' @param strata (file) A tab delimited file with a minimum of
#' 2 columns \code{INDIVIDUALS, STRATA} for VCF files and 3 columns for DArT files
#' \code{TARGET_ID, INDIVIDUALS, STRATA}.
#' \itemize{
#' \item \code{TARGET_ID:} it's the column header of the DArT file.
#' \item \code{STRATA:} is the grouping column, here the sex info.
#' 3 values: \code{M} for male, \code{F} for female and \code{U} for unknown.
#' Anything else is converted to \code{U}.
#' \item \code{INDIVIDUALS:} Is how you want your samples to be named.
#' }
#' Default: \code{strata = NULL}, the function will look for sex info in the
#' tidy data or GDS data (individuals.meta section).
#' You can easily build the strata file by starting with the output of these
#' functions: \code{\link{extract_dart_target_id}} and \code{\link{extract_individuals_vcf}}

#' @param coverage.thresholds (optional, integer) The minimum coverage required
#' to account for the genotype.
#' Default: \code{coverage.thresholds = 1}.
#'
#'
#' @param filters (optional, logical) When \code{filters = TRUE}, the data will
#' be filtered for
#' monomorphic loci, missingness of individuals and heterozygosity of individuals.
#' CAUTION: we advice to use these filter, since not filtering or filtering too
#' stingently will results in
#' false positve or false negative detections.
#' Default: \code{filters = TRUE}.
#'
#' @param interactive.filter (optional, logical)
#' Default: \code{interactive.filter = TRUE}.
#'
#' @inheritParams radiator_common_arguments
#'
#'
#' @details
#' This function takes DArT and RAD data to find markers that have a specific
#'  pattern that is linked to sex.
#' The function hypothesizes the presence of sex-chromosomes in you
#' species/population. The tests are designed to identify markers that are
#' located on putative heterogametic (Y or W) or homogametic (X or Z) chromosomes.
#' \emph{Note:} Violating Assumptions or Prerequisites (see below) can lead to
#' false positive or the absence of detection of sex-linked markers.

#' @section Assumptions:
#' \enumerate{
#' \item \strong{Genetic Sex Determination System} over a e.g. environmental-sex-determination system.
#' \item \strong{Genome coverage:} restriction sites randomly spread throughout the whole genome.
#' \item \strong{Mutations:} Processes such as sex-specific mutations
#' in the restriction sites could lead to false positive results.
#' \item \strong{Deletions/duplications:} Processes such as
#' sex-specific deletions or duplications could lead to false positive results.
#' \item \strong{Homology:} The existence of homologous sequences between the
#' homogametic and heterogametic chromosomes could lead to false negative
#' results.
#' \item \strong{Absence of population signal}
#' }

#' @section Prerequisites:
#' \enumerate{
#' \item \strong{Sample size:} Ideally, the data must have enough individuals (n > 100).
#' \item \strong{Batch effect:} Sex should be randomized on lanes/chips during sequencing.
#' \item \strong{Sex ratio:} Dataset with equal ratio work best.
#' \item \strong{Genotyping rate}: for DArT data, if the minimum call rate is
#' > 0.5 ask DArT to lower their filtering threshold.
#' RADseq data, lower markers missingness thresholds during filtering
#' (e.g. stacks \code{r} and \code{p}).
#' \item \strong{Identity-by-Missingness:} Absence of artifactual pattern
#' of missingness (\href{https://github.com/thierrygosselin/grur}{see missing visualization})
#' \item \strong{Low genotyping error rate:} see \code{\link{detect_het_outliers}}
#' and \href{https://github.com/eriqande/whoa}{whoa}.
#' \item \strong{Low heterozygosity miscall rate:} see \code{\link{detect_het_outliers}} and
#' \href{https://github.com/eriqande/whoa}{whoa}.
#' \item \strong{Absence of pattern of heterozygosity driven by missingness:}
#' see \code{\link{detect_mixed_genomes}}.
#' \item \strong{Absence of paralogous sequences in the data}.
#' }


#' @section Sex methods:
#' \strong{Heterogametic sex-markers:}
#' \itemize{
#' \item \strong{Presence/Absence method}: To identify markers on Y or W chromosomes,
#' we look at the presence or
#' absence of a marker between females and males. More specifically, if a marker
#' is always present in males but never in females, they are putatively located
#' in the Y-chromosome; and vice versa for the W-linked markers.
#' }
#' \strong{Homogametic sex-markers:}
#' We have two different methods to identify markers on X or Z chromosomes:
#' \itemize{
#' \item \strong{Heterozygosity method:} By looking at the heterozygosity of a marker between sexes, we can
#' identify markers that are always homozygous in one sex (e.g. males for an
#' XY system), while exhibiting an intermediate range of heterozygosity in the
#' other sex (0.1 - 0.5).
#' \item \strong{Coverage method:} If the data includes count information, this function will look for
#' markers that have double the number of counts for either of the sexes.
#' For example if an XY/XX system is present, females are expected to have
#' double the number of counts for markers on the X chromosome.
#' }

#' @section Life cycle:
#'
#' Machine Learning approaches (Random Forest and Extreme Gradient Boosting Trees)
#' are currently been tested. They usually show a lower discovery rate but tend to
#' perform better with new samples.
#'
#'

#' @seealso Eric Anderson's \href{https://github.com/eriqande/whoa}{whoa} package.

#' @export
#' @rdname sexy_markers

#' @return A list with ...

#' @examples
#' \dontrun{

#' # The minumum
#' sex.markers <- radiator::sexy_markers(
#'     data = "shark.dart.data.csv",
#'     strata = "shark.strata.tsv")

#' # This will use the default: interactive version and a list is created and to view the sex markers
#' }


#' @author Floriaan Devloo-Delva \email{Floriaan.Devloo-Delva@@csiro.au} and with help from
#' Thierry Gosselin \email{thierrygosselin@@icloud.com}

sexy_markers <- function(data,
                         silicodata = NULL,
                         boost.analysis = FALSE,
                         strata = NULL,
                         coverage.thresholds = 1,
                         filters = TRUE,
                         interactive.filter = TRUE,
                         parallel.core = parallel::detectCores() - 1,
                         ...
) {

  # Test
  # library(radiator)
  # setwd("C:/Users/fdevloo/OneDrive - University of Tasmania/PhD/4.Sex markers/2.Code/Gosselin Function_Radiator")
  # silicodata = NULL
  # coverage.thresholds = 1 ##NEEDS TO BE SET TO 1 if genotype data
  # boost.analysis = FALSE
  # filters = TRUE
  # interactive.filter = TRUE
  # parallel.core = 1
  # species = "shark"
  # population = "Australia"
  # tau = 0.03
  # threshold.y.markers = NULL
  # threshold.y.silico.markers = NULL
  # sex.id.input = NULL
  # threshold.x.markers.qr = NULL

  # parallel.core = parallel::detectCores() - 1

  # setwd("/Users/thierry/Dropbox/partage/Sex-marker/1.Data/G.galeus")
  # silicodata = "SchoolShark_silico_counts.csv"
  # data = "SchoolShark_SNP_counts.csv"
  # strata = "School_strata_Counts.tsv"

  # data = "SchoolShark_SNP_genotype_2Row.csv"
  # strata = "School_strata_2ROW.tsv"
  # silicodata = "SchoolShark_silico_counts.csv"
  # data = "Ggar_Report_DGl17-2561_4_moreOrders_SNP_singlerow_2.csv"
  # strata = "Ggar_strata.tsv"
  # data = "Silvertip.vcf"
  # strata = "silvertip_strata.tsv"
  # data = "Dusky_SNP_genotype_1row.csv"
  # strata = "Dusky_strata.tsv"
  # data = "example_vcf2dadi_ferchaud_2015.vcf"
  # strata = "strata.stickleback.tsv"


  if (boost.analysis) message("Under construction: come back next week... ")
  # Cleanup---------------------------------------------------------------------
  verbose <- TRUE
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose)
    message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose)
    message("\nComputation time, overall: ", round(timing[[3]]), " sec"),
    add = TRUE)
  on.exit(
    if (verbose)
      cat(
        "############################ sex markers completed #############################\n"
      ),
    add = TRUE
  )

  # required packages ----------------------------------------------------------
  if (!"quantreg" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install quantreg for this option:\n
                 install.packages("quantreg")
                 ')
  }
  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("species", "population", "tau",
                "threshold.y.markers", "threshold.y.silico.markers",
                "sex.id.input", "threshold.x.markers.qr"),
    verbose = FALSE
  )

  # Folders---------------------------------------------------------------------
  wd <- path.folder <- radiator::generate_folder(
    f = NULL,
    rad.folder = "sexy_markers",
    internal = FALSE,
    prefix_int = FALSE,
    file.date = file.date,
    verbose = TRUE
  )

  # write the dots file
  write_rad(
    data = rad.dots,
    path = path.folder,
    filename = stringi::stri_join("radiator_sexy_markers_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = FALSE,
    verbose = TRUE
  )


  # Detect format --------------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  if (!data.type %in% c("SeqVarGDSClass", "gds.file", "dart", "vcf.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }

  # GDS file and object --------------------------------------------------------
  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (!"SeqVarTools" %in% utils::installed.packages()[, "Package"]) {
      rlang::abort(
        'Please install SeqVarTools for this option:\n
        install.packages("BiocManager")
        BiocManager::install("SeqVarTools")'
      )
    }

    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = verbose)
      data.type <- "SeqVarGDSClass"
    }

    # data <- gds2tidy(gds = data, parallel.core = parallel.core)
    # data.type <- "tbl_df"
  }

  # VCF files ------------------------------------------------------------------
  if (data.type == "vcf.file") {
    data <- radiator::read_vcf(
      data = data,
      strata = strata,
      path.folder = path.folder,
      internal = TRUE,
      parallel.core = parallel.core,
      verbose = FALSE
    )
  }


  # DArT files ------------------------------------------------------------------
  if (data.type == "dart") {
    data <- radiator::read_dart(
      data = data,
      strata = strata,
      tidy.dart = FALSE,
      parallel.core = parallel.core,
      path.folder = path.folder,
      internal = TRUE,
      verbose = FALSE
    )
  }


  # clean the strata -----------------------------------------------------------
  strata <- radiator::extract_individuals_metadata(gds = data)
  # if (rlang::has_name(strata, "POP_ID")) data %<>% dplyr::rename(STRATA = POP_ID)


  # Detect source --------------------------------------------------------------
  data.source <- radiator::extract_data_source(gds = data)
  if (!data.type %in% c("SeqVarGDSClass", "gds.file", "dart", "vcf.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }


  # clean strata----------------------------------------------------------------
  strata %<>%
    dplyr::mutate(
      STRATA = stringi::stri_trans_toupper(str = STRATA),
      STRATA = stringi::stri_sub(str = STRATA, from = 1, to = 1),
      STRATA = replace(x = STRATA, !STRATA %in% c("F", "M"), "U")
    )
  #checks
  strata.groups <- unique(sort(strata$STRATA))
  if (length(strata.groups) > 3 || length(strata.groups) < 2) {
    rlang::abort("The strata requires a minimum of 2 groups: F and M")
  }
  if (length(strata.groups) == 2 &&
      all(c("F", "M") %in% strata.groups) == FALSE) {
    rlang::abort("The strata requires a minimum of 2 groups: F and M")
  }

  # Generate new strata and write to disk
  # strata <- generate_strata(data)
  readr::write_tsv(x = strata,
                   path = file.path(path.folder, "strata.cleaned.tsv"))

  #update strata in GDS



  # Filter monomorphic ---------------------------------------------------------
  data <- radiator::filter_monomorphic(
    data = data,
    parallel.core = parallel.core,
    verbose = FALSE,
    internal = FALSE,
    path.folder = path.folder
  )


  # Filters----------------------------------------------------------------------
  if (filters) {
    data <- radiator::filter_individuals(
      data = data,
      interactive.filter = interactive.filter,
      filter.individuals.missing = "outliers",
      filter.individuals.heterozygosity = "outliers",
      internal = FALSE,
      path.folder = path.folder,
      parallel.core = parallel.core,
      verbose = FALSE
    ) %>%
      radiator::filter_ld(
        data = .,
        interactive.filter = interactive.filter,
        filter.short.ld = "mac",
        parallel.core = parallel.core,
        verbose = FALSE,
        internal = FALSE,
        path.folder = path.folder
      )

  }

  # New strata after filter in case id were removed
  strata <- radiator::extract_individuals_metadata(
    gds = data,
    ind.field.select = c("TARGET_ID", "INDIVIDUALS", "STRATA"),
    whitelist = TRUE
  )

  # check Sex Ratio ------------------------------------------------------------
  sex.ratio <- dplyr::filter(strata, STRATA != "U") %>%
    dplyr::count(STRATA, name = "SEX_RATIO")
  message("\n\nSex-ratio (F/M): ", round((sex.ratio$SEX_RATIO[sex.ratio$STRATA == "F"] / sex.ratio$SEX_RATIO[sex.ratio$STRATA == "M"]), 2))


  gds.bk <- data
  # gds.bk -> data


  # Extract from the GDS the tidy data... ###ISSUE FOR VCF & 1row
  # data <-

  # temp <- radiator::extract_genotypes_metadata(gds = data, index.only = TRUE)
  # if (is.null(temp)) {
  #   temp <- radiator::extract_genotypes_metadata(gds = data)
  # }

  data <- radiator::extract_genotypes_metadata(
    gds = data,
    genotypes.meta.select = c("MARKERS", "INDIVIDUALS", "GT_BIN", "READ_DEPTH"),
    whitelist = TRUE
  ) %>%
    radiator::join_strata(data = .,
                          strata = strata,
                          verbose = FALSE)


  # SILICO files ----------------------------------------------------------------

  if (!is.null(silicodata)) {
    data.type <- radiator::detect_genomic_format(silicodata)
    silicodata <- radiator::read_dart(
      data = silicodata,
      strata = strata,
      tidy.dart = FALSE,
      parallel.core = parallel.core,
      path.folder = path.folder,
      internal = TRUE,
      verbose = TRUE
    ) %>%
      dplyr::rename(., MARKERS = CLONE_ID)

    if (max(silicodata$VALUE, na.rm = TRUE) > 1) {
      data.source <- c("counts", data.type, data.source)
    } else{
      data.source <- c("genotype", data.type, data.source)
    }
  }

  ### add warning about not using count data
  if (!(any(c("counts", "vcf.file") %in% data.source))) {
    message(
      "Your data does not have Read Depth information; the analysis based on Read Depth will not be performed"
    )
  }


  # Sex markers-----------------------------------------------------------------
  # 1. read depth per strata and markers markers: xz
  # 2. presence/absence per strata and markers: yw
  # 3. heterozygosity per strata and markers: xz

  # Summarise DArT counts/silico and VCFs---------------------------------------
  if (is.null(tau)) tau <- 0.03

  res$sum <- summarize_sex(
    data = data,
    silicodata = silicodata,
    coverage.thresholds = coverage.thresholds,
    data.source = data.source,
    tau = tau
  )
  data.sum <- res$sum$data.sum
  silico.sum <- res$sum$silico.sum

  # Figures -------------------------------------------------------------------

  # Here we start to generate the figures for each sex determination mechanism
  # The function for figures was moved outside the function below
  # you need to load the function first so that it's accessible
  # e.g.

  SexID <- "visually"

  #### P/A #####

  # perhaps for count data, set coverage threshold, based on plot

  ###* For data.sum ####
  ### SCATTER
  plot.filename <- "sex_markers_PA_scatter_plot"
  if (!is.null(species) && !is.null(population)) {
    plot.filename <-
      stringi::stri_join(plot.filename, species, population, sep = "_")
  }
  plot.filename <-
    stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

  scat.fig <- sex_markers_plot(
    data = data.sum,
    x = "PRESENCE_ABSENCE_F",
    y = "PRESENCE_ABSENCE_M",
    scat = TRUE,
    qreg = FALSE,
    tuckey = FALSE,
    x.title = "Proportion of females",
    y.title = "Proportion of males",
    subtitle = if (is.null(population)) {
      paste0("Sex is ", SexID, " assigned")
    } else {
      paste0(population, ": Sex is ", SexID, " assigned")
    },
    title = "Absence of each SNP marker between females and males",
    plot.filename = plot.filename
  )
  print(scat.fig)

  ### TUCKEY
  plot.filename <- "sex_markers_PA_tuckey_plot"
  if (!is.null(species) && !is.null(population)) {
    plot.filename <-
      stringi::stri_join(plot.filename, species, population, sep = "_")
  }
  plot.filename <-
    stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

  scat.fig <- sex_markers_plot(
    data = data.sum,
    x = "MEAN",
    y = "DIFF",
    scat = TRUE,
    qreg = FALSE,
    tuckey = TRUE,
    x.title = "Mean of males and females",
    y.title = "Difference between males and females (<- Y/W ->)",
    subtitle = if (is.null(population)) {
      paste0("Sex is ", SexID, " assigned")
    } else {
      paste0(population, ": Sex is ", SexID, " assigned")
    },
    title = "Tukey mean-difference plot of each SNP marker between females \nand males",
    plot.filename = plot.filename
  )
  print(scat.fig)

  # Interacive selection of threshold
  if (interactive.filter) {
    filter.y.markers <- radiator::radiator_question(x = "Have a look at the figure: do you want to filter the data (y/n): ", answer.opt = c("y", "n"))

    if (filter.y.markers == "y") {
      threshold.y.markers <-
        radiator::radiator_question(x = "Choose the threshold for y sex markers, note that the Y-axis is inverted (-1 to 1): ",
                                    minmax = c(-1, 1))
    } else {
      threshold.y.markers <- NULL
    }
  }

  # Select and print sex markers
  if (!is.null(threshold.y.markers)) {
    if (threshold.y.markers < 0) {
      y.markers <- dplyr::filter(data.sum, DIFF < threshold.y.markers)$MARKERS
    } else {
      y.markers <- dplyr::filter(data.sum, DIFF > threshold.y.markers)$MARKERS
    }
    res$heterogametic.markers <- y.markers
  } else{
    y.markers <- NULL
  }


  ###* For silico.sum ####
  if ("silico.dart" %in% data.source) {
    ### SCATTER
    plot.filename <- "sex_markers_SILICO_PA_scatter_plot"
    if (!is.null(species) && !is.null(population)) {
      plot.filename <-
        stringi::stri_join(plot.filename, species, population, sep = "_")
    }
    plot.filename <-
      stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

    scat.fig <- sex_markers_plot(
      data = silico.sum,
      x = "PRESENCE_ABSENCE_F",
      y = "PRESENCE_ABSENCE_M",
      scat = TRUE,
      qreg = FALSE,
      tuckey = FALSE,
      x.title = "Proportion of females",
      y.title = "Proportion of males",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Absence of each SILICO marker between females and males",
      plot.filename = plot.filename
    )
    print(scat.fig)

    ### TUCKEY
    plot.filename <- "sex_markers_SILICO_PA_tuckey_plot"
    if (!is.null(species) && !is.null(population)) {
      plot.filename <-
        stringi::stri_join(plot.filename, species, population, sep = "_")
    }
    plot.filename <-
      stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

    scat.fig <- sex_markers_plot(
      data = silico.sum,
      x = "MEAN",
      y = "DIFF",
      scat = TRUE,
      qreg = FALSE,
      tuckey = TRUE,
      x.title = "Mean of males and females",
      y.title = "Difference between males and females (<- Y/W ->)",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Tukey mean-difference plot of each SILICO marker between females \nand males",
      plot.filename = plot.filename
    )
    print(scat.fig)


    # Interacive selection of threshold
    if (interactive.filter) {
      filter.y.markers <- radiator::radiator_question(x = "Have a look at the figure: do you want to filter the data (y/n): ", answer.opt = c("y", "n"))

      if (filter.y.markers == "y") {
        threshold.y.silico.markers <-
          radiator::radiator_question(x = "Choose the threshold for y SILICO sex markers, note that the Y-axis is inverted (-1 to 1): ",
                                      minmax = c(-1, 1))
      } else {
        threshold.y.silico.markers <- NULL
        y.silico.markers <- NULL
      }
    }

    # Select and print sex markers
    if (!is.null(threshold.y.silico.markers)) {
      if (threshold.y.silico.markers < 0) {
        y.silico.markers <- dplyr::filter(silico.sum, DIFF < threshold.y.silico.markers)$MARKERS
      } else {
        y.silico.markers <- dplyr::filter(silico.sum, DIFF > threshold.y.silico.markers)$MARKERS
      }
      res$heterogametic.silico.markers <- y.silico.markers
    }
  } else{
    y.silico.markers <- NULL
  }

  ####* SET GENETIC SEX STRATA ####
  if (exists("y.markers") && !rlang::is_empty(y.markers)) {
    y.data <- dplyr::filter(data, MARKERS %in% y.markers) %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::summarise(MEAN_RD = mean(READ_DEPTH)) %>%
      dplyr::ungroup(.) %>%
      dplyr::left_join(strata, by = "INDIVIDUALS") %>%
      dplyr::mutate(VISUAL_SEX = STRATA) %>%
      dplyr::mutate(
        GENETIC_SEX =
          dplyr::case_when(
            is.na(MEAN_RD) | MEAN_RD < coverage.thresholds ~ "F",
            !(is.na(MEAN_RD) |
                MEAN_RD < coverage.thresholds) ~ "M",
            VISUAL_SEX == "F" &
              (is.na(MEAN_RD) |
                 MEAN_RD < coverage.thresholds) ~ "U"
          )
      ) %>%
      dplyr::mutate(STRATA = GENETIC_SEX)



    # print visual and gentic sex table
    SumTable <-
      tibble::tibble(
        Visual_Sex = c("M", "M", "F", "F", "F", "U", "U"),
        Genetic_Sex = c("M", "F", "F", "M", "U", "M", "F"),
        Individuals = c(
          length(which(
            y.data$VISUAL_SEX == "M" & y.data$GENETIC_SEX == "M"
          )),
          length(which(
            y.data$VISUAL_SEX == "M" & y.data$GENETIC_SEX == "F"
          )),
          length(which(
            y.data$VISUAL_SEX == "F" & y.data$GENETIC_SEX == "F"
          )),
          length(which(
            y.data$VISUAL_SEX == "F" & y.data$GENETIC_SEX == "M"
          )),
          length(which(
            y.data$VISUAL_SEX == "F" & y.data$GENETIC_SEX == "U"
          )),
          length(which(
            y.data$VISUAL_SEX == "U" & y.data$GENETIC_SEX == "M"
          )),
          length(which(
            y.data$VISUAL_SEX == "U" & y.data$GENETIC_SEX == "F"
          ))
        )
      )
    print(SumTable)
    readr::write_tsv(x = SumTable, path = file.path(path.folder, "summary_table_visual_genetic_sex.tsv"))

    if (exists("y.silico.markers") && !rlang::is_empty(y.silico.markers)) {
      y.silico.data <-
        dplyr::filter(silicodata, MARKERS %in% y.silico.markers) %>%
        dplyr::group_by(INDIVIDUALS) %>%
        dplyr::summarise(MEAN_RD = mean(VALUE)) %>%
        dplyr::ungroup(.) %>%
        dplyr::left_join(strata, by = "INDIVIDUALS") %>%
        dplyr::mutate(VISUAL_SEX = STRATA) %>%
        dplyr::mutate(GENETIC_SEX =
                        dplyr::if_else(condition =
                                         is.na(MEAN_RD) |
                                         MEAN_RD < coverage.thresholds, "F", "M"
                        )) %>%
        dplyr::mutate(STRATA = GENETIC_SEX)

      # print visual and gentic sex table
      SumTable <-
        tibble::tibble(
          Visual_Sex = c("M", "M", "F", "F", "U", "U"),
          Genetic_Sex = c("M", "F", "F", "M", "M", "F"),
          Individuals = c(
            length(which(y.silico.data$VISUAL_SEX == "M" & y.silico.data$GENETIC_SEX == "M"
            )),
            length(which(y.silico.data$VISUAL_SEX == "M" & y.silico.data$GENETIC_SEX == "F"
            )),
            length(which(y.silico.data$VISUAL_SEX == "F" & y.silico.data$GENETIC_SEX == "F"
            )),
            length(which(y.silico.data$VISUAL_SEX == "F" & y.silico.data$GENETIC_SEX == "M"
            )),
            # length(which(y.silico.data$VISUAL_SEX == "F" & y.silico.data$GENETIC_SEX == "U"
            # )),
            length(which(y.silico.data$VISUAL_SEX == "U" & y.silico.data$GENETIC_SEX == "M"
            )),
            length(which(y.silico.data$VISUAL_SEX == "U" & y.silico.data$GENETIC_SEX == "F"
            )))
        )
      print(SumTable)
      readr::write_tsv(x = SumTable, path = file.path(path.folder, "summary_table_visual_genetic_sex_silico.tsv"))


      # print genetic SNP and gentic SILICO sex table
      SumTable <-
        tibble::tibble(
          GENETIC_SEX_SNP = c("M", "M", "F", "F", "U", "U"),
          Genetic_Sex_SILICO = c("M", "F", "F", "M", "M", "F"),
          Individuals = c(
            length(which(y.data$GENETIC_SEX == "M" & y.silico.data$GENETIC_SEX == "M")
            ),
            length(
              which(y.data$GENETIC_SEX == "M" & y.silico.data$GENETIC_SEX == "F")
            ),
            length(
              which(y.data$GENETIC_SEX == "F" & y.silico.data$GENETIC_SEX == "F")
            ),
            length(
              which(y.data$GENETIC_SEX == "F" & y.silico.data$GENETIC_SEX == "M")
            ),
            length(
              which(y.data$GENETIC_SEX == "U" & y.silico.data$GENETIC_SEX == "M")
            ),
            length(
              which(y.data$GENETIC_SEX == "U" & y.silico.data$GENETIC_SEX == "F")
            ))
        )
      print(SumTable)
      readr::write_tsv(x = SumTable, path = file.path(path.folder, "summary_table_genetic_sex_silico.tsv"))


      # Interacive selection which sex info
      if (is.null(sex.id.input)) sex.id.input <- "1"
      if (interactive.filter) {
        sex.id.input <- radiator::radiator_question(
          x = "Do you want to continue based on (1) visual, (2) genetic dart or (3) genetic silico sex?\nWe advise (3) for better results",
          answer.opt = c("1","2", "3") # Not working with numbers...
        )
        sex.id.input <- as.integer(sex.id.input)
      }
    } else if (rlang::is_empty(y.silico.markers)) {
      # Interacive selection which sex info
      if (is.null(sex.id.input))
        sex.id.input <- "1"
      if (interactive.filter) {
        sex.id.input <- radiator::radiator_question(x = "Do you want to continue based on (1) visual or (2) genetic dart?\nWe advise (2) for better results",
                                                    answer.opt = c("1", "2"))
        sex.id.input <- as.integer(sex.id.input)
      }
    }

    #set new sexID for Het analysis
    if (sex.id.input == 2) {
      SexID <- "genetically (SNP)"
      ### recalculate data based on new sexID??
      data.genetic <- dplyr::rename(data, VISUAL_STRATA = STRATA) %>%
        dplyr::left_join(y.data, by = "INDIVIDUALS") %>%
        dplyr::mutate(GENETIC_STRATA = STRATA)
        radiator::write_rad(data = data.genetic, path = file.path(wd, "sexy_markers_temp"))

        if ("silico.dart" %in% data.source) {
          data.silico.genetic <-
            dplyr::rename(silicodata, VISUAL_STRATA = STRATA) %>%
            dplyr::left_join(y.data, by = "INDIVIDUALS") %>%
            dplyr::mutate(GENETIC_STRATA = STRATA)
          radiator::write_rad(data = data.silico.genetic,
                              path = file.path(wd, "sexy_markers_silico_temp"))
        }else{data.silico.genetic <- NULL}

    } else if (sex.id.input == 3) {
      SexID <- "genetically (SILICO)"
      ### recalculate data based on new sexID??
      data.genetic <- dplyr::rename(data, VISUAL_STRATA = STRATA) %>%
        dplyr::left_join(y.silico.data, by = "INDIVIDUALS") %>%
        dplyr::mutate(GENETIC_STRATA = STRATA)
        radiator::write_rad(data = data.genetic, path = file.path(wd, "sexy_markers_temp"))

        if ("silico.dart" %in% data.source) {
          data.silico.genetic <-
            dplyr::rename(silicodata, VISUAL_STRATA = STRATA) %>%
            dplyr::left_join(y.silico.data, by = "INDIVIDUALS") %>%
            dplyr::mutate(GENETIC_STRATA = STRATA)
          radiator::write_rad(data = data.silico.genetic,
                              path = file.path(wd, "sexy_markers_silico_temp"))
        } else{data.silico.genetic <- NULL}
    }


    # silico
    if (sex.id.input != 1){
      res$sum <- summarize_sex(
        data = data.genetic,
        silicodata = data.silico.genetic,
        coverage.thresholds = coverage.thresholds,
        data.source = data.source,
        tau = tau
      )
      data.sum <- res$sum$data.sum
      silico.sum <- res$sum$silico.sum
    }
  }

  # Remove markers that are already extracted and have high missingness
  ##SHOULD MAKE THE MISSINGNESS INTERACTIVE
  data.sum <- dplyr::filter(data.sum, !(MARKERS %in% y.markers | MISSINGNESS > 0.2))
  if(!is.null(silicodata)){
    silico.sum <- dplyr::filter(silico.sum, !(MARKERS %in% y.silico.markers| MISSINGNESS > 0.2))
  }


  #### HET ####
  # You want to remove the Y-linked markers and do a filter on missingness

  ### SCATTER
  plot.filename <- "sex_markers_HET_scat_plot"
  if (!is.null(species) && !is.null(population)) {
    plot.filename <-
      stringi::stri_join(plot.filename, species, population, sep = "_")
  }
  plot.filename <-
    stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

  scat.fig <- sex_markers_plot(
    data = data.sum,
    x = "MEAN_HET_F",
    y = "MEAN_HET_M",
    scat = TRUE,
    qreg = FALSE,
    x.title = "Proportion of females",
    y.title = "Proportion of males",
    subtitle = if (is.null(population)) {
      paste0("Sex is ", SexID, " assigned")
    } else {
      paste0(population, ": Sex is ", SexID, " assigned")
    },
    title = "Heterozygosity of each SNP marker between females and males",
    plot.filename = plot.filename
  )
  print(scat.fig)

  ###QREG -> Tau value is important
  plot.filename <- "sex_markers_HET_qr_plot"
  if (!is.null(species) && !is.null(population)) {
    plot.filename <-
      stringi::stri_join(plot.filename, species, population, sep = "_")
  }
  plot.filename <-
    stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

  qr.fig <- sex_markers_plot(
    data = data.sum,
    x = "QR_RESIDUALS",
    y = "MEAN_HET_M",
    qreg = TRUE,
    tuckey = FALSE,
    x.title = "Quantile residuals for M~F (<- X/Z ->)",
    y.title = "Proportion of heterozygosity in males",
    subtitle = if (is.null(population)) {
      paste0("Sex is ", SexID, " assigned; tau = ", tau)
    } else {
      paste0(population, ": Sex is ", SexID, " assigned; tau = ", tau)
    },
    title = "Quantile residual plot of each SNP marker between females and males",
    plot.filename = plot.filename
  )
  print(qr.fig)

  # Interacive selection of threshold
  if (interactive.filter) {
    filter.x.markers <- radiator::radiator_question(x = "Have a look at the figure: do you want to filter the data (y/n): ", answer.opt = c("y", "n"))

    if (filter.x.markers == "y") {
      threshold.x.markers.qr <-
        radiator::radiator_question(x = "Choose the threshold for x sex markers(-1 to 1): ",
                                    minmax = c(-1, 1))
    } else {
      threshold.x.markers.qr <- NULL
    }
  }

  if (!is.null(threshold.x.markers.qr)) {
    if (threshold.x.markers.qr < 0) {
      x.markers <- dplyr::filter(data.sum, QR_RESIDUALS < threshold.x.markers.qr)$MARKERS
    } else{
      x.markers <- dplyr::filter(data.sum, QR_RESIDUALS > threshold.x.markers.qr)$MARKERS
    }
    res$homogametic.het.markers <- x.markers
  }



  #### READ DEPTH ####
  ####* For data.sum ####
  if (all(c("dart", "counts") %in% data.source)) {
    ### SCATTER
    plot.filename <- "sex_markers_RD_scat_plot"
    if (!is.null(species) && !is.null(population)) {
      plot.filename <-
        stringi::stri_join(plot.filename, species, population, sep = "_")
    }

    plot.filename <-
      stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

    scat.fig <- sex_markers_plot(
      data = data.sum,
      x = "MEAN_READ_DEPTH_F",
      y = "MEAN_READ_DEPTH_M",
      scat = TRUE,
      qreg = FALSE,
      tuckey = FALSE,
      RD = TRUE,
      x.title = "Mean counts for females",
      y.title = "Mean counts for males",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Average counts of each marker between females and males",
      plot.filename = plot.filename
    )
    print(scat.fig)


    ##HIST
    plot.filename <- "sex_markers_RD_hist_plot"
    if (!is.null(species) && !is.null(population)) {
      plot.filename <-
        stringi::stri_join(plot.filename, species, population, sep = "_")
    }

    plot.filename <-
      stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

    scat.fig <- sex_markers_plot(
      data = data.sum,
      x = "RATIO",
      qreg = FALSE,
      tuckey = FALSE,
      hist = TRUE,
      x.title = "Female - male count ratio",
      y.title = "",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Histogram of females counts over males counts for each marker",
      plot.filename = plot.filename
    )
    print(scat.fig)

    ##NEEDS AN INTERACTIVE ZOOM

    ##HIST2
    plot.filename <- "sex_markers_RD_hist_subsetted_plot"
    if (!is.null(species) && !is.null(population)) {
      plot.filename <-
        stringi::stri_join(plot.filename, species, population, sep = "_")
    }

    plot.filename <-
      stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

    scat.fig <- sex_markers_plot(
      data = dplyr::filter(data.sum, RATIO > 1.2),
      x = "RATIO",
      qreg = FALSE,
      tuckey = FALSE,
      hist = TRUE,
      x.title = "Female - male count ratio",
      y.title = "",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Histogram of females counts over males counts for each marker",
      plot.filename = plot.filename
    )
    print(scat.fig)


    # Interacive selection of threshold
    if (interactive.filter) {
      filter.x.markers <- radiator::radiator_question(x = "Have a look at the figure: do you want to filter the data (y/n): ", answer.opt = c("y", "n"))

      if (filter.x.markers == "y") {
        threshold.x.markers.RD <-
          radiator::radiator_question(x = "Choose the RATIO threshold for x sex markers: ", minmax = c(-Inf,Inf) )
      } else {
        threshold.x.markers.RD <- NULL
      }
    }

    if (!is.null(threshold.x.markers.RD)) {
      if (threshold.x.markers.RD > 1) {
        x.markers <- dplyr::filter(data.sum, RATIO > threshold.x.markers.RD)$MARKERS
      } else{
        x.markers <- dplyr::filter(data.sum, RATIO < threshold.x.markers.RD)$MARKERS
      }
      res$homogametic.RD.markers <- x.markers
    }
  }





  ####* For silico.sum ####
  if (all(c("silico.dart", "counts") %in% data.source)) {
    ### SCATTER
    plot.filename <- "sex_markers_SILICO_RD_scat_plot"
    if (!is.null(species) && !is.null(population)) {
      plot.filename <-
        stringi::stri_join(plot.filename, species, population, sep = "_")
    }

    plot.filename <-
      stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

    scat.fig <- sex_markers_plot(
      data = silico.sum,
      x = "MEAN_READ_DEPTH_F",
      y = "MEAN_READ_DEPTH_M",
      scat = TRUE,
      qreg = FALSE,
      tuckey = FALSE,
      RD = TRUE,
      x.title = "Mean counts for females",
      y.title = "Mean counts for males",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Average counts of each silico marker between females and males",
      plot.filename = plot.filename
    )
    print(scat.fig)


    ##HIST
    plot.filename <- "sex_markers_SILICO_RD_hist_plot"
    if (!is.null(species) && !is.null(population)) {
      plot.filename <-
        stringi::stri_join(plot.filename, species, population, sep = "_")
    }

    plot.filename <-
      stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

    scat.fig <- sex_markers_plot(
      data = silico.sum,
      x = "RATIO",
      qreg = FALSE,
      tuckey = FALSE,
      hist = TRUE,
      x.title = "Female - male count ratio",
      y.title = "",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Histogram of females counts over males counts for each SILICO",
      plot.filename = plot.filename
    )
    print(scat.fig)


    ##HIST2
    plot.filename <- "sex_markers_SILICO_RD_hist_subsetted_plot"
    if (!is.null(species) && !is.null(population)) {
      plot.filename <-
        stringi::stri_join(plot.filename, species, population, sep = "_")
    }

    plot.filename <-
      stringi::stri_join(path.folder, "/", plot.filename, "_", file.date, ".pdf")

    scat.fig <- sex_markers_plot(
      data = dplyr::filter(silico.sum, RATIO > 1.2),
      x = "RATIO",
      qreg = FALSE,
      tuckey = FALSE,
      hist = TRUE,
      x.title = "Female - male count ratio",
      y.title = "",
      subtitle = if (is.null(population)) {
        paste0("Sex is ", SexID, " assigned")
      } else {
        paste0(population, ": Sex is ", SexID, " assigned")
      },
      title = "Histogram of females counts over males counts for each SILICO",
      plot.filename = plot.filename
    )
    print(scat.fig)


    # Interacive selection of threshold
    if (interactive.filter) {
      filter.x.markers <-
        radiator::radiator_question(x = "Have a look at the figure: do you want to filter the data (y/n): ", answer.opt = c("y", "n"))

      if (filter.x.markers == "y") {
        threshold.x.markers.RD.silico <-
          radiator::radiator_question(x = "Choose the RATIO threshold for x SILICO sex markers: ", minmax = c(-Inf, Inf))
      } else {
        threshold.x.markers.RD.silico <- NULL
      }
    }

    if (!is.null(threshold.x.markers.RD.silico)) {
      if (threshold.x.markers.RD.silico > 1) {
        x.markers <-
          dplyr::filter(silico.sum, RATIO > threshold.x.markers.RD.silico)$MARKERS
      } else{
        x.markers <-
          dplyr::filter(silico.sum, RATIO < threshold.x.markers.RD.silico)$MARKERS
      }
      res$homogametic.RD.silico.markers <- x.markers
    }
  }




  ## TODO:
  # counts & vcf
  # silico
  # res$homogametic.rd.silico.markers <- x.markers


  # Export -------------------------------------------------------------------

  # summary of all sex-markers per method (PA, HET, RD)

  # common markers plot
  #

  # FASTA file with sex markers for all methods
  if(class(data) == "SeqVarGDSClass"){
    # Set sex-marker to whitelist and allign the sex-marker method with the markers
    meta <- radiator::extract_markers_metadata(
      gds = data,
      markers.meta.select = c("MARKERS", "SEQUENCE"),
      radiator.node = TRUE,
      whitelist = FALSE,
      blacklist = FALSE,
      verbose = FALSE
    )

  }



  return(res)
}#End sexy_markers


# INTERNAL NESTED FUNCTION------------------------------------------------------
#' @title sex_markers_plot
#' @description Function to generate the different figures required.
#' @keywords internal
#' @export

sex_markers_plot <- function(
  data, x, y, x.title = "x title", y.title = "y title", subtitle = "subtitle", title = "Big title",
  width = 15, height = 15, scat = FALSE,tuckey = FALSE, qreg = FALSE, hist = FALSE, RD = FALSE, plot.filename = NULL) {
  x <- dplyr::sym(x)

  element.text <- ggplot2::element_text(size = 10, face = "bold")

  if (qreg) data %<>% dplyr::filter(!is.na(QR_RESIDUALS))
  if (hist) {
    sex.plot <- ggplot2::ggplot(
      data = data,
      ggplot2::aes(x = !!x)
    ) +
      ggplot2::geom_histogram(
        col = "black",
        fill = "grey",
        alpha = .2,
        binwidth = 0.1)
  } else{
    y <- dplyr::sym(y)
    sex.plot <- ggplot2::ggplot(
      data = data,
      ggplot2::aes(x = !!x, y = !!y)) +
      ggplot2::geom_point(size = 2, alpha = 0.3)

    if (tuckey) {
      sex.plot <- sex.plot + ggplot2::scale_y_reverse()
    } else if (qreg) {
      sex.plot <- sex.plot
    } else if (RD){
      sex.plot <- sex.plot + ggplot2::geom_abline(slope = c(2,1,1/2))
    } else if (scat) {
      sex.plot <- sex.plot + ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::scale_y_continuous(limits = c(0, 1))
    } else {
      sex.plot <- sex.plot + ggplot2::geom_abline()
    }
  }

  sex.plot <- sex.plot +
    ggplot2::labs(x = x.title, y = y.title, subtitle = subtitle, title = title) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
      axis.title.y = element.text,
      axis.title.x = element.text,
      axis.text.x = element.text)
  # print(sex.plot)

  ggplot2::ggsave(
    filename = plot.filename,
    plot = sex.plot,
    width = width,
    height = height,
    dpi = 300,
    units = "cm",
    useDingbats = FALSE
  )
  return(sex.plot)
}

#' @title summarize_sex
#' @description Function to generate the different figures required.
#' @keywords internal
#' @export

summarize_sex <- function (data, silicodata, data.source, coverage.thresholds = NULL, tau = 0.3) {
  if (tibble::has_name(data, "READ_DEPTH")) {
    # With DArT count data and VCFs

    mis <-
      dplyr::select(data, MARKERS, INDIVIDUALS, STRATA, READ_DEPTH, GT_BIN) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(MISSINGNESS = length(INDIVIDUALS[is.na(GT_BIN)]) / length(INDIVIDUALS)) %>%
      dplyr::ungroup(.)

    data.sum <-
      dplyr::select(data, MARKERS, INDIVIDUALS, STRATA, READ_DEPTH, GT_BIN) %>%
      dplyr::group_by(STRATA, MARKERS) %>%
      dplyr::summarise(
        PRESENCE_ABSENCE = length(INDIVIDUALS[READ_DEPTH < coverage.thresholds]) / length(INDIVIDUALS),
        MEAN_HET = length(INDIVIDUALS[GT_BIN == 1]) / length(INDIVIDUALS),
        MEAN_READ_DEPTH = mean(READ_DEPTH, na.rm = TRUE)
      ) %>%
      dplyr::ungroup(.) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = MARKERS ~ STRATA,

        value.var = c("PRESENCE_ABSENCE", "MEAN_READ_DEPTH", "MEAN_HET")
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(
        MEAN = (PRESENCE_ABSENCE_M + PRESENCE_ABSENCE_F) / 2,
        DIFF = PRESENCE_ABSENCE_M - PRESENCE_ABSENCE_F,
        RATIO = MEAN_READ_DEPTH_F / MEAN_READ_DEPTH_M
      ) %>%
      dplyr::filter(RATIO != 0) %>%
      dplyr::ungroup(.) %>%
      dplyr::left_join(mis, by = "MARKERS")


    data.sum[is.na(data.sum)] <- NA
    mis <- NULL

    if (is.null(tau)) tau <- 0.03
    #Quantile regression plot
    data.sum <- dplyr::bind_rows(
      dplyr::filter(data.sum, is.na(MEAN_HET_M) | is.na(MEAN_HET_F)),
      dplyr::filter(data.sum, !is.na(MEAN_HET_M) &
                      !is.na(MEAN_HET_F)) %>%
        dplyr::mutate(
          QR_RESIDUALS = quantreg::rq(MEAN_HET_M ~ MEAN_HET_F,
                                      tau = tau,
                                      data = .) %$%
            residuals
        )
    )
  } else if (tibble::has_name(data, "GT_BIN") ){  #genotype data
    mis <-
      dplyr::select(data, MARKERS, INDIVIDUALS, STRATA, GT_BIN) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(MISSINGNESS = length(INDIVIDUALS[is.na(GT_BIN)]) / length(INDIVIDUALS)) %>%
      dplyr::ungroup(.)

    data.sum <-
      dplyr::select(data, MARKERS, INDIVIDUALS, STRATA, GT_BIN) %>%
      dplyr::group_by(STRATA, MARKERS) %>%
      dplyr::summarise(
        PRESENCE_ABSENCE = length(INDIVIDUALS[is.na(GT_BIN)]) / length(INDIVIDUALS),##ISSUE
        MEAN_HET = length(INDIVIDUALS[GT_BIN == 1]) / length(INDIVIDUALS)
      ) %>%
      dplyr::ungroup(.) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = MARKERS ~ STRATA,
        value.var = c("PRESENCE_ABSENCE", "MEAN_HET")
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(
        MEAN = (PRESENCE_ABSENCE_M + PRESENCE_ABSENCE_F) / 2,
        DIFF = PRESENCE_ABSENCE_M - PRESENCE_ABSENCE_F
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::left_join(mis, by = "MARKERS")


    data.sum[is.na(data.sum)] <- NA
    mis <- NULL

    #Quantile regression plot
    data.sum <- dplyr::bind_rows(
      dplyr::filter(data.sum, is.na(MEAN_HET_M) | is.na(MEAN_HET_F)),
      dplyr::filter(data.sum, !is.na(MEAN_HET_M) &
                      !is.na(MEAN_HET_F)) %>%
        dplyr::mutate(
          QR_RESIDUALS = quantreg::rq(MEAN_HET_M ~ MEAN_HET_F,
                                      tau = tau,
                                      data = .) %$%
            residuals
        )
    )

  }
    else {
    data.sum <- NULL
  }


  if ("silico.dart" %in% data.source) {
    mis <-
      dplyr::select(silicodata, MARKERS, INDIVIDUALS, STRATA, VALUE) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(MISSINGNESS = length(INDIVIDUALS[VALUE < coverage.thresholds])
                       / length(INDIVIDUALS)) %>%
      dplyr::ungroup(.)

    ### silico counts
    if ("counts" %in% data.source) {
      silico.sum <- silicodata %>%
        # dplyr::rename(., READ_DEPTH = VALUE) %>%
        dplyr::group_by(STRATA, MARKERS) %>%
        dplyr::summarise(
          PRESENCE_ABSENCE = length(INDIVIDUALS[VALUE < coverage.thresholds])
          / length(INDIVIDUALS[!is.na(VALUE)]),
          MEAN_READ_DEPTH = mean(VALUE, na.rm = TRUE)
        ) %>%
        dplyr::ungroup(.) %>%
        data.table::as.data.table(.) %>%
        data.table::dcast.data.table(
          data = .,
          formula = MARKERS ~ STRATA ,
          value.var = c("PRESENCE_ABSENCE", "MEAN_READ_DEPTH"),
          fill = NA
        ) %>%
        tibble::as_tibble(.) %>%
        dplyr::mutate(
          RATIO = MEAN_READ_DEPTH_F / MEAN_READ_DEPTH_M,
          MEAN = (PRESENCE_ABSENCE_M + PRESENCE_ABSENCE_F) / 2,
          DIFF = PRESENCE_ABSENCE_M - PRESENCE_ABSENCE_F
        )%>%
        dplyr::left_join(mis, by = "MARKERS")
    } else {
      silico.sum <- silicodata %>%
        dplyr::group_by(STRATA, MARKERS) %>%
        dplyr::summarise(PRESENCE_ABSENCE = length(INDIVIDUALS[VALUE < coverage.thresholds])
                         / length(INDIVIDUALS[!is.na(VALUE)])) %>%
        dplyr::ungroup(.) %>%
        data.table::as.data.table(.) %>%
        data.table::dcast.data.table(
          data = .,
          formula = MARKERS ~ STRATA,
          value.var = "PRESENCE_ABSENCE",
          fill = NA
        ) %>%
        tibble::as_tibble(.) %>%
        dplyr::rename(
          PRESENCE_ABSENCE_M = M,
          PRESENCE_ABSENCE_F = F
        ) %>%
        dplyr::mutate(
          MEAN = (PRESENCE_ABSENCE_M + PRESENCE_ABSENCE_F) / 2,
          DIFF = PRESENCE_ABSENCE_M - PRESENCE_ABSENCE_F
        )%>%
        dplyr::left_join(mis, by = "MARKERS")
    }
  } else {
    silico.sum <- NULL
  }
  return(list(data.sum = data.sum, silico.sum = silico.sum))
}


#' @title Export FASTA
#' @description Function to write fasta files for sex markers
#' @keywords internal
#' @export
# FASTA from gds


# FASTA file (different for silico)
write_fasta <- function (sexmarkdf, filename) {
  afile <- file(filename, open = 'w')
  if (!is.null(species) && !is.null(population)) {
    for (i in 1:length(sexmarkdf$MARKERS)) {
      cat(
        paste(
          '>',Species,"|",population,'|',sexmarkdf$method[i],
          sexmarkdf$MARKERS[i],'\n',
          sexmarkdf$SEQUENCE[i],'\n',
          sep = ''
        ),
        sep = '',
        file = afile
      )
    }
  } else {
    for (i in 1:length(sexmarkdf$MARKERS)) {
      cat(
        paste(
          '>', sexmarkdf$method[i],"|",sexmarkdf$MARKERS[i],'\n',
          sexmarkdf$SEQUENCE[i],'\n',
          sep = ''
        ),
        sep = '',
        file = afile
      )
    }
  }
  close(afile)
}#write_fasta
