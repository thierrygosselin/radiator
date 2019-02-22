# Make genomic input file tidy

#' @name tidy_genomic_data
#' @title Transform common genomic dataset format in a tidy data frame
#' @description Transform genomic data set produced by massive parallel
#' sequencing pipeline (e.g.GBS/RADseq,
#' SNP chip, DArT, etc) into a tidy format.
#' The use of blacklist and whitelist along
#' several filtering options are available to prune the dataset.
#' Several arguments are available to make your data population-wise and easily
#' rename the pop id.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data 12 options for input: VCFs (SNPs or Haplotypes,
#' to make the vcf population ready),
#' plink, stacks haplotype file, genind (library(adegenet)),
#' genlight (library(adegenet)), gtypes (library(strataG)), genepop, DArT,
#' and a data frame in long/tidy or wide format. To verify that radiator detect
#' your file format use \code{\link{detect_genomic_format}} (see example below).
#' Documented in \strong{Input genomic datasets} of \code{\link{tidy_genomic_data}}.
#'

#' @param strata (optional)
#' The strata file is a tab delimited file with a minimum of 2 columns headers:
#' \code{INDIVIDUALS} and \code{STRATA}. Documented in \code{\link{read_strata}}.
#' Default: \code{strata = NULL}.


#' @param filename (optional) The function uses \code{\link[fst]{write.fst}},
#' to write the tidy data frame in
#' the working directory. The file extension appended to
#' the \code{filename} provided is \code{.rad}.
#' With default: \code{filename = NULL}, the tidy data frame is
#' in the global environment only (i.e. not written in the working directory...).


#' @inheritParams radiator_common_arguments
#' @param ... (optional) To pass further arguments for fine-tuning the function.

#' @return The output in your global environment is a tidy data frame.
#' If \code{filename} is provided, the tidy data frame is also
#' written in the working directory with file extension \code{.rad}.
#' The file is written with the
#' \href{https://github.com/fstpackage/fst}{Lightning Fast Serialization of Data Frames for R} package.
#' To read the file back in R use \code{\link[fst]{read.fst}}.

#' @section Input genomic datasets:
#' \enumerate{
#' \item VCF files must end with \code{.vcf}: documented in \code{\link{tidy_vcf}}
#'
#' \item PLINK files must end with \code{.tped}: documented in \code{\link{tidy_plink}}
#'
#' \item genind object from
#' \href{https://github.com/thibautjombart/adegenet}{adegenet}:
#' documented in \code{\link{tidy_genind}}.
#'
#' \item genlight object from
#' \href{https://github.com/thibautjombart/adegenet}{adegenet}:
#' documented in \code{\link{tidy_genlight}}.
#'
#' \item gtypes object from
#' \href{https://github.com/EricArcher/strataG}{strataG}:
#' documented in \code{\link{tidy_gtypes}}.
#'
#' \item dart data from \href{http://www.diversityarrays.com}{DArT}:
#' documented in \code{\link{tidy_dart}}.
#'
#' \item genepop file must end with \code{.gen}, documented in \code{\link{tidy_genepop}}.
#'
#' \item fstat file must end with \code{.dat}, documented in \code{\link{tidy_fstat}}.
#'
#' \item haplotype file created in STACKS (e.g. \code{data = "batch_1.haplotypes.tsv"}).
#' To make the haplotype file population ready, you need the \code{strata} argument.
#'
#' \item Data frames: documented in \code{\link{tidy_wide}}
#' }

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#'
#'
#' \item \code{vcf.metadata} (optional, logical or string).
#' Default: \code{vcf.metadata = TRUE}. Documented in \code{\link{tidy_vcf}}.
#'
#'
#' \item \code{vcf.stats} (optional, logical).
#' Default: \code{vcf.stats = TRUE}.
#' Documented in \code{\link{tidy_vcf}}.
#'
#'
#' \item \code{whitelist.markers} (optional, path or object) To keep only markers in a whitelist.
#' Default \code{whitelist.markers = NULL}.
#' Documented in \code{\link{read_whitelist}}.
#'
#' \item \code{blacklist.id} (optional) Default: \code{blacklist.id = NULL}.
#' Ideally, managed in the strata file.
#' Documented in \code{\link{read_strata}} and \code{\link{read_blacklist_id}}.
#'
#' \item \code{filter.common.markers} (optional, logical).
#' Default: \code{filter.common.markers = TRUE},
#' Documented in \code{\link{filter_common_markers}}.
#'
#' \item \code{filter.monomorphic} (logical, optional) Should the monomorphic
#' markers present in the dataset be filtered out ?
#' Default: \code{filter.monomorphic = TRUE}.
#' Documented in \code{\link{filter_monomorphic}}.
#' }




#' @export
#' @rdname tidy_genomic_data
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom stringi stri_join stri_replace_all_fixed stri_extract_all_fixed stri_replace_all_regex stri_sub stri_pad_left stri_count_fixed stri_replace_na
#' @importFrom stats var median quantile
#' @importFrom purrr map flatten keep discard
#' @importFrom data.table fread as.data.table
#' @importFrom tidyr spread gather unite separate
#' @importFrom utils count.fields
#' @importFrom readr write_tsv read_tsv
#' @importFrom tibble as_tibble
#' @importFrom rlang has_name
#' @importFrom parallel detectCores

#' @seealso \code{\link{detect_genomic_format}} and \code{\link{genomic_converter}}

#' @examples
#' \dontrun{
#' #To verify your file is detected by radiator as the correct format:
#' radiator::detect_genomic_format(data = "populations.snps.vcf")
#'
#'
#' # using VCF file as input
#' require(SeqVarTools)
#' tidy.vcf <- tidy_genomic_data(
#'    data = "populations.snps.vcf", strata = "strata.treefrog.tsv",
#'    whitelist.markers = "whitelist.vcf.txt")
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_genomic_data <- function(
  data,
  strata = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {
  if (verbose) {
    cat("################################################################################\n")
    cat("######################### radiator::tidy_genomic_data ##########################\n")
    cat("################################################################################\n")
  }
  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date/time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("######################### tidy_genomic_data completed ##########################\n"), add = TRUE)

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("path.folder", "parameters","keep.allele.names", "blacklist.id",
                "whitelist.markers", "filter.common.markers",
                "filter.monomorphic", "vcf.metadata", "vcf.stats",
                "blacklist.genotypes", "internal"),
    verbose = verbose
  )

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data is missing")

  if (!gt.vcf.nuc && !gt) {
    rlang::abort("At least one of gt.vcf.nuc or gt must be TRUE")
  }

  # Folders---------------------------------------------------------------------
  path.folder <- generate_folder(
    f = path.folder,
    rad.folder = "radiator_tidy_genomic",
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  # write the dots file
  write_rad(
    data = rad.dots,
    path = path.folder,
    filename = stringi::stri_join("radiator_tidy_genomic_data_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = internal,
    verbose = verbose
  )

  # File type detection----------------------------------------------------------
  skip.tidy.wide <- FALSE # initiate for data frame below
  data.type <- radiator::detect_genomic_format(data)

  # Import whitelist of markers-------------------------------------------------
  whitelist.markers <- read_whitelist(whitelist.markers, verbose)

  # Import blacklist id --------------------------------------------------------
  blacklist.id <- read_blacklist_id(blacklist.id, verbose)

  # Strata----------------------------------------------------------------------
  strata.df <- read_strata(
    strata = strata,
    pop.id = TRUE,
    blacklist.id = blacklist.id,
    verbose = verbose) %$% strata

  # GDS file -------------------------------------------------------------------
  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install SeqVarTools for this option:\n
                   install.packages("BiocManager")
                   BiocManager::install("SeqVarTools")')
    }

    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = verbose)
    }
    data <- gds2tidy(gds = data, parallel.core = parallel.core)
    data.type <- "tbl_df"
  }
  # Import VCF------------------------------------------------------------------
  if (data.type == "vcf.file") { # VCF
    if (verbose) message("Importing and tidying the VCF...")
    input <- radiator::tidy_vcf(
      data = data,
      strata = strata.df,
      vcf.metadata = vcf.metadata,
      parallel.core = parallel.core,
      verbose = verbose,
      whitelist.markers = whitelist.markers,
      filename = NULL,
      vcf.stats = vcf.stats,
      gt.vcf.nuc = TRUE,
      gt.vcf = TRUE,
      gt = TRUE,
      gt.bin = TRUE,
      path.folder = path.folder
    )
    biallelic <- radiator::detect_biallelic_markers(input)
    # biallelic <- input$biallelic
  } # End import VCF

  # Import PLINK ---------------------------------------------------------------
  if (data.type == "plink.file") { # PLINK
    if (verbose) message("Importing the PLINK files...")

    input <- tidy_plink(
      data = data,
      strata = strata.df,
      verbose = verbose,
      whitelist.markers = whitelist.markers,
      blacklist.id = blacklist.id
    )

    biallelic <- input$biallelic
    input <- input$input
  } # End import PLINK

  # Import stacks haplotypes----------------------------------------------------
  if (data.type == "haplo.file") { # Haplotype file
    if (verbose) message("Importing STACKS haplotype file")

    strata.df <- strata_haplo(
      strata = strata.df,
      data = data,
      blacklist.id = blacklist.id)

    # import header row
    want <- tibble::tibble(
      INFO = "CATALOG",
      COL_TYPE = "c") %>%
      dplyr::bind_rows(
        dplyr::select(strata.df, INFO = INDIVIDUALS) %>%
          dplyr::mutate(
            COL_TYPE = rep("c", n()),
            INFO = clean_ind_names(INFO)
          ))

    haplo.col.type <- readr::read_tsv(
      file = data,
      n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character())) %>%
      tidyr::gather(data = .,key = DELETE, value = INFO) %>%
      dplyr::mutate(INFO = clean_ind_names(INFO)) %>%
      dplyr::select(-DELETE) %>%
      dplyr::mutate(INFO = clean_ind_names(INFO)) %>%
      dplyr::left_join(want, by = "INFO") %>%
      dplyr::mutate(COL_TYPE = stringi::stri_replace_na(str = COL_TYPE, replacement = "_")) %>%
      dplyr::select(COL_TYPE)

    haplo.col.type[1,1] <- "c"

    haplo.col.type <- purrr::flatten_chr(haplo.col.type) %>% stringi::stri_join(collapse = "")

    # readr now faster/easier than fread...
    input <- readr::read_tsv(
      file = data, col_names = TRUE, na = "-",
      col_types = haplo.col.type)

    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = c("# Catalog ID", "Catalog ID", "# Catalog Locus ID"),
      replacement = c("LOCUS", "LOCUS", "LOCUS"), vectorize_all = FALSE)

    if (rlang::has_name(input, "Seg Dist")) {
      input <- dplyr::select(.data = input, -`Seg Dist`)
    }

    n.catalog.locus <- dplyr::n_distinct(input$LOCUS)
    n.individuals <- ncol(input) - 1

    message("\nNumber of loci in catalog: ", n.catalog.locus)
    message("Number of individuals: ", n.individuals)
    input <- tidyr::gather(
      data = input,
      key = "INDIVIDUALS",
      value = "GT_VCF_NUC", # previously using "GT_HAPLO"
      -LOCUS)

    input$INDIVIDUALS <- radiator::clean_ind_names(input$INDIVIDUALS)

    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      input <- filter_whitelist(
        data = input, whitelist.markers = whitelist.markers)
    }

    # remove consensus markers
    if (verbose) message("\nScanning for consensus markers...")
    consensus.markers <- dplyr::filter(input, GT_VCF_NUC == "consensus") %>%
      dplyr::distinct(LOCUS)

    if (length(consensus.markers$LOCUS) > 0) {
      input <- suppressWarnings(dplyr::anti_join(input, consensus.markers, by = "LOCUS"))
      readr::write_tsv(consensus.markers, "radiator.tidy.genomic.data.consensus.markers.tsv")
    }
    if (verbose) message("    number of consensus markers removed: ", dplyr::n_distinct(consensus.markers$LOCUS))
    consensus.markers <- NULL

    # population levels and strata
    if (!is.null(strata)) {
      if (rlang::has_name(input, "POP_ID")) {
        input %<>% dplyr::select(-POP_ID)
      }
      input %<>% dplyr::left_join(strata.df, by = "INDIVIDUALS")
      check.ref <- TRUE
    }
    # # using pop.levels and pop.labels info if present
    # input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)

    # # Pop select
    # if (!is.null(pop.select)) {
    #   if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
    #   input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    # }

    # removing errors and potential paralogs (GT with > 2 alleles)
    if (verbose) message("Scanning for artifactual genotypes...")
    input <- input %>%
      dplyr::mutate(POLYMORPHISM = stringi::stri_count_fixed(GT_VCF_NUC, "/"))

    blacklist.paralogs <- input %>%
      dplyr::filter(POLYMORPHISM > 1) %>%
      dplyr::select(LOCUS, INDIVIDUALS)

    if (verbose) message("    number of genotypes with more than 2 alleles: ", length(blacklist.paralogs$LOCUS))
    if (length(blacklist.paralogs$LOCUS) > 0) {
      input <- input %>%
        dplyr::mutate(GT_VCF_NUC = replace(GT_VCF_NUC, which(POLYMORPHISM > 1), NA)) %>%
        dplyr::select(-POLYMORPHISM)

      readr::write_tsv(blacklist.paralogs, "blacklist.genotypes.paralogs.tsv")
    }
    blacklist.paralogs <- NULL

    if (verbose) message("Calculating REF/ALT alleles...")
    # Prep for REF/ALT alleles and new genotype coding
    # part below could be parallelized if necessary, test with larger dataset for bottleneck...
    input <- input %>%
      dplyr::mutate(
        GT_VCF_NUC = dplyr::if_else(
          POLYMORPHISM == 0,
          stringi::stri_join(GT_VCF_NUC, "/", GT_VCF_NUC), GT_VCF_NUC,
          missing = "./."),
        GT_VCF_NUC = dplyr::if_else(stringi::stri_detect_fixed(GT_VCF_NUC, "N"),
                                    "./.", GT_VCF_NUC)
      ) %>%
      dplyr::select(-POLYMORPHISM)

    input.temp <- radiator::calibrate_alleles(
      data = input,
      biallelic = FALSE,
      parallel.core = parallel.core,
      verbose = verbose)
    input <- input.temp$input
    input.temp <- NULL
    biallelic <- FALSE
    input <- dplyr::rename(input, LOCUS = MARKERS)
  } # End import haplotypes file

  # Import genepop--------------------------------------------------------------
  if (data.type == "genepop.file") {
    if (verbose) message("Tidying the genepop file ...")
    input <- radiator::tidy_genepop(data = data, tidy = TRUE)
    skip.tidy.wide <- TRUE
  }

  # Import DArT ----------------------------------------------------------------
  if (data.type == "dart") {
    if (verbose) message("Tidying DArT data...")
    input <- radiator::tidy_dart(
      data = data,
      strata = strata,
      verbose = FALSE,
      parallel.core = parallel.core)
    skip.tidy.wide <- TRUE

    # if (rlang::has_name(strata.df, "NEW_ID")) {
    #   strata.df <- strata.df %>%
    #     dplyr::select(-INDIVIDUALS) %>%
    #     dplyr::rename(INDIVIDUALS = NEW_ID)
    # }
  }# End dart

  # Import fst.file ------------------------------------------------------------
  if (data.type == "fst.file") {
    if (verbose) message("Importing the fst.file as a data frame...")
    input <- read_rad(data = data)
    skip.tidy.wide <- TRUE
  } # End fst.file

  # Import GENIND--------------------------------------------------------------
  if (data.type == "genind") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the genind object ...")
    input <- radiator::tidy_genind(data = data, gds = FALSE,
                                   keep.allele.names = keep.allele.names)
    data <- NULL
    # remove unwanted sep in id and pop.id names
    input <- input %>%
      dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                       .funs = clean_ind_names) %>%
      dplyr::mutate_at(.tbl = ., .vars = "POP_ID",
                       .funs = clean_pop_names)
    skip.tidy.wide <- TRUE
  } # End tidy genind

  # Import GENLIGHT ------------------------------------------------------------
  if (data.type == "genlight") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the genlight object ...")
    input <- radiator::tidy_genlight(data = data, gds = FALSE)
    data <- NULL
    # remove unwanted sep in id and pop.id names
    input <- input %>%
      dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                       .funs = clean_ind_names) %>%
      dplyr::mutate_at(.tbl = ., .vars = "POP_ID",
                       .funs = clean_pop_names)
    biallelic <- TRUE
    skip.tidy.wide <- TRUE
  } # End tidy genlight

  # Import STRATAG gtypes ------------------------------------------------------
  if (data.type == "gtypes") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the gtypes object ...")
    input <- tidy_gtypes(data) %>%
      dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                       .funs = clean_ind_names) %>%
      dplyr::mutate_at(.tbl = ., .vars = "POP_ID",
                       .funs = clean_pop_names)
    data <- NULL
    skip.tidy.wide <- TRUE
  } # End tidy gtypes

  # Import DF-------------------------------------------------------------------
  if (data.type == "tbl_df" || skip.tidy.wide) { # DATA FRAME OF GENOTYPES
    if (verbose) message("Importing the data frame ...")
    if (!skip.tidy.wide) {
      input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
      data <- NULL
    }

    if (!is.null(whitelist.markers)) {
      input <- filter_whitelist(
        data = input, whitelist.markers = whitelist.markers)
    }

    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      if (verbose) message("Filtering with blacklist of individuals")
      input  %<>% dplyr::filter(!INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
      check.ref <- TRUE
    }


    # population levels and strata
    if (!is.null(strata)) {
      if (rlang::has_name(input, "POP_ID")) {
        input %<>% dplyr::select(-POP_ID)
      }
      input %<>% dplyr::left_join(strata.df, by = "INDIVIDUALS")
      check.ref <- TRUE
    }

    # initiate ref check
    if (rlang::has_name(input, "REF")) {
      check.ref <- FALSE
    } else {
      check.ref <- TRUE
    }

    if (check.ref) {
      input.temp <- radiator::calibrate_alleles(data = input)
      input <- input.temp$input
      biallelic <- input.temp$biallelic
      input.temp <- NULL
    } else {
      biallelic <- radiator::detect_biallelic_markers(data = input)
    }

  } # End import data frame of genotypes

  # END IMPORT DATA-------------------------------------------------------------

  # strata integration ---------------------------------------------------------
  if (!is.null(strata)) {
    strata.df <- dplyr::ungroup(input) %>%
      dplyr::distinct(POP_ID, INDIVIDUALS)
  }

  # Blacklist genotypes --------------------------------------------------------
  if (is.null(blacklist.genotypes)) { # no Whitelist
    if (verbose) message("Erasing genotype: no")
  } else {
    input <- filter_blacklist_genotypes(
      data = input,
      blacklist.genotypes = blacklist.genotypes,
      verbose = verbose)
  } # End erase genotypes

  # dump unused object
  blacklist.id <- whitelist.markers <- whitelist.markers.ind <- NULL
  want <- blacklist.genotypes <- NULL


  # Unique markers id ----------------------------------------------------------
  # we want to keep LOCUS in the vcf, but not in the other type of input file
  # if (rlang::has_name(input, "LOCUS") && !rlang::has_name(input, "MARKERS")) {
  #   colnames(input) <- stringi::stri_replace_all_fixed(
  #     str = colnames(input),
  #     pattern = "LOCUS",
  #     replacement = "MARKERS",
  #     vectorize_all = FALSE
  #   )
  # }

  # filter.common.markers ------------------------------------------------------
  if (is.null(strata)) filter.common.markers <- FALSE
  if (filter.common.markers) {
    input <- filter_common_markers(
      data = input,
      verbose = TRUE,
      path.folder = path.folder,
      parameters = parameters)
  } # End common markers

  # filter_monomorphic----------------------------------------------------------
  if (filter.monomorphic) {
    input <- filter_monomorphic(
      data = input,
      verbose = TRUE,
      path.folder = path.folder,
      parameters = parameters,
      internal = TRUE)
  } # End filter.monomorphic


  # Results --------------------------------------------------------------------
  if (!is.null(strata)) {
    input %<>% dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  } else {
    input %<>% dplyr::arrange(INDIVIDUALS, MARKERS)
  }

  # Write to working directory -------------------------------------------------
  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".rad")
    message("\nWriting tidy data set:\n", tidy.name)
    write_rad(data = input, path = file.path(path.folder, tidy.name))
  }
  # tidy.name <- stringi::stri_join(filename, ".rad")
  # message("\nWriting tidy data set:\n", tidy.name)
  # write_rad(data = input, path = tidy.name)

  n.markers <- length(unique(input$MARKERS))
  if (rlang::has_name(input, "CHROM")) {
    n.chromosome <- length(unique(input$CHROM))
  } else {
    n.chromosome <- "no chromosome info"
  }
  n.individuals <- length(unique(input$INDIVIDUALS))
  if(!is.null(strata)) n.pop <- length(unique(input$POP_ID))

  if (verbose) {
    cat("################################### RESULTS ####################################\n")
    if (!is.null(filename)) {
      message("Tidy data written in global environment and working directory")
    } else {
      message("Tidy data written in global environment")
    }
    message("Data format: ", data.type)
    if (biallelic) {
      message("Biallelic data")
    } else{
      message("Multiallelic data")
    }

    message("\nTidy genomic data:")
    message("    Number of markers: ", n.markers)
    message("    Number of chromosome/contig/scaffold: ", n.chromosome)
    message("    Number of individuals: ", n.individuals)
    if (!is.null(strata)) message("    Number of populations: ", n.pop)
  }
  return(input)
} # tidy genomic data


# Internal nested Function -----------------------------------------------------

#' @title strata_haplo
#' @description Manage strata
#' @rdname strata_haplo
#' @keywords internal
#' @export
strata_haplo <- function(strata = NULL, data = NULL, blacklist.id = NULL) {

  if (is.null(strata)) {
    message("No strata file provided")
    message("    generating a strata with 1 grouping")
    if (is.null(data)) rlang::abort("data required to generate strata")
    strata.df <- readr::read_tsv(
      file = data,
      n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character())) %>%
      tidyr::gather(data = .,key = DELETE, value = INDIVIDUALS) %>%
      dplyr::mutate(INDIVIDUALS = clean_ind_names(INDIVIDUALS)) %>%
      dplyr::select(-DELETE) %>%
      dplyr::filter(!INDIVIDUALS %in% c("Catalog ID", "Cnt")) %>%
      dplyr::distinct(INDIVIDUALS) %>%
      dplyr::mutate(STRATA = rep("pop1", n()))
  } else {
    if (is.vector(strata)) {
      suppressMessages(
        strata.df <- readr::read_tsv(
          file = strata, col_names = TRUE,
          # col_types = col.types
          col_types = readr::cols(.default = readr::col_character())
        ))
    } else {
      strata.df <- strata
    }
  }

  colnames(strata.df) <- stringi::stri_replace_all_fixed(
    str = colnames(strata.df),
    pattern = "STRATA",
    replacement = "POP_ID",
    vectorize_all = FALSE
  )
  # Remove potential whitespace in pop_id
  strata.df$POP_ID <- clean_pop_names(strata.df$POP_ID)
  colnames.strata <- colnames(strata.df)

  # clean ids
  strata.df$INDIVIDUALS <- clean_ind_names(strata.df$INDIVIDUALS)

  # filtering the strata if blacklist id available
  if (!is.null(blacklist.id)) {
    strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
  }

  strata.df <- dplyr::distinct(strata.df, POP_ID, INDIVIDUALS)
  return(strata.df)
}#End strata_haplo
