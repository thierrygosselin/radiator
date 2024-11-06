# tidy_wide --------------------------------------------------------------------
# read data frames in long/tidy and wide format

#' @name tidy_wide

#' @title Read/Import a tidy genomic data frames.

#' @description Read/Import and tidy genomic data frames. If data is in
#' wide format, the functions will gather the data.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A file in the working directory or object in the global environment
#' in wide or long (tidy) formats. See details for more info.
#'
#' \emph{How to get a tidy data frame ?}
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' \code{\link{tidy_genomic_data}}.

#' @param import.metadata (optional, logical) With \code{import.metadata = TRUE}
#' the metadata (anything else than the genotype) will be imported for the long
#' format exclusively. Default: \code{import.metadata = FALSE}, no metadata.

#' @return A tidy data frame in the global environment.
#' @export
#' @rdname tidy_wide

#' @details \strong{Input data:}
#'
#' To discriminate the long from the wide format,
#' the function \pkg{radiator} \code{\link[radiator]{tidy_wide}} searches
#' for \code{MARKERS} in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns:
#' \code{INDIVIDUALS}, \code{STRATA} (that refers to any grouping of individuals),
#' the remaining columns are the markers in separate columns storing genotypes.
#'
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info.
#' (e.g. from a VCF see \pkg{radiator} \code{\link{tidy_genomic_data}}).
#' A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{STRATA},
#' \code{MARKERS} and \code{GT} for the genotypes.
#' The remaining columns are considered metadata info.
#'
#' \strong{Genotypes with separators:}
#' ALL separators will be removed.
#' Genotypes should be coded with 3 integers for each alleles.
#' 6 integers in total for the genotypes.
#' e.g. \code{001002 or 111333} (for heterozygote individual).
#' 6 integers WITH separator: e.g. \code{001/002 or 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}, and will
#' be removed.
#'
#'
#' \strong{separators in STRATA, INDIVIDUALS and MARKERS:}
#' Some separators can interfere with packages or codes and are cleaned by radiator.
#' \itemize{
#' \item MARKERS: \code{/}, \code{:}, \code{-} and \code{.} are changed to an
#' underscore
#' \code{_}.
#' \item STRATA: white spaces in population names are replaced by underscore.
#' \item INDIVIDUALS: \code{_} and \code{:} are changed to a dash \code{-}
#' }
#'
#' \emph{How to get a tidy data frame ?}
#' \pkg{radiator} \code{\link{tidy_genomic_data}} can transform 6 genomic data formats
#' in a tidy data frame.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_wide <- function(data, import.metadata = FALSE) {

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) rlang::abort("Input file argument is missing")

  if (is.vector(data)) {# for file in the working directory
    if (stringi::stri_detect_fixed(
      str = stringi::stri_sub(str = data, from = -4, to = -1),
      pattern = ".tsv")) {
      data %<>% readr::read_tsv(file = ., col_types = readr::cols(.default = readr::col_character()))
    } else if (radiator::detect_genomic_format(data) == "fst.file") {
      data %<>% radiator::read_rad(data = .)
    }
  }

  # No longer supported, remove POP_ID
  data %<>% dplyr::rename(STRATA = tidyselect::any_of("POP_ID"))

  # Determine long (tidy) or wide dataset
  if (rlang::has_name(data, "LOCUS") && !rlang::has_name(data, "MARKERS")) {
    data %<>% dplyr::rename(MARKERS = LOCUS)
  }

  if (!"MARKERS" %in% colnames(data) && !"LOCUS" %in% colnames(data)) {

    grouping.cols <- "INDIVIDUALS"
    if (rlang::has_name(data, "STRATA")) grouping.cols <- c("STRATA", "INDIVIDUALS")
    data %<>%
      radiator::rad_long(
        x = .,
        cols = grouping.cols,
        names_to = "MARKERS",
        values_to = "GT",
        variable_factor = FALSE
      )
  }

  if (!import.metadata) {
    want <- c("STRATA", "INDIVIDUALS", "MARKERS", "CHROM", "LOCUS", "POS", "GT",
              "GT_VCF_NUC", "GT_VCF", "GT_BIN")
    data %<>% dplyr::select(tidyselect::any_of(want))
  }

  # cleaning...
  if (rlang::has_name(data, "MARKERS")) data$MARKERS %<>% clean_markers_names(.)
  data$INDIVIDUALS %<>% clean_ind_names(.)
  if (rlang::has_name(data, "STRATA")) data$STRATA %<>% clean_pop_names(.)

  # Make sure no data groupings exists
  data %<>% dplyr::ungroup(.)
  return(data)
}#End tidy_wide


# tidy2wide --------------------------------------------------------------------
#' @title tidy2wide
#' @description tidy2wide
#' @keywords internal
#' @export
tidy2wide <- function(
  x = NULL,
  gds = NULL,
  individuals = NULL,
  markers = NULL,
  tidy = TRUE,
  wide = TRUE,
  wide.markers = TRUE
) {
  res <- list()
  if (is.null(markers) && !is.null(gds)) {
    markers <- radiator::extract_markers_metadata(
      gds = gds,
      markers.meta.select = "MARKERS",
      whitelist = TRUE
    ) %>%
      dplyr::pull()
  }
  if (is.null(individuals) && !is.null(gds)) {
    individuals <- radiator::extract_individuals_metadata(
      gds = gds,
      ind.field.select = "INDIVIDUALS",
      whitelist = TRUE
    ) %>%
      dplyr::pull()
  }

  n.markers <- length(markers)
  n.ind <- length(individuals)


  # sample and markers in long format...
  res$data.tidy <- tibble::tibble(
    INDIVIDUALS = rep(individuals, n.markers),
    MARKERS = sort(rep(markers, n.ind)),
  )

  if (wide) {
    if (wide.markers) {
      res$data.wide <- radiator::rad_wide(
        x = res$data.tidy,
        formula = "INDIVIDUALS ~ MARKERS",
        values_from = names(x)
      )
    } else {
      res$data.wide <- radiator::rad_wide(
        x = res$data.tidy,
        formula = "MARKERS ~ INDIVIDUALS",
        values_from = names(x)
      )
    }
  }
  if (!tidy) res$data.tidy <- NULL
  return(res)
}# End tidy2wide




# tidy_genomic_data ------------------------------------------------------------

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

#' @param data 14 options for input (\strong{diploid data only}): VCFs (SNPs or Haplotypes,
#' to make the vcf population ready),
#' plink (tped, bed), stacks haplotype file, genind (library(adegenet)),
#' genlight (library(adegenet)), gtypes (library(strataG)), genepop, DArT,
#' and a data frame in long/tidy or wide format. To verify that radiator detect
#' your file format use \code{\link{detect_genomic_format}} (see example below).
#' Documented in \strong{Input genomic datasets} of \code{\link{tidy_genomic_data}}.
#'
#' \strong{DArT and VCF data}: \pkg{radiator} was not meant to generate alleles
#' and genotypes if you are using a VCF file with no genotype
#' (only genotype likelihood: GL or PL).
#' Neither is \pkg{radiator} able to magically generate a genind object
#' from a SilicoDArT dataset. Please look at the first few lines of your dataset
#' to understand it's limit before asking raditor to convert or filter your dataset.

#' @param strata (optional)
#' The strata file is a tab delimited file with a minimum of 2 columns headers:
#' \code{INDIVIDUALS} and \code{STRATA}. Documented in \code{\link{read_strata}}.
#' DArT data: a third column \code{TARGET_ID} is required.
#' Documented on \code{\link{read_dart}}. Also use the strata read function to
#' blacklist individuals.
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
#' \item PLINK files must end with \code{.tped} or \code{.bed}: documented in \code{\link{tidy_plink}}
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
#' documented in \code{\link{read_dart}}.
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

#' @seealso \code{\link{detect_genomic_format}} and \code{\link{genomic_converter}}

#' @examples
#' \dontrun{
#' #To verify your file is detected by radiator as the correct format:
#' radiator::detect_genomic_format(data = "populations.snps.vcf")
#'
#'
#' # using VCF file as input
#' require(SeqArray)
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

  # test
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE
  # path.folder = NULL
  # parameters = NULL
  # keep.allele.names = NULL
  # blacklist.id = NULL
  # whitelist.markers = NULL
  # filter.common.markers = TRUE
  # filter.monomorphic = TRUE
  # vcf.metadata = TRUE
  # vcf.stats = TRUE
  # blacklist.genotypes = NULL
  # internal = FALSE

  # Cleanup-------------------------------------------------------------------
  radiator_function_header(f.name = "tidy_genomic_data", verbose = verbose)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- radiator_tic()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing, verbose = verbose), add = TRUE)
  on.exit(radiator_function_header(f.name = "tidy_genomic_data", start = FALSE, verbose = verbose), add = TRUE)
  res <- list()

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
    verbose = FALSE
  )

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data is missing")

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
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  # File type detection----------------------------------------------------------
  skip.tidy.wide <- FALSE # initiate for data frame below
  data.type <- radiator::detect_genomic_format(data)

  # Import whitelist of markers-------------------------------------------------
  whitelist.markers %<>% read_whitelist(whitelist.markers = ., verbose)

  # Import blacklist id --------------------------------------------------------
  blacklist.id %<>% read_blacklist_id(blacklist.id = ., verbose)

  # Strata----------------------------------------------------------------------
  strata.df <- read_strata(
    strata = strata,
    pop.id = FALSE,
    # pop.id = TRUE,
    blacklist.id = blacklist.id,
    verbose = verbose) %$% strata


  # Import and tidy files ----------------------------------------------------

  # GDS file -------------------------------------------------------------------
  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

    if (data.type == "gds.file") data %<>% radiator::read_rad(data = ., verbose = verbose)
    data %<>% gds2tidy(gds = ., parallel.core = parallel.core)
    data.type <- "tbl_df"
  }

  # Import VCF------------------------------------------------------------------
  if (data.type == "vcf.file") { # VCF
    if (verbose) message("Importing and tidying the VCF...")
    input <- radiator::tidy_vcf(
      data = data,
      strata = strata.df,
      parallel.core = parallel.core,
      verbose = verbose,
      whitelist.markers = whitelist.markers,
      filter.monomorphic = filter.monomorphic,
      filter.common.markers = filter.common.markers,
      filename = NULL,
      vcf.metadata = vcf.metadata,
      vcf.stats = vcf.stats,
      gt.vcf.nuc = TRUE,
      gt.vcf = TRUE,
      gt = TRUE,
      gt.bin = TRUE,
      path.folder = path.folder,
      internal = TRUE,
      tidy.check = FALSE
    )
    biallelic <- radiator::detect_biallelic_markers(input)
  } # End import VCF

  # Import PLINK ---------------------------------------------------------------
  if (data.type %in% c("plink.tped.file", "plink.bed.file")) { # PLINK
    if (verbose) message("Importing the PLINK file...")

    input <- radiator::tidy_plink(
      data = data,
      parallel.core = parallel.core,
      verbose = verbose
    )

    biallelic <- input$biallelic
    input <- input$input
    strata <- "not.null" # to trigger

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
      col_types = readr::cols(.default = readr::col_character())
    ) %>%
      radiator::rad_long(
        x = .,
        cols = tidyselect::everything(),
        names_to = "DELETE",
        values_to = "DELETE",
        variable_factor = FALSE
      ) %>%
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
      pattern = c("# Catalog ID", "Catalog ID", "# Catalog Locus ID", "Catalog.ID"),
      replacement = c("LOCUS", "LOCUS", "LOCUS", "LOCUS"),
      vectorize_all = FALSE
    )

    if (rlang::has_name(input, "Seg Dist")) input %<>% dplyr::select(-`Seg Dist`)

    n.catalog.locus <- dplyr::n_distinct(input$LOCUS)
    n.individuals <- ncol(input) - 1

    message("\nNumber of loci in catalog: ", n.catalog.locus)
    message("Number of individuals: ", n.individuals)

    input %<>%
      radiator::rad_long(
        x = .,
        cols = "LOCUS",
        names_to = "INDIVIDUALS",
        values_to =  "GT_VCF_NUC",
        variable_factor = FALSE
      )

    input$INDIVIDUALS %<>% radiator::clean_ind_names(.)

    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      input %<>%
        filter_whitelist(data = ., whitelist.markers = whitelist.markers)
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
      input %<>% join_strata(strata = strata.df)
      check.ref <- TRUE
    }

    # removing errors and potential paralogs (GT with > 2 alleles)
    if (verbose) message("Scanning for artifactual genotypes...")
    input %<>%
      dplyr::mutate(POLYMORPHISM = stringi::stri_count_fixed(GT_VCF_NUC, "/"))

    blacklist.paralogs <- input %>%
      dplyr::filter(POLYMORPHISM > 1) %>%
      dplyr::select(LOCUS, INDIVIDUALS)

    if (verbose) message("    number of genotypes with more than 2 alleles: ", length(blacklist.paralogs$LOCUS))
    if (length(blacklist.paralogs$LOCUS) > 0) {
      input %<>%
        dplyr::mutate(GT_VCF_NUC = replace(GT_VCF_NUC, which(POLYMORPHISM > 1), NA))
      readr::write_tsv(blacklist.paralogs, "blacklist.genotypes.paralogs.tsv")
    }
    blacklist.paralogs <- NULL

    if (verbose) message("Calculating REF/ALT alleles...")
    # Prep for REF/ALT alleles and new genotype coding
    # part below could be parallelized if necessary, test with larger dataset for bottleneck...
    input %<>%
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
      verbose = verbose)
    input <- input.temp$input
    input.temp <- NULL
    biallelic <- FALSE

    # test
    # Add a column MARKERS
    # input <- dplyr::rename(input, LOCUS = MARKERS)
    input %<>% dplyr::mutate(LOCUS = MARKERS)



  } # End import haplotypes file

  # Import genepop--------------------------------------------------------------
  if (data.type == "genepop.file") {
    if (verbose) message("Tidying the genepop file ...")
    input <- radiator::tidy_genepop(data = data, tidy = TRUE)
    skip.tidy.wide <- TRUE
    if (is.null(strata) && is.null(strata.df)) {
      strata.df <- strata <- radiator::generate_strata(data = input)
    }
  }

  # Import DArT ----------------------------------------------------------------
  if (data.type == "dart") {
    if (verbose) message("Tidying DArT data...")
    input <- radiator::read_dart(
      data = data,
      strata = strata,
      verbose = FALSE,
      parallel.core = parallel.core,
      tidy.dart = TRUE)
    skip.tidy.wide <- TRUE
  }# End dart

  # Import GENIND--------------------------------------------------------------
  if (data.type == "genind") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the genind object ...")
    input <- radiator::tidy_genind(data = data, gds = FALSE)
    data <- NULL
    # remove unwanted sep in id and pop.id names
    input %<>%
      dplyr::mutate(dplyr::across(.cols = "INDIVIDUALS", .fns = clean_ind_names)) %>%
      dplyr::mutate(dplyr::across(.cols = "STRATA", .fns = clean_pop_names))
    skip.tidy.wide <- TRUE
  } # End tidy genind

  # Import GENLIGHT ------------------------------------------------------------
  if (data.type == "genlight") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the genlight object ...")
    input <- radiator::tidy_genlight(data = data, gds = FALSE)
    data <- NULL
    # remove unwanted sep in id and pop.id names
    input %<>%
      dplyr::mutate(dplyr::across(.cols = "INDIVIDUALS", .fns = clean_ind_names)) %>%
      dplyr::mutate(dplyr::across(.cols = "STRATA", .fns = clean_pop_names))
    biallelic <- TRUE
    skip.tidy.wide <- TRUE
  } # End tidy genlight

  # Import STRATAG gtypes ------------------------------------------------------
  if (data.type == "gtypes") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the gtypes object ...")
    input <- tidy_gtypes(data) %>%
      dplyr::mutate(dplyr::across(.cols = "INDIVIDUALS", .fns = clean_ind_names)) %>%
      dplyr::mutate(dplyr::across(.cols = "STRATA", .fns = clean_pop_names))
    data <- NULL
    skip.tidy.wide <- TRUE
  } # End tidy gtypes

  # Import fst.file ------------------------------------------------------------
  if (data.type == "fst.file") {
    if (verbose) message("Importing the fst.file...")
    input <- read_rad(data = data)
    data.type <- "tbl_df"
    skip.tidy.wide <- TRUE
  } # End fst.file

  # Import DF-------------------------------------------------------------------
  if (data.type == "tbl_df" || skip.tidy.wide) { # DATA FRAME OF GENOTYPES
    if (!skip.tidy.wide) {
      if (verbose) message("Importing the data frame ...")
      input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
      data <- NULL
    }

    if (!is.null(whitelist.markers)) {
      input %<>% filter_whitelist(data = ., whitelist.markers = whitelist.markers)
    }

    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      if (verbose) message("Filtering with blacklist of individuals")
      input  %<>% dplyr::filter(!INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
      check.ref <- TRUE
    }


    # population levels and strata
    if (!is.null(strata)) {
      input  %<>% join_strata(strata = strata.df)
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
  if (is.null(strata)) filter.common.markers <- FALSE

  # Blacklist genotypes --------------------------------------------------------
  if (is.null(blacklist.genotypes)) { # no Whitelist
    if (verbose) message("Erasing genotype: no")
  } else {
    input %<>%
      filter_blacklist_genotypes(
        data = .,
        blacklist.genotypes = blacklist.genotypes,
        verbose = verbose
      )
  } # End erase genotypes

  # dump unused object
  blacklist.id <- whitelist.markers <- whitelist.markers.ind <- NULL
  want <- blacklist.genotypes <- NULL

  # radiator_parameters-----------------------------------------------------------
  filters.parameters <- radiator_parameters(
    generate = TRUE,
    initiate = TRUE,
    update = FALSE,
    parameter.obj = parameters,
    data = input,
    path.folder = path.folder,
    file.date = file.date,
    internal = FALSE,
    verbose = verbose)

  # filter_common_markers ------------------------------------------------------
  input %<>%
    filter_common_markers(
      data = .,
      filter.common.markers = filter.common.markers,
      verbose = verbose,
      path.folder = path.folder,
      parameters = filters.parameters,
      internal = TRUE)

  # filter_monomorphic----------------------------------------------------------
  input %<>%
    filter_monomorphic(
      data = .,
      filter.monomorphic = filter.monomorphic,
      verbose = verbose,
      path.folder = path.folder,
      parameters = filters.parameters,
      internal = TRUE)


  # Results --------------------------------------------------------------------
  # is this really necessary ? Very long with huge dataset... (bottleneck even)
  # turned off 20201212
  # if (!is.null(strata)) {
  #   input %<>% dplyr::arrange(STRATA, INDIVIDUALS, MARKERS)
  # } else {
  #   input %<>% dplyr::arrange(INDIVIDUALS, MARKERS)
  # }

  # Write to working directory -------------------------------------------------
  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".rad")
    message("\nWriting tidy data set:\n", tidy.name)
    write_rad(data = input, path = file.path(path.folder, tidy.name))
  }
  n.markers <- length(unique(input$MARKERS))
  if (rlang::has_name(input, "CHROM")) {
    n.chromosome <- length(unique(input$CHROM))
  } else {
    n.chromosome <- "no chromosome info"
  }
  n.individuals <- length(unique(input$INDIVIDUALS))
  if (!is.null(strata)) n.pop <- length(unique(input$STRATA))

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
    if (!is.null(strata)) message("    Number of strata: ", n.pop)
    message("    Number of individuals: ", n.individuals)
  }
  return(input)
} # tidy genomic data
