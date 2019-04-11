# extract_dart_target_id---------------------------------------------------------
#' @name extract_dart_target_id

#' @title Extract \href{http://www.diversityarrays.com}{DArT} target id

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users. The function allows to extract DArT
#' target id from a DArT file. To help prepare the appropriate STRATA file.

#' @inheritParams read_dart

#' @param write With default \code{write = TRUE}, the dart target id column is
#' written in a file in the working directory.

#' @return A tidy dataframe with a \code{TARGET_ID} column. Spaces are remove and
#' UPPER case is used.

#' @export
#' @rdname extract_dart_target_id

#' @examples
#' \dontrun{
#' # Built a strata file:
#' strata <- radiator::extract_dart_target_id("mt.dart.file.csv") %>%
#'     dplyr::mutate(
#'         INDIVIDUALS = "new id you want to give",
#'         STRATA = "fill this"
#'     ) %>%
#'     readr::write_tsv(x = ., path = "my.new.dart.strata.tsv")
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and Peter Grewe \email{peter.grewe@csiro.au}

extract_dart_target_id <- function(data, write = TRUE) {
  ##TEST
  # write = FALSE
  # write = TRUE

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # Check DArT format file -----------------------------------------------------
  dart.check <- check_dart(data)

  # Import data ----------------------------------------------------------------
  dart.target.id <- readr::read_delim(
    file = data,
    delim = dart.check$tokenizer.dart,
    skip = dart.check$skip.number,
    n_max = 1,
    na = "-",
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character())) %>%
    t %>%
    magrittr::set_colnames(x = ., value = "TARGET_ID") %>%
    tibble::as_tibble(.)


  if (dart.check$star.number > 0) {
    dart.target.id %<>%
      dplyr::filter(dplyr::row_number() > dart.check$star.number)
  } else {
    # This is a string of known DArT col header not wanted
    discard <- c(
      "TARGET_ID", "ALLELE_ID", "CLONE_ID", "CLUSTER_TEMP_INDEX", "ALLELE_SEQUENCE",
      "CLUSTER_CONSENSUS_SEQUENCE", "CLUSTER_SIZE", "ALLELE_SEQ_DIST", "SNP",
      "SNP_POSITION", "CALL_RATE", "ONE_RATIO_REF", "ONE_RATIO_SNP", "FREQ_HOM_REF",
      "FREQ_HOM_SNP", "FREQ_HETS", "PIC_REF", "PIC_SNP", "AVG_PIC", "AVG_COUNT_REF",
      "AVG_COUNT_SNP", "RATIO_AVG_COUNT_REF_AVG_COUNT_SNP",
      "FREQ_HETS_MINUS_FREQ_MIN_HOM",
      "ALLELE_COUNTS_CORRELATION", "AGGREGATE_TAGS_TOTAL", "DERIVED_CORR_MINUS_SEED_CORR",
      "REP_REF", "REP_SNP", "REP_AVG", "PIC_REP_REF", "PIC_REP_SNP", "TOTAL_PIC_REP_REF_TEST",
      "TOTAL_PIC_REP_SNP_TEST", "BIN_ID", "BIN_SIZE", "ALLELE_SEQUENCE_REF",
      "ALLELE_SEQUENCE_SNP", "TRIMMED_SEQUENCE_REF", "TRIMMED_SEQUENCE", "ONE_RATIO",
      "PIC", "AVG_READ_DEPTH", "STDEV_READ_DEPTH", "Q_PMR", "REPRODUCIBILITY", "MAF",
      "TOT_COUNTS", "TOTAL_PIC_REP_TEST", "PIC_REP")

    discard.genome <- c("CHROM_|CHROM_POS_|ALN_CNT_|ALN_EVALUE_")

    dart.target.id %<>%
      dplyr::filter(!radiator_snakecase(x = TARGET_ID) %in% discard) %>%
      dplyr::filter(!stringi::stri_detect_regex(
        str = stringi::stri_trans_toupper(TARGET_ID),
        pattern = discard.genome, negate = FALSE))
  }

  dart.target.id %<>%
    dplyr::mutate(
      TARGET_ID = clean_ind_names(x = TARGET_ID),
      TARGET_ID = stringi::stri_trans_toupper(TARGET_ID)
      )
  if (write) readr::write_tsv(x = dart.target.id, path = "dart.target.id.tsv")

  # Check that DArT file as good target id written -----------------------------
  if (nrow(dart.target.id) != length(unique(dart.target.id$TARGET_ID))) {
    rlang::warn(
      "Non unique TARGET_ID or sample names used in the DArT file.
Solution: edit manually")
  }
  return(dart.target.id)
}#End extract_dart_target_id


# read_dart -------------------------------------------------------------------------
# Import, filter and transform a dart output file to different formats

#' @name read_dart

#' @title Read and tidy \href{http://www.diversityarrays.com}{DArT} output files.

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users. The function generate a GDS object/file
#' and optionally, a tidy dataset using
#' \href{http://www.diversityarrays.com}{DArT} files.

#' @param data One of the DArT output files. 6 formats used by DArT are recognized
#' by radiator.
#' recognised:
#' \enumerate{
#' \item \code{1row}: Genotypes are in 1 row and coded (0, 1, 2, -).
#' \code{0 for 2 reference alleles REF/REF}, \code{1 for 2 alternate alleles ALT/ALT},
#' \code{2 for heterozygote REF/ALT}, \code{- for missing}.
#' \item \code{2rows}: No genotypes. It's absence/presence, 0/1, of the REF and ALT alleles.
#' Sometimes called binary format.
#' \item \code{counts}: No genotypes, It's counts/read depth for the REF and ALT alleles.
#' Sometimes just called count data.
#' \item \code{silico.dart}: SilicoDArT data. No genotypes, no REF or ALT alleles.
#' It's a file coded as absence/presence, 0/1, for the presence of sequence in
#' the clone id.
#' \item \code{silico.dart.counts}: SilicoDArT data. No genotypes, no REF or ALT alleles.
#' It's a file coded as absence/presence, with counts for the presence of sequence in
#' the clone id.
#' \item \code{dart.vcf}: For DArT VCFs, please use \code{\link{read_vcf}}.
#' }
#'
#' Depending on the number of markers, these format will be recoded similarly to
#' VCF files (dosage of alternate allele, see details).
#'
#' The function can import \code{.csv} or \code{.tsv} files.
#'
#'
#' If you encounter a problem, sent me your data so that I can update
#' the function.

#' @param strata A tab delimited file or object with 3 columns.
#' Columns header is:
#' \code{TARGET_ID}, \code{INDIVIDUALS} and \code{STRATA}.
#' Note: the column \code{STRATA} refers to any grouping of individuals.
#' You need to make sure that
#' the column \code{TARGET_ID} match the id used by DArT. With the \code{counts}
#' format the \code{TARGET_ID} is a series of integer.
#' With \code{1row} and \code{2rows} the \code{TARGET_ID} is actually the sample
#' name submitted to DArT.
#' The column \code{INDIVIDUALS} and \code{STRATA} will be kept in the tidy data.
#' Only individuals in the strata file are kept in the tidy, i.e. that the strata
#' is also used as a whitelist of individuals/strata.
#' Silico DArT data is currently used to detect sex markers, so the \code{STRATA}
#' column should be filed with sex information: \code{M} or \code{F}.
#'
#' See example on how to extract the TARGET_ID of your DArT file.
#'
#' \href{https://www.dropbox.com/s/utq2h6o00v55kep/example.dart.strata.tsv?dl=0)}{example.dart.strata.tsv}.

#' @param tidy.dart (logical, optional) Generate a tidy dataset.
#' Default:\code{tidy.dart = FALSE}.


#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments
#' @inheritParams filter_whitelist



#' @param ... (optional) To pass further argument for fine-tuning the function.

#' @return A radiator GDS file and tidy dataframe with several columns depending on DArT file:
#' \code{silico.dart:} A tibble with 5 columns: \code{CLONEID, SEQUENCE, VALUE, INDIVIDUALS, STRATA}.
#' This object is also saved in the directory (file ending with .rad).
#'
#' Common to \code{1row, 2rows and counts}: A GDS file is automatically generated.
#' To have a tidy tibble, the argument \code{tidy.dart = TRUE} must be used.
#'
#' \enumerate{
#' \item VARIANT_ID: generated by radiator and correspond the markers in integer.
#' \item MARKERS: generated by radiator and correspond to CHROM + LOCUS + POS separated by 2 underscores.
#' \item CHROM: the chromosome info, for de novo: CHROM_1.
#' \item LOCUS: the locus info.
#' \item POS: the SNP id on the LOCUS.
#' \item COL: the position of the SNP on the short read.
#' \item REF: the reference allele.
#' \item ALT: the alternate allele.
#' \item INDIVIDUALS: the sample name.
#' \item STRATA/POP_ID: populations id of the sample.
#' \item GT_BIN: the genotype based on the number of alternate allele in the genotype
#' (the count/dosage of the alternate allele). \code{0, 1, 2, NA}.
#' \item REP_AVG: the reproducibility average, output specific of DArT.
#' }
#' Other columns potentially in the tidy tibble:
#' \enumerate{
#' \item GT: the genotype in 6 digit format \emph{à la genepop}.
#' \item GT_VCF: the genotype in VCF format \code{0/0, 0/1, 1/1, ./.}.
#' \item GT_VCF_NUC: the genotype in VCF format, but keeping the nucleotide information.
#' \code{A/A, A/T, T/T, ./.}
#' \item AVG_COUNT_REF: the coverage for the reference allele, output specific of DArT.
#' \item AVG_COUNT_SNP: the coverage for the alternate allele, output specific of DArT.
#' \item READ_DEPTH: the number of reads used for the genotype (count data).
#' \item ALLELE_REF_DEPTH: the number of reads of the reference allele (count data).
#' \item ALLELE_ALT_DEPTH: the number of reads of the alternate allele (count data).
#' }
#'
#'
#' Written in the working directory:
#' \itemize{
#' \item The radiator GDS file
#' \item The DArT metadata information
#' \item The tidy DArT data
#' \item The strata file associated with this tidy dataset
#' \item The allele dictionary is a tibble with columns:
#' \code{MARKERS, CHROM, LOCUS, POS, REF, ALT}.
#' }

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#' \item \code{whitelist.markers}: detailed in \code{\link[radiator]{filter_whitelist}}.
#' Defautl: \code{whitelist.markers = NULL}.
#' \item \code{missing.memory} (option, path)
#' This argument allows to erase genotypes that have bad statistics.
#' It's the path to a file \code{.rad} file that contains 3 columns:
#' \code{MARKERS, INDIVIDUALS, ERASE}. The file is produced by several radiator
#' functions. For DArT data, \code{\link[radiator]{filter_rad}} generate the file.
#' Defautl: \code{missing.memory = NULL}. Currently not used.
#' \item \code{path.folder}: (optional, path) To write output in a specific folder.
#' Default: \code{path.folder = NULL}. The working directory is used.
#' \item \code{pop.levels}: detailed in \code{\link[radiator]{tidy_genomic_data}}.
#' }


#' @export
#' @rdname read_dart

#' @examples
#' \dontrun{
#' clownfish.dart.tidy <- radiator::read_dart(
#'     data = "clownfish.dart.csv",
#'     strata = "clownfish.strata.tsv"
#'     )
#' }
#' @seealso \code{\link{extract_dart_target_id}}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

read_dart <- function(
  data,
  strata,
  filename = NULL,
  tidy.dart = FALSE,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1,
  ...
) {

  # # for testing
  # filename = NULL
  # verbose = TRUE
  # parallel.core = parallel::detectCores() - 1
  # whitelist.markers = NULL
  # missing.memory <- NULL
  # path.folder = NULL
  # radiator.pipeline = NULL
  # internal <- FALSE
  # tidy.dart = FALSE
  # tidy.check = FALSE
  # gt = NULL
  # gt.bin = NULL
  # gt.vcf = NULL
  # gt.vcf.nuc = NULL
  # pop.levels = NULL

  if (verbose) {
    cat("################################################################################\n")
    cat("############################## radiator::read_dart #############################\n")
    cat("################################################################################\n")
  }

  # Cleanup-------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(message("\nDArT conversion timing: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("############################## completed read_dart #############################\n"), add = TRUE)

  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("whitelist.markers", "missing.memory",
                "path.folder", "internal", "pop.levels",
                "gt", "gt.bin", "gt.vcf", "gt.vcf.nuc",
                "tidy.check"
    ),
    verbose = FALSE
  )
  if (is.null(tidy.check)) tidy.check <- FALSE
  if (tidy.check) tidy.dart <- TRUE
  if (!tidy.dart) {
    gt <- gt.vcf <- gt.vcf.nuc <- FALSE
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data is missing")

  # Folders---------------------------------------------------------------------
  wf <- path.folder <- generate_folder(
    f = path.folder,
    rad.folder = "read_dart",
    prefix_int = FALSE,
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  # radiator.folder <- generate_folder(
  #   f = path.folder,
  #   rad.folder = "import_dart",
  #   prefix_int = TRUE,
  #   internal = FALSE,
  #   file.date = file.date,
  #   verbose = verbose)

  # write the dots file
  write_rad(
    data = rad.dots,
    path = path.folder,
    filename = stringi::stri_join("radiator_tidy_dart_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = internal,
    verbose = verbose
  )

  # Generate filenames ---------------------------------------------------------
  filename.gds <- generate_filename(
    name.shortcut = filename,
    path.folder = path.folder,
    date = TRUE,
    extension = "gds.rad")

  filename <- generate_filename(
    name.shortcut = filename,
    path.folder = path.folder,
    date = TRUE,
    extension = "rad")

  meta.filename <- generate_filename(
    name.shortcut = "radiator_tidy_dart_metadata",
    path.folder = path.folder,
    date = TRUE,
    extension = "rad")

  dic.filename <- generate_filename(
    name.shortcut = "radiator_tidy_dart_allele_dictionary",
    path.folder = path.folder,
    date = TRUE,
    extension = "tsv")

  # Import data ---------------------------------------------------------------
  data <- import_dart(
    data = data,
    strata = strata,
    pop.levels = pop.levels,
    parallel.core = parallel.core,
    verbose = verbose
  )
  strata <- data$strata
  dart.format <- data$dart.format
  data <- data$data

  # Whitelist ------------------------------------------------------------------
  if (!is.null(whitelist.markers)) {
    data %<>% filter_whitelist(data = ., whitelist.markers = whitelist.markers)
  }

  # Silico DArT ----------------------------------------------------------------
  if ("silico.dart" %in% dart.format) {
    data <- data.table::as.data.table(data) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("CLONEID", "SEQUENCE"),
        variable.name = "TARGET_ID",
        variable.factor = FALSE,
        value.name = "VALUE"
      ) %>%
      tibble::as_tibble(.)

    # TODO : I don't think this is necessary with new technique
    strata %<>% dplyr::filter(TARGET_ID %in% data$TARGET_ID)
    if (nrow(strata) == 0) {
      rlang::abort("No more individuals in your data, check data and strata ID names...")
    }

    data %<>% dplyr::filter(TARGET_ID %in% strata$TARGET_ID)
    if (nrow(data) == 0) {
      rlang::abort("No more individuals in your data, check data and strata ID names...")
    }

    suppressWarnings(
      data %<>% dplyr::left_join(strata, by = "TARGET_ID") %>%
        dplyr::select(-TARGET_ID) %>%
        dplyr::mutate(VALUE = as.integer(VALUE))
    )
    strata <- generate_strata(data, pop.id = FALSE)
    n.clone <- length(unique(data$CLONEID))

    filename <- generate_filename(
      name.shortcut = "radiator.silico.dart",
      path.folder = path.folder,
      date = TRUE,
      extension = "rad")
    write_rad(
      data = data,
      path = filename$filename,
      filename = filename$filename.short,
      write.message = "standard",
      verbose = verbose
    )

    if (verbose) cat("################################### SUMMARY ####################################\n")
    message("\nNumber of clones: ", n.clone)
    summary_strata(strata)
    return(data)
  }


  # STRATEGY tidy dart  --------------------------------------------------------
  # Depending on the number of markers ...
  # All this can be overwritten in ... argument
  # gt.bin is the dosage of ALT allele: 0, 1, 2 NA
  ERASE <- NULL
  n.markers <- length(unique(data$MARKERS))
  gt.bin <- TRUE

  # gt.vcf is genotype coding in the VCF: 0/0, 0/1, 1/1, ./.
  if (is.null(gt.vcf)) {
    if (n.markers < 5000) gt.vcf <- TRUE
    if (n.markers >= 5000 && n.markers < 30000) gt.vcf <- TRUE
    if (n.markers >= 30000) gt.vcf <- FALSE
  }
  # gt.vcf.nuc is genotype coding in the VCF but with nucleotides: A/C, ./.
  if (is.null(gt.vcf.nuc)) {
    if (n.markers < 5000) gt.vcf.nuc <- TRUE
    if (n.markers >= 5000 && n.markers < 30000) gt.vcf.nuc <- FALSE
    if (n.markers >= 30000) gt.vcf.nuc <- FALSE
  }
  # gt is genotype coding a la genepop: 001002, 000000
  if (is.null(gt)) {
    if (n.markers < 5000) gt <- TRUE
    if (n.markers >= 5000 && n.markers < 30000) gt <- FALSE
    if (n.markers >= 30000) gt <- FALSE
  }

  # dart2gds -------------------------------------------------------------------
  # check if some rows are missing...
  if (anyNA(unique(data$REF))) {
    if (length(unique(data$MARKERS[!is.na(data$REF)])) != length(unique(data$MARKERS[is.na(data$REF)]))) {
      rlang::abort("The DArT file is missing rows...")
    }
  }

  # Markers meta
  markers.meta <- suppressWarnings(
    dplyr::ungroup(data) %>%
      dplyr::select(-dplyr::one_of(strata$TARGET_ID)) %>%
      dplyr::filter(!is.na(REF) | !is.na(ALT)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::mutate(FILTERS = "whitelist") %>%
      dplyr::select(FILTERS, VARIANT_ID, MARKERS, CHROM, LOCUS, POS, COL, REF, ALT, dplyr::everything(.))
  )

  write_rad(
    data = markers.meta,
    path = meta.filename$filename,
    filename = meta.filename$filename.short,
    write.message = "standard",
    verbose = verbose
  )

  # if (verbose) message("File written: ", meta.filename$)
  want <- c("VARIANT_ID", "MARKERS", "REF", "ALT", strata$INDIVIDUALS)
  suppressWarnings(
    data %<>%
      dplyr::select(dplyr::one_of(want)) %>%
      dplyr::arrange(MARKERS, REF))

  #GDS -----------------------------------------------------------------------
  data <- suppressWarnings(
    dart2gds(
      data = data,
      strata = strata,
      markers.meta = markers.meta,
      filename = filename.gds$filename,
      dart.format = dart.format,
      gt.vcf = gt.vcf,
      gt.vcf.nuc = gt.vcf.nuc,
      gt = gt,
      parallel.core = parallel.core,
      verbose = verbose
    )
  )

  # Tidy data -----------------------------------------------------------------
  n.chrom <- length(unique(markers.meta$CHROM))
  n.locus <- length(unique(markers.meta$LOCUS))
  n.snp <- length(unique(markers.meta$MARKERS))


  # message("Closing the GDS connection")
  # SeqArray::seqClose(object = data.gds)

  # Tidy check ---------------------------------------------------------------
  if (tidy.dart && tidy.check && n.markers > 20000) {
    cat("\n\n################################## IMPORTANT ###################################\n")
    message("Tidying DArT data with ", n.markers, " SNPs is not optimal:")
    message("    1. a computer with lots of RAM is required")
    message("    2. it's very slow to generate")
    message("    3. it's very slow to run codes after")
    message("    4. for most non model species this number of markers is not realistic...")
    message("\nRecommendation:")
    message("    1. stop here and just use the GDS")
    message("    2. filter your dataset. e.g. with filter_rad")
    message("\nIdeally target a maximum of ~ 10 000 - 20 0000 unlinked SNPs\n")

    if (n.markers > 20000) tidy.dart <- FALSE
    tidy.dart <- radiator_question(
      x = "\nContinue tidying the DArT data (y/n) ?",
      answer.opt = c("Y", "N", "Yes", "No", "YES", "NO", "yes", "no", "y", "n"))
    if (any(c("y", "Y", "Yes", "YES", "yes") %in% tidy.dart)) {
      tidy.dart <- TRUE
      message("Tidying the large DArT file...")
    } else {
      message("\nKeeping only the GDS object/file")
      tidy.dart <- FALSE
    }
  }

  # Print genotypes tidying
  if (verbose) {
    message("\nGenotypes formats generated with ", n.markers, " SNPs: ")
    message("    GT_BIN (the dosage of ALT allele: 0, 1, 2 NA): ", gt.bin)
    message("    GT_VCF (the genotype coding VCFs: 0/0, 0/1, 1/1, ./.): ", gt.vcf)
    message("    GT_VCF_NUC (the genotype coding in VCFs, but with nucleotides: A/C, ./.): ", gt.vcf.nuc)
    message("    GT (the genotype coding 'a la genepop': 001002, 001001, 000000): ", gt)
  }

  if (tidy.dart) {
    notwanted <- c("REF", "ALT")
    tidy.data <- markers.meta %>%
      dplyr::left_join(
        extract_genotypes_metadata(gds = data, whitelist = TRUE) %>%
          dplyr::select(-dplyr::one_of(notwanted))
        , by = "MARKERS"
      )

    # Check that merging was successful
    if (!identical(tidy.data$VARIANT_ID.x, tidy.data$VARIANT_ID.y)) {
      rlang::abort("Tidying problem: raise an issue on GitHub...")
    } else {
      tidy.data %<>%
        dplyr::select(-VARIANT_ID.y) %>%
        dplyr::rename(VARIANT_ID = VARIANT_ID.x)
    }

    # Calibration of ref/alt alleles ------------------------------------------
    # Now done during coding of genotypes

    # Final strata ---------------------------------------------------------
    strata.filename <- generate_filename(
      name.shortcut = "radiator_tidy_dart_strata",
      path.folder = path.folder,
      date = TRUE,
      extension = "tsv")
    strata <- extract_individuals_metadata(gds = data, whitelist = TRUE)
    readr::write_tsv(x = strata, path = strata.filename$filename)

    if (!is.null(strata)) {
      if (rlang::has_name(tidy.data, "TARGET_ID")) {
        tidy.data %<>%
          dplyr::left_join(strata, by = "TARGET_ID") %>%
          dplyr::select(-TARGET_ID) %>%
          dplyr::rename(POP_ID = STRATA)
      } else {
        tidy.data %<>%
          dplyr::left_join(strata, by = "INDIVIDUALS") %>%
          dplyr::rename(POP_ID = STRATA)
      }
    } else {
      tidy.data %<>% dplyr::mutate(POP_ID = 1L)
    }


    # Erase genotypes
    if (!is.null(missing.memory)) {
      message("Using missing.memory file to erase genotypes")
      missing <- radiator::read_rad(data = missing.memory) %>%
        dplyr::arrange(MARKERS, INDIVIDUALS)

      data %<>% dplyr::arrange(MARKERS, INDIVIDUALS)

      #check identical markers
      same.markers <- identical(unique(missing$MARKERS), unique(data$MARKERS))
      same.individuals<- identical(unique(missing$INDIVIDUALS), unique(data$INDIVIDUALS))
      if (!same.markers || !same.individuals) {
        message("note: data and missing memory don't share all the same markers and/or individuals")
        data %<>% dplyr::left_join(missing, by = c("MARKERS", "INDIVIDUALS"))
        #%>% dplyr::mutate(ERASE = replace(ERASE, which(is.na(ERASE)), FALSE))
        data$ERASE[is.na(data$ERASE)] <- FALSE # faster
        which.missing <- which(data$ERASE)
        data %<>% dplyr::select(-ERASE)
      } else {
        which.missing <- which(missing$ERASE)
      }

      message("Erasing genotypes and genotypes metadata...")
      if (rlang::has_name(data, "GT_BIN")) data$GT_BIN[which.missing] <- NA
      if (rlang::has_name(data, "GT")) data$GT[which.missing] <- "000000"
      if (rlang::has_name(data, "GT_VCF")) data$GT_VCF[which.missing] <- "./."
      if (rlang::has_name(data, "GT_VCF_NUC")) data$GT_VCF_NUC[which.missing] <- "./."
      if (rlang::has_name(data, "READ_DEPTH")) data$READ_DEPTH[which.missing] <- NA
      if (rlang::has_name(data, "ALLELE_REF_DEPTH")) data$ALLELE_REF_DEPTH[which.missing] <- NA
      if (rlang::has_name(data, "ALLELE_ALT_DEPTH")) data$ALLELE_ALT_DEPTH[which.missing] <- NA
    }#End missing.memory


    # write tidy
    write_rad(data = tidy.data, path = filename$filename)

    # Whitelist
    # readr::write_tsv(x = markers.meta, file.path(wf, "whitelist.markers.tidy.dart.tsv"))
    write_rad(data = markers.meta,
              path = wf,
              filename = "whitelist.markers.tidy.dart.tsv",
              tsv = TRUE,
              write.message = "standard",
              verbose = verbose)
    data <- tidy.data
  }#tidy

  # Results --------------------------------------------------------------------
  if (verbose) cat("################################### SUMMARY ####################################\n")
  message("\nNumber of chrom: ", n.chrom)
  message("Number of locus: ", n.locus)
  message("Number of SNPs: ", n.snp)
  summary_strata(strata)
  return(data)

}#End read_dart

# INTERNAL FUNCTIONS------------------------------------------------------------
#' @title import_dart
#' @description Read DArT file
#' @rdname import_dart
#' @keywords internal
#' @export
import_dart <- function(
  data,
  strata,
  pop.levels = NULL,
  path.folder = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {

  # pop.levels = NULL
  # path.folder = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE

  message("Reading DArT file...")

  # Check that DArT file as good target id written -----------------------------
  target.id <- extract_dart_target_id(data, write = FALSE)
  n.ind.dart <- nrow(target.id)
  if (verbose) message("    Number of individuals: ", n.ind.dart)

  # Check DArT format file -----------------------------------------------------
  dart.check <- check_dart(data)
  silico.dart <- "silico.dart" %in% dart.check$data.type

  # Strata file ------------------------------------------------------------------
  strata.df <- radiator::read_strata(
    strata = strata,
    pop.levels = pop.levels, pop.labels = NULL,
    pop.select = NULL, blacklist.id = NULL,
    keep.two = FALSE, verbose = verbose)
  pop.levels <- strata.df$pop.levels
  pop.labels <- strata.df$pop.labels
  pop.select <- strata.df$pop.select
  blacklist.id <- strata.df$blacklist.id
  strata.df <- strata.df$strata
  n.ind.strata <- nrow(strata.df)

  # Check that TARGET_ID in strata match TARGET_ID in the DArT file ------------
  if (verbose) message("Using individuals in strata file to filter individuals in DArT file")
  if (n.ind.dart != n.ind.strata) {
    strata.id.check <- strata.df %>%
      dplyr::mutate(IN_DART = stringi::stri_trans_toupper(strata.df$TARGET_ID)
                    %in% stringi::stri_trans_toupper(target.id$TARGET_ID))
    strata.id.pass <- !FALSE %in% (unique(strata.id.check$IN_DART))
    if (!strata.id.pass) {
      problem.filename <- generate_filename(
        name.shortcut = "radiator_tidy_dart_strata_problem",
        path.folder = path.folder,
        date = TRUE,
        extension = "tsv")

      readr::write_tsv(
        x = strata.id.check,
        path = problem.filename$filename)
      rlang::abort("\nSome of the samples in the strata are not found in the DArT file.
                     For more info: ", problem.filename$filename.short)
    }
    if (verbose) message("\nNote: Careful if using DArT statistics generated for all samples...")
    strata.id.check <- NULL
  } else {
    if (!identical(sort(target.id$TARGET_ID), sort(strata.df$TARGET_ID))) {
      rlang::abort("\nThe DArT and strata files don't have the same TARGET_IDs")
    }
  }
  blacklist.id <- c(
    dplyr::filter(strata.df, !TARGET_ID %in% target.id$TARGET_ID) %$% TARGET_ID,
    dplyr::filter(target.id, !TARGET_ID %in% strata.df$TARGET_ID) %$% TARGET_ID
  ) %>% unique
  if (length(blacklist.id) == 0) blacklist.id <- NULL

  message("Number of blacklisted samples: ", length(blacklist.id))
  target.id <- NULL

  # need to check for duplicate names... yes happening all the time
  duplicate.id.strata <- length(strata.df$INDIVIDUALS) - length(unique(strata.df$INDIVIDUALS))

  if (duplicate.id.strata > 0) {
    message("Duplicated individuals names found in the strata.\n   number of duplicate names = ", duplicate.id.strata, "\n")
    rlang::abort("\nFix the strata with unique names and\nverify the DArT file for the same issue, adjust accordingly...")
  }

  # use.readr <- FALSE # fread is now a lot faster...switching to fread...
  # if (use.readr) {
  #   dart.col.type <- readr::read_delim(
  #     file = data,
  #     delim = dart.check$tokenizer.dart,
  #     skip = dart.check$skip.number,
  #     n_max = 1,
  #     na = "-",
  #     col_names = FALSE,
  #     col_types = readr::cols(.default = readr::col_character()))
  #
  #   if (silico.dart) {
  #     info <- c("CLONEID", "ALLELESEQUENCE", "SEQUENCE", "TRIMMEDSEQUENCE")
  #     info.type <- c("c", "c", "c", "c")
  #   } else {
  #     info <- c("ALLELEID", "SNP", "SNPPOSITION", "CALLRATE",
  #               "AVGCOUNTREF", "AVGCOUNTSNP", "REPAVG", "CLONEID", "AVGREADDEPTH",
  #               "REPRODUCIBILITY", "CLUSTERCONSENSUSSEQUENCE", "ONERATIOREF", "ONERATIOSNP")
  #     info.type <- c("c", "c", "i", "d", "d", "d", "d", "c", "d", "d", "c", "d", "d")
  #   }
  #
  #   want <- tibble::tibble(INFO = info, COL_TYPE = info.type) %>%
  #     dplyr::bind_rows(
  #       dplyr::select(strata.df, INFO = TARGET_ID) %>%
  #         dplyr::mutate(
  #           COL_TYPE = rep("c", n()),
  #           INFO = stringi::stri_trans_toupper(INFO)))
  #
  #   dart.col.type %<>%
  #     t %>%
  #     magrittr::set_colnames(x = ., value = "INFO") %>%
  #     tibble::as_tibble(.) %>%
  #     dplyr::mutate(
  #       INFO = stringi::stri_trans_toupper(INFO),
  #       INFO = clean_ind_names(INFO)
  #     ) %>%
  #     dplyr::left_join(want, by = "INFO") %>%
  #     dplyr::mutate(COL_TYPE = stringi::stri_replace_na(str = COL_TYPE, replacement = "_")) %>%
  #     dplyr::select(COL_TYPE) %>%
  #     purrr::flatten_chr(.) %>%
  #     stringi::stri_join(collapse = "")
  #   want <- NULL
  #
  #   data <- readr::read_delim(
  #     file = data,
  #     delim = dart.check$tokenizer.dart,
  #     col_names = TRUE,
  #     col_types = dart.col.type,
  #     na = c("-", " ", "", "NA"),
  #     skip = dart.check$skip.number
  #   )
  #   dart.col.type <- NULL
  # } else {#fread
    # There's no big difference really here if we import everything and filter after...
    data <- suppressWarnings(
      data.table::fread(
        file = data,
        header = TRUE,
        strip.white = TRUE,
        na.strings = c("-", "NA"),
        stringsAsFactors = FALSE,
        skip = "CallRate",
        drop = blacklist.id,
        select = NULL,
        showProgress = TRUE,
        nThread = parallel.core,
        verbose = FALSE) %>%
        tibble::as_tibble(.) %>%
        # We want snakecase not camelcase
        # Change the TARGET_ID by INDIVIDUALS...
        clean_dart_colnames(data = ., strata = strata.df)
    )

  # }
  # check <- tibble::tibble(BLACKLIST_ID = colnames(data)) %>%
  #   dplyr::filter(BLACKLIST_ID %in% blacklist.id)

  # Check sequence
  # check.seq <- data %>%
  #   dplyr::select(ALLELE_SEQUENCE, TRIMMED_SEQUENCE) %>%
  #   dplyr::filter(ALLELE_SEQUENCE != TRIMMED_SEQUENCE)

  if (silico.dart) {
    colnames(data) <- stringi::stri_replace_all_fixed(
      str = colnames(data),
      pattern = c("ALLELE_SEQUENCE", "TRIMMED_SEQUENCE"),
      replacement = c("SEQUENCE", "SEQUENCE"),
      vectorize_all = FALSE)
    # }
  } else {# non silico dart data
    colnames(data) %<>%
      stringi::stri_replace_all_fixed(
        str = .,
        pattern = c("ALLELE_ID","SNP_POSITION"),
        replacement = c("LOCUS", "POS"),
        vectorize_all = FALSE)


    # keep consensus sequence if found
    # or rename TRIMMED_SEQUENCE
    if (rlang::has_name(data, "CLUSTER_CONSENSUS_SEQUENCE")) {
      data %<>% dplyr::rename(SEQUENCE = CLUSTER_CONSENSUS_SEQUENCE)
    } else if (rlang::has_name(data, "TRIMMED_SEQUENCE")) {
      data %<>% dplyr::rename(SEQUENCE = TRIMMED_SEQUENCE)
    } else {
      if (rlang::has_name(data, "ALLELE_SEQUENCE")) {
        data %<>% dplyr::rename(SEQUENCE = ALLELE_SEQUENCE)
      }
    }


    if (rlang::has_name(data, "CLONE_ID")) {
      if (rlang::has_name(data, "LOCUS")) {
        data %<>% dplyr::select(-CLONE_ID)
      } else {
        colnames(data) %<>%
          stringi::stri_replace_all_fixed(
            str = .,
            pattern = "CLONE_ID",
            replacement = "LOCUS",
            vectorize_all = FALSE)
      }
    }

    if (!rlang::has_name(data, "LOCUS")) {
      rlang::abort("\nProblem tidying DArT dataset: contact author")
    }

    # necessary steps...observed with DArT file using ref genome ---------------
    data %<>% dplyr::filter(!is.na(LOCUS))
    # Check for duplicate rows (sometimes people combine DArT data...)
    if (rlang::has_name(data, "POS")) {
      data %<>% dplyr::arrange(LOCUS, POS)
      data.dup <- nrow(dplyr::distinct(data, LOCUS, SNP, POS, CALL_RATE, .keep_all = FALSE))
    } else {
      data %<>% dplyr::arrange(LOCUS)
      data.dup <- nrow(dplyr::distinct(data, LOCUS, SNP, CALL_RATE, .keep_all = FALSE))
    }

    # make sure no duplicates
    if (nrow(data) != data.dup) {
      message("Duplicate rows were identified")
      message("    using distinct rows")
      message("    check data if downstream problems")
      data %<>% dplyr::distinct(LOCUS, SNP, POS, CALL_RATE, .keep_all = TRUE)
    }
    data.dup <- NULL

    # Screen for duplicate names -------------------------------------------------
    id <- purrr::discard(
      .x = colnames(data),
      .p = colnames(data) %in% c("LOCUS", "SNP", "POS", "CALL_RATE", "AVG_COUNT_REF",
                                 "AVG_COUNT_SNP", "REP_AVG")
    )
    dup.id <- length(id) - length(unique(id))
    if (dup.id > 0) {
      rlang::abort(stringi::stri_join("Duplicated individual names in the data: ", dup.id))
    }
    # removing unused object
    remove.list <- id <- dup.id <- NULL

    # clean locus, generate MARKERS and VARIANT_ID
    data %<>% clean_dart_locus(.)
  }

  # DArT characteristics--------------------------------------------------------
  if (verbose)  message("\nDArT characteristics:")
  # dart.format:
  # "1row" !binary
  # "2rows" binary
  # "counts" binary
  # "silico.dart"
  dart.format <- detect_dart_format(
    x = data,
    target.id = strata.df$TARGET_ID,
    verbose = TRUE)

  return(res = list(data = data, strata = strata.df, dart.format = dart.format))
}#End import_dart

# clean_dart_locus--------------------------------------------------------------
#' @title clean_dart_locus
#' @description Clean LOCUS and generate VARIANT_ID and MARKERS
#' @rdname clean_dart_locus
#' @keywords internal
#' @export
clean_dart_locus <- function(x, fast = TRUE) {

  if (fast) {
    # x <- data
    if (!rlang::has_name(x, "CHROM")) x %<>% dplyr::mutate(CHROM = "CHROM_1")
    want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL",
              "REF", "ALT", "SEQUENCE",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
    x %<>%
      dplyr::mutate(
        COL = stringi::stri_extract_first_regex(
          str = LOCUS,
          pattern = "[-][0-9]+[\\:]"),
        COL = stringi::stri_replace_all_fixed(
          str = COL,
          pattern = c("-", ":"),
          replacement = c("", ""),
          vectorize_all = FALSE),
        COL = as.integer(COL),
        LOCUS = stringi::stri_extract_first_regex(str = LOCUS, pattern = "^[0-9]+"),
        REF = stringi::stri_extract_first_regex(str = SNP, pattern = "[A-Z]"),
        ALT = stringi::stri_extract_last_regex(str = SNP, pattern = "[A-Z]"),
        MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
        VARIANT_ID = as.integer(factor(MARKERS)),
        SNP = NULL
      ) %>%
      dplyr::select(dplyr::one_of(want), dplyr::everything()) %>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("MARKERS", "CHROM", "LOCUS", "POS"),
        .funs = as.character
      ) %>%
      dplyr::arrange(CHROM, LOCUS, POS, REF)

  } else {
    want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL",
              "REF", "ALT", "SEQUENCE",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
    suppressWarnings(
      x %<>%
        tidyr::separate(col = LOCUS,
                        into = c("LOCUS", "NOT_USEFUL"),
                        sep = "\\|",
                        extra = "drop"
        ) %>%
        dplyr::select(-NOT_USEFUL) %>%
        tidyr::separate(col = SNP,
                        into = c("NOT_USEFUL", "KEEPER"),
                        sep = ":",
                        extra = "drop") %>%
        dplyr::select(-NOT_USEFUL) %>%
        tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">") %>%
        dplyr::mutate(
          CHROM = rep("CHROM_1", n()),
          MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
          VARIANT_ID = as.integer(factor(MARKERS)),
          COL = POS
        ) %>%
        dplyr::select(dplyr::one_of(want), dplyr::everything()) %>%
        dplyr::mutate_at(
          .tbl = .,
          .vars = c("MARKERS", "CHROM", "LOCUS", "POS"),
          .funs = as.character
        ) %>%
        dplyr::arrange(CHROM, LOCUS, POS, REF)
    )
  }
}#End clean_dart_locus

# detect_dart_format-------------------------------------------------------------
#' @title detect_dart_format
#' @description Detect the dart genotype format: 1row, 2rows or counts
#' @rdname detect_dart_format
#' @keywords internal
#' @export
detect_dart_format <- function(x = NULL, target.id = NULL, verbose = TRUE) {
  # dart.format:
  # silico.dart
  # 1row
  # 2rows
  # counts

  if (rlang::has_name(x, "CLONEID") && rlang::has_name(x, "SEQUENCE")) {
    if (verbose) message("DArT SNP format: silico DArT")
    dart.format <- "silico.dart"
  } else {

    binary <- anyNA(x$REF)

    if (!binary) {
      if (verbose) message("DArT SNP format: genotypes in 1 Row")
      dart.format <- "1row"
    } else {
      count.data <- x %>%
        dplyr::select(
          dplyr::one_of(
            sample(x = target.id, size = min(10, floor(0.1 * length(target.id))))
          )
        ) %>%
        dplyr::mutate_all(.tbl = ., .funs = as.numeric) %>%
        purrr::flatten_dbl(.) %>%
        unique(.)

      count.data <- any(count.data > 1, na.rm = TRUE)

      if (count.data) {
        dart.format <- "counts"
        if (verbose) message("DArT SNP format: alleles coverage in 2 Rows counts")
      } else {
        dart.format <- "2rows"
        if (verbose) message("DArT SNP format: alleles absence/presence in 2 Rows")
      }
    }
  }
  return(dart.format)
}#End detect_dart_format

# dart2gds----------------------------------------------------------------
#' @title dart2gds
#' @description Transform dart to GDS format
#' @rdname dart2gds
#' @keywords internal
#' @export

dart2gds <- function(
  data,
  strata = NULL,
  markers.meta,
  filename,
  dart.format,
  gt.vcf = NULL,
  gt.vcf.nuc = NULL,
  gt = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {
  # gt.vcf <- gt.vcf.nuc <- gt <- NULL

  message("Generating genotypes...")
  if (dart.format == "1row") data.source <- c("dart", "1row")
  if (dart.format == "2rows") data.source <- c("dart", "2rows")
  if (dart.format == "counts") data.source <- c("dart", "counts")
  parallel.core.temp <- max(1, 2 * floor(parallel.core / 2))
  geno.coding <- "alt.dos"
  dp <- NULL

  genotypes.meta <- dplyr::left_join(
    data,
    dplyr::distinct(data, MARKERS) %>%
      dplyr::mutate(
        SPLIT_VEC = split_vec_row(., 4, parallel.core = parallel.core.temp)
      )
    , by = "MARKERS") %>%
    dplyr::group_split(SPLIT_VEC, keep = FALSE) %>%
    .radiator_parallel_mc(
      X = .,
      FUN = generate_geno_v1,
      mc.cores = parallel.core.temp,
      data.source = data.source,
      gt.vcf = gt.vcf,
      gt.vcf.nuc = gt.vcf.nuc,
      gt = gt
    ) %>%
    dplyr::bind_rows(.) %>%
    tibble::add_column(.data = ., FILTERS = "whitelist", .before = 1) %>%
    dplyr::arrange(MARKERS, INDIVIDUALS)

  data <- genotypes.meta %>%
    dplyr::select(VARIANT_ID, INDIVIDUALS, GT_BIN) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = VARIANT_ID ~ INDIVIDUALS,
      value.var = "GT_BIN"
    ) %>%
    tibble::as_tibble(.) %>%
    dplyr::arrange(VARIANT_ID) %>%
    tibble::column_to_rownames(.data = ., var = "VARIANT_ID")

  # Generate GDS ---------------------------------------------------------------
  data <- radiator_gds(
    genotypes.df = data,
    geno.coding = geno.coding,
    strata = strata,
    biallelic = TRUE,
    markers.meta = markers.meta,
    genotypes.meta = genotypes.meta,
    dp = dp,
    filename = filename,
    data.source = data.source,
    open = TRUE,
    verbose = verbose)
  if (verbose) message("done!")
  return(data)
}# End dart2gds

# generate_geno----------------------------------------------------------------
#' @title generate_geno_v1------------------------------------------------------
#' @description Generate the genotypes and prep for GDS
#' @rdname generate_geno_v1
#' @keywords internal
#' @export
generate_geno_v1 <- function(
  x,
  data.source,
  gt = FALSE,
  gt.vcf.nuc = FALSE,
  gt.vcf = FALSE
) {
  message("Generating genotypes and calibrating REF/ALT alleles...")
  # required func:
  switch_allele_count <- function(x, data.source) {
    if ("1row" %in% data.source) {
      x <- dplyr::case_when(
        x == 1 ~ 2,
        x == 2 ~ 1,
        x == 0 ~ 0
      )
    }
    if ("2rows" %in% data.source) {
      # here we want count of alternate allele instead...
      # x <- as.integer(dplyr::recode(.x = as.character(x), "0" = "1", "1" = "0"))
      # case_when is much faster than recode...
      x <- dplyr::case_when(
        x == 0 ~ 1,
        x == 1 ~ 0
      )
    }
    return(x)
  }# End switch_allele_count

  #1-row
  if ("1row" %in% data.source) {
    # x <- genotypes.meta[[1]]
    # x <- data
    res <- x %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("VARIANT_ID", "MARKERS", "REF", "ALT"),
        variable.name = "INDIVIDUALS",
        value.name = "GT_BIN",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(GT_BIN = as.integer(GT_BIN))
    x <- NULL


    res %<>% dplyr::mutate(GT_BIN = switch_allele_count(GT_BIN, data.source = data.source))

    # Counts for allele calibration
    switch <- dplyr::select(res, MARKERS, GT_BIN) %>%
      dplyr::filter(!is.na(GT_BIN)) %>%
      dplyr::count(GT_BIN, MARKERS) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        REF_COUNT = sum((2 * n[GT_BIN == 0]), n[GT_BIN == 1], na.rm = TRUE),
        ALT_COUNT = sum((2 * n[GT_BIN == 2]), n[GT_BIN == 1], na.rm = TRUE)
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(dplyr::if_else(REF_COUNT < ALT_COUNT, TRUE, FALSE)) %$%
      MARKERS

    n.switch <- length(switch)
    if (n.switch > 0) {
      message("Calibration REF/ALT based on counts of alleles: ", n.switch)
      res <- dplyr::filter(res, !MARKERS %in% switch) %>%
        dplyr::bind_rows(
          dplyr::filter(res, MARKERS %in% switch) %>%
            dplyr::rename(ALT = REF, REF = ALT)
        )
    }
    switch <- NULL
  }#1row genotypes
  #2-rows
  if ("2rows" %in% data.source) {
    # x <- genotypes.meta[[1]]
    res <- dplyr::filter(x, !is.na(REF)) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("VARIANT_ID", "MARKERS", "REF", "ALT"),
        variable.name = "INDIVIDUALS",
        value.name = "A2",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::bind_cols(
        dplyr::filter(x, is.na(REF)) %>%
          dplyr::select(-REF, -ALT, -VARIANT_ID) %>%
          data.table::as.data.table(.) %>%
          data.table::melt.data.table(
            data = .,
            id.vars = "MARKERS",
            variable.name = "INDIVIDUALS",
            value.name = "A1",
            variable.factor = FALSE) %>%
          tibble::as_tibble(.)
      ) %>%
      dplyr::mutate_at(.tbl = ., .vars = c("A1", "A2"), .funs = as.integer)
    x <- NULL

    if (!identical(res$MARKERS, res$MARKERS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-MARKERS1)
    }

    if (!identical(res$INDIVIDUALS, res$INDIVIDUALS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-INDIVIDUALS1)
    }

    res %<>% dplyr::mutate(GT_BIN = switch_allele_count(A1, data.source) + A2, A1 = NULL, A2 = NULL)

    # Counts for allele calibration
    switch <- dplyr::select(res, MARKERS, GT_BIN) %>%
      dplyr::filter(!is.na(GT_BIN)) %>%
      dplyr::count(GT_BIN, MARKERS) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        REF_COUNT = sum((2 * n[GT_BIN == 0]), n[GT_BIN == 1], na.rm = TRUE),
        ALT_COUNT = sum((2 * n[GT_BIN == 2]), n[GT_BIN == 1], na.rm = TRUE)
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(dplyr::if_else(REF_COUNT < ALT_COUNT, TRUE, FALSE)) %$%
      MARKERS

    n.switch <- length(switch)
    if (n.switch > 0) {
      message("Calibration REF/ALT based on counts of alleles: ", n.switch)
      res <- dplyr::filter(res, !MARKERS %in% switch) %>%
        dplyr::bind_rows(
          dplyr::filter(res, MARKERS %in% switch) %>%
            dplyr::rename(ALT = REF, REF = ALT)
        )
    }
    switch <- NULL
  }#2rows genotypes

  # counts
  # x <- genotypes.meta[[1]]
  if ("counts" %in% data.source) {
    res <- dplyr::filter(x, !is.na(REF)) %>%
      dplyr::arrange(MARKERS) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("VARIANT_ID", "MARKERS", "REF", "ALT"),
        variable.name = "INDIVIDUALS",
        value.name = "ALLELE_ALT_DEPTH",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::bind_cols(
        dplyr::filter(x, is.na(REF)) %>%
          dplyr::select(-REF, -ALT, -VARIANT_ID) %>%
          data.table::as.data.table(.) %>%
          data.table::melt.data.table(
            data = .,
            id.vars = "MARKERS",
            variable.name = "INDIVIDUALS",
            value.name = "ALLELE_REF_DEPTH",
            variable.factor = FALSE) %>%
          tibble::as_tibble(.)
      )
    x <- NULL

    if (!identical(res$MARKERS, res$MARKERS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-MARKERS1)
    }

    if (!identical(res$INDIVIDUALS, res$INDIVIDUALS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-INDIVIDUALS1)
    }

    res %<>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
        .funs = as.numeric
      ) %>%
      dplyr::mutate(READ_DEPTH = ALLELE_REF_DEPTH + ALLELE_ALT_DEPTH)

    # Coverage for allele calibration
    want <- c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH")
    switch <- dplyr::group_by(res, MARKERS) %>%
      dplyr::summarise_at(.tbl = ., .vars = want, .funs = sum, na.rm = TRUE) %>%
      dplyr::mutate_at(.tbl = ., .vars = want, .funs = round, digits = 0) %>%
      dplyr::mutate_at(.tbl = ., .vars = want, .funs = as.integer) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(dplyr::if_else(ALLELE_ALT_DEPTH > ALLELE_REF_DEPTH, TRUE, FALSE)) %$%
      MARKERS
    n.switch <- length(switch)
    if (n.switch > 0) {
      message("Calibration REF/ALT based on read depth of alleles: ", n.switch)
      res <- dplyr::filter(res, !MARKERS %in% switch) %>%
        dplyr::bind_rows(
          dplyr::filter(res, MARKERS %in% switch) %>%
            dplyr::rename(
              ALT = REF,
              REF = ALT,
              ALLELE_REF_DEPTH = ALLELE_ALT_DEPTH,
              ALLELE_ALT_DEPTH = ALLELE_REF_DEPTH
            )
        )
    }
    switch <- NULL

    # GT_BIN
    res %<>%
      dplyr::mutate(
        GT_BIN = dplyr::case_when(
          ALLELE_REF_DEPTH > 0 & ALLELE_ALT_DEPTH == 0 ~ 0,
          ALLELE_REF_DEPTH > 0 & ALLELE_ALT_DEPTH > 0 ~ 1,
          ALLELE_REF_DEPTH == 0 & ALLELE_ALT_DEPTH > 0 ~ 2
        )
      ) %>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
        .funs = replace_by_na, what = 0
      )
  }# Counts
  if (gt.vcf) {
    res %<>%
      dplyr::mutate(
        GT_VCF = dplyr::case_when(
          GT_BIN == 0 ~ "0/0", GT_BIN == 1 ~ "0/1", GT_BIN == 2 ~ "1/1",
          is.na(GT_BIN) ~ "./.")
      )
  }
  if (gt.vcf.nuc) {
    res %<>%
      dplyr::mutate(
        GT_VCF_NUC = dplyr::case_when(
          GT_BIN == 0 ~ stringi::stri_join(REF, REF, sep = "/"),
          GT_BIN == 1 ~ stringi::stri_join(REF, ALT, sep = "/"),
          GT_BIN == 2 ~ stringi::stri_join(ALT, ALT, sep = "/"),
          is.na(GT_BIN) ~ "./.")
      )
  }
  if (gt) {
    res %<>%
      dplyr::mutate(
        GT = stringi::stri_replace_all_fixed(
          str = GT_VCF_NUC,
          pattern = c("A", "C", "G", "T", "/", ".."),
          replacement = c("001", "002", "003", "004", "", "000000"),
          vectorize_all = FALSE)
      )
  }

  return(res)
}# End generate_geno_v1

#' @title generate_geno
#' @description Generate the genotypes and prep for GDS
#' @rdname generate_geno_v2
#' @keywords internal
#' @export
generate_geno_v2 <- function(
  x,
  data.source,
  gt = FALSE,
  gt.vcf.nuc = FALSE,
  gt.vcf = FALSE
) {
  message("Generating genotypes and calibrating REF/ALT alleles...")
  # required func:
  switch_allele_count <- function(x, data.source) {
    if ("1row" %in% data.source) {
      x <- dplyr::case_when(
        x == 1 ~ 2,
        x == 2 ~ 1,
        x == 0 ~ 0
      )
    }
    if ("2rows" %in% data.source) {
      # here we want count of alternate allele instead...
      # x <- as.integer(dplyr::recode(.x = as.character(x), "0" = "1", "1" = "0"))
      # case_when is much faster than recode...
      x <- dplyr::case_when(
        x == 0 ~ 1,
        x == 1 ~ 0
      )
    }
    return(x)
  }# End switch_allele_count

  #1-row
  if ("1row" %in% data.source) {
    # x <- genotypes.meta[[1]]
    # x <- data
    res <- x %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("VARIANT_ID", "MARKERS", "REF", "ALT"),
        variable.name = "INDIVIDUALS",
        value.name = "GT_BIN",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(GT_BIN = as.integer(GT_BIN))
    x <- NULL


    res %<>% dplyr::mutate(GT_BIN = switch_allele_count(GT_BIN, data.source = data.source))

    # Counts for allele calibration
    switch <- dplyr::select(res, MARKERS, GT_BIN) %>%
      dplyr::filter(!is.na(GT_BIN)) %>%
      dplyr::count(GT_BIN, MARKERS) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        REF_COUNT = sum((2 * n[GT_BIN == 0]), n[GT_BIN == 1], na.rm = TRUE),
        ALT_COUNT = sum((2 * n[GT_BIN == 2]), n[GT_BIN == 1], na.rm = TRUE)
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(dplyr::if_else(REF_COUNT < ALT_COUNT, TRUE, FALSE)) %$%
      MARKERS

    n.switch <- length(switch)
    if (n.switch > 0) {
      message("Calibration REF/ALT based on counts of alleles: ", n.switch)
      res <- dplyr::filter(res, !MARKERS %in% switch) %>%
        dplyr::bind_rows(
          dplyr::filter(res, MARKERS %in% switch) %>%
            dplyr::rename(ALT = REF, REF = ALT)
        )
    }
    switch <- NULL
  }#1row genotypes
  #2-rows
  if ("2rows" %in% data.source) {
    # x <- genotypes.meta[[1]]
    res <- dplyr::filter(x, !is.na(REF)) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("VARIANT_ID", "MARKERS", "REF", "ALT"),
        variable.name = "INDIVIDUALS",
        value.name = "A2",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::bind_cols(
        dplyr::filter(x, is.na(REF)) %>%
          dplyr::select(-REF, -ALT, -VARIANT_ID) %>%
          data.table::as.data.table(.) %>%
          data.table::melt.data.table(
            data = .,
            id.vars = "MARKERS",
            variable.name = "INDIVIDUALS",
            value.name = "A1",
            variable.factor = FALSE) %>%
          tibble::as_tibble(.)
      ) %>%
      dplyr::mutate_at(.tbl = ., .vars = c("A1", "A2"), .funs = as.integer)
    x <- NULL

    if (!identical(res$MARKERS, res$MARKERS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-MARKERS1)
    }

    if (!identical(res$INDIVIDUALS, res$INDIVIDUALS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-INDIVIDUALS1)
    }

    res %<>% dplyr::mutate(GT_BIN = switch_allele_count(A1, data.source) + A2, A1 = NULL, A2 = NULL)

    # Counts for allele calibration
    switch <- dplyr::select(res, MARKERS, GT_BIN) %>%
      dplyr::filter(!is.na(GT_BIN)) %>%
      dplyr::count(GT_BIN, MARKERS) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        REF_COUNT = sum((2 * n[GT_BIN == 0]), n[GT_BIN == 1], na.rm = TRUE),
        ALT_COUNT = sum((2 * n[GT_BIN == 2]), n[GT_BIN == 1], na.rm = TRUE)
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(dplyr::if_else(REF_COUNT < ALT_COUNT, TRUE, FALSE)) %$%
      MARKERS

    n.switch <- length(switch)
    if (n.switch > 0) {
      message("Calibration REF/ALT based on counts of alleles: ", n.switch)
      res <- dplyr::filter(res, !MARKERS %in% switch) %>%
        dplyr::bind_rows(
          dplyr::filter(res, MARKERS %in% switch) %>%
            dplyr::rename(ALT = REF, REF = ALT)
        )
    }
    switch <- NULL
  }#2rows genotypes

  # counts
  # x <- genotypes.meta[[1]]
  if ("counts" %in% data.source) {
    res <- dplyr::filter(x, !is.na(REF)) %>%
      dplyr::arrange(MARKERS) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("VARIANT_ID", "MARKERS", "REF", "ALT"),
        variable.name = "INDIVIDUALS",
        value.name = "ALLELE_ALT_DEPTH",
        variable.factor = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::bind_cols(
        dplyr::filter(x, is.na(REF)) %>%
          dplyr::select(-REF, -ALT, -VARIANT_ID) %>%
          data.table::as.data.table(.) %>%
          data.table::melt.data.table(
            data = .,
            id.vars = "MARKERS",
            variable.name = "INDIVIDUALS",
            value.name = "ALLELE_REF_DEPTH",
            variable.factor = FALSE) %>%
          tibble::as_tibble(.)
      )
    x <- NULL

    if (!identical(res$MARKERS, res$MARKERS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-MARKERS1)
    }

    if (!identical(res$INDIVIDUALS, res$INDIVIDUALS1)) {
      rlang::abort("Contact author, DArT tiding problem")
    } else {
      res %<>% dplyr::select(-INDIVIDUALS1)
    }

    res %<>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
        .funs = as.numeric
      ) %>%
      dplyr::mutate(READ_DEPTH = ALLELE_REF_DEPTH + ALLELE_ALT_DEPTH)

    # Coverage for allele calibration
    want <- c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH")
    switch <- dplyr::group_by(res, MARKERS) %>%
      dplyr::summarise_at(.tbl = ., .vars = want, .funs = sum, na.rm = TRUE) %>%
      dplyr::mutate_at(.tbl = ., .vars = want, .funs = round, digits = 0) %>%
      dplyr::mutate_at(.tbl = ., .vars = want, .funs = as.integer) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(dplyr::if_else(ALLELE_ALT_DEPTH > ALLELE_REF_DEPTH, TRUE, FALSE)) %$%
      MARKERS
    n.switch <- length(switch)
    if (n.switch > 0) {
      message("Calibration REF/ALT based on read depth of alleles: ", n.switch)
      res <- dplyr::filter(res, !MARKERS %in% switch) %>%
        dplyr::bind_rows(
          dplyr::filter(res, MARKERS %in% switch) %>%
            dplyr::rename(
              ALT = REF,
              REF = ALT,
              ALLELE_REF_DEPTH = ALLELE_ALT_DEPTH,
              ALLELE_ALT_DEPTH = ALLELE_REF_DEPTH
            )
        )
    }
    switch <- NULL

    # GT_BIN
    res %<>%
      dplyr::mutate(
        GT_BIN = dplyr::case_when(
          ALLELE_REF_DEPTH > 0 & ALLELE_ALT_DEPTH == 0 ~ 0,
          ALLELE_REF_DEPTH > 0 & ALLELE_ALT_DEPTH > 0 ~ 1,
          ALLELE_REF_DEPTH == 0 & ALLELE_ALT_DEPTH > 0 ~ 2
        )
      ) %>%
      dplyr::mutate_at(
        .tbl = .,
        .vars = c("READ_DEPTH", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
        .funs = replace_by_na, what = 0
      )
  }# Counts
  if (gt.vcf) {
    res %<>%
      dplyr::mutate(
        GT_VCF = dplyr::case_when(
          GT_BIN == 0 ~ "0/0", GT_BIN == 1 ~ "0/1", GT_BIN == 2 ~ "1/1",
          is.na(GT_BIN) ~ "./.")
      )
  }
  if (gt.vcf.nuc) {
    res %<>%
      dplyr::mutate(
        GT_VCF_NUC = dplyr::case_when(
          GT_BIN == 0 ~ stringi::stri_join(REF, REF, sep = "/"),
          GT_BIN == 1 ~ stringi::stri_join(REF, ALT, sep = "/"),
          GT_BIN == 2 ~ stringi::stri_join(ALT, ALT, sep = "/"),
          is.na(GT_BIN) ~ "./.")
      )
  }
  if (gt) {
    res %<>%
      dplyr::mutate(
        GT = stringi::stri_replace_all_fixed(
          str = GT_VCF_NUC,
          pattern = c("A", "C", "G", "T", "/", ".."),
          replacement = c("001", "002", "003", "004", "", "000000"),
          vectorize_all = FALSE)
      )
  }

  return(res)
}# End generate_geno_v2

## Merge dart-------------------------------------------------------------------

#' @title Merge DArT files
#' @description This function allows to merge 2 DArT files (filtered of not).

#' @param dart1 Full path of the first DArT file.
#' @param strata1 Full path of the first strata file for dart1.
#' @param dart2 Full path of the second DArT file.
#' @param strata2 Full path of the second strata file for dart2.
#' @param keep.rad Unless the default is changed, the function removes the
#' temporary tidy dart files generated for each DArT datasets during import.
#' Default: \code{keep.rad = FALSE}.
#' @param filename Name of the merged DArT file.
#' By default, the function gives the merged data the filename:\code{merge_dart}
#' with date and time appended. The function will also append the date and time
#' to the filename provided.
#' The data is written in the working directory.
#' Default: \code{filename = NULL}.
#' @inheritParams tidy_genomic_data
#' @inheritParams read_dart
#' @inheritParams radiator_common_arguments
#' @details
#' The function average across markers the columns: CALL_RATE, REP_AVG,
#' AVG_COUNT_REF and AVG_COUNT_SNP, when found in the data.
#' For DArT, theses columns represent:
#' \itemize{
#' \item CALL_RATE: is the proportion of samples for which the genotype was called.
#' \item REP_AVG: is the proportion of technical replicate assay pairs for which
#' the marker score is consistent.
#' \item AVG_COUND_REF and AVG_COUND_SNP: the mean coverage for the reference and
#' alternate alleles, respectively.
#' }
#'
#' The function removes markers with starting with 1000 that are not immortalized by DArT
#'
#'
#' When the argument \strong{common.markers} is kept to TRUE, the function
#' produces an
#' \href{https://github.com/hms-dbmi/UpSetR}{UpSet plot} to visualize the number
#' of markers common or not between populations. The plot is not saved automatically,
#' this as to be done manually by the user.

#' @return The function returns a list in the global environment and 2 data frames
#' in the working directory. The dataframes are the tidy dataset and a strata file
#' of the 2 merged DArT files.

#' @examples
#' \dontrun{
#' # The simplest way to run the function:
#' sum <- radiator::merge_dart(
#' dart1 = "bluefin_larvae.tsv", strata1 = "strata1_bft_larvae.tsv",
#' dart1 = "bluefin_adults.csv", strata2 = "strata2_bft_adults.tsv",
#' filename = "bluefin_combined")
#' }
#' @rdname merge_dart
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

merge_dart <- function(
  dart1, strata1,
  dart2, strata2,
  #pop.select = NULL,
  # filter.monomorphic = TRUE,
  # filter.common.markers = TRUE,
  keep.rad = FALSE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  cat("#######################################################################\n")
  cat("########################## radiator::merge_dart #######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  opt.change <- getOption("width")
  options(width = 70)
  # manage missing arguments -----------------------------------------------------
  if (missing(dart1)) rlang::abort("dart1 file missing")
  if (missing(dart2)) rlang::abort("dart2 file missing")
  if (missing(strata1)) rlang::abort("strata1 file missing")
  if (missing(strata2)) rlang::abort("strata2 file missing")

  # Filename -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  if (is.null(filename)) {
    filename <- stringi::stri_join("merged_dart_", format(Sys.time(), "%Y%m%d@%H%M"), ".rad")
    strata.filename <- stringi::stri_join("merged_dart_strata_", format(Sys.time(), "%Y%m%d@%H%M"), ".tsv")
  } else {
    filename <- stringi::stri_join(filename, "_",format(Sys.time(), "%Y%m%d@%H%M"), ".rad")
    strata.filename <- stringi::stri_join(filename, "_strata_",format(Sys.time(), "%Y%m%d@%H%M"), ".tsv")
  }

  # File type detection---------------------------------------------------------
  data.type <- radiator::detect_genomic_format(dart1)

  # import data ----------------------------------------------------------------
  if (!is.null(pop.select)) {
    pop.select <- clean_pop_names(pop.select)
    message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
  }

  # Columns we want
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "GT_VCF_NUC",
            "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")

  # the DArT data
  if (data.type == "dart") {
    message("Importing and tidying dart1...")
    input <- suppressWarnings(radiator::read_dart(
      data = dart1, strata = strata1,
      filename = "temp.tidy.dart1", verbose = FALSE) %>%
        dplyr::select(dplyr::one_of(want)))
    message("Importing and tidying dart2...")
    dart2 <- suppressWarnings(radiator::read_dart(
      data = dart2, strata = strata2,
      filename = "temp.tidy.dart2", verbose = FALSE) %>%
        dplyr::select(dplyr::one_of(want)))
    if (!keep.rad) {
      message("Removing temporary tidy DArT files...")
      file.remove("temp.tidy.dart1.rad")
      file.remove("temp.tidy.dart2.rad")
    }
  }

  # Import the filtered DArT data
  if (data.type == "fst.file") {
    message("Importing the filtered dart1...")
    input <- suppressWarnings(radiator::read_rad(dart1) %>% dplyr::select(dplyr::one_of(want)))
    message("Importing the filtered dart2...")
    dart2 <- suppressWarnings(radiator::read_rad(dart2) %>% dplyr::select(dplyr::one_of(want)))
  }

  # Keeping selected pop -------------------------------------------------------
  if (!is.null(pop.select)) {
    input <- suppressWarnings(dplyr::filter(input, POP_ID %in% pop.select))
    dart2 <- suppressWarnings(dplyr::filter(dart2, POP_ID %in% pop.select))
  }

  # cleaning up non-immortalized markers ---------------------------------------
  message("Removing markers starting with 1000 (non-immortalized DArT markers)")
  markers.before <- length(unique(input$MARKERS)) + length(unique(dart2$MARKERS))

  input <- suppressWarnings(dplyr::filter(input, !stringi::stri_detect_regex(str = LOCUS, pattern = "^1000")))
  dart2 <- suppressWarnings(dplyr::filter(dart2, !stringi::stri_detect_regex(str = LOCUS, pattern = "^1000")))
  markers.after <- length(unique(input$MARKERS)) + length(unique(dart2$MARKERS))
  message("    Markers removed: ", markers.before - markers.after)

  # merging DArT tidy data -----------------------------------------------------
  # check alt alleles
  # test1 <- dplyr::filter(input, (stringi::stri_count_fixed(str = ALT, pattern = ",") + 1) > 1)
  # test2 <- dplyr::filter(dart2, (stringi::stri_count_fixed(str = ALT, pattern = ",") + 1) > 1)

  # we keep common column
  dart.col <- dplyr::intersect(colnames(input), colnames(dart2))

  input <- suppressWarnings(
    dplyr::select(input, dplyr::one_of(dart.col)) %>%
      dplyr::bind_rows(dplyr::select(dart2, dplyr::one_of(dart.col)))
  )
  dart2 <- NULL

  # Pop select cleanup ---------------------------------------------------------
  if (!is.null(pop.select)) {
    if (is.factor(input$POP_ID)) input$POP_ID <- droplevels(input$POP_ID)
  }

  # Remove weird markers
  input <- radiator::detect_biallelic_problems(data = input, verbose = FALSE, parallel.core = parallel.core)$biallelic.data

  # Averaging across markers the call rate and other DArT markers metadata statistics
  # Note to myself: Might be easier/faster to use mutate_if
  if (rlang::has_name(input, "CALL_RATE")) {
    message("Averaging across markers the call rate")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(CALL_RATE = mean(CALL_RATE)) %>%
      dplyr::ungroup(.)
  }

  if (rlang::has_name(input, "REP_AVG")) {
    message("Averaging across markers the REP_AVG")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(REP_AVG = mean(REP_AVG)) %>%
      dplyr::ungroup(.)
  }

  message("Adjusting REF/ALT alleles...")
  input <- radiator::calibrate_alleles(
    data = input, biallelic = NULL,
    parallel.core = parallel.core,
    verbose = TRUE)$input

  if (rlang::has_name(input, "POLYMORPHIC")) {
    input <- dplyr::select(input, -POLYMORPHIC)
  }

  if (rlang::has_name(input, "AVG_COUNT_REF")) {
    message("Averaging across markers the coverage for the REF allele")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(AVG_COUNT_REF = mean(AVG_COUNT_REF)) %>%
      dplyr::ungroup(.)
  }

  if (rlang::has_name(input, "AVG_COUNT_ALT")) {
    message("Averaging across markers the coverage for the ALT allele")
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(AVG_COUNT_REF = mean(AVG_COUNT_REF)) %>%
      dplyr::ungroup(.)
  }

  if (filter.monomorphic) {
    test <- radiator::filter_monomorphic(data = input, verbose = TRUE)
  }

  if (filter.common.markers) {
    input <- radiator::filter_common_markers(data = input, fig = TRUE, verbose = TRUE) %$% input
    message("\n    **** Manually save the figure ****")
  }

  # Write tidy in the working directory
  radiator::write_rad(data = input, path = filename)

  strata <- dplyr::select(input, INDIVIDUALS, STRATA = POP_ID) %>%
    dplyr::distinct(INDIVIDUALS, STRATA) %>%
    readr::write_tsv(x = ., path = strata.filename)

  # results --------------------------------------------------------------------
  message("\nMerged DArT file and strata file in the working directory:")
  message("    ", filename)
  message("    ", strata.filename)
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  options(width = opt.change)
  return(res = list(merged.dart = input, strata = strata))
}#End merge_dart

# tidy_dart_metadata--------------------------------------------------------------
#' @name tidy_dart_metadata

#' @title Import and tidy \href{http://www.diversityarrays.com}{DArT} metadata.

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' The function generate a tidy dataset of
#' \href{http://www.diversityarrays.com}{DArT} markers and associated metadata.
#' Usefull to filter before importing the actual dataset.

#' @param data DArT output file. Note that most popular formats used by DArT are
#' recognised (1- and 2- row format, also called binary, and count data.).
#' If you encounter a problem, sent me your data so that I can update
#' the function. The function can import \code{.csv} or \code{.tsv} files.


#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments

#' @return A tidy dataframe with these columns:
#' \enumerate{
#' \item MARKERS: generated by radiator and correspond to CHROM + LOCUS + POS
#' separated by 2 underscores.
#' \item CHROM: the chromosome, for de novo: CHROM_1.
#' \item LOCUS: the locus.
#' \item POS: the SNP id on the LOCUS.
#' \item REF: the reference allele.
#' \item ALT: the alternate allele.
#' \item CALL_RATE: call rate output specific of DArT.
#' \item AVG_COUNT_REF: the coverage for the reference allele, output specific of DArT.
#' \item AVG_COUNT_SNP: the coverage for the alternate allele, output specific of DArT.
#' \item REP_AVG: the reproducibility average, output specific of DArT.
#' }

#' @export
#' @rdname tidy_dart_metadata

#' @examples
#' \dontrun{
#' clownfish.dart.tidy <- radiator::tidy_dart_metadata(
#' data = "clownfish.dart.tsv",
#' verbose = TRUE)
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_dart_metadata <- function(
  data,
  filename = NULL,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1
) {
  # Cleanup-------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  # res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nTiming: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("############################## completed ##############################\n"), add = TRUE)


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")
  # Date and Time --------------------------------------------------------------
  if (!is.null(filename)) {
    meta.filename <- stringi::stri_join(filename, "_metadata_", file.date,".rad")
  }

  data.type <- radiator::detect_genomic_format(data)

  if (!data.type %in% c("dart", "fst.file")) {
    rlang::abort("Contact author to show your DArT data, problem duting import")
  }

  if (verbose) message("Importing DArT markers metadata")

  # Import metadata-------------------------------------------------------------
  if (data.type == "dart") {
    if (stringi::stri_detect_fixed(
      str = stringi::stri_sub(str = data, from = -4, to = -1),
      pattern = ".csv")) {
      csv <- TRUE
    } else {
      csv <- FALSE
    }

    dart.check <- check_dart(data)
    if (!dart.check$data.type %in% c("dart", "silico.dart")) {
      rlang::abort("Contact author to show your DArT data, problem during import")
    } else {
      skip.number <- dart.check$skip.number
    }


    if (csv) {
      dart.col.type <- readr::read_csv(
        file = data,
        skip = skip.number, n_max = 1,
        na = "-",
        col_names = FALSE,
        col_types = readr::cols(.default = readr::col_character()))
    } else {
      dart.col.type <- readr::read_tsv(
        file = data,
        skip = skip.number, n_max = 1,
        na = "-",
        col_names = FALSE,
        col_types = readr::cols(.default = readr::col_character()))
    }

    want <- tibble::tibble(
      INFO = c("ALLELEID", "SNP", "SNPPOSITION", "CALLRATE",
               "AVGCOUNTREF", "AVGCOUNTSNP", "REPAVG"),
      COL_TYPE = c("c", "c", "i", "d", "d", "d", "d"))

    dart.col.type <- dart.col.type %>%
      tidyr::gather(data = .,key = DELETE, value = INFO) %>%
      dplyr::select(-DELETE) %>%
      dplyr::mutate(INFO = stringi::stri_trans_toupper(INFO)) %>%
      dplyr::left_join(want, by = "INFO") %>%
      dplyr::mutate(COL_TYPE = stringi::stri_replace_na(str = COL_TYPE, replacement = "_")) %>%
      dplyr::select(COL_TYPE) %>%
      purrr::flatten_chr(.) %>% stringi::stri_join(collapse = "")

    if (csv) {
      input <- suppressMessages(suppressWarnings(
        readr::read_csv(
          file = data,
          skip = skip.number,
          na = "-",
          col_names = TRUE,
          col_types = dart.col.type)
      ))
    } else {
      input <- suppressMessages(suppressWarnings(
        readr::read_tsv(
          file = data,
          skip = skip.number,
          na = "-",
          col_names = TRUE,
          col_types = dart.col.type)
      ))
    }

    colnames(input) <- stringi::stri_trans_toupper(colnames(input))
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = c("AVGCOUNTREF", "AVGCOUNTSNP", "REPAVG", "ALLELEID", "SNPPOSITION", "CALLRATE"),
      replacement = c("AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG", "LOCUS", "POS", "CALL_RATE"),
      vectorize_all = FALSE)

    # Check for duplicate rows (sometimes people combine DArT data...)----------
    input.dup <- dplyr::distinct(input, LOCUS, SNP, POS, CALL_RATE, .keep_all = FALSE)

    # make sure no duplicates
    if (nrow(input) != nrow(input.dup)) {
      input.dup <- NULL
      message("Duplicate rows were identified")
      message("    using distinct rows")
      message("    check input data if downstream problems")
      input <- dplyr::distinct(input, LOCUS, SNP, POS, CALL_RATE, .keep_all = TRUE)
    }
    input.dup <- NULL

    # Tidying data ---------------------------------------------------------------
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT", "CALL_RATE",
              "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")

    input <- suppressWarnings(
      input %>%
        tidyr::separate(col = LOCUS, into = c("LOCUS", "NOT_USEFUL"), sep = "\\|", extra = "drop") %>%
        dplyr::select(-NOT_USEFUL) %>%
        tidyr::separate(col = SNP, into = c("NOT_USEFUL", "KEEPER"), sep = ":", extra = "drop") %>%
        dplyr::select(-NOT_USEFUL) %>%
        tidyr::separate(col = KEEPER, into = c("REF", "ALT"), sep = ">") %>%
        dplyr::mutate(
          CHROM = rep("CHROM_1", n()),
          MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__")) %>%
        dplyr::select(dplyr::one_of(want), dplyr::everything()) %>%
        dplyr::mutate_at(.tbl = ., .vars = c("MARKERS", "CHROM", "LOCUS", "POS"), .funs = as.character) %>%
        dplyr::arrange(CHROM, LOCUS, POS, REF) %>%
        dplyr::filter(!is.na(REF) | !is.na(ALT)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::arrange(MARKERS))
  }

  if (data.type == "fst.file") {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT", "CALL_RATE",
              "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG")
    input <- suppressWarnings(
      radiator::read_rad(data) %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE))
  }

  if (!is.null(filename)) {
    radiator::write_rad(data = input, path = meta.filename)
    short.name <- list.files(path = ".", pattern = "metadata")
    if (length(short.name) > 1) {
      short.name <- file.info(short.name) %>%
        tibble::rownames_to_column(df = ., var = "FILE") %>%
        dplyr::filter(mtime == max(mtime))
      short.name <- short.name$FILE
    }
    message("Marker's metadata file written:\n    ", short.name)
  }
  # Results --------------------------------------------------------------------
  if (verbose) {
    n.chrom <- length(unique(input$CHROM))
    n.locus <- length(unique(input$LOCUS))
    n.snp <- length(unique(input$MARKERS))
    message("\nNumber of chrom: ", n.chrom)
    message("Number of locus: ", n.locus)
    message("Number of SNPs: ", n.snp)
  }
  return(input)
}#End dart_markers_metadata


# generate_geno----------------------------------------------------------------
#' @title clean_dart_colnames
#' @description clean_dart_colnames (only the DArT columns = snakecase...)
#' @rdname clean_dart_colnames
#' @keywords internal
#' @export
clean_dart_colnames <- function(data, strata) {
  keeper <- length(colnames(data)) - length(strata$TARGET_ID)
  colnames(data) <- c(
    radiator::radiator_snakecase(x = colnames(data)[1:keeper]),
    strata$TARGET_ID
  )

  colnames(data) <- tibble::tibble(TARGET_ID = colnames(data)) %>%
    dplyr::left_join(strata, by = "TARGET_ID") %>%
    dplyr::mutate(
      INDIVIDUALS = dplyr::if_else(
        is.na(INDIVIDUALS), TARGET_ID, INDIVIDUALS)
      ) %$% INDIVIDUALS

  # Below generate errors when some id are very close... ID-10 and ID-1
  # colnames(data) <- stringi::stri_replace_all_fixed(
  #   str = colnames(data),
  #   pattern = strata$TARGET_ID,
  #   replacement = strata$INDIVIDUALS,
  #   vectorize_all = FALSE
  #   )
  return(data)
}#End clean_dart_colnames
