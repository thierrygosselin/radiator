# Read PLINK -------------------------------------------------------------------
#' @name read_plink
#' @title Reads PLINK tped and bed files
#' @description The function reads PLINK tped and bed files.
#' radiator prefers the use of BED file. These files are converted to
#' a connection SeqArray \href{https://github.com/zhengxwen/SeqArray}{SeqArray}
#' GDS object/file of class \code{SeqVarGDSClass} (Zheng et al. 2017).
#' The Genomic Data Structure (GDS) file format is detailed in
#' \href{https://github.com/zhengxwen/gdsfmt}{gdsfmt}.
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data The PLINK file.
#' \itemize{
#' \item {bi-allelic data only}. For haplotypes use VCF.
#' \item \code{tped} file format: the corresponding \code{tfam} file must be in the directory.
#' \item \code{bed} file format: IS THE PREFERRED format, the corresponding
#' \code{fam} and \code{bim} files must be in the directory.
#' }

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' radiator will append \code{.gds.rad} to the filename.
#' If the filename chosen exists in the working directory,
#' the default \code{radiator_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments


#' @details
#' Large PLINK files will require the use of BED plink format. Look below
#' in the example for conversion with PLINK.
#'
#' Large PLINK bed files will take longer to import and transform in GDS, but
#' after the file is generated, you can close your computer and
#' come back to it a month later and it's now a matter of sec to open a connection!


#### To do ....

# @section Advance mode:
#
# \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
# \enumerate{
# \item \code{path.folder}: to write ouput in a specific path
# (used internally in radiator). Default: \code{path.folder = getwd()}.
# If the supplied directory doesn't exist, it's created.
# \item \code{random.seed}: (integer, optional) For reproducibility, set an integer
# that will be used inside codes that uses randomness. With default,
# a random number is generated, printed and written in the directory.
# Default: \code{random.seed = NULL}.
# \item \code{subsample.markers.stats}: By default, when no filters are
# requested and that the number of markers is > 200K,
# 0.2 of markers are randomly selected to generate the
# statistics (individuals and markers). This is an all-around and
# reliable number.
# In doubt, overwrite this value by using 1 (all markers selected) and
# expect a small computational cost.
# }

#' @return
#' For \code{tped} the function returns a list object with the non-modified \code{tped} and
#' the strata corresponding to the \code{tfam}.
#' With \code{bed}, the function returns a GDS object.

#' @export
#' @rdname read_plink

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.
#'
#' @references
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795


#' @examples
#' \dontrun{
#' data <- radiator::read_plink(data = "my_plink_file.bed")
#' # when conversion is required from TPED to BED, in Terminal:
#' # plink --tfile my_plink_file --make-bed --allow-no-sex --allow-extra-chr --chr-set 95
#' }

#' @seealso
#' \href{https://www.cog-genomics.org/plink/1.9/}{PLINK}


read_plink <- function(
  data,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  # #Test
  # filename <- NULL
  # parallel.core <- parallel::detectCores() - 1
  # verbose <- TRUE

  # Cleanup---------------------------------------------------------------------
  radiator_function_header(f.name = "read_plink", verbose = verbose)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(future.globals.maxSize = Inf)
  options(width = 70)
  timing <- radiator_tic()
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(if (verbose) radiator_toc(timing), add = TRUE)
  on.exit(radiator_function_header(f.name = "read_plink", start = FALSE, verbose = verbose), add = TRUE)

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("PLINK file missing")

  # Warning concerning plink files----------------------------------------------
  plink.format <- radiator::detect_genomic_format(data = data)

  # PLINK TPED -----------------------------------------------------------------
  if (plink.format == "plink.tped.file") {
    message("Reading PLINK tped...")

    # Get file size
    plink.size <- file.size(data)
    if (!plink.size <= 100000000) {
      message("\n\nNote: huge PLINK tped files are no longer recommended inside radiator")
      message("Convert to PLINK BED with: --tfile [] --make-bed --allow-no-sex --allow-extra-chr --chr-set [] with PLINK command")
      rlang::abort("Read the documentation of radiator::read_plink")
    }

    # Importing file -----------------------------------------------------------
    fam.file <- stringi::stri_replace_all_fixed(
      str = data,
      pattern = ".tped",
      replacement = ".tfam",
      vectorize_all = FALSE
    )
    res$strata <- readr::read_delim(
      file = fam.file,
      delim = " ",
      col_names = c("STRATA", "INDIVIDUALS"),
      col_types = "cc____"
    ) %>%
      radiator::read_strata(strata = .) %$%
      strata
    fam.file <- NULL

    # preparing header for tped file
    tped.header.prep <- res$strata %>%
      dplyr::select(INDIVIDUALS) %>%
      dplyr::mutate(
        NUMBER = seq(1, n()),
        ALLELE1 = rep("A1", n()), ALLELE2 = rep("A2", n())
      ) %>%
      radiator::rad_long(
        x = .,
        cols = c("INDIVIDUALS", "NUMBER"),
        names_to = "ALLELES_GROUP",
        values_to = "ALLELES"
        ) %>%
      dplyr::arrange(NUMBER) %>%
      dplyr::select(-ALLELES_GROUP) %>%
      tidyr::unite(INDIVIDUALS_ALLELES, c(INDIVIDUALS, ALLELES), sep = "_", remove = FALSE) %>%
      dplyr::arrange(NUMBER) %>%
      dplyr::mutate(NUMBER = seq(from = (1 + 4), to = n() + 4)) %>%
      dplyr::select(-ALLELES)
    tped.header.names <- c("CHROM", "LOCUS", "POS", tped.header.prep$INDIVIDUALS_ALLELES)
    tped.header.integer <- c(1, 2, 4, tped.header.prep$NUMBER)
    tped.header.prep <- NULL

    # import PLINK
    res$data <- data.table::fread(
      input = data,
      sep = " ",
      header = FALSE,
      stringsAsFactors = FALSE,
      verbose = FALSE,
      select = tped.header.integer,
      col.names = tped.header.names,
      showProgress = TRUE,
      data.table = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(
        CHROM = as.character(CHROM),
        LOCUS = as.character(LOCUS),
        POS = as.character(POS)
        ) %>%
      tidyr::unite(data = ., col = "MARKERS", c("CHROM", "LOCUS", "POS"), remove = FALSE, sep = "__")

    # Unused objects
    tped.header.integer <- tped.header.names <- NULL
    return(res)
  }

  # PLINK BED ------------------------------------------------------------------
  if (plink.format == "plink.bed.file") {
    message("Reading PLINK bed file...")

    # Required package -----------------------------------------------------------
    radiator_packages_dep(package = "SeqVarTools", cran = FALSE, bioc = TRUE)

    # Function call and dotslist -------------------------------------------------
    rad.dots <- radiator_dots(
      func.name = as.list(sys.call())[[1]],
      fd = rlang::fn_fmls_names(),
      args.list = as.list(environment()),
      dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
      keepers = c("internal", "path.folder", "parameters"),
      verbose = FALSE
    )

    # Folders---------------------------------------------------------------------
    wf <- path.folder <- generate_folder(
      f = path.folder,
      rad.folder = "read_plink",
      prefix_int = FALSE,
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    radiator.folder <- generate_folder(
      f = path.folder,
      rad.folder = "import_gds",
      prefix_int = TRUE,
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    # write the dots file
    write_rad(
      data = rad.dots,
      path = radiator.folder,
      filename = stringi::stri_join("radiator_read_plink_args_", file.date, ".tsv"),
      tsv = TRUE,
      internal = internal,
      write.message = "Function call and arguments stored in: ",
      verbose = verbose
    )

    filename <- generate_filename(
      name.shortcut = filename,
      path.folder = radiator.folder,
      extension = "gds")

    filename.short <- filename$filename.short
    filename <- filename$filename

    timing.plink <- proc.time()

    # read_plink -----------------------------------------------------------------
    bim.file <- stringi::stri_replace_all_fixed(
      str = data,
      pattern = ".bed",
      replacement = ".bim",
      vectorize_all = FALSE
    )
    fam.file <- stringi::stri_replace_all_fixed(
      str = data,
      pattern = ".bed",
      replacement = ".fam",
      vectorize_all = FALSE
    )

    if (verbose) message("Generating GDS...")
    gds <- SeqArray::seqBED2GDS(
      bed.fn = data,
      fam.fn = fam.file,
      bim.fn = bim.file,
      out.gdsfn = filename,
      compress.geno = "ZIP_RA",
      compress.annotation = "ZIP_RA",
      verbose = verbose
    ) %>%
      SeqArray::seqOpen(gds.fn = ., readonly = FALSE)
    if (verbose) message("done! timing: ", round((proc.time() - timing.plink)[[3]]), " sec\n")

    # PLINK: Summary ----------------------------------------------------------------
    summary_gds(gds = gds, verbose = TRUE)
    if (verbose) message("\nFile written: ", filename.short)

    # Generate radiator skeleton -------------------------------------------------
    radiator.gds <- radiator_gds_skeleton(gds)

    # data.source ----------------------------------------------------------------
    update_radiator_gds(gds = gds, node.name = "data.source", value = "plink")

    # bi- or multi-alllelic--------------------------------------------------
    biallelic <- detect_biallelic_markers(data = gds, verbose = verbose)

    # Clean sample id---------------------------------------------------------
    individuals.plink <- tibble::tibble(
      INDIVIDUALS_PLINK = SeqArray::seqGetData(gds, "sample.id")) %>%
      dplyr::mutate(INDIVIDUALS_CLEAN = radiator::clean_ind_names(INDIVIDUALS_PLINK))

    if (!identical(individuals.plink$INDIVIDUALS_PLINK, individuals.plink$INDIVIDUALS_CLEAN)) {
      if (verbose) message("Cleaning PLINK's sample names")
      clean.id.filename <- stringi::stri_join("cleaned.plink.id.info_", file.date, ".tsv")
      readr::write_tsv(x = individuals.plink,
                       file = file.path(radiator.folder, clean.id.filename))

      update_radiator_gds(gds = gds, node.name = "id.clean", value = individuals.plink)
    }
    # replace id
    update_radiator_gds(
      gds = gds,
      radiator.gds = FALSE,
      node.name = "sample.id",
      value = individuals.plink$INDIVIDUALS_CLEAN,
      replace = TRUE)

    individuals <- dplyr::select(individuals.plink, INDIVIDUALS = INDIVIDUALS_CLEAN)
    individuals.plink <- NULL

    # sync id with STRATA---------------------------------------------------------
    if (verbose) message("Using .fam file for strata...")
    strata <- readr::read_delim(
      file = fam.file,
      delim = " ",
      col_names = c("STRATA", "INDIVIDUALS"),
      col_types = "cc____"
    ) %>%
      radiator::read_strata(strata = .) %$%
      strata

    id.levels <- individuals$INDIVIDUALS
    individuals %<>%
      dplyr::left_join(
        join_strata(individuals, strata, verbose = verbose) %>%
          dplyr::mutate(FILTERS = "whitelist")
        , by = "INDIVIDUALS"
      ) %>%
      dplyr::mutate(FILTERS = tidyr::replace_na(data = FILTERS, replace = "filter.stata"))

    strata <- generate_strata(data = dplyr::filter(individuals, FILTERS == "whitelist"), pop.id = FALSE)

    #Update GDS node
    update_radiator_gds(gds = gds, node.name = "individuals.meta", value = individuals, sync = TRUE)

    # PLINK: Markers metadata  ------------------------------------------------------
    markers.meta <- extract_markers_metadata(gds = gds)

    # PLINK: reference genome or de novo -------------------------------------------
    ref.genome <- detect_ref_genome(data = gds, verbose = verbose)

    # Generate MARKERS column and fix types --------------------------------------
    markers.meta %<>%
      dplyr::mutate(
        dplyr::across(
          .cols = c(CHROM, LOCUS, POS),
          .fns = radiator::clean_markers_names
        )
      ) %>%
      dplyr::mutate(
        MARKERS = stringi::stri_join(CHROM, LOCUS, POS, sep = "__"),
        REF = SeqArray::seqGetData(gdsfile = gds, var.name = "$ref"),
        ALT = SeqArray::seqGetData(gdsfile = gds, var.name = "$alt")
      )

    # PLINK file with duplicate markers... sometimes tagged ISOFORMS...
    dup.markers <- length(markers.meta$MARKERS) - length(unique(markers.meta$MARKERS))
    if (dup.markers > 0) {
      message("\nNumber of duplicate MARKERS id: ", dup.markers)
      message("Adding integer to differentiate...")
      markers.meta %<>%
        dplyr::arrange(MARKERS) %>%
        dplyr::mutate(MARKERS_NEW = MARKERS) %>%
        dplyr::group_by(MARKERS_NEW) %>%
        dplyr::mutate(
          MARKERS = stringi::stri_join(MARKERS, seq(1, n(), by = 1), sep = "_")
        ) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-MARKERS_NEW)
    }

    # ADD MARKERS META to GDS
    update_radiator_gds(gds = gds, node.name = "markers.meta", value = markers.meta)

    # replace chromosome info in GDS
    # Why ? well snp ld e.g. will otherwise be performed by chromosome and with de novo data = by locus...
    update_radiator_gds(
      gds = gds,
      radiator.gds = FALSE,
      node.name = "chromosome",
      value = markers.meta$CHROM,
      replace = TRUE
    )
    # # radiator_parameters: generate --------------------------------------------
    filters.parameters <- radiator_parameters(
      generate = TRUE,
      initiate = FALSE,
      update = FALSE,
      parameter.obj = parameters,
      path.folder = radiator.folder,
      file.date = file.date,
      verbose = verbose,
      internal = internal)

    # radiator_parameters: initiate --------------------------------------------
    # with original PLINK values
    filters.parameters <- radiator_parameters(
      generate = FALSE,
      initiate = TRUE,
      update = TRUE,
      parameter.obj = filters.parameters,
      data = gds,
      filter.name = "plink",
      param.name = "original values in PLINK + strata",
      values = "",
      path.folder = path.folder,
      file.date = file.date,
      internal = internal,
      verbose = verbose
    )
    return(gds)
  } # end plink's BED format
} # End read_plink


# tidy_plink -------------------------------------------------------------------
#' @name tidy_plink
#' @title Tidy PLINK tped and bed files

#' @description Transform bi-allelic PLINK files in \code{.tped} or {.bed} formats
#' into a tidy dataset.
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data The PLINK file.
#' \itemize{
#' \item {bi-allelic data only}. For haplotypes use VCF.
#' \item \code{tped} file format: the corresponding \code{tfam} file must be in the directory.
#' \item \code{bed} file format: IS THE PREFERRED format, the corresponding
#' \code{fam} and \code{bim} files must be in the directory.
#' }

#' @inheritParams radiator_common_arguments

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#' \item \code{calibrate.alleles}: logical. For \code{tped} files, if
#' \code{calibrate.alleles = FALSE} the function runs faster
#' but REF/ALT alleles may not be calibrated. The default assumes the users or
#' sotware producing the PLINK file calibrated the alleles.
#' Default: \code{calibrate.alleles = FALSE}.
#' }


#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.
#'
#' @references
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795

#' @examples
#' \dontrun{
#' data <- radiator::tidy_plink(data = "my_plink_file.bed", verbose = TRUE)
#'
#'
#' # when conversion is required from TPED to BED, in Terminal:
#' # plink --tfile my_plink_file --make-bed --allow-no-sex --allow-extra-chr --chr-set 95
#' }

#' @seealso
#' \href{https://www.cog-genomics.org/plink/1.9/}{PLINK}
#'
#' \code{\link[radiator]{read_plink}}

#' @return
#' A tidy tibble of the PLINK file.
#'


#' @export
#' @rdname tidy_plink
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_plink <- function(
  data,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...
) {
  # Test
  # verbose = TRUE
  # strata = NULL
  # calibrate.alleles = TRUE
  # parallel.core = parallel::detectCores() - 1

  # Cleanup---------------------------------------------------------------------
  radiator_function_header(f.name = "tidy_plink", verbose = verbose)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(future.globals.maxSize = Inf)
  options(width = 70)
  timing <- radiator_tic()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(radiator_toc(timing), add = TRUE)
  on.exit(radiator_function_header(f.name = "tidy_plink", start = FALSE, verbose = verbose), add = TRUE)

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("PLINK file missing")

  # Function call and dotslist -------------------------------------------------
  path.folder <- filename <- calibrate.alleles <- NULL
  rad.dots <- radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    deprecated = c("blacklist.id", "pop.select", "pop.levels", "pop.labels"),
    keepers = c("calibrate.alleles", "filename", "internal", "path.folder", "parameters"),
    verbose = FALSE
  )

  # Warning concerning plink files----------------------------------------------
  plink.format <- radiator::detect_genomic_format(data = data)

  data <- radiator::read_plink(
    data = data,
    filename = filename,
    parallel.core = parallel.core,
    verbose = FALSE,
    internal = TRUE,
    path.folder = path.folder,
    parameters = parameters
  )

  if (plink.format == "plink.tped.file") {
    # Make tidy ------------------------------------------------------------------
    if (verbose) message("Tidying the PLINK tped file ...")
    # Filling GT and new separating INDIVIDUALS from ALLELES
    # combining alleles
    strata <- data$strata
    data <- radiator::rad_long(
      x = data$data,
      cols = c("MARKERS", "CHROM", "LOCUS", "POS"),
      names_to = "INDIVIDUALS_ALLELES",
      values_to = "GT"
    )

    # detect GT coding
    detect.gt.coding <- unique(sample(x = data$GT, size = 100, replace = FALSE))
    gt.letters <- c("A", "C", "G", "T")

    if (TRUE %in% unique(gt.letters %in% detect.gt.coding)) {
      if (verbose) message("Genotypes coded with letters")
      gt.letters.df <- tibble::tibble(
        GT = c("A", "C", "G", "T", "0"),
        NEW_GT = c("001", "002", "003", "004", "000")
      )
      data  %<>%
        dplyr::left_join(
          gt.letters.df, by = "GT") %>%
        dplyr::select(-GT) %>%
        dplyr::rename(GT = NEW_GT)
      gt.letters.df <- NULL
    } else {
      if (verbose) message("Genotypes coded with integers")
      data  %<>%
        dplyr::mutate(GT = stringi::stri_pad_left(str = GT, width = 3, pad = "0"))
    }
    detect.gt.coding <- gt.letters <- NULL

    data %<>%
      tidyr::separate(
        data = .,
        col = INDIVIDUALS_ALLELES,
        into = c("INDIVIDUALS", "ALLELES"),
        sep = "_") %>%
      radiator::rad_wide(
        x = .,
        formula = "MARKERS + CHROM + LOCUS + POS + INDIVIDUALS ~ ALLELES",
        values_from = "GT"
        ) %>%
      dplyr::ungroup(.) %>%
      tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>%
      dplyr::select(MARKERS, CHROM, LOCUS, POS, INDIVIDUALS, GT)

    # population levels and strata
    if (verbose) message("Integrating the tfam/strata file...")

    data %<>% dplyr::left_join(strata, by = "INDIVIDUALS")
    strata <- NULL

    # removing untyped markers across all-pop
    remove.missing.gt <- data %>%
      dplyr::select(LOCUS, GT) %>%
      dplyr::filter(GT != "000000")

    untyped.markers <- dplyr::n_distinct(data$LOCUS) - dplyr::n_distinct(remove.missing.gt$LOCUS)
    if (untyped.markers > 0) {
      if (verbose) message("Number of marker with 100 % missing genotypes: ", untyped.markers)
      data <- suppressWarnings(
        dplyr::semi_join(data,
                         remove.missing.gt %>%
                           dplyr::distinct(LOCUS, .keep_all = TRUE),
                         by = "LOCUS")
      )
    }

    # Unused objects
    remove.missing.gt <- NULL

    # detect if biallelic give vcf style genotypes
    # biallelic <- radiator::detect_biallelic_markers(input)
    # filename <- internal <- parameters <- path.folder <- calibrate.alleles <- NULL
    # rm(filename, internal, parameters, path.folder, calibrate.alleles)

    if (calibrate.alleles) {
      data %<>% radiator::calibrate_alleles(data = ., verbose = verbose)
      return(res = list(input = data$input, biallelic = data$biallelic))
    } else {
      return(res = list(input = data, biallelic = radiator::detect_biallelic_markers(data)))
    }
  } #End tidy tped

  if (plink.format == "plink.bed.file") {
    # tidy_plink folder ------------------------------------------------------------
    # tidy.folder <- generate_folder(
    #   f = path.folder,
    #   rad.folder = "tidy_plink",
    #   prefix_int = TRUE,
    #   internal = FALSE,
    #   file.date = file.date,
    #   verbose = verbose)

    # write the dots file: after the GDS import...
    # write_rad(
    #   data = rad.dots,
    #   path = tidy.folder,
    #   filename = stringi::stri_join("radiator_tidy_plink_args_", file.date, ".tsv"),
    #   tsv = TRUE,
    #   internal = internal,
    #   verbose = verbose
    # )
    # Get info markers and individuals -----------------------------------------
    gds.info <- radiator::summary_gds(gds = data, verbose = FALSE)
    # markers.meta <- extract_markers_metadata(gds = data, whitelist = TRUE)
    # individuals <- extract_individuals_metadata(
    #   gds = data,
    #   ind.field.select = "INDIVIDUALS",
    #   whitelist = TRUE
    # )
    n.markers <- gds.info$n.markers
    n.individuals <- gds.info$n.ind

    cat("\n\n################################## IMPORTANT ###################################\n")
    message("Tidying PLINK file with ", n.markers, " SNPs is not optimal:")
    message("    1. a computer with lots of RAM is required")
    message("    2. it's very slow to generate")
    message("    3. it's very slow to run codes after")
    message("    4. for most non model species this number of markers is not realistic...")
    message("\nRecommendation:")
    message("    1. filter your dataset. e.g. with filter_rad")
    message("\nIdeally target a maximum of ~ 10 000 - 20 0000 unlinked SNPs\n")

    message("\nGenotypes formats generated with ", n.markers, " SNPs: ")
    message("    GT_BIN (the dosage of ALT allele: 0, 1, 2 NA)")

    tidy.data <- gds2tidy(
      gds = data,
      markers.meta = NULL,
      calibrate.alleles = FALSE
    )
    # checks
    # dplyr::n_distinct(tidy.data$INDIVIDUALS)
    # dplyr::n_distinct(tidy.data$MARKERS)
    # filename <- internal <- parameters <- path.folder <- calibrate.alleles <- NULL
    # rm(filename, internal, parameters, path.folder, calibrate.alleles)
    return(res = list(input = tidy.data, biallelic = "biallelic"))
  }#End tidy bed


} # End import PLINK


# write_plink ------------------------------------------------------------------

#' @name write_plink
#' @title Write a plink tped/tfam file from a tidy data frame

#' @description Write a plink file from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.

#' @param filename (optional) The file name prefix for tped/tfam files
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_plink_}.

#' @export
#' @rdname write_plink
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
#' Bender D, et al.
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_plink <- function(data, filename = NULL) {

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) data %<>% radiator::tidy_wide(data = ., import.metadata = TRUE)

  tped <- data %>%
    dplyr::arrange(INDIVIDUALS) %>%
    dplyr::mutate(
      COL1 = rep("0", n()),
      COL3 = rep("0", n()),
      COL4 = rep("0", n())
    ) %>%
    dplyr::select(COL1, MARKERS, COL3, COL4, INDIVIDUALS, GT) %>%
    dplyr::mutate(
      A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
      A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
    ) %>%
    dplyr::select(-GT) %>%
    radiator::rad_long(
      x = .,
      cols = c("COL1", "MARKERS", "COL3", "COL4", "INDIVIDUALS"),
      names_to = "ALLELES",
      values_to = "GENOTYPE"
      ) %>%
    dplyr::mutate(
      GENOTYPE = as.character(as.numeric(GENOTYPE)),
      GENOTYPE = stringi::stri_pad_left(GENOTYPE, width = 2, pad = "0")
    ) %>%
    dplyr::arrange(INDIVIDUALS, ALLELES) %>%
    tidyr::unite(INDIVIDUALS_ALLELES, INDIVIDUALS, ALLELES, sep = "_") %>%
    radiator::rad_wide(x = ., formula = "COL1 + MARKERS +COL3 + COL4 ~ INDIVIDUALS_ALLELES", values_from = "GENOTYPE") %>%
    dplyr::arrange(MARKERS)

  tfam <- dplyr::distinct(.data = data, POP_ID, INDIVIDUALS) %>%
    dplyr::arrange(INDIVIDUALS) %>%
    dplyr::mutate(
      COL3 = rep("0",n()),
      COL4 = rep("0",n()),
      COL5 = rep("0",n()),
      COL6 = rep("-9",n())
    )

  # Create a filename to save the output files ********************************
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
    filename.tped <- stringi::stri_join("radiator_plink_", file.date, ".tped")
    filename.tfam <- stringi::stri_join("radiator_plink_", file.date, ".tfam")
  } else {
    filename.tped <- stringi::stri_join(filename, ".tped")
    filename.tfam <- stringi::stri_join(filename, ".tfam")
  }
  readr::write_delim(x = tped, file = filename.tped, col_names = FALSE, delim = " ")
  readr::write_delim(x = tfam, file = filename.tfam, col_names = FALSE, delim = " ")
} # end write_plink
