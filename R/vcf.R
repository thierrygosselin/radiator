# extract_individuals_vcf-------------------------------------------------------
#' @title Extract individuals from vcf file
#' @description Function that returns the individuals present in a vcf file.
#' Useful to create a strata file or
#' to make sure you have the right individuals in your VCF.
#' @param data (character) The path to the vcf file.
#' @rdname extract_individuals_vcf
#' @export
#' @return A tibble with a column: \code{INDIVIDUALS}.
#' @seealso \pkg{radiator} \code{\link{read_strata}}
#' @author Thierry Gosselin \email{thierrygosselin@icloud.com}
extract_individuals_vcf <- function(data) {
  temp.file <-
    suppressWarnings(suppressMessages(
      readr::read_table(file = data, n_max = 200, col_names = "HEADER")
      ))
  skip.number <- which(stringi::stri_detect_fixed(str = temp.file$HEADER,
                                                  pattern = "#CHROM")) - 1
  temp.file <- NULL
  remove <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  id <- tibble::data_frame(INDIVIDUALS = colnames(readr::read_tsv(
    file = data,
    n_max = 1,
    skip = skip.number,
    col_types = readr::cols(.default = readr::col_character())) %>%
      dplyr::select(-dplyr::one_of(remove))))
  return(id)
}#End extract_individuals_vcf

# extract_info_vcf-------------------------------------------------------

#' @title extract_info_vcf
#' @description Extract vcf information
#' @rdname extract_info_vcf
#' @keywords internal
#' @export
extract_info_vcf <- function(vcf) {
  res <- list()
  # print(vcf, all=TRUE, attribute=TRUE)
  # tic()
  vcf.info <- SeqArray::seqVCF_Header(vcf.fn = vcf, getnum = TRUE)
  # toc()
  res$vcf.source <- vcf.info$header$value[2]
  res$n.ind <- vcf.info$num.sample
  res$n.markers <- vcf.info$num.variant
  res$sample.id <- vcf.info$sample.id
  return(res)
}#End extract_info_vcf


# check_header_source_vcf
#' @title Check the vcf header and detect vcf source
#' @description Check the vcf header and detect vcf source
#' @rdname check_header_source_vcf
#' @keywords internal
#' @export
check_header_source_vcf <- function(vcf) {

  check.header <- SeqArray::seqVCF_Header(vcf)
  problematic.id <- c("AD", "AO", "QA", "GL")
  problematic.id <-
    purrr::keep(
      .x = problematic.id,
      .p = problematic.id %in% check.header$format$ID
      )
  for (p in problematic.id) {
    check.header$format[check.header$format$ID == p, "Number"] <- "."
  }
  # check.header$format

  check.source <- check.header$header$value[check.header$header$id == "source"]
  is.stacks <- stringi::stri_detect_fixed(str = check.source, pattern = "Stacks")
  if (is.stacks) {
    stacks.2 <- keep.stacks.gl <- stringi::stri_detect_fixed(
      str = check.source,
      pattern = "Stacks v2")
    keep.stacks.gl <- TRUE
    if (!keep.stacks.gl) {
      check.header$format <- dplyr::filter(check.header$format, ID != "GL")
    }
    markers.info <- NULL
    overwrite.metadata <- NULL
  } else {
    stacks.2 <- FALSE
    markers.info <- NULL
    overwrite.metadata <- NULL
  }
  return(
    res = list(
      source = stacks.2,
      check.header = check.header,
      markers.info = markers.info,
      overwrite.metadata = overwrite.metadata
    )
  )
}#End check_header_source_vcf


# vcf_strata -------------------------------------------------------------------
#' @name vcf_strata
#' @title Join stratification metadata to a VCF (population-aware VCF)
#' @description Include stratification metadata, e.g. population-level information,
#' to the \code{FORMAT} field of a VCF file.
#' @param data A VCF file

#' @param strata (optional) A tab delimited file at least 2 columns
#' with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} and any other columns can be any hierarchical grouping.
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.

#' @param filename (optional) The file name for the modifed VCF,
#' written to the working directory. Default: \code{filename = NULL} will make a
#' custom filename with data and time.

#' @export
#' @rdname vcf_strata
#' @importFrom purrr detect_index
#' @importFrom readr read_delim write_tsv read_lines
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom tibble as_data_frame add_row
#' @importFrom stringi stri_replace_all_fixed stri_detect_fixed
#' @importFrom utils write.table count.fields
#' @importFrom dplyr left_join mutate rename ungroup group_by
#' @importFrom tidyr unite_ spread

#' @return A VCF file in the working directory with new \code{FORMAT} field(s)
#' correponding to the strata column(s).

#' @seealso
#' \href{https://vcftools.github.io}{VCF web page}
#'
#' \href{VCF specification page}{https://vcftools.github.io/specs.html}

#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


vcf_strata <- function(data, strata, filename = NULL) {
  # data <- "batch_1.vcf"
  # strata <- "strata.sturgeon.12pop.tsv"
  # filename <- NULL
  # data <- "example_vcf2dadi_ferchaud_2015.vcf"
  # strata <- "strata.stickleback.tsv"

  cat("#######################################################################\n")
  cat("######################### radiator: vcf_strata ##########################\n")
  cat("#######################################################################\n")

  # Checking for missing and/or default arguments ******************************
  if (missing(data)) rlang::abort("Input file missing")
  if (missing(strata)) rlang::abort("Strata file missing")

  # import the first 50 lines
  quick.scan <- readr::read_lines(file = data, n_max = 75)

  # Function to detect where CHROM line starts
  detect_vcf_header <- function(x) {
    stringi::stri_detect_fixed(str = x, pattern = "CHROM", negate = FALSE)
  }

  # Detect the index
  max.vcf.header <- purrr::detect_index(.x = quick.scan, .p = detect_vcf_header) - 1

  # import VCF header and add a row containing the new format field
  vcf.header <- readr::read_delim(file = data, n_max = max.vcf.header, col_names = "VCF_HEADER", delim = "\n")

  # import the vcf file, no filters etc.
  message("Importing the VCF file")
  input <- data.table::fread(
    input = data,
    sep = "\t",
    stringsAsFactors = FALSE,
    header = TRUE,
    skip = "CHROM",
    showProgress = TRUE,
    verbose = FALSE
  ) %>%
    tibble::as_data_frame()

  # transform in long format
  input <- data.table::melt.data.table(
    data = data.table::as.data.table(input),
    id.vars = c("#CHROM", "POS", "ID",  "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),
    variable.name = "INDIVIDUALS",
    variable.factor = FALSE,
    value.name = "FORMAT_ID"
  ) %>%
    tibble::as_data_frame() %>%
    dplyr::mutate(
      INDIVIDUALS = stringi::stri_replace_all_fixed(
        str = INDIVIDUALS,
        pattern = c("_", ":"),
        replacement = c("-", "-"),
        vectorize_all = FALSE)
    )

  # population levels and strata  ---------------------------------------------
  message("Importing the strata file")
  if (is.vector(strata)) {
    # message("strata file: yes")
    number.columns.strata <- max(utils::count.fields(strata, sep = "\t"))
    col.types <- stringi::stri_join(rep("c", number.columns.strata), collapse = "")
    strata.df <- readr::read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>%
      dplyr::rename(POP_ID = STRATA)
  } else {
    # message("strata object: yes")
    colnames(strata) <- stringi::stri_replace_all_fixed(
      str = colnames(strata),
      pattern = "STRATA",
      replacement = "POP_ID",
      vectorize_all = FALSE
    )
    strata.df <- strata
  }

  strata.number <- length(strata.df) - 1
  strata.colnames <- purrr::discard(.x = colnames(strata.df), .p = colnames(strata.df) %in% "INDIVIDUALS")

  # Replace unwanted whitespace pattern in the strata
  strata.df <- strata.df %>%
    dplyr::mutate(
      INDIVIDUALS = stringi::stri_replace_all_fixed(
        str = INDIVIDUALS,
        pattern = c("_", ":"),
        replacement = c("-", "-"),
        vectorize_all = FALSE
      ),
      POP_ID = stringi::stri_replace_all_fixed(
        str = POP_ID,
        pattern = " ",
        replacement = "_",
        vectorize_all = FALSE
      )
    )


  # Join strata to input and merge strata column to FORMAT field
  message("Joining the strata to the VCF into new field format...")
  input <- input %>%
    dplyr::left_join(strata.df, by = "INDIVIDUALS") %>%
    tidyr::unite_(
      data = .,
      col = "FORMAT_ID",
      from = c("FORMAT_ID", strata.colnames),
      sep = ":",
      remove = TRUE
    ) %>%
    dplyr::group_by(`#CHROM`, POS, ID,  REF, ALT, QUAL, FILTER, INFO, FORMAT) %>%
    tidyr::spread(data = ., key = INDIVIDUALS, value = FORMAT_ID) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
      FORMAT = stringi::stri_join(
        FORMAT,
        stringi::stri_join(strata.colnames, collapse = ":"),
        sep = ":",
        collapse = NULL
      )
    )

  # Filename ------------------------------------------------------------------
  message("Writing to the working directory...")
  if (is.null(filename)) {
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "", vectorize_all = FALSE)
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    filename <- stringi::stri_join("radiator_vcf_file_", file.date, ".vcf")
  } else {
    filename <- stringi::stri_join(filename, ".vcf")
  }
  # File format ----------------------------------------------------------------
  # write_delim(x = data_frame("##fileformat=VCFv4.2"), path = filename, delim = " ", append = FALSE, col_names = FALSE)
  vcf.header[1,] <- "##fileformat=VCFv4.3"

  # File date ------------------------------------------------------------------
  file.date <- stringi::stri_replace_all_fixed(Sys.Date(), pattern = "-", replacement = "")
  file.date <- stringi::stri_join("##fileDate=", file.date, sep = "")
  # write_delim(x = data_frame(file.date), path = filename, delim = " ", append = TRUE, col_names = FALSE)
  vcf.header[2,] <- file.date

  # Source ---------------------------------------------------------------------
  # write_delim(x = data_frame(stringi::stri_join("##source=radiator_v.", utils::packageVersion("radiator"))), path = filename, delim = " ", append = TRUE, col_names = FALSE)
  # vcf.header[3,] <- stringi::stri_replace_all_fixed(str = vcf.header[3,], pattern = '"', replacement = "", vectorize_all = FALSE)
  # vcf.header[3,]<- stringi::stri_join(vcf.header[3,], "and radiator v.", utils::packageVersion("radiator"))

  # New FORMAT -----------------------------------------------------------------
  for (i in strata.colnames) {
    vcf.header <- tibble::add_row(
      .data = vcf.header,
      VCF_HEADER = stringi::stri_join(
        "##FORMAT=<ID=", i, ',Number=1,Type=Character,Description="New strata",Source="radiator",Version="', utils::packageVersion("radiator"), '">')
    )
  }
  # VCF HEADER  ------------------------------------------------------------------
  utils::write.table(x = vcf.header, file = filename, sep = " ", append = FALSE, col.names = FALSE, quote = FALSE, row.names = FALSE)

  # Write the data   -------------------------------------------------------------
  suppressWarnings(readr::write_tsv(x = input, path = filename, append = TRUE, col_names = TRUE))

  cat("############################## completed ##############################\n")
}#vcf_strata


## Split vcf--------------------------------------------------------------------

#' @title Split a VCF file
#' @description This function allows to split a VCF file in several VCFs,
#' based on individuals or populations.

#' @param strata A file identical to the strata file usually used in radiator,
#' with an additional column named: \code{SPLIT}.
#' This new column contains numerical values
#' (e.g. 1, 1, 1, ..., 2, 2, 2, 2, ..., 3, 3, ...),
#' that indicate for each INDIVIDUALS/STRATA, how to split.
#' The number of VCF to split to is based on the max value found in the column
#' \code{SPLIT}, above this would result in 3 VCF files created).

#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_common_arguments

#' @importFrom stringi stri_replace_all_fixed stri_replace_na stri_join stri_count_fixed
#' @importFrom tibble as_data_frame data_frame add_column add_row
#' @importFrom dplyr select rename n_distinct distinct mutate summarise group_by ungroup arrange left_join full_join semi_join anti_join bind_rows bind_cols if_else
#' @importFrom readr write_tsv read_tsv
#' @importFrom tidyr separate gather
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid stat_smooth
#' @importFrom purrr flatten_chr map_df flatten_dbl

#' @return The function returns in the global environment a list with
#' the different tidy dataset from the split vcf. In the working directory,
#' the splitted VCF files with \code{"_1", "_2"} in the name.

#' @examples
#' \dontrun{
#' split.data <- radiator::split_vcf(
#' data = "batch_1.vcf",
#' strata = "strata.split.tsv",
#' blacklist.id = "blacklisted.id.txt",
#' whitelist.markers = "whitelist.loci.txt")
#' }


#' @rdname split_vcf
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

split_vcf <- function(
  data,
  strata,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  cat("#######################################################################\n")
  cat("########################### radiator::split_vcf #########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # manage missing arguments -----------------------------------------------------
  if (missing(data)) rlang::abort("data/vcf file missing")
  if (missing(strata)) rlang::abort("strata file missing")

  # if (!is.null(pop.levels) & is.null(pop.labels)) {
  #   pop.levels <- stringi::stri_replace_all_fixed(
  #     pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  #   pop.labels <- pop.levels
  # }
  # if (!is.null(pop.labels) & is.null(pop.levels)) rlang::abort("pop.levels is required if you use pop.labels")
  # if (!is.null(pop.labels)) {
  #   if (length(pop.labels) != length(pop.levels)) rlang::abort("pop.labels and pop.levels must have the same length (number of groups)")
  #   pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  # }

  # Filename -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  file.date <- stringi::stri_replace_all_fixed(
    Sys.time(),
    pattern = " EDT",
    replacement = "",
    vectorize_all = FALSE
  )
  file.date <- stringi::stri_replace_all_fixed(
    file.date,
    pattern = c("-", " ", ":"),
    replacement = c("", "@", ""),
    vectorize_all = FALSE
  )
  file.date <- stringi::stri_sub(file.date, from = 1, to = 13)

  filename <- stringi::stri_join("radiator_split_vcf_", file.date)

  # import data ----------------------------------------------------------------
  # if (!is.null(blacklist.id)) {# With blacklist of ID
  #   if (is.vector(blacklist.id)) {
  #     suppressMessages(blacklist <- readr::read_tsv(blacklist.id, col_names = TRUE))
  #   } else {
  #     if (!tibble::has_name(blacklist.id, "INDIVIDUALS")) {
  #       rlang::abort("Blacklist of individuals should have 1 column named: INDIVIDUALS")
  #     }
  #     blacklist <- blacklist.id
  #   }
  #   blacklist$INDIVIDUALS <- stringi::stri_replace_all_fixed(
  #     str = blacklist$INDIVIDUALS,
  #     pattern = c("_", ":"),
  #     replacement = c("-", "-"),
  #     vectorize_all = FALSE
  #   )
  #
  #   # remove potential duplicate id
  #   blacklist <- dplyr::distinct(.data = blacklist, INDIVIDUALS)
  # }


  split <- suppressMessages(readr::read_tsv(file = strata))
  strata <- dplyr::select(split, -SPLIT)
  split <- dplyr::select(split, -STRATA)

  # if (!is.null(blacklist.id)) {
  #   split <- dplyr::anti_join(x = split, y = blacklist, by = "INDIVIDUALS")
  # }

  split$INDIVIDUALS <- stringi::stri_replace_all_fixed(
    str = split$INDIVIDUALS,
    pattern = c("_", ":"),
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )

  # Function required
  split_vcf <- function(data, filename) {
    split.id <- unique(data$SPLIT)
    filename <- stringi::stri_join(filename, "_", split.id)
    radiator::write_vcf(data = dplyr::select(data, -SPLIT),
                        pop.info = FALSE, filename = filename)
  }


  message("Importing and tidying the vcf...")
  input <- suppressMessages(
    radiator::tidy_genomic_data(
      data = data,
      strata = strata,
      vcf.metadata = FALSE,
      whitelist.markers = whitelist.markers,
      filename = NULL,
      verbose = FALSE) %>%
      dplyr::full_join(split, by = "INDIVIDUALS") %>%
      split(x = ., f = .$SPLIT)
  )

  split <- strata <- blacklist <- NULL

  .radiator_parallel(
    X = input, FUN = split_vcf, mc.cores = parallel.core, filename = filename)


  # results --------------------------------------------------------------------
  message("Split VCFs were written in the working directory")
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(input)
}#End split_vcf




## Merge vcf

#' @title Merge VCF files
#' @description This function allows to merge 2 VCF files.

#' @param vcf1 First VCF file.
#' @param strata1 strata file for vcf1.
#' @param vcf2 Second VCF file.
#' @param strata2 strata file for vcf2.
#' @param filename Name of the merged VCF file.
#' With the default, the function gives a filename based on date and time.
#' Default: \code{filename = NULL}.
#' @inheritParams tidy_genomic_data

#' @importFrom stringi stri_replace_all_fixed stri_replace_na stri_join stri_count_fixed
#' @importFrom tibble as_data_frame data_frame add_column add_row
#' @importFrom dplyr select rename n_distinct distinct mutate summarise group_by ungroup arrange left_join full_join semi_join anti_join bind_rows bind_cols if_else
#' @importFrom readr write_tsv read_tsv
#' @importFrom tidyr separate gather
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid stat_smooth
#' @importFrom purrr flatten_chr map_df flatten_dbl

#' @return The function returns in the global environment a tidy dataset with
#' the merged VCF files and the merged VCF in the working directory.

#' @examples
#' \dontrun{
#' # The simplest way to run the function:
#' sum <- radiator::merge_vcf(
#' vcf1 = "batch_1.vcf", strata1 = "strata1_brook_charr.tsv",
#' vcf1 = "batch_2.vcf", strata2 = "strata2_brook_charr.tsv",
#' pop.select = c("QC", "ON", "NE"),
#' maf.thresholds = c(0.002, 0.001),
#' maf.pop.num.threshold = 1,
#' maf.approach = "SNP",maf.operator = "OR",
#' filename = "my_new_VCF.vcf"
#' }

# merge_vcf----------------------------------------------------------------
# @rdname merge_vcf
# @export
# @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# merge_vcf <- function(
#   vcf1, strata1,
#   vcf2, strata2,
#   whitelist.markers = NULL,
#   filename = NULL,
#   parallel.core = parallel::detectCores() - 1
# ) {
#   cat("#######################################################################\n")
#   cat("########################### radiator::merge_vcf #########################\n")
#   cat("#######################################################################\n")
#   timing <- proc.time()
#
#   # manage missing arguments -----------------------------------------------------
#   if (missing(vcf1)) rlang::abort("vcf1 file missing")
#   if (missing(vcf2)) rlang::abort("vcf2 file missing")
#   if (missing(strata1)) rlang::abort("strata1 file missing")
#   if (missing(strata2)) rlang::abort("strata2 file missing")
#
#
#   # if (!is.null(pop.levels) & is.null(pop.labels)) {
#   #   pop.levels <- stringi::stri_replace_all_fixed(
#   #     pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
#   #   pop.labels <- pop.levels
#   # }
#   # if (!is.null(pop.labels) & is.null(pop.levels)) rlang::abort("pop.levels is required if you use pop.labels")
#   # if (!is.null(pop.labels)) {
#   #   if (length(pop.labels) != length(pop.levels)) rlang::abort("pop.labels and pop.levels must have the same length (number of groups)")
#   #   pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
#   # }
#
#   # Filename -------------------------------------------------------------------
#   # Get date and time to have unique filenaming
#   if (is.null(filename)) {
#     file.date <- format(Sys.time(), "%Y%m%d@%H%M")
#   }
#
#   # import data ----------------------------------------------------------------
#   message("Importing and tidying the vcf1...")
#   input <- suppressMessages(radiator::tidy_genomic_data(
#     data = vcf1,
#     strata = strata1,
#     vcf.metadata = FALSE,
#     whitelist.markers = whitelist.markers,
#     filename = NULL,
#     verbose = FALSE))
#
#   message("Importing and tidying the vcf2...")
#   # Also Using pop.levels and pop.labels info if present
#   input <- suppressWarnings(
#     dplyr::bind_rows(
#       input,
#       suppressMessages(
#         radiator::tidy_genomic_data(
#           data = vcf2,
#           strata = strata2,
#           vcf.metadata = FALSE,
#           whitelist.markers = whitelist.markers,
#           filename = NULL,
#           verbose = FALSE))))
#
#   message("Adjusting REF/ALT alleles...")
#   input <- radiator::calibrate_alleles(
#     data = input,
#     parallel.core = parallel.core,
#     verbose = TRUE)$input
#
#   # if (filter.monomorphic) {
#   #   input <- radiator::filter_monomorphic(data = input, verbose = TRUE)
#   # }
#   #
#   # if (filter.common.markers) {
#   #   input <- radiator::filter_common_markers(data = input, verbose = TRUE)$input
#   # }
#   #
#   # if (!is.null(filter.mac)) {
#   #   input <- filter_maf(
#   #     data = input,
#   #     interactive.filter = FALSE,
#   #     filter.mac = filter.mac,
#   #     parallel.core = parallel.core,
#   #     verbose = FALSE)$tidy.filtered.mac
#   # }
#
#   # Write VCF in the working directory
#   radiator::write_vcf(data = input, pop.info = FALSE, filename = filename)
#
#   # results --------------------------------------------------------------------
#   message("Merged VCF in the working directory: ", filename, ".vcf")
#   timing <- proc.time() - timing
#   message("\nComputation time: ", round(timing[[3]]), " sec")
#   cat("############################## completed ##############################\n")
#   return(input)
# }#End merge_vcf
