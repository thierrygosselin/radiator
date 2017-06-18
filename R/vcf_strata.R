# vcf_strata
#' @name vcf_strata
#' @title Join stratification metadata to a VCF (population-aware VCF)
#' @description Include stratification metadata, e.g. population-level information,
#' to the \code{FORMAT} field of a VCF file.
#' @param data A VCF file

#' @param strata (optional for data frame and PLINK files, 
#' required for VCF and haplotypes files) A tab delimited file at least 2 columns 
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
  if (missing(data)) stop("Input file missing")
  if (missing(strata)) stop("Strata file missing")
  
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
}
