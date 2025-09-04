# Detect detect_genomic_format

#' @name detect_genomic_format
#' @title Used internally in radiator to detect the file format
#' @description Detect the file format of genomic data set.
#' @param data 15 options for input: VCFs (SNPs or Haplotypes,
#' to make the vcf population ready),
#' plink (tped and bed), stacks haplotype file, genind (library(adegenet)),
#' genlight (library(adegenet)), gtypes (library(strataG)), genepop, DArT,
#' and a data frame in long/tidy or wide format. To verify that radiator detect
#' your file format use \code{\link{detect_genomic_format}} (see example below).
#' New addition to radiator: the Apache Parquet columnar storage file format that
#' will replace the fst (Lightning Fast Serialiation) format.
#' Documented in \strong{Input genomic datasets} of \code{\link{tidy_genomic_data}}.

#' @param guess (character) In development, guess faster the type of file.
#' Default: \code{guess = NULL}.

#' @return One of these file format:
#' \itemize{
#' \item tbl_df: for a data frame
#' \item genind: for a genind object
#' \item genlight: for a genlight object
#' \item gtypes: for a gtypes object
#' \item vcf.file: for a vcf file
#' \item plink.tped.file: for a plink tped file
#' \item plink.bed.file: for a plink bed file
#' \item genepop.file: for a genepop file
#' \item haplo.file: for a stacks haplotypes file
#' \item fstat.file: for a fstat file
#' \item dart: for a DArT file
#' \item fst.file: for a file ending with .rad
#' \item SeqVarGDSClass: for SeqArray GDS file.
#' \item arrow parquet: for a Apache Parquet columnar storage file format,
#' used by arrow R package.
#' }

#' @rdname detect_genomic_format
# @keywords internal
#' @export
#' @examples
#' \dontrun{
#' #To verify your file is detected by radiator as the correct format:
#' radiator::detect_genomic_format(data = "populations.snps.vcf")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_genomic_format <- function(data, guess = NULL){
  data.type <- NULL # in case we don't find the format

  if (is.null(guess)) {
    if (!is.vector(data)) {
      if (tibble::has_name(data, "INDIVIDUALS")) {
        data.type <- "tbl_df" #"df.file"
      } else {
        data.type <- class(data)[1]
        if (!data.type %in% c("genind", "genlight", "gtypes", "SeqVarGDSClass")) rlang::abort("Input file not recognised")
      }
      return(data.type)
    } # not vector
    if (is.vector(data)) {
      data.type <- suppressWarnings(readLines(con = data, n = 1L))
      file.ending <- stringi::stri_sub(str = data, from = -4, to = -1)

      # Specifics with gz file ...
      if (stringi::stri_detect_fixed(str = file.ending, pattern = "gz")) {
        file.ending <- stringi::stri_sub(str = data, from = -7, to = -4)
      }

      # VCF
      # if (identical(data.type, "##fileformat=VCF") || file.ending == ".vcf") {
      if (stringi::stri_detect_fixed(str = data.type, pattern = "##fileformat=VCF") || file.ending == ".vcf") {
        data.type <- "vcf.file"
        # message("File type: VCF")
        return(data.type)
      } #End VCF


      # PLINK
      if (file.ending == "tped") {
        data.type <- "plink.tped.file"
        # message("File type: PLINK tped")
        if (!file.exists(stringi::stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
          rlang::abort("Missing tfam file with the same prefix as your tped")
        }
        return(data.type)
      }#End plink

      if (file.ending == ".bed") {
        data.type <- "plink.bed.file"

        # BED file requires bim and fam files...
        bim.file <- file.exists(stringi::stri_replace_all_fixed(str = data, pattern = ".bed", replacement = ".bim", vectorize_all = FALSE))
        fam.file <- file.exists(stringi::stri_replace_all_fixed(str = data, pattern = ".bed", replacement = ".fam", vectorize_all = FALSE))
        if (FALSE %in% c(bim.file, fam.file)) {
          rlang::abort("Missing fam or bim file(s) with the same prefix as your PLINK bed file")
        }
        return(data.type)
      }#end plink


      # TIBBLE
      if (stringi::stri_detect_fixed(str = data.type, pattern = "POP_ID") | stringi::stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") | stringi::stri_detect_fixed(str = data.type, pattern = "MARKERS") | stringi::stri_detect_fixed(str = data.type, pattern = "LOCUS")) {
        data.type <- "tbl_df" #"df.file"
        # message("File type: data frame of genotypes")
        return(data.type)
      }
      if (stringi::stri_detect_fixed(str = data.type, pattern = "Catalog")) {
        data.type <- "haplo.file"
        return(data.type)
      }
      if (file.ending == ".gen") {
        data.type <- "genepop.file"
        return(data.type)
      }
      if (file.ending == ".dat") {
        data.type <- "fstat.file"
        return(data.type)
      }

      # .rad, fst file or gds
      if (TRUE %in% (c(".rad", ".gds") %in% file.ending)) {
        if (stringi::stri_detect_fixed(str = data.type, pattern = "COREARRAY")) {
          data.type <- "gds.file"
        } else {
          data.type <- "fst.file"
          message("This file format will be deprecated soon:")
          message("    Save to another format (e.g. parquet) or ")
          message("    Use the package fst to open and save in the future")
        }
        return(data.type)
      }

      # Arrow Parquet
      if (file.ending == "quet") {
        data.type <- "arrow.parquet"
        return(data.type)
      }


      # DArT data
      data.type <- detect_dart_format(data = data, verbose = FALSE) %$% data.type
      if ("dart" %in% data.type) {
      return(data.type)
      }


      if (file.ending == ".rds") {
        temp <- readRDS(file = "rat8rec_reassign.rds")
        data.type <- class(temp)[1]
        temp <- NULL
        if (!data.type %in% c("genind", "genlight", "gtypes")) {
          rlang::abort("Input file not recognised")
        }
        return(data.type)
      }
    } # end vector detection
  } #End no guess

  if (guess == "dart") {
    data.type <- detect_dart_format(data = data, verbose = FALSE)
    data.type <- data.type$data.type
    return(data.type)
  } # guess DART
  # End guess
  if (is.null(data.type)) rlang::abort("Genomic format not detected: send email to developer")
  return(data.type)
} # End detect_genomic_format

# INTERNAL functions -----------------------------------------------------------
# detect_dart_format-------------------------------------------------------------
#' @title detect_dart_format
#' @description Detect the dart genotype format: 1row, 2rows or counts
#' @rdname detect_dart_format
#' @keywords internal
#' @export
detect_dart_format <- function(data, verbose = TRUE) {
  # dart.format:
  # silico.dart
  # 1row
  # 2rows
  # counts

  # Detect how to read the file
  if (stringi::stri_detect_fixed(
    str = stringi::stri_sub(str = data, from = -4, to = -1),
    pattern = ".csv")) {
    tokenizer.dart <- ","
  } else {
    tokenizer.dart <- "\t"
  }

  # read the 1st line
  data.type <- suppressWarnings(readLines(con = data, n = 1L))

  # presence of * in the headers
  dart.headers <- TRUE %in% (stringi::stri_detect_fixed(str = data.type, pattern = c("*\t", "*,")))
  n.col <- stringi::stri_count_fixed(str = data.type, pattern = tokenizer.dart) + 1

  # deconstruct
  skip.number <- 0
  star.number <- 0

  if (dart.headers) {
    temp.file <- suppressWarnings(suppressMessages(readr::read_table(file = data, n_max = 20, col_names = "HEADER")))

    # skip.number <- which(
    #   stringi::stri_detect_regex(str = temp.file$HEADER,pattern = "^[:Letter:]")
    # ) - 1
    skip.number <- which(
      stringi::stri_detect_regex(str = temp.file$HEADER,pattern = "^[[:alnum:]]")
    )[1] - 1

    star.number <- stringi::stri_count_fixed(str = data.type, pattern = "*")
    data.type <- readr::read_lines(file = data, skip = skip.number, n_max = skip.number + 1)[1]
  }

  dart.clone.id <- stringi::stri_detect_fixed(str = data.type, pattern = "CloneID")
  dart.allele.id <- stringi::stri_detect_fixed(str = data.type, pattern = "AlleleID")
  dart.snp.position <- stringi::stri_detect_fixed(str = data.type, pattern = "SnpPosition")
  dart.allele.sequece <- stringi::stri_detect_fixed(str = data.type, pattern = "AlleleSequence")

  # DArT or no
  if (dart.snp.position || dart.allele.id || dart.clone.id || dart.allele.sequece) {
    data.type <- "dart"
  } else {
    return(res = list(data.type = NULL))
  }

  # What DArT format
  if (dart.clone.id && !dart.allele.id) {
    dart.format <- "silico.dart"
    if (verbose) message("DArT format: silico DArT")
  } else {
    temp.line <- readr::read_lines(file = data, skip = skip.number, n_max = 1)
    dart.col <- c(
      if (stringi::stri_detect_regex(str = temp.line, pattern = "SNP")) "SNP",
      if (stringi::stri_detect_regex(str = temp.line, pattern = "REF")) "REF",
      if (stringi::stri_detect_regex(str = temp.line, pattern = "Variant")) "Variant"
    )

    binary <- vroom::vroom(
      file = data,
      delim = tokenizer.dart,
      # col_select = tidyselect::all_of(1:star.number),
      col_selec = tidyselect::all_of(dart.col),
      skip = skip.number,
      na = c("-", "NA",""),
      guess_max = 20,
      altrep = TRUE,
      skip_empty_rows = TRUE,
      trim_ws = FALSE,
      progress = FALSE,
      show_col_types = FALSE,
      .name_repair = "minimal"
    ) %>%
     anyNA(.)

    if (!binary) {
      if (verbose) message("DArT format: genotypes in 1 Row")
      dart.format <- "1row"
    } else {
      count.data <- vroom::vroom(
        file = data,
        delim = tokenizer.dart,
        col_select = tidyselect::last_col(3):tidyselect::last_col(),
        skip = skip.number,
        n_max = 50,
        na = c("-", "NA",""),
        guess_max = 20,
        altrep = TRUE,
        skip_empty_rows = TRUE,
        trim_ws = FALSE,
        progress = FALSE,
        show_col_types = FALSE,
        .name_repair = "minimal"
      )
      count.data <- any(count.data > 1, na.rm = TRUE)

      if (count.data) {
        dart.format <- "counts"
        if (verbose) message("DArT format: alleles coverage in 2 Rows counts")
      } else {
        dart.format <- "2rows"
        if (verbose) message("DArT format: alleles absence/presence in 2 Rows")
      }
    }
  }

  # New 20250904
  # denovo or ref genome used

  # default
  ref.genome <- FALSE
  ref.genome <- any(stringi::stri_detect_regex(str = temp.file$HEADER, pattern = "scaffold"))

  return(
    res = list(
      data.type = data.type,
      dart.format = dart.format,
      ref.genome = ref.genome,
      skip.number = skip.number,
      star.number = star.number,
      tokenizer.dart = tokenizer.dart,
      n.col = n.col
    )
  )
}#End detect_dart_format
