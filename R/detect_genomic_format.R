# Detect detect_genomic_format

#' @name detect_genomic_format
#' @title Used internally in radiator to detect the file format
#' @description Detect file format of genomic data set.
#' @param data 14 options for input: VCFs (SNPs or Haplotypes,
#' to make the vcf population ready),
#' plink (tped and bed), stacks haplotype file, genind (library(adegenet)),
#' genlight (library(adegenet)), gtypes (library(strataG)), genepop, DArT,
#' and a data frame in long/tidy or wide format. To verify that radiator detect
#' your file format use \code{\link{detect_genomic_format}} (see example below).
#' Documented in \strong{Input genomic datasets} of \code{\link{tidy_genomic_data}}.

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

detect_genomic_format <- function(data){

  if (!is.vector(data)) {
    if (tibble::has_name(data, "INDIVIDUALS")) {
      data.type <- "tbl_df" #"df.file"
    } else {
      data.type <- class(data)[1]
      if (!data.type %in% c("genind", "genlight", "gtypes", "SeqVarGDSClass")) rlang::abort("Input file not recognised")
    }
  } else {
    data.type <- suppressWarnings(readLines(con = data, n = 1L))
    # data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
    file.ending <- stringi::stri_sub(str = data, from = -4, to = -1)

    # problem with gz file ...
    if (stringi::stri_detect_fixed(str = file.ending, pattern = "gz")) {
      file.ending <- stringi::stri_sub(str = data, from = -7, to = -4)
    }

    # if (identical(data.type, "##fileformat=VCF") || file.ending == ".vcf") {
    if (stringi::stri_detect_fixed(str = data.type, pattern = "##fileformat=VCF") || file.ending == ".vcf") {
      data.type <- "vcf.file"
      # message("File type: VCF")
    }

    if (file.ending == "tped") {
      data.type <- "plink.tped.file"
      # message("File type: PLINK tped")
      if (!file.exists(stringi::stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
        rlang::abort("Missing tfam file with the same prefix as your tped")
      }
      return(data.type)
    }


    if (file.ending == ".bed") {
      data.type <- "plink.bed.file"
      # BED file requires bim and fam files...
      bim.file <- file.exists(stringi::stri_replace_all_fixed(str = data, pattern = ".bed", replacement = ".bim", vectorize_all = FALSE))
      fam.file <- file.exists(stringi::stri_replace_all_fixed(str = data, pattern = ".bed", replacement = ".fam", vectorize_all = FALSE))
      if (FALSE %in% c(bim.file, fam.file)) {
        rlang::abort("Missing fam or bim file(s) with the same prefix as your PLINK bed file")
      }
      return(data.type)
    }

    if (stringi::stri_detect_fixed(str = data.type, pattern = "POP_ID") | stringi::stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") | stringi::stri_detect_fixed(str = data.type, pattern = "MARKERS") | stringi::stri_detect_fixed(str = data.type, pattern = "LOCUS")) {
      data.type <- "tbl_df" #"df.file"
      # message("File type: data frame of genotypes")
      return(data.type)
    }
    if (stringi::stri_detect_fixed(str = data.type, pattern = "Catalog")) {
      data.type <- "haplo.file"
    }
    if (file.ending == ".gen") {
      data.type <- "genepop.file"
    }

    if (file.ending == ".dat") {
      data.type <- "fstat.file"
    }

    # .rad, fst file or gds
    if (TRUE %in% (c(".rad", ".gds") %in% file.ending)) {
      if (stringi::stri_detect_fixed(str = data.type, pattern = "COREARRAY")) {
        data.type <- "gds.file"
      } else {
        data.type <- "fst.file"
      }
    }

    # DArT data
    dart.temp <- check_dart(data)
    if (dart.temp$data.type %in% c("dart", "silico.dart")) {
      data.type <- dart.temp$data.type
    }

    if (file.ending == ".rds") {
      temp <- readRDS(file = "rat8rec_reassign.rds")
      data.type <- class(temp)[1]
      temp <- NULL
      if (!data.type %in% c("genind", "genlight", "gtypes")) {
        rlang::abort("Input file not recognised")
      }
    }
  } # end file type detection
  return(data.type)
} # End detect_genomic_format

# INTERNAL functions -----------------------------------------------------------
#' @title check_dart
#' @description check dart
#' @rdname check_dart
#' @export
#' @keywords internal
check_dart <- function(data) {
  if (stringi::stri_detect_fixed(
    str = stringi::stri_sub(str = data, from = -4, to = -1),
    pattern = ".csv")) {
    tokenizer.dart <- ","
  } else {
    tokenizer.dart <- "\t"
  }
  # data.type <- readChar(con = data, nchars = 200L, useBytes = TRUE)
  data.type <- suppressWarnings(readLines(con = data, n = 1L))

  dart.with.header <- TRUE %in% (stringi::stri_detect_fixed(str = data.type, pattern = c("*\t", "*,")))

  if (dart.with.header) {
    temp.file <- suppressWarnings(suppressMessages(readr::read_table(file = data, n_max = 20, col_names = "HEADER")))
    skip.number <- which(
      stringi::stri_detect_regex(str = temp.file$HEADER,pattern = "^[:Letter:]")
    ) -1
    star.number <- stringi::stri_count_fixed(str = data.type, pattern = "*")
    data.type <- readr::read_lines(file = data, skip = skip.number, n_max = skip.number + 1)[1]
  } else {
    skip.number <- 0
    star.number <- 0
  }
  dart.clone.id <- stringi::stri_detect_fixed(str = data.type, pattern = "CloneID")
  dart.allele.id <- stringi::stri_detect_fixed(str = data.type, pattern = "AlleleID")
  dart.snp.position <- stringi::stri_detect_fixed(str = data.type, pattern = "SnpPosition")

  if (dart.snp.position || dart.allele.id) {
    data.type <- "dart"
  }

  if (dart.clone.id && !dart.allele.id) {
    data.type <- "silico.dart"
  }
  # if (!data.type %in% c("dart", "silico.dart")) {
  #   rlang::abort("Contact author to show your DArT data, problem during import")
  # }
  return(
    res = list(
      data.type = data.type,
      skip.number = skip.number,
      star.number = star.number,
      tokenizer.dart = tokenizer.dart
    )
  )
}#End check_dart
