# Detect detect_genomic_format

#' @name detect_genomic_format
#' @title Used internally in radiator to detect the file format
#' @description Detect file format of genomic data set.
#' @param data A genomic data set in the global environment

#' @return One of these file format:
#' \itemize{
#' \item tbl_df: for a data frame
#' \item genind: for a genind object
#' \item genlight: for a genlight object
#' \item gtypes: for a gtypes object
#' \item vcf.file: for a vcf file
#' \item plink.file: for a plink file
#' \item genepop.file: for a genepop file
#' \item haplo.file: for a stacks haplotypes file
#' \item fstat.file: for a fstat file
#' \item dart: for a DArT file
#' \item fst.file: for a file ending with .rad
#' \item SeqVarGDSClass: for SeqArray GDS file.
#' }

#' @rdname detect_genomic_format
#' @importFrom stringi stri_detect_fixed stri_replace_all_fixed
#' @importFrom tibble has_name
# @keywords internal
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_genomic_format <- function(data){

  if (!is.vector(data)) {
    if (tibble::has_name(data, "POP_ID") && tibble::has_name(data, "INDIVIDUALS")) {
      data.type <- "tbl_df" #"df.file"
    } else {
      data.type <- class(data)[1]
      if (!data.type %in% c("genind", "genlight", "gtypes", "SeqVarGDSClass")) stop("Input file not recognised")
    }
  } else {
    data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
    file.ending <- stringi::stri_sub(str = data, from = -4, to = -1)

    if (identical(data.type, "##fileformat=VCF") || file.ending == ".vcf") {
      data.type <- "vcf.file"
      # message("File type: VCF")
    }

    if (file.ending == ".tped") {
      data.type <- "plink.file"
      # message("File type: PLINK")
      if (!file.exists(stringi::stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
        stop("Missing tfam file with the same prefix as your tped")
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
  data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
  dart.with.header <- TRUE %in% (stringi::stri_detect_fixed(str = data.type, pattern = c("*\t", "*,")))
  if (dart.with.header) {
    temp.file <- suppressWarnings(suppressMessages(readr::read_table(file = data, n_max = 20, col_names = "HEADER")))
    # skip.number <- which(stringi::stri_detect_fixed(str = temp.file$HEADER,
    #                                                 pattern = "AlleleID")) -1
    skip.number <- which(stringi::stri_detect_regex(str = temp.file$HEADER,
                                                    pattern = "^[:Letter:]")) -1
    data.type <- readr::read_lines(file = data, skip = skip.number, n_max = skip.number + 1)[1] #%>% stringi::stri_sub(str = ., from = 1, to = 16)
  } else {
    skip.number <- 0
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
  return(res = list(data.type = data.type, skip.number = skip.number))
}#End check_dart
