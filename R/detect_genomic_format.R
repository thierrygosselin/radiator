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
#' }

#' @rdname detect_genomic_format
#' @importFrom stringi stri_detect_fixed stri_replace_all_fixed
#' @importFrom tibble has_name
# @keywords internal
#' @export
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

detect_genomic_format <- function(data){

  if (!is.vector(data)) {
    if (tibble::has_name(data, "POP_ID") & tibble::has_name(data, "INDIVIDUALS") & tibble::has_name(data, "MARKERS")) {
      data.type <- "tbl_df" #"df.file"
    } else {
      data.type <- class(data)[1]

      if (!data.type %in% c("genind", "genlight", "gtypes")) stop("Input file not recognised")

      # old code
      # if (adegenet::is.genind(data)) {
      #   data.type <- "genind.file"
      #   # message("File type: genind object")
      # } else if (class(data)[1] == "gtypes") {
      #   data.type <- "gtypes"
      #   } else {
      #   stop("Input file not recognised")
      # }
    }
  } else {
    data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)

    if (identical(data.type, "##fileformat=VCF") | stringi::stri_detect_fixed(str = data, pattern = ".vcf")) {
      data.type <- "vcf.file"
      # message("File type: VCF")
    }

    if (stringi::stri_detect_fixed(str = data, pattern = ".tped")) {
      data.type <- "plink.file"
      # message("File type: PLINK")
      if (!file.exists(stringi::stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
        stop("Missing tfam file with the same prefix as your tped")
      }
    }

    if (stringi::stri_detect_fixed(str = data.type, pattern = "POP_ID") | stringi::stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") | stringi::stri_detect_fixed(str = data.type, pattern = "MARKERS") | stringi::stri_detect_fixed(str = data.type, pattern = "LOCUS")) {
      data.type <- "tbl_df" #"df.file"
      # message("File type: data frame of genotypes")
    }
    if (stringi::stri_detect_fixed(str = data.type, pattern = "Catalog")) {
      data.type <- "haplo.file"
    }
    if (stringi::stri_detect_fixed(
      str = stringi::stri_sub(str = data, from = -4, to = -1),
      pattern = ".gen")) {
      data.type <- "genepop.file"
    }
    if (stringi::stri_detect_fixed(
      str = stringi::stri_sub(str = data, from = -4, to = -1),
      pattern = ".dat")) {
      data.type <- "fstat.file"
    }

    # fst file
    if (stringi::stri_detect_fixed(
      str = stringi::stri_sub(str = data, from = -4, to = -1),
      pattern = ".rad")) {
      data.type <- "fst.file"
    }

    # DArT data
    dart.with.header <- stringi::stri_detect_fixed(str = data.type, pattern = "*\t")
    if (dart.with.header) {
      temp.file <- suppressWarnings(suppressMessages(readr::read_table(file = data, n_max = 20, col_names = "HEADER")))
      skip.number <- which(stringi::stri_detect_fixed(str = temp.file$HEADER,
                                                      pattern = "AlleleID"))
      data.type <- readr::read_lines(file = data, skip = 5, n_max = 6)[1] %>%
        stringi::stri_sub(str = ., from = 1, to = 16)
    }
    dart.clone.id <- stringi::stri_detect_fixed(str = data.type, pattern = "CloneID")
    dart.allele.id <- stringi::stri_detect_fixed(str = data.type, pattern = "AlleleID")

    if (dart.clone.id || dart.allele.id) {
      data.type <- "dart"
    }
  } # end file type detection
  return(data.type)
} # End detect_genomic_format
