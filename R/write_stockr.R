# write a stockR data set
#' @name write_stockr
#' @title Write a stockR dataset from a tidy data frame or GDS file or object.

#' @description Write a stockR dataset (Fost et al. submitted).
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @inheritParams radiator_common_arguments
#' @param filename (optional) The stockr object is written in the working
#' directory. The file is written with \code{radiator_stockr_DATE@TIME.RData} and
#' can be open with readRDS.
#' Default: \code{filename = NULL}.

#' @export
#' @rdname write_stockr
#' @references Foster et al. submitted

#' @return The object generated is a matrix with
#' dimension: MARKERS x INDIVIDUALS. The genotypes are coded like PLINK:
#' 0, 1 or 2 alternate allele. 0: homozygote for the reference allele,
#' 1: heterozygote, 2: homozygote for the alternate allele.
#' Missing genotypes have NA. The object also as 2 attributes.
#' \code{attributes(data)$grps} with \code{STRATA/POP_ID} of the individuals and
#' \code{attributes(data)$sample.grps} filled with \code{INDIVIDUALS}.
#' Both attributes can be used inside \emph{stockR}.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

write_stockr <- function(data, filename = NULL, verbose = TRUE) {
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = verbose)
    }
    data <- gds2tidy(gds = data, parallel.core = parallel::detectCores() - 1)
    data.type <- "tbl_df"
  } else {
    if (is.vector(data)) {
      data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
    }
  }

  if (!rlang::has_name(data, "POP_ID") && rlang::has_name(data, "STRATA")) {
    data %<>% dplyr::rename(POP_ID = STRATA)
  }


  data <- dplyr::select(data, MARKERS, POP_ID, INDIVIDUALS, GT_BIN) %>%
    dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
  strata <- dplyr::distinct(data, INDIVIDUALS, POP_ID)

  data <- suppressWarnings(
    dplyr::select(data, MARKERS, INDIVIDUALS, GT_BIN) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = MARKERS ~ INDIVIDUALS,
        value.var = "GT_BIN"
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::select(MARKERS, strata$INDIVIDUALS) %>%
      tibble::remove_rownames(.) %>%
      tibble::column_to_rownames(.data = ., var = "MARKERS")) %>%
    as.matrix(.)

  attr(data,"grps") <- strata$POP_ID
  attr(data,"sample.grps") <- factor(strata$INDIVIDUALS)


  if (is.null(filename)) {
    filename.temp <- generate_filename(extension = "stockr")
  } else {
    filename.temp <- generate_filename(name.shortcut = filename, extension = "stockr")
  }
  filename.short <- filename.temp$filename.short
  filename <- filename.temp$filename
  saveRDS(object = data, file = filename)
  if (verbose) message("File written: ", filename.short)
  return(data)
} # end write_stockr
