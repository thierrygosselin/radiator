## Split vcf

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
  monomorphic.out = FALSE,
  common.markers = FALSE,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  blacklist.id = NULL,
  whitelist.markers = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  cat("#######################################################################\n")
  cat("########################### radiator::split_vcf #########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()

  # manage missing arguments -----------------------------------------------------
  if (missing(data)) stop("data/vcf file missing")
  if (missing(strata)) stop("strata file missing")

  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(
      pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }

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
  if (!is.null(blacklist.id)) {# With blacklist of ID
    if (is.vector(blacklist.id)) {
      suppressMessages(blacklist <- readr::read_tsv(blacklist.id, col_names = TRUE))
    } else {
      if (!tibble::has_name(blacklist.id, "INDIVIDUALS")) {
        stop("Blacklist of individuals should have 1 column named: INDIVIDUALS")
      }
      blacklist <- blacklist.id
    }
    blacklist$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = blacklist$INDIVIDUALS,
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )

    # remove potential duplicate id
    blacklist <- dplyr::distinct(.data = blacklist, INDIVIDUALS)
  }


  split <- suppressMessages(readr::read_tsv(file = strata))
  strata <- dplyr::select(split, -SPLIT)
  split <- dplyr::select(split, -STRATA)

  if (!is.null(blacklist.id)) {
    split <- dplyr::anti_join(x = split, y = blacklist, by = "INDIVIDUALS")
  }

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
      blacklist.id = blacklist.id,
      whitelist.markers = whitelist.markers,
      monomorphic.out = monomorphic.out,
      common.markers = common.markers,
      pop.levels = pop.levels,
      pop.labels = pop.labels,
      pop.select = pop.select,
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
