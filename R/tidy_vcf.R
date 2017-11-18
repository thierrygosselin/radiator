#' @name tidy_vcf

#' @title Tidy a vcf file (bi/multi-allelic).

#' @description Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.
#' Highly recommended to use \code{\link[radiator]{tidy_genomic_data}} or
#' \code{\link[radiator]{genomic_converter}}. Those two functions allows
#' to manipulate and prune the dataset with blacklists and whitelists along
#' several other filtering options.

#' @inheritParams tidy_genomic_data

#' @param ... (optional) To pass further argument for fine-tuning the tidying
#' (details below).

#' @export
#' @rdname tidy_vcf
#' @importFrom rlang UQ



#' @return The output in your global environment is a tidy data frame.

#' @details
#' \strong{... :dot dot dot arguments}
#' 5 arguments are available:
#' \enumerate{
#' \item whitelist.markers
#' \item blacklist.id
#' \item pop.select
#' \item pop.levels
#' \item pop.labels
#' }
#' Documentation for these arguments is detailed
#' in \code{\link[radiator]{tidy_genomic_data}}, the only difference here, the
#' arguments are objects are in the global environment, not as file in the directory.
#' If files are essential in your pipeline, please use \code{\link[radiator]{tidy_genomic_data}}.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_vcf <- function(
  data,
  strata,
  vcf.metadata = FALSE,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...) {

  if (is.null(strata)) stop("strata argument is required")

  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)
  want <- c("whitelist.markers", "blacklist.id", "pop.select", "pop.levels", "pop.labels")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  vcf.dots <- dotslist[names(dotslist) %in% want]
  whitelist.markers <- vcf.dots[["whitelist.markers"]]
  blacklist.id <- vcf.dots[["blacklist.id"]]
  pop.select <- vcf.dots[["pop.select"]]
  pop.levels <- vcf.dots[["pop.levels"]]
  pop.labels <- vcf.dots[["pop.labels"]]

  # managing vcf.metadata and what approach to use to import faster
  if (is.logical(vcf.metadata) && !vcf.metadata) {
    import.pegas <- TRUE # very fast but no metadata
  } else {
    import.pegas <- FALSE #vcfR as metadata but much slower
  }

  # detect stacks (to manage ID that changed purposes over versions)
  stacks.vcf <- readr::read_lines(file = data, skip = 2, n_max = 1) %>%
    stringi::stri_detect_fixed(str = ., pattern = "Stacks")

  # import vcf with pegas (fastest, but only GT no metadata)
  if (import.pegas) {
    # change names of columns and CHROM column modif
    input <- pegas::VCFloci(file = data, quiet = verbose) %>%
      dplyr::select(-QUAL, -INFO, -FORMAT) %>%
      dplyr::rename(LOCUS = ID) %>%
      dplyr::mutate(
        CHROM = stringi::stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1")
      ) %>%
      dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = as.character) %>%
      tibble::rownames_to_column(df = ., var = "KEEP") %>%
      dplyr::mutate(KEEP = as.integer(KEEP))
  } else {# import with vcfR
    suppressMessages(vcf.data <- vcfR::read.vcfR(file = data, verbose = FALSE))
    input <- tibble::as_data_frame(vcf.data@fix) %>%
      dplyr::select(-QUAL, -INFO) %>%
      dplyr::rename(LOCUS = ID) %>%
      dplyr::mutate(
        CHROM = stringi::stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1")
      ) %>%
      tibble::rownames_to_column(df = ., var = "KEEP") %>%
      dplyr::mutate(KEEP = as.integer(KEEP))
  }

  # bi- or multi-alllelic VCF
  alt.num <- max(unique(
    stringi::stri_count_fixed(str = unique(input$ALT), pattern = ","))) + 1

  if (alt.num > 1) {
    biallelic <- FALSE
    message("VCF is multi-allelic")
  } else {
    biallelic <- TRUE
    message("VCF is biallelic")
  }


  # Check for duplicate indentifiers
  duplicate.markers <- nrow(input) - nrow(dplyr::distinct(input, CHROM, LOCUS, POS))

  if (duplicate.markers > 0) {
    duplicate.markers <- input %>%
      dplyr::group_by(CHROM, LOCUS, POS) %>%
      dplyr::tally(.) %>%
      dplyr::filter(n > 1) %>%
      dplyr::arrange(CHROM, LOCUS, POS)
    readr::write_tsv(x = duplicate.markers, path = "duplicated.markers.tsv")
    stop("VCF file contains duplicated CHROM, ID and POS\n
         A file named: duplicated.markers.tsv was written in the directory")
  }

  # Scan and filter with FILTER column
  filter.check <- dplyr::distinct(input, FILTER)

  if (nrow(filter.check) > 1) {
    message("Filtering markers based on VCF FILTER column")
    nrow.before <- nrow(input)
    input <- dplyr::filter(input, FILTER %in% "PASS")
    nrow.after <- nrow(input)
    message("    Number of markers before = ", nrow.before)
    message("    Number of markers removed = ", nrow.before - nrow.after)
    message("    Number of markers after = ", nrow.after)
  }
  input <- dplyr::select(input, -FILTER)
  filter.check <- NULL

  # GATK VCF file sometimes have "." in the LOCUS column: replace by CHROM
  # platypus VCF file sometimes have NA in LOCUS column: replace by POS
  weird.locus <- unique(input$LOCUS)
  if (length(weird.locus) <= 1) {
    # if (is.na(weird.locus)) {
    input$LOCUS <- input$POS
    # } else {
    #   input$LOCUS <- input$CHROM
    # }
  }
  weird.locus <- NULL #unused object

  # Unique MARKERS column --------------------------------------------------------

  # Since stacks v.1.44 ID as LOCUS + COL (from sumstats) the position of the SNP on the locus.
  # Choose the first 100 markers to scan
  detect.snp.col <- sample(x = unique(input$LOCUS), size = 100, replace = FALSE) %>%
    stringi::stri_detect_fixed(str = ., pattern = "_") %>%
    unique

  if (detect.snp.col && stacks.vcf) {
    if (nrow(input) > 30000) {
      input <- input %>%
        dplyr::mutate(
          SPLIT_VEC = split_vec_row(x = ., cpu.rounds = 3,
                                    parallel.core = parallel.core)) %>%
        split(x = ., f = .$SPLIT_VEC) %>%
        .radiator_parallel_mc(
          X = .,
          FUN = split_vcf_id,
          mc.cores = parallel.core
        ) %>%
        dplyr::bind_rows(.) %>%
        dplyr::select(-SPLIT_VEC)
    } else {
      input <- dplyr::rename(input, ID = LOCUS) %>%
        tidyr::separate(data = ., col = ID, into = c("LOCUS", "COL"),
                        sep = "_", extra = "drop", remove = FALSE) %>%
        dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = as.character) %>%
        dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = clean_markers_names) %>%
        tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE)
    }
  } else {
    input <- dplyr::mutate_at(
      .tbl = input,
      .vars = c("CHROM", "POS", "LOCUS"), .funs = clean_markers_names) %>%
      tidyr::unite(
        data = .,
        MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE)
  }

  # Filter with whitelist of markers and FILTER column -------------------------
  if (!is.null(whitelist.markers)) {

    if (!biallelic) {
      if (ncol(whitelist.markers) >= 3) {
        message("Note: whitelist with CHROM LOCUS POS columns and VCF haplotype:
                If the whitelist was not created from this VCF,
                the filtering could result in loosing all the markers.
                The POS column is different in biallelic and multiallelic file...\n")

        message("Discarding the POS column in the whitelist")
        whitelist.markers <- dplyr::select(whitelist.markers, -POS)
      }

      if (ncol(whitelist.markers) == 1 && tibble::has_name(whitelist.markers, "MARKERS")) {
        message("Note: whitelist MARKERS column and VCF haplotype:
                If the whitelist was not created from this VCF,
                the filtering could result in loosing all the markers.
                The POS column used in the MARKERS column is different in biallelic and multiallelic file...\n")
      }
    }
    message("Filtering: ", nrow(whitelist.markers), " markers in whitelist")
    columns.names.whitelist <- colnames(whitelist.markers)
    input <- suppressWarnings(
      dplyr::semi_join(
        input, whitelist.markers, by = columns.names.whitelist))
    if (nrow(input) == 0) stop("No markers left in the dataset, check whitelist...")
  }

  # keep vector
  keep.markers <- dplyr::select(input, KEEP) %>%
    dplyr::arrange(KEEP) %>%
    purrr::flatten_int(.)

  # import genotypes
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "ID", "COL", "REF", "ALT", "INDIVIDUALS", "GT")
  # bk <- input
  if (import.pegas) {
    if (verbose) message("Working on the vcf...")
    input.gt <- suppressWarnings(
      pegas::read.vcf(
        file = data,
        which.loci = keep.markers,
        quiet = verbose) %>%
        `colnames<-`(input$MARKERS
        ) %>%
        tibble::as_data_frame(.) %>%
        tibble::rownames_to_column(., var = "INDIVIDUALS") %>%
        # dplyr::mutate_all(.tbl = ., .funs = as.character) %>% # very long
        tidyr::gather(data = ., key = MARKERS, value = GT, -INDIVIDUALS))

    input <- suppressWarnings(
      dplyr::full_join(input.gt, input, by = "MARKERS") %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                         .funs = clean_ind_names) %>%
        dplyr::mutate(GT = stringi::stri_replace_na(
          str = GT, replacement = "./."))
    )

    input.gt <- NULL
  } else {
    # filter the vcf.data
    vcf.data <- vcf.data[keep.markers,]
    keep.markers <- NULL

    input <- suppressWarnings(dplyr::select(
      .data = input,
      dplyr::one_of(want)))#MARKERS, CHROM, LOCUS, POS, REF, ALT)

    filter.check <- NULL

    input <- dplyr::bind_cols(
      input,
      parse_genomic(x = "GT", data = vcf.data, return.alleles = TRUE,
                    verbose = verbose))

    input <- suppressWarnings(tidyr::gather(
      data = input,
      key = INDIVIDUALS,
      value = GT,
      -dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS", "ID", "COL", "REF", "ALT"))) %>%
        dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                         .funs = clean_ind_names))

    # metadata
    if (is.logical(vcf.metadata)) {
      overwrite.metadata <- NULL
    } else {
      overwrite.metadata <- vcf.metadata
      vcf.metadata <- TRUE
    }

    if (vcf.metadata) {
      if (verbose) message("Keeping vcf metadata: yes")
      # detect FORMAT fields available
      have <- suppressWarnings(vcfR::vcf_field_names(vcf.data, tag = "FORMAT")$ID)
      # current version doesn't deal well with PL with 3 fields separated with ","
      want <- c("DP", "AD", "GL", "PL", "HQ", "GQ", "GOF", "NR", "NV")
      if (!is.null(overwrite.metadata)) want <- overwrite.metadata
      parse.format.list <- purrr::keep(.x = have, .p = have %in% want)

      # work on parallelization of this part
      input <- dplyr::bind_cols(
        input,
        purrr::map(parse.format.list, parse_genomic, data = vcf.data,
                   gather.data = TRUE, verbose = verbose) %>%
          dplyr::bind_cols(.))

      # some VCF are too big to fit in memory multiple times for parallel processing...
      # work in progress
      # } else {
      # if (length(parse.format.list) < parallel.core) {
      # parallel.core.parse <- length(parse.format.list)
      # } else {
      # parallel.core.parse <- parallel.core
      # }
      #   if (verbose) message("Parsing and tidying: ", stringi::stri_join(parse.format.list, collapse = ", "))
      #   parsed.format <- list()
      #   parsed.format <- .radiator_parallel(
      #   # parsed.format <- parallel::mclapply(
      #     X = parse.format.list,
      #     FUN = parse_genomic,
      #     mc.cores = parallel.core.parse,
      #     data = vcf.data,
      #     gather.data = TRUE
      #   ) %>% dplyr::bind_cols(.)
      #   input <- dplyr::bind_cols(input, parsed.format)
      #   parsed.format <- NULL
      # }
    } else {
      if (verbose) message("Keeping vcf metadata: no")
    }
    vcf.data <- NULL
  }#End vcfR GT and metadata import

  # Import blacklist id --------------------------------------------------------
  if (!is.null(blacklist.id)) {# With blacklist of ID
    if (is.vector(blacklist.id)) {
      suppressMessages(blacklist.id <- readr::read_tsv(blacklist.id, col_names = TRUE))
    } else {
      if (!tibble::has_name(blacklist.id, "INDIVIDUALS")) {
        stop("Blacklist of individuals should have 1 column named: INDIVIDUALS")
      }
    }
    blacklist.id$INDIVIDUALS <- radiator::clean_ind_names(blacklist.id$INDIVIDUALS)

    # remove potential duplicate id
    dup <- dplyr::distinct(.data = blacklist.id, INDIVIDUALS)
    blacklist.id.dup <- nrow(blacklist.id) - nrow(dup)
    if (blacklist.id.dup >1) {
      message("Duplicate id's in blacklist: ", blacklist.id.dup)
      blacklist.id <- dup
    }
    dup <- blacklist.id.dup <- NULL
    message("Number of individuals in blacklist: ", nrow(blacklist.id))

    input <- dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS")
  }

  # Population levels and strata------------------------------------------------
  if (verbose) message("Making the vcf population-wise...")
  if (is.vector(strata)) {
    suppressMessages(
      strata.df <- readr::read_tsv(
        file = strata, col_names = TRUE,
        # col_types = col.types
        col_types = readr::cols(.default = readr::col_character())
      ))
  } else {
    strata.df <- strata
  }
  colnames(strata.df) <- stringi::stri_replace_all_fixed(
    str = colnames(strata.df),
    pattern = "STRATA",
    replacement = "POP_ID",
    vectorize_all = FALSE
  )
  # Remove potential whitespace in pop_id
  strata.df$POP_ID <- radiator::clean_pop_names(strata.df$POP_ID)
  colnames.strata <- colnames(strata.df)

  # clean ids
  strata.df$INDIVIDUALS <- radiator::clean_ind_names(strata.df$INDIVIDUALS)

  # filtering the strata if blacklist id available
  if (!is.null(blacklist.id)) {
    strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
  }


# check that names match between strata and input before going further
if (!identical(sort(unique(input$INDIVIDUALS)), sort(unique(strata.df$INDIVIDUALS)))) {
  stop("The individuals in the strata file don't match the individuals in the vcf file")
}

input <- dplyr::left_join(x = input, y = strata.df, by = "INDIVIDUALS")

# Using pop.levels and pop.labels info if present-------------------------------
input <- radiator::change_pop_names(
  data = input, pop.levels = pop.levels, pop.labels = pop.labels)

# Pop select--------------------------------------------------------------------
if (!is.null(pop.select)) {
  pop.select <- clean_pop_names(pop.select)

  if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
  input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
  input$POP_ID <- droplevels(input$POP_ID)
}

# Haplotypes or biallelic VCF---------------------------------------------------
# recoding genotype
if (biallelic) {# biallelic VCF
  if (verbose) message("Recoding bi-allelic VCF...")
  input <- dplyr::rename(input, GT_VCF_NUC = GT)
} else {#multi-allelic vcf
  if (verbose) message("Recoding VCF haplotype...")
  input <- dplyr::rename(input, GT_HAPLO = GT)
}

if (verbose) message("Calculating REF/ALT alleles...")
input <- radiator::change_alleles(
  data = input,
  biallelic = biallelic,
  parallel.core = parallel.core,
  verbose = verbose)$input

# Re ordering columns
want <- c("MARKERS", "CHROM", "LOCUS", "POS", "ID", "COL", "INDIVIDUALS", "POP_ID",
          "REF", "ALT", "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN",
          "POLYMORPHIC")

input <- suppressWarnings(
  dplyr::select(input, dplyr::one_of(want), dplyr::everything()))

# More VCF cleaning here -----------------------------------------------------
# check, parse and clean FORMAT columns
# Some software that produce vcf do strange thing and don't follow convention
# You need the GT field to clean correctly the remaining fields...
if (vcf.metadata) {

  input <-  dplyr::mutate_at(
    .tbl = input, .vars = parse.format.list, .funs =  replace_by_na)

  split.vec <- split_vec_row(input, 3, parallel.core = parallel.core)

  # Cleaning AD (ALLELES_DEPTH)
  if (tibble::has_name(input, "AD")) {
    if (verbose) message("AD column: splitting coverage info into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH")
    input <- clean_ad(x = input, split.vec = split.vec,
                      parallel.core = parallel.core)
  }#End cleaning AD column

  # Cleaning DP and changing name to READ_DEPTH
  if (tibble::has_name(input, "DP")) {
    if (verbose) message("DP column: cleaning and renaming to READ_DEPTH")
    input <- dplyr::rename(.data = input, READ_DEPTH = DP) %>%
      dplyr::mutate(
        READ_DEPTH = dplyr::if_else(GT_VCF == "./.", as.numeric(NA_character_),
                                    as.numeric(READ_DEPTH))
      )
  }#End cleaning DP column

  # PL for biallelic as 3 values:
  if (tibble::has_name(input, "PL")) {
    if (verbose) message("PL column (normalized, phred-scaled likelihoods for genotypes): separating into PROB_HOM_REF, PROB_HET and PROB_HOM_ALT")
    # Value 1: probability that the site is homozgyous REF
    # Value 2: probability that the sample is heterzygous
    # Value 2: probability that it is homozygous ALT
    input <- clean_pl(x = input, split.vec = split.vec,
                      parallel.core = parallel.core)
  }#End cleaning PL column

  # GL cleaning
  if (tibble::has_name(input, "GL")) {
    if (verbose) message("GL column: cleaning Genotype Likelihood column")
    input <- clean_gl(x = input,
                      split.vec = split.vec,
                      parallel.core = parallel.core)
  }#End cleaning GL column

  # Cleaning GQ: Genotype quality as phred score
  if (tibble::has_name(input, "GQ")) {
    if (verbose) message("GQ column: Genotype Quality")
    input <- dplyr::mutate(
      input,
      GQ = dplyr::if_else(GT_VCF == "./.", as.numeric(NA_character_), as.numeric(GQ))
    )
  }#End cleaning GQ column

  # Cleaning GOF: Goodness of fit value
  if (tibble::has_name(input, "GOF")) {
    if (verbose) message("GOF column: Goodness of fit value")
    input <- dplyr::mutate(
      input,
      GOF = dplyr::if_else(GT_VCF == "./.", as.numeric(NA_character_), as.numeric(GOF))
    )
  }#End cleaning GOF column

  # Cleaning NR: Number of reads covering variant location in this sample
  if (tibble::has_name(input, "NR")) {
    if (verbose) message("NR column: splitting column into the number of variant")
    input <- clean_nr(x = input,
                      split.vec = split.vec,
                      parallel.core = parallel.core)
  }#End cleaning NR column

  # Cleaning NV: Number of reads containing variant in this sample
  if (tibble::has_name(input, "NV")) {
    if (verbose) message("NV column: splitting column into the number of variant")
    input <- clean_nv(x = input,
                      split.vec = split.vec,
                      parallel.core = parallel.core)
  }#End cleaning NV column

  split.vec <- NULL

  # Re ordering columns
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID",
            "REF", "ALT", "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN")

  input <- suppressWarnings(
    dplyr::select(input, dplyr::one_of(want), dplyr::everything()))

}# end cleaning columns

# Sort id
input <- dplyr::arrange(input, POP_ID, INDIVIDUALS)

return(input)
}#End tidy_vcf


# Internal nested Function -----------------------------------------------------
#' @title parse_genomic
#' @description function to parse the format field and tidy the results of VCF
#' @rdname parse_genomic
#' @keywords internal
#' @export
parse_genomic <- function(
  x, data = NULL, mask = FALSE, gather.data = FALSE, return.alleles = FALSE,
  verbose = TRUE) {
  format.name <- x

  if (verbose) message("Parsing and tidying: ", format.name)
  x <- tibble::as_data_frame(vcfR::extract.gt(
    x = data, element = format.name, mask = mask, return.alleles = return.alleles,
    IDtoRowNames = FALSE, convertNA = FALSE))

  if (format.name == "GT") {
    colnames(x) <- stringi::stri_replace_all_fixed(
      str = colnames(x),
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )
  }

  if (gather.data) {
    x <- dplyr::mutate(x, ID = seq(1, n()))
    x <- tidyr::gather(data = x, key = INDIVIDUALS, value = rlang::UQ(format.name), -ID) %>%
      dplyr::select(-ID, -INDIVIDUALS)
  }
  return(x)
}#End parse_genomic

#' @title split_vcf_id
#' @description split VCF ID in parallel
#' @rdname split_vcf_id
#' @keywords internal
#' @export

split_vcf_id <- function(x) {
  res <- dplyr::rename(x, ID = LOCUS) %>%
    tidyr::separate(data = ., col = ID, into = c("LOCUS", "COL"),
                    sep = "_", extra = "drop", remove = FALSE) %>%
    dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = as.character) %>%
    dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = clean_markers_names) %>%
    tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE)
  return(res)
}#End split_vcf_id



#' @title clean_ad
#' @description Clean allele depth field in VCF.
#' AD column: splitting coverage info into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH
#' @rdname clean_ad
#' @keywords internal
#' @export
clean_ad <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {
  clean <- function(x) {
    res <- suppressWarnings(
      x %>%
        dplyr::mutate(
          AD = dplyr::if_else(GT_VCF == "./.", NA_character_, AD)) %>%
        tidyr::separate(AD, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
                        sep = ",", extra = "drop") %>%
        dplyr::mutate(
          ALLELE_REF_DEPTH = as.numeric(
            stringi::stri_replace_all_regex(
              ALLELE_REF_DEPTH, "^0$", "NA", vectorize_all = TRUE)),
          ALLELE_ALT_DEPTH = as.numeric(
            stringi::stri_replace_all_regex(
              ALLELE_ALT_DEPTH, "^0$", "NA", vectorize_all = TRUE))
        ) %>%
        dplyr::select(-GT_VCF)
    )
    return(res)
  }
  x <- dplyr::bind_cols(
    x,
    dplyr::ungroup(x) %>%
      dplyr::select(GT_VCF, AD) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel(
        X = ., FUN = clean, mc.cores = parallel.core) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_ad


#' @title clean_pl
#' @description Clean PL column.
#' PL column (normalized, phred-scaled likelihoods for genotypes):
#' separating into PROB_HOM_REF, PROB_HET and PROB_HOM_ALT")
#' Value 1: probability that the site is homozgyous REF
#' Value 2: probability that the sample is heterzygous
#' Value 2: probability that it is homozygous ALT
#' @rdname clean_pl
#' @keywords internal
#' @export
clean_pl <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {
  clean <- function(x) {
    res <- x %>%
      dplyr::mutate(
        PL = dplyr::if_else(GT_VCF == "./.", NA_character_, PL)) %>%
      tidyr::separate(
        data = ., PL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
        sep = ",", extra = "drop", remove = FALSE) %>%
      dplyr::mutate_at(
        .tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
        .funs = as.numeric) %>%
      dplyr::select(-GT_VCF)
    return(res)
  }#End clean
  x <- dplyr::bind_cols(
    dplyr::select(x, -PL),
    dplyr::ungroup(x) %>%
      dplyr::select(GT_VCF, PL) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel(
        X = ., FUN = clean, mc.cores = parallel.core) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_pl

#' @title clean_gl
#' @description Clean GL column.
#' GL column: cleaning Genotype Likelihood column
#' @rdname clean_gl
#' @keywords internal
#' @export

clean_gl <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {
  x <- x %>%
    dplyr::mutate(
      GL = dplyr::if_else(GT_VCF == "./.", NA_character_, GL),
      GL = suppressWarnings(
        stringi::stri_replace_all_fixed(
          GL, c(".,.,.", ".,", ",."), c("NA", "", ""), vectorize_all = FALSE))
    )

  # check GL and new stacks version with no GL
  all.missing <- all(is.na(x$GL))

  if (!all.missing) {
    gl.clean <- max(
      unique(stringi::stri_count_fixed(
        str = unique(sample(x = x$GL, size = 100, replace = FALSE)),
        pattern = ",")
      ), na.rm = TRUE
    )

    if (gl.clean == 2) {
      message("GL column: separating into PROB_HOM_REF, PROB_HET and PROB_HOM_ALT")
      # Value 1: probability that the site is homozgyous REF
      # Value 2: probability that the sample is heterzygous
      # Value 2: probability that it is homozygous ALT
      # system.time(input2 <- input %>%
      #   tidyr::separate(data = ., GL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"), sep = ",", extra = "drop", remove = FALSE) %>%
      #   dplyr::mutate_at(.tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"), .funs = as.numeric)
      # )
      clean <- function(x) {
        res <- x %>%
          tidyr::separate(
            data = ., GL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
            sep = ",", extra = "drop", remove = FALSE) %>%
          dplyr::mutate_at(
            .tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
            .funs = as.numeric)
        return(res)
      }

      x <- dplyr::bind_cols(
        dplyr::select(x, -GL),
        dplyr::ungroup(x) %>%
          dplyr::select(GL) %>%
          split(x = ., f = split.vec) %>%
          .radiator_parallel(
            X = ., FUN = clean, mc.cores = parallel.core) %>%
          dplyr::bind_rows(.))

    } else {
      x$GL <- suppressWarnings(as.numeric(x$GL))
    }
  } else {
    message("    GL values are all missing: removing column")
    x <- dplyr::select(x, -GL)
  }
  return(x)
}#End clean_gl

#' @title clean_nr
#' @description Cleaning NR: Number of reads covering variant location in
#' this sample
#' @rdname clean_nr
#' @keywords internal
#' @export
clean_nr <- function(x, split.vec, parallel.core = parallel::detectCores() - 1){
  nr.col <- max(unique(stringi::stri_count_fixed(str = unique(x$NR), pattern = ","))) + 1
  nr.col.names <- stringi::stri_join(rep("NR_", nr.col), seq(1:nr.col))
  nr.col <- NULL

  clean <- function(x, nr.col.names = NULL) {
    res <- tidyr::separate(data = x, col = NR, into = nr.col.names,
                           sep = ",", extra = "drop", remove = FALSE)
    return(res)
  }#End clean

  x <- dplyr::bind_cols(
    dplyr::select(x, -NR),
    dplyr::ungroup(x) %>%
      dplyr::mutate(NR = dplyr::if_else(GT_VCF == "./.", NA_character_, NR)) %>%
      dplyr::select(NR) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel_mc(
        X = ., FUN = clean, mc.cores = parallel.core,
        nr.col.names = nr.col.names) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_nr

#' @title clean_nv
#' @description Cleaning NV: Number of reads containing variant in this sample
#' @rdname clean_nv
#' @keywords internal
#' @export
clean_nv <- function(x, split.vec, parallel.core = parallel::detectCores() - 1) {

  nv.col <- max(unique(stringi::stri_count_fixed(str = unique(x$NV), pattern = ","))) + 1
  nv.col.names <- stringi::stri_join(rep("NV_", nv.col), seq(1:nv.col))
  nv.col <- NULL

  clean <- function(x, nv.col.names = NULL) {
    res <- tidyr::separate(
      data = x, col = NV, into = nv.col.names,
      sep = ",", extra = "drop", remove = FALSE)
    return(res)
  }

  x <- dplyr::bind_cols(
    dplyr::select(x, -NV),
    dplyr::ungroup(x) %>%
      dplyr::mutate(NV = dplyr::if_else(GT_VCF == "./.", NA_character_, NV)) %>%
      dplyr::select(NV) %>%
      split(x = ., f = split.vec) %>%
      .radiator_parallel_mc(
        X = ., FUN = clean, mc.cores = parallel.core,
        nv.col.names = nv.col.names) %>%
      dplyr::bind_rows(.))
  return(x)
}#End clean_nv
