#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL



# .onUnload <- function(libpath) {
#   library.dynam.unload("radiator", libpath)
# }


#' @title split_vec_row
#' @description Split input into chunk for parallel processing
#' @rdname split_vec_row
#' @keywords internal
#' @export
split_vec_row <- function(x, cpu.rounds, parallel.core = parallel::detectCores() - 1) {
  n.row <- nrow(x)
  split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.row - 1) / n.row) + 1))
  return(split.vec)
}#End split_vec_row


#' @title replace_by_na
#' @description Fast removal of NA
#' @rdname replace_by_na
#' @keywords internal
#' @export
replace_by_na <- function(data, what = ".") {
  replace(data, which(data == what), NA)
}#End replace_by_na

#' @title separate_gt
#' @description Separate genotype field
#' @rdname separate_gt
#' @keywords internal
#' @export
separate_gt <- function(
  x,
  sep = "/",
  gt = "GT_VCF_NUC",
  exclude = c("LOCUS", "INDIVIDUALS", "POP_ID"),
  cpu.rounds = 10,
  parallel.core = parallel::detectCores() - 1
) {
  # sep <-  "/"
  # x <- input2
  n.row <- nrow(x)
  split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.row - 1) / n.row) + 1))

  separate_genotype <- function(x, sep, gt, exclude){
    res <- tidyr::separate(
      data = x,
      col = gt, into = c("ALLELE1", "ALLELE2"),
      sep = sep,
      extra = "drop", remove = TRUE
    ) %>%
      tidyr::gather(data = ., key = ALLELE_GROUP, value = HAPLOTYPES, -dplyr::one_of(exclude))
    return(res)
  }

  res <- split(x = x, f = split.vec) %>%
    .radiator_parallel_mc(
      X = .,
      FUN = separate_genotype,
      mc.cores = parallel.core,
      sep = sep, gt = gt, exclude = exclude) %>%
    dplyr::bind_rows(.)
  return(res)
}#End separate_gt


#' @title compute_maf
#' @description Compute MAF
#' @rdname compute_maf
#' @keywords internal
#' @export
compute_maf <- function(x, biallelic) {
  if (tibble::has_name(x, "GT_BIN") && biallelic) {
    maf.data <- x %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::summarise(
        NN = as.numeric(2 * n()),
        PP = as.numeric(2 * length(GT_BIN[GT_BIN == 0])),
        PQ = as.numeric(length(GT_BIN[GT_BIN == 1])),
        QQ = as.numeric(2 * length(GT_BIN[GT_BIN == 2]))
      ) %>%
      # need this step because seen cases where the minor allele is not minor
      dplyr::mutate(
        PP = PP + PQ,
        QQ = QQ + PQ,
        PQ = NULL) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        NN_G = sum(NN),
        PP_G = sum(PP),
        QQ_G = sum(QQ)) %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::mutate(
        ALT = dplyr::if_else(PP_G < QQ_G, PP, QQ),
        MAF_LOCAL = (ALT / NN),
        PP = NULL,
        QQ = NULL,
        NN = NULL) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        ALT = dplyr::if_else(PP_G < QQ_G, PP_G, QQ_G),
        MAF_GLOBAL = (ALT / NN_G),
        ALT = NULL,
        PP_G = NULL,
        QQ_G = NULL,
        NN_G = NULL) %>%
      dplyr::ungroup(.)
  } else {
    if (!tibble::has_name(x, "GT_VCF_NUC")) {
      maf.data <- x %>%
        dplyr::select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(GT, 1, 3),
          A2 = stringi::stri_sub(GT, 4,6)
        ) %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID))

      maf.local <- maf.data %>%
        dplyr::group_by(MARKERS, POP_ID, GT) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_LOCAL = n / n.al.tot) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL)

      maf.data <- maf.data %>%
        dplyr::group_by(MARKERS, GT) %>%
        dplyr::tally(.) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(n.al.tot = sum(n)) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::summarise(MAF_GLOBAL = n / n.al.tot) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(MARKERS, MAF_GLOBAL) %>%
        dplyr::left_join(maf.local, by = c("MARKERS")) %>%
        dplyr::select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
    } else {
      maf.data <- x %>%
        dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_VCF_NUC) %>%
        tidyr::separate(
          data = .,
          col = GT_VCF_NUC, into = c("A1", "A2"),
          sep = "/",
          extra = "drop", remove = TRUE
        ) %>%
        tidyr::gather(
          data = ., key = ALLELE_GROUP, value = HAPLOTYPES,
          -dplyr::one_of(c("MARKERS", "INDIVIDUALS", "POP_ID"))) %>%
        dplyr::select(-ALLELE_GROUP) %>%
        dplyr::group_by(MARKERS, HAPLOTYPES, POP_ID) %>%
        dplyr::tally(.) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, HAPLOTYPES), fill = list(n = 0)) %>%
        dplyr::group_by(MARKERS, POP_ID) %>%
        dplyr::mutate(N_LOCAL = sum(n)) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(N_GLOBAL = sum(n)) %>%
        dplyr::arrange(MARKERS, POP_ID) %>%
        dplyr::group_by(MARKERS, POP_ID, HAPLOTYPES) %>%
        dplyr::mutate(MAF_LOCAL = n / N_LOCAL) %>%
        dplyr::group_by(MARKERS, HAPLOTYPES) %>%
        dplyr::mutate(
          MAF_GLOBAL = sum(n) / N_GLOBAL,
          n = NULL,
          N_LOCAL = NULL,
          N_GLOBAL = NULL
        ) %>%
        dplyr::ungroup(.)

      ref.info <- dplyr::distinct(maf.data, MARKERS, HAPLOTYPES, MAF_GLOBAL) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::filter(MAF_GLOBAL == max(MAF_GLOBAL)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(REF = rep("REF", n()), MAF_GLOBAL = NULL) %>%
        dplyr::bind_rows(
          dplyr::distinct(maf.data, MARKERS, HAPLOTYPES, MAF_GLOBAL) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::filter(MAF_GLOBAL == min(MAF_GLOBAL)) %>%
            dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
            dplyr::ungroup(.) %>%
            dplyr::mutate(REF = rep("MAF", n()), MAF_GLOBAL = NULL)
        )

      maf.data <- dplyr::left_join(maf.data, ref.info, by = c("MARKERS", "HAPLOTYPES")) %>%
        dplyr::mutate(REF = stringi::stri_replace_na(REF, replacement = "ALT"))
      ref.info <- NULL
    }
  }
  maf.local <- NULL
  return(maf.data)
}#End compute_maf
