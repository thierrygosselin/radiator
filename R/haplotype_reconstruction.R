# haplotype_reconstruction
#' @title haplotype_reconstruction
#' @description Reconstruct haplotypes
#' @rdname haplotype_reconstruction
#' @inheritParams radiator_common_arguments
#' @keywords internal
#' @export
haplotype_reconstruction <- function(
  data,
  parallel.core = parallel::detectCores() - 1
) {
  # data <- haplo.reconstruction
  reconstruct <- carrier::crate(function(m, data) {
    `%>%` <- magrittr::`%>%`
    `%<>%` <- magrittr::`%<>%`
    # m <- "102632"
    # data <- data
    data <- dplyr::filter(data, MARKERS %in% m)
    n.snp <- unique(data$SNP_N)
    data <- tidyr::separate(
      data = data,
      col = HAPLOTYPES,
      into = as.character(seq(1, n.snp, 1)),
      sep = 1:(n.snp - 1),
      remove = FALSE
    ) %>%
      rad_long(
        x = .,
        cols = c("MARKERS", "HAPLOTYPES", "SNP_N"),
        names_to = "SNP",
        values_to = "NUC",
        variable_factor = FALSE
      ) %>%
      dplyr::mutate(SNP = as.integer(SNP)) %>%
      dplyr::group_by(SNP) %>%
      dplyr::mutate(
        POLYMORPHIC = dplyr::if_else(length(unique(NUC)) > 1,
                                     "polymorphic", "monomorphic")) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(POLYMORPHIC == "polymorphic") %>%
      dplyr::select(-POLYMORPHIC) %>%
      dplyr::arrange(SNP, HAPLOTYPES) %>%
      rad_wide(
        x = .,
        formula = "MARKERS + HAPLOTYPES + SNP_N ~ SNP",
        values_from = "NUC"
        ) %>%
      tidyr::unite(
        data = ., col = HAPLOTYPES_NEW,
        -c(MARKERS, HAPLOTYPES, SNP_N),
        sep = "") %>%
      dplyr::ungroup(.) %>%
      dplyr::select(MARKERS, HAPLOTYPES, HAPLOTYPES_NEW)
    return(data)
  })#End reconstruct

  res <- radiator_future(
    .x = dplyr::ungroup(data),
    .f = reconstruct,
    flat.future = "dfr",
    split.vec = FALSE,
    split.with = "MARKERS",
    split.chunks = 4L,
    parallel.core = parallel.core,
  )
  return(res)
}#End haplotype_reconstruction
