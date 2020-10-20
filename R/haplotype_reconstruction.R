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
      into = as.character(seq(1, n.snp, 1)), sep = 1:(n.snp - 1), remove = FALSE) %>%
      data.table::as.data.table(.) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = c("MARKERS", "HAPLOTYPES", "SNP_N"),
        variable.name = "SNP",
        value.name = "NUC",
        variable.factor = FALSE
      ) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(SNP = as.integer(SNP)) %>%
      dplyr::group_by(SNP) %>%
      dplyr::mutate(
        POLYMORPHIC = dplyr::if_else(length(unique(NUC)) > 1,
                                     "polymorphic", "monomorphic")) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(POLYMORPHIC == "polymorphic") %>%
      dplyr::select(-POLYMORPHIC) %>%
      dplyr::arrange(SNP, HAPLOTYPES) %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = MARKERS + HAPLOTYPES + SNP_N ~ SNP,
        value.var = "NUC"
      ) %>%
      tibble::as_tibble(.) %>%
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
