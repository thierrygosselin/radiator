# haplotype_reconstruction
#' @title haplotype_reconstruction
#' @description Reconstruct haplotypes
#' @rdname haplotype_reconstruction
#' @keywords internal
#' @export
haplotype_reconstruction <- function(
  data,
  parallel.core = parallel::detectCores() - 1
  ) {
  # data <- haplo.reconstruction
  data <- dplyr::ungroup(data)
  markers <- dplyr::distinct(data, MARKERS) %>% purrr::flatten_chr(.)

  reconstruct <- function(m, data) {
    # m <- "102632"
    # data <- data
    data <- dplyr::filter(data, MARKERS %in% m)
    n.snp <- unique(data$SNP_N)
    data <- tidyr::separate(
      data = data,
      col = HAPLOTYPES,
      into = as.character(seq(1, n.snp, 1)), sep = 1:(n.snp-1), remove = FALSE) %>%
      tidyr::gather(
        data = ., key = "SNP",
        value = "NUC",
        -c(MARKERS, HAPLOTYPES, SNP_N)) %>%
      dplyr::mutate(SNP = as.integer(SNP)) %>%
      dplyr::group_by(SNP) %>%
      dplyr::mutate(
        POLYMORPHIC = dplyr::if_else(length(unique(NUC)) > 1,
                                     "polymorphic", "monomorphic")) %>%
      dplyr::ungroup(.) %>%
      dplyr::filter(POLYMORPHIC == "polymorphic") %>%
      dplyr::select(-POLYMORPHIC) %>%
      dplyr::arrange(SNP, HAPLOTYPES) %>%
      dplyr::group_by(MARKERS, HAPLOTYPES, SNP_N) %>%
      tidyr::spread(SNP, NUC, convert = TRUE) %>%
      tidyr::unite(
        data = ., col = HAPLOTYPES_NEW,
        -c(MARKERS, HAPLOTYPES, SNP_N),
        sep = "") %>%
      dplyr::ungroup(.) %>%
      dplyr::select(MARKERS, HAPLOTYPES, HAPLOTYPES_NEW)
    return(data)
  }

  res <- .radiator_parallel(
    X = markers,
    FUN = reconstruct,
    mc.cores = parallel.core,
    data = data
  ) %>%
    dplyr::bind_rows(.)
  return(res)
}#End haplotype_reconstruction
