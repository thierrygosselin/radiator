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
separate_gt <- function(x, sep = "/", gt = "GT_VCF_NUC", exclude = c("LOCUS", "INDIVIDUALS", "POP_ID"), cpu.rounds = 10, parallel.core = parallel::detectCores() - 1) {
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
