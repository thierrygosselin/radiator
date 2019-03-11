# replace_by_na ----------------------------------------------------------------
#' @title replace_by_na
#' @description Fast removal of NA
#' @rdname replace_by_na
#' @keywords internal
#' @export
replace_by_na <- function(data, what = ".") {
  replace(data, which(data == what), NA)
}#End replace_by_na

# distance2tibble---------------------------------------------------------------
#' @title distance2tibble
#' @description melt the distance matrix into a tibble
#' @rdname distance2tibble
#' @export
#' @keywords internal
distance2tibble <- function(
  x,
  remove.diag = TRUE,
  na.diag = FALSE,
  remove.lower = TRUE,
  relative = TRUE,
  pop.levels = NULL,
  distance.class.double = TRUE
) {
  # x <- dist.computation
  x <- as.matrix(x)
  if (remove.diag || na.diag) diag(x) <- NA
  if (remove.lower) x[lower.tri(x)] <- NA
  x <- dplyr::bind_cols(tibble::tibble(ID1 = rownames(x)),
                        tibble::as_tibble(x)) %>%
    data.table::as.data.table(.) %>%
    data.table::melt.data.table(
      data = ., id.vars = "ID1", variable.name = "ID2", value.name = "DISTANCE",
      variable.factor = FALSE) %>%
    tibble::as_tibble(.)

  if (!na.diag) x  %<>% dplyr::filter(!is.na(DISTANCE))

  if (distance.class.double) {
    x %<>% dplyr::mutate(DISTANCE = as.double(as.character(DISTANCE)))
  }

  x %<>% dplyr::arrange(DISTANCE)

  if (relative) {
    x  %<>% dplyr::mutate(DISTANCE_RELATIVE = DISTANCE/max(DISTANCE))
  }
  if (!is.null(pop.levels)) {
    x  %<>% dplyr::mutate(
      ID1 = factor(x = ID1, levels = pop.levels, ordered = TRUE),
      ID2 = factor(x = ID2, levels = pop.levels, ordered = TRUE)
    )
  }

  return(x)
}#End distance2tibble


