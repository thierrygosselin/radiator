# using future and future.apply -------------------------------------------------
#' @name radiator_future
#' @title radiator parallel function
#' @description Updating radiator to use future
# @inheritParams future::plan
# @inheritParams future::availableCores
#' @inheritParams future.apply::future_apply
#' @rdname radiator_future
#' @keywords internal
radiator_future <- function(
  X,
  FUN,
  parallel.core = parallel::detectCores() - 1,
  split.tibble = NULL,
  bind.rows = FALSE,
  flatten.int = FALSE,
  ...
) {
  if (parallel.core == 1L) {
    future::plan(strategy = "sequential")
  } else {
    future::plan(strategy = "multisession", workers = parallel.core)
  }
  if (!is.null(split.tibble)) X %<>% split(x = ., f = split.tibble)
  res <- future.apply::future_apply(X = X, FUN = FUN, ...)
  if (parallel.core > 1L) future::plan(strategy = "sequential")
  if (bind.rows) res %<>% dplyr::bind_rows(.)
  if (flatten.int) res %<>% purrr::flatten_int(.)
  return(res)
}#End radiator_parallel
