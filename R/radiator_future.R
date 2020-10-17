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
    .x,
    .f,
    parallel.core = parallel::detectCores() - 1,
    split.with = NULL,
    flat.future = c("int", "chr", "dfr", "dfc", "walk", "drop"),
    ...
  ) {

    # argument for flattening the results
    flat.future <- match.arg(
      arg = flat.future,
      choices = c("int", "chr", "dfr", "dfc", "walk", "drop"),
      several.ok = FALSE
    )

    if (!is.null(split.with)) {
      if (length(split.with) == 1 && is.character(split.with)) {
        .x %<>% dplyr::group_split(.tbl = ., .data[[split.with]])
      } else {
        .x %<>% split(x = ., f = split.with)
      }
    }

    if (parallel.core == 1L) {
      future::plan(strategy = "sequential")
    } else {
      parallel.core <- parallel_core_opt(parallel.core = parallel.core)
      lx <- length(.x)
      if (lx < parallel.core) {
        future::plan(strategy = "multisession", workers = lx)
      } else {
        future::plan(strategy = "multisession", workers = parallel.core)
      }
    }

    # .x <- future.apply::future_apply(X = .x, FUN = .f, ...)
    # capture dots
    # d <- rlang::dots_list(..., .ignore_empty = "all", .preserve_empty = TRUE, .homonyms = "first")
    # if (bind.rows) .x %<>% dplyr::bind_rows(.)



    # Run the function in parallel and account for dots-dots-dots argument
    rad_map <- switch(flat.future,
                      int = {furrr::future_map_int},
                      chr = {furrr::future_map_chr},
                      dfr = {furrr::future_map_dfr},
                      dfc = {furrr::future_map_dfc},
                      walk = {furrr::future_walk},
                      drop = {furrr::future_map}
    )

    opts <- furrr::furrr_options(globals = FALSE)
    if (length(list(...)) == 0) {
      .x %<>% rad_map(.x = ., .f = .f, .options = opts)
    } else {
      .x %<>% rad_map(.x = ., .f = .f, ..., .options = opts)
    }

    if (parallel.core > 1L) future::plan(strategy = "sequential")
    return(.x)
  }#End radiator_parallel
