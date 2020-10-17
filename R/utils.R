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

  if (na.diag || remove.diag) x  %<>% dplyr::filter(!is.na(DISTANCE))

  if (distance.class.double) {
    x %<>% dplyr::mutate(DISTANCE = as.double(as.character(DISTANCE)))
  }

  x %<>% dplyr::arrange(DISTANCE)

  if (relative) {
    x  %<>% dplyr::mutate(DISTANCE_RELATIVE = DISTANCE / max(DISTANCE))
  }
  if (!is.null(pop.levels)) {
    x  %<>% dplyr::mutate(
      ID1 = factor(x = ID1, levels = pop.levels, ordered = TRUE),
      ID2 = factor(x = ID2, levels = pop.levels, ordered = TRUE)
    )
  } else {
    x  %<>% dplyr::mutate(
      ID1 = factor(x = ID1),
      ID2 = factor(x = ID2)
    )
  }

  return(x)
}#End distance2tibble


# split_vec_row ----------------------------------------------------------------
#' @title split_vec_row
#' @description Split input into chunk for parallel processing
#' @rdname split_vec_row
#' @keywords internal
#' @export
split_vec_row <- function(x, cpu.rounds, parallel.core = parallel::detectCores() - 1) {
  if (!is.integer(x)) {
    n.row <- nrow(x)
  } else {
    n.row <- x
  }
  split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.row - 1) / n.row) + 1))
  return(split.vec)
}#End split_vec_row

#' @title split_tibble_rows
#' @description Split rows of tibble for parallel processing
#' @rdname split_tibble_rows
#' @keywords internal
#' @export
split_tibble_rows <- function(
  x,
  lines.cpu = 1000, #lines per CPU rounds
  parallel.core = parallel::detectCores() - 1,
  group.split = TRUE # does dplyr: group_by and group_split
) {
  n.row <- nrow(x)
  n.cores <- parallel::detectCores()
  if (parallel.core > n.cores) parallel.core <- n.cores
  if (n.row < parallel.core) parallel.core <- n.row
  if (lines.cpu > n.row) lines.cpu <- n.row
  lines.rounds <- parallel.core * lines.cpu
  x$SPLIT_VEC <- sort(rep_len(x = 1:floor(n.row / lines.rounds), length.out = n.row))
  if (group.split) {
    x %<>%
      dplyr::group_by(SPLIT_VEC) %>%
      dplyr::group_split(.tbl = ., .keep = FALSE)
  }
  return(x)
}#End split_tibble_rows


# parallel_core_opt ------------------------------------------------------------
#' @title parallel_core_opt
#' @description Optimization of parallel core argument for radiator
#' @keywords internal
#' @export
parallel_core_opt <- function(parallel.core = NULL, max.core = NULL) {
  # strategy:
  # minimum of 1 core and a maximum of all the core available -2
  # even number of core
  # test
  # parallel.core <- 1
  # parallel.core <- 2
  # parallel.core <- 3
  # parallel.core <- 11
  # parallel.core <- 12
  # parallel.core <- 16
  # max.core <- 5
  # max.core <- 50
  # max.core <- NULL

  # Add-ons options
  # to control the max and min number to use...

  if (is.null(parallel.core)) {
    parallel.core <- parallel::detectCores() - 2
  } else {
    parallel.core <- floor(parallel.core / 2) * 2
    parallel.core <- max(1, min(parallel.core, parallel::detectCores() - 2))
  }

  if (is.null(max.core)) {
    parallel.core.opt <- parallel.core
  } else {
    parallel.core.opt <- min(parallel.core, floor(max.core / 2) * 2)
  }
  return(parallel.core.opt)
}#End parallel_core_opt
