# mclappply_win ----------------------------------------------------------------
#' @name mclapply_win
#' @title hack to switch function for parallel computation based on OS
#' @description \code{parallel::mclapply} doesn't work on Windows,
#' because forking is not supported.
#' This function defines a socket version of mclapply for windows computer
#' An implementation that switch automatically the parallel process when detecting
#' the os.
#' The code below was inspired from
#' \pkg{parallel} \code{\link{mclapply}},
#' \href{https://github.com/nathanvan}{Nathan VanHoudnos},
#' \href{https://github.com/kvnkuang/pbmcapply}{Kevin Kuang},
#' \href{https://github.com/psolymos/pbapply}{Peter Solymos} and
#' \href{https://github.com/EricArcher/}{Eric Archer}.
# @inheritParams parallel::mclapply
# Doesnt work and throws an error for bad markup so have to do it manually until
# parallel fix this bug
#' @param X see \pkg{parallel} \code{\link{mclapply}}
#' @param FUN see \pkg{parallel} \code{\link{mclapply}}
#' @param ... see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.preschedule see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.set.seed see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.silent see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.cores see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.cleanup see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.allow.recursive see \pkg{parallel} \code{\link{mclapply}}
#' @param affinity.list see \pkg{parallel} \code{\link{mclapply}}

# @return For mclapply, a list of the same length as X and named by X.
#' @rdname radiator_parallel
#' @export
#' @keywords internal

mclapply_win <- function(
  X,
  FUN,
  ...,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE,
  mc.silent = FALSE,
  mc.cores = getOption("mc.cores", 2L),
  mc.cleanup = TRUE,
  mc.allow.recursive = TRUE,
  affinity.list = NULL
) {

  # Create a cluster
  if (is.null(mc.cores)) {
    mc.cores <- parallel::detectCores() - 1
  }
  cl <- parallel::makeCluster(mc.cores)

  # We need to find the names of the loaded packages and export them to cluster
  tryCatch(
    {
      loaded.packages <- c(
        utils::sessionInfo()$basePkgs, #Base packages
        names(utils::sessionInfo()$otherPkgs) #Additional packages
      )

      #Export the packages to the clusters
      parallel::clusterExport(cl, 'loaded.packages', envir = environment())

      # parLapply --------------------------------------------------------------
      # Load the libraries on all the clusters
      parallel::parLapply(
        cl,
        1:length(cl),
        function(xx){
          lapply(loaded.packages, function(yy) {
            require(yy , character.only = TRUE)})
        }
      )

      # We want the enclosing environment, not the calling environment
      cluster_export <- function(cl, FUN) {
        env <- environment(FUN)
        while (!identical(env, globalenv())) {
          env <- parent.env(env)
          parallel::clusterExport(cl, ls(all.names = TRUE, envir = env), envir = env)
        }
        parallel::clusterExport(cl, ls(all.names = TRUE, envir = env), envir = env)
      } # End cluster_export

      cluster_export(cl, FUN)

      # Run the lapply in parallel, with a special case for the ... arguments
      if (length(list(...)) == 0) {
        return(parallel::parLapply(cl = cl, X = X, fun = FUN))
      } else {
        return(parallel::parLapply(cl = cl, X = X, fun = FUN, ...))
      }
    }, finally = {
      parallel::stopCluster(cl) #Stop the cluster
    }
  )#End tryCatch
}#End mclapply_win


# radiator_parallel_mc--------------------------------------------------------------
# Overwrite the serial version of mclapply on Windows only
# @title Enable parallel execution on Windows
# @description Internal hack to enable parallel execution of \pkg{assigner}
#' functions on Windows.
# @inheritParams parallel::mclapply
#' @return For mclapply, a list of the same length as X and named by X.
#' @rdname radiator_parallel
#' @param X see \pkg{parallel} \code{\link{mclapply}}
#' @param FUN see \pkg{parallel} \code{\link{mclapply}}
#' @param ... see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.preschedule see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.set.seed see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.silent see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.cores see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.cleanup see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.allow.recursive see \pkg{parallel} \code{\link{mclapply}}
#' @param affinity.list see \pkg{parallel} \code{\link{mclapply}}
#' @keywords internal
#' @export
radiator_parallel_mc <- switch(
  Sys.info()[['sysname']],
  Windows = {mclapply_win},
  Linux   = {parallel::mclapply},
  Darwin  = {parallel::mclapply}
)

# radiator_parallel with progress bar -------------------------------------------
# Overwrite the serial version of mclapply on Windows only
# @title Enable parallel execution on Windows
# @description Internal hack to enable parallel execution of \pkg{assigner}
#' functions on Windows.
# @inheritParams parallel::mclapply
#' @return For mclapply, a list of the same length as X and named by X.
#' @rdname radiator_parallel
#' @param X see \pkg{parallel} \code{\link{mclapply}}
#' @param FUN see \pkg{parallel} \code{\link{mclapply}}
#' @param ... see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.preschedule see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.set.seed see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.silent see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.cores see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.cleanup see \pkg{parallel} \code{\link{mclapply}}
#' @param mc.allow.recursive see \pkg{parallel} \code{\link{mclapply}}
#' @param affinity.list see \pkg{parallel} \code{\link{mclapply}}
#' @keywords internal
#' @export
radiator_parallel <- switch(
  Sys.info()[['sysname']],
  Windows = {mclapply_win},
  Linux   = {pbmcapply::pbmclapply},
  Darwin  = {pbmcapply::pbmclapply}
)
#End radiator_parallel with progress bar

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
  parallel.core.opt
  return(parallel.core.opt)
}#End parallel_core_opt


