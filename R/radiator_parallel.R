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

# @return For mclapply, a list of the same length as X and named by X.
#' @importFrom utils sessionInfo
#' @importFrom parallel detectCores makeCluster clusterExport mclapply parLapply stopCluster
#' @rdname radiator_parallel
#' @export
#' @keywords internal

mclapply_win <- function(
  X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
  mc.silent = FALSE, mc.cores = NULL, mc.cleanup = TRUE, mc.allow.recursive = TRUE
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

      # Load the libraries on all the clusters
      parallel::parLapply(
        cl, 1:length(cl), function(xx){
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


# .radiator_parallel_mc--------------------------------------------------------------
# Overwrite the serial version of mclapply on Windows only
# @name .radiator_parallel
# @title Enable parallel execution on Windows
# @description Internal hack to enable parallel execution of \pkg{assigner}
#' functions on Windows.
# @inheritParams parallel::mclapply
#' @return For mclapply, a list of the same length as X and named by X.
# @importFrom parallel detectCores makeCluster clusterExport mclapply parLapply stopCluster
# @importFrom pbmcapply pbmclapply
#' @rdname radiator_parallel
#' @keywords internal
#' @export
.radiator_parallel_mc <- switch(
  Sys.info()[['sysname']],
  Windows = {mclapply_win},
  # Linux   = {mclapply_progress_bar},
  Linux   = {parallel::mclapply},
  # Linux   = {pbmcapply::pbmclapply},
  # Darwin  = {mclapply_progress_bar}
  Darwin  = {parallel::mclapply}
  # Darwin  = {pbmcapply::pbmclapply}
)

# .radiator_parallel with progress bar -------------------------------------------
# Overwrite the serial version of mclapply on Windows only
# @name .radiator_parallel
# @title Enable parallel execution on Windows
# @description Internal hack to enable parallel execution of \pkg{assigner}
#' functions on Windows.
# @inheritParams parallel::mclapply
#' @return For mclapply, a list of the same length as X and named by X.
# @importFrom parallel detectCores makeCluster clusterExport mclapply parLapply stopCluster
#' @importFrom pbmcapply pbmclapply
#' @rdname radiator_parallel
#' @keywords internal
#' @export
.radiator_parallel <- switch(
  Sys.info()[['sysname']],
  Windows = {mclapply_win},
  # Linux   = {parallel::mclapply},
  # Darwin  = {parallel::mclapply}
  Linux   = {pbmcapply::pbmclapply},
  Darwin  = {pbmcapply::pbmclapply}
)
