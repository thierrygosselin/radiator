# bcftools helpers -------------------------------------------------------------

#' @title Check that bcftools is available
#' @description
#' Internal helper that checks whether \code{bcftools} can be executed.
#' It runs \code{bcftools.path --version} and errors if the command fails.
#'
#' @param bcftools.path (character)
#' Path or name of the \code{bcftools} executable.
#' Default: \code{bcftools.path = "bcftools"}.
#'
#' @export
#' @keywords internal
bcftools_require <- function(bcftools.path = "bcftools") {
  ok <- tryCatch(
    {
      res <- sys::exec_internal(bcftools.path, "--version")
      identical(res$status, 0L)
    },
    error = function(e) FALSE
  )

  if (!ok) {
    rlang::abort(
      stringi::stri_join(
        "bcftools not found or not executable. Check `bcftools.path` (\"",
        bcftools.path, "\")."
      )
    )
  }

  invisible(TRUE)
} # end bcftools_require


#' @title Run a bcftools command and log stderr
#' @description
#' Internal helper around \code{sys::exec_internal()} that:
#' \itemize{
#'   \item runs \code{bcftools.path} with the supplied \code{args},
#'   \item appends stderr to a log file (if provided),
#'   \item optionally aborts on non-zero exit status.
#' }
#'
#' @param bcftools.path (character) Path or name of the \code{bcftools} executable.
#' @param args (character) Vector of arguments passed to bcftools.
#' @param log.file (optional, character) Path to a log file. If not \code{NULL},
#' stderr is appended to this file with a small header.
#' @param label (optional, character) Short label used in the log header.
#' @param verbose (logical) If \code{TRUE}, the constructed command line is printed.
#' Default: \code{verbose = TRUE}.
#' @param fail_on_status (logical) If \code{TRUE} (default), a non-zero exit
#' status triggers an error. If \code{FALSE}, the status is returned invisibly.
#'
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item \code{status} – exit status (integer),
#'   \item \code{stdout} – character scalar with STDOUT,
#'   \item \code{stderr} – character scalar with STDERR.
#' }
#' @export
#' @keywords internal
bcftools_exec <- function(
    bcftools.path,
    args,
    log.file       = NULL,
    label          = NULL,
    verbose        = TRUE,
    fail_on_status = TRUE
) {
  if (verbose) {
    cmd.pretty <- paste(c(bcftools.path, args), collapse = " ")
    message("bcftools cmd", if (!is.null(label)) paste0(" [", label, "]"), ":\n  ", cmd.pretty)
    if (!is.null(log.file)) message("Log \u2192 ", log.file)
  }

  res <- sys::exec_internal(bcftools.path, args)

  stdout_txt <- rawToChar(res$stdout)
  stderr_txt <- rawToChar(res$stderr)

  if (!is.null(log.file) && nzchar(stderr_txt)) {
    cat(
      "## ", if (!is.null(label)) label else "bcftools",
      " @ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      stderr_txt, "\n",
      file   = log.file,
      append = TRUE,
      sep    = ""
    )
  }

  if (fail_on_status && !identical(res$status, 0L)) {
    msg <- paste0(
      "bcftools command failed",
      if (!is.null(label)) paste0(" [", label, "]") else "",
      ".",
      if (!is.null(log.file)) paste0(" See log: ", log.file) else "",
      if (nzchar(stderr_txt)) paste0("\n\nbcftools stderr:\n", stderr_txt) else ""
    )
    rlang::abort(msg)
  }

  invisible(list(
    status = res$status,
    stdout = stdout_txt,
    stderr = stderr_txt
  ))
} # end bcftools_exec
