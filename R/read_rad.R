# read .rad file

#' @name read_rad
#' @title Read radiator file ending \code{.gds, .rad, .gds.rad, .arrow.parquet}.
#' @description Fast read of \code{.gds, .rad, .gds.rad, .arrow.parquet} files.

#' The function uses \code{\link[arrow]{read_parquet}} or
#' CoreArray Genomic Data Structure (\href{https://github.com/zhengxwen/gdsfmt}{GDS})
#' file system.
#'
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.


#' @param data A file in the working directory ending with .arrow.parquet or .gds,
#' and produced by radiator, assigner or grur.
#' @param columns (optional) For arrow.parquet file. Column names to read. The default is to read all all columns.
#' Default: \code{columns = NULL}.
#' @param allow.dup (optional, logical) To allow the opening of a GDS file with
#' read-only mode when it has been opened in the same R session.
#' Default: \code{allow.dup = FALSE}.
#' @param check (optional, logical) Verify that GDS number of samples and markers
#' match.
#' Default: \code{check = TRUE}.
#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @details For GDS file system, \strong{read_rad} will open the GDS connection file
#' set the filters (variants and samples) based on the info found in the file.

#' @return A radiator tidy data frame or
#' GDS object (with read/write permissions) in the global environment.
#' @export
#' @rdname read_rad
#' @seealso
#' \href{https://github.com/apache/arrow/}{arrow}
#' \href{https://github.com/zhengxwen/gdsfmt}{GDS}
#' \code{\link{read_rad}}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @examples
#' \dontrun{
#' require(SeqArray)
#' shark <- radiator::read_rad(data = "data.shark.gds")
#' turtle <- radiator::read_rad(data = "data.turtle.arrow.parquet")
#' }

read_rad <- function(
    data,
    columns = NULL,
    allow.dup = FALSE,
    check = TRUE,
    verbose = FALSE
) {

  ## TEST
  # columns = NULL
  # from = 1
  # to = NULL
  # as.data.table = FALSE
  # old.format = FALSE

  # detect format---------------------------------------------------------------
  data.type <- detect_genomic_format(data)

  # arrow parquet file ---------------------------------------------------------
  if ("arrow.parquet" %in% data.type) {
    data <- arrow::read_parquet(file = data, col_select = columns)
    message("If the dataset is not readable, use:")
    message("obj <- arrow::read_parquet(file = 'your.file.arrow.parquet')")
    return(data)
  }

  # FST FILE -------------------------------------------------------------------
  if ("fst.file" %in% data.type) {

    # since version 0.8.4 there is a distinction between old and new format...
    # Catch error while reading
    read_rad_fst <- function(
    data) {
      message("This file format will be deprecated soon:")
      message("    Save to another format (e.g. parquet) or ")
      message("    Use the package fst to open and save in the future")
      fst::read_fst(
        path = data,
        columns = NULL,
        from = 1,
        to = 2,
        as.data.table = FALSE,
        old_format = FALSE
      ) %>%
        tibble::as_tibble(.)
    }


    safe_rad_fst <- purrr::safely(.f = read_rad_fst)
    data.safe <- safe_rad_fst(data)
    if (is.null(data.safe$error)) {
      # return(data.safe$result)
      data <- fst::read_fst(
        path = data,
        columns = columns,
        from = from,
        to = to,
        as.data.table = as.data.table
      ) %>%
        tibble::as_tibble(.)
      return(data)
    } else {
      data.old <- suppressWarnings(fst::read_fst(path = data, old_format = TRUE)) %>%
        radiator::write_rad(data = ., path = data)
      data.old <- NULL
      data <- fst::read_fst(
        path = data,
        columns = columns,
        old_format = old.format) %>%
        tibble::as_tibble(.)
      message("\nThis .rad file was created with an earlier version of the fst package")
      message("A new version with the same name was written")
    }
    return(data)
  }#End fst.file

  # GDS file -------------------------------------------------------------------
  if ("gds.file" %in% data.type) {
    if (verbose) message("Opening GDS file connection")

    # if (!is.character(gds)) data <- data$filename


    seq_open_temp <- function(data, allow.dup) {
      SeqArray::seqOpen(gds.fn = data, readonly = allow.dup, allow.duplicate = allow.dup)
    }#End seq_open_temp

    safe_seq_open <- purrr::safely(.f = seq_open_temp)

    data.safe <- safe_seq_open(data, allow.dup)

    if (is.null(data.safe$error)) {
      data <- data.safe$result
    } else {
      temp <- gdsfmt::openfn.gds(filename = data, readonly = FALSE, allow.fork = FALSE, allow.duplicate = TRUE)
      if (gdsfmt::get.attr.gdsn(temp$root)$FileFormat == "SNP_ARRAY") {
        gdsfmt::closefn.gds(gdsfile = temp)
        radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)
        temp.name <- stringi::stri_join("radiator_temp_", data)
        message("The input file is a SNP GDS file, conversion using SeqArray::seqSNP2GDS")
        SeqArray::seqSNP2GDS(data, temp.name)
        data <- SeqArray::seqOpen(gds.fn = temp.name, readonly = allow.dup, allow.duplicate = allow.dup)
        message("SeqArray GDS file generated: ", temp.name)
      }
    }

    # if (allow.dup) {
    #   data <- SeqArray::seqOpen(gds.fn = data, readonly = TRUE, allow.duplicate = TRUE)
    # } else {
    #   data <- SeqArray::seqOpen(gds.fn = data, readonly = FALSE)
    # }

    s <- extract_individuals_metadata(
      gds = data,
      ind.field.select = "INDIVIDUALS",
      whitelist = TRUE
    ) %$%
      INDIVIDUALS
    m <- extract_markers_metadata(
      gds = data,
      markers.meta.select = "VARIANT_ID",
      whitelist = TRUE
    ) %$%
      VARIANT_ID
    if (verbose) message("Setting filters to:")
    if (verbose) message("    number of samples: ", length(s))
    if (verbose) message("    number of markers: ", length(m))

    SeqArray::seqSetFilter(object = data,
                           variant.id = m,
                           sample.id = as.character(s),
                           verbose = FALSE)

    # Checks--------------------------------------------------------------------
    if (check) {
      check <- SeqArray::seqGetFilter(data)
      if (length(check$sample.sel[check$sample.sel]) != length(s)) {
        rlang::abort("Number of samples don't match, contact author")
      }
      if (length(check$variant.sel[check$variant.sel]) != length(m)) {
        rlang::abort("Number of markers don't match, contact author")
      }
    }
    return(data)
  }#End gds.file
}#End read_rad


#' @name write_rad
#' @title Write tidy genomic data file or close GDS file
#' @description When datasets are a tidy data frame, the function provides a
#' fast way to write a \code{.arrow.parquet} file. It superseded the
#' \code{.rad} file format that is essentially the \code{.fst} format provided by the
#' the package \href{https://github.com/fstpackage/fst}{fst}.
#' The \code{.fst} end was replaced by \code{.rad} to remove the confusion
#' with population genetics statistic fst ...
#' The decision to use arrow parquet format from Apache was taken because the
#' package is easier to install than \code{fst}, write and read files faster and
#' files sizes are are also smaller.
#' When the object is a CoreArray Genomic Data Structure
#' (\href{https://github.com/zhengxwen/gdsfmt}{GDS}) file system, the function
#' set filters (variants and samples) based on the info found in the file and
#' close the connection with the GDS file.
#'
#'
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data An object in the global environment: tidy data frame or GDS connection file

#' @param path (optional) For tidy data frame, the path to write the data on disk.
#' Default: \code{path = NULL}.

#' @param filename (optional) Name of the file (when using tsv files).
#' Default: \code{filename = NULL}.

#' @param tsv (optinal, logical) To trigger saving using \code{readr::write_tsv}.
#' Default: \code{tsv = FALSE}.

#' @param internal (optional, logical) This is used inside radiator internal code and it stops
#' from writting the file.
#' Default: \code{internal = FALSE}.

#' @param append (optional, logical) If \code{FALSE}, will overwrite existing file.
#' If \code{TRUE}, will append to existing file.
#' In both cases, if file does not exist a new file is created.
#' Default: \code{append = FALSE}.

#' @param col.names (optional, logical) Write columns names at the top of the file?
#' Must be either TRUE or FALSE.
#' Default: \code{col.names = TRUE}.

#' @param write.message (optional, character) Print a message in the console
#' after writting file.
#' With \code{write.message = NULL}, nothing is printed in the console.
#' Default: \code{write.message = "standard"}. This will print
#' \code{message("File written: ", folder_short(filename))}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @return A file written in the working directory or nothing if it's a GDS connection file.
#' @export
#' @rdname write_rad
#'
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @seealso
#' \href{https://github.com/fstpackage/fst}{fst}
#' \href{https://github.com/zhengxwen/gdsfmt}{GDS}
#' \code{\link{read_rad}}
#'
#' @examples
#' \dontrun{
#' require(SeqArray)
#' radiator::write_rad(data = tidy.data, path = "data.shark.rad")
#' radiator::write_rad(data = gds.object)
#' }


write_rad <- function(
    data,
    path = NULL,
    filename = NULL,
    tsv = FALSE,
    internal = FALSE,
    append = FALSE,
    col.names = TRUE,
    write.message = "standard",
    verbose = FALSE
) {

  if (!internal) {
    if (tsv) {
      # write the dots file
      if (!is.null(path)) {
        path.filename <- file.path(path, filename)
      } else {
        path.filename <- filename
      }
      readr::write_tsv(x = data, file = path.filename, append = append, col_names = col.names)

      if (!is.null(write.message) && verbose) {
        if (write.message == "standard") {
          message("File written: ", folder_short(filename))
        } else {
          message(write.message,  folder_short(filename))
        }
      }

    } else {
      # detect format---------------------------------------------------------------
      data.type <- class(data)

      if ("SeqVarGDSClass" %in% data.type) {
        s <- extract_individuals_metadata(
          gds = data,
          ind.field.select = "INDIVIDUALS",
          whitelist = TRUE
        ) %$%
          INDIVIDUALS
        m <- extract_markers_metadata(
          gds = data,
          markers.meta.select = "VARIANT_ID",
          whitelist = TRUE
        ) %$%
          VARIANT_ID
        if (verbose) message("Setting filters to:")
        if (verbose) message("    number of samples: ", length(s))
        if (verbose) message("    number of markers: ", length(m))
        SeqArray::seqSetFilter(object = data,
                               variant.id = m,
                               sample.id = as.character(s),
                               verbose = FALSE)
        if (verbose) message("Closing connection with GDS file:\n", data$filename)
        gds.filename <- data$filename
        SeqArray::seqClose(data)
        return(gds.filename)
      } else {
        if (is.null(path)) rlang::abort("The function requires the path of the file")
        if (verbose) cli::cli_progress_step(msg = "Writing arrow parquet tidy dataset...")
        tibble::as_tibble(data) %>%
          arrow::write_parquet(x = ., sink = path)
        if (verbose) cli::cli_progress_done()
        # tibble::as_tibble(data) %>%
        #   fst::write_fst(x = ., path = path, compress = 85)
        if (!is.null(write.message) && verbose) {
          if (write.message == "standard") {
            message("File written: ", folder_short(filename))
          } else {
            write.message
          }
        }
      }
    }
  }
}#End write_rad
