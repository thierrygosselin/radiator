# write a LDna object from a tidy data frame

#' @name ld_long_distance
#'
#' @title Generate SNP long distance Linkage Disequilibrium statistics
#'
#' @description Generate SNP long distance Linkage Disequilibrium statistics
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.
#' \strong{The genotypes are biallelic.}


#' @param filename (optional) TTO DO
#' he file name of the LDna lower matrix file.
#' Radiator will append \code{.ldna.rds} to the filename.
#' If filename chosen is already present in the
#' working directory, the default \code{radiator_datetime.ldna.rds} is chosen.
#' With default, \code{filename = NULL}, no file is generated, only an object in
#' the Global Environment.
#' To read the data back into R, use readRDS("filename.ldna.rds").

#' @inheritParams tidy_genomic_data

#' @param ldna (optional, logical) Output an LDna lower matrix object.
#' Default: \code{ldna = FALSE}.
#'
#' @param ld.threshold (optional, double) The threshold to prune SNP based on
#' Long Distance Threshold.
#' Default: \code{ld.threshold = NULL}.

#' @param ... (optional) To pass further argument for fine-tuning the
#' function (see details).

#' @export
#' @rdname ld_long_distance

#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom tibble has_name
#' @importFrom tidyr spread

#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' (2012) A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606

#' @details The function requires \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}
#'
#' To install SNPRelate:
#' source("https://bioconductor.org/biocLite.R")
#' biocLite("SNPRelate").
#'#'
#'
#' \strong{Further arguments passed via the \emph{dots-dots-dots}:}
#' \itemize{
#' \item keep.gds Default \code{keep.gds = FALSE}, the SNPRelate GDS object
#' generated is removed after completion of the LDna object.
#' }
#'

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

ld_long_distance <- function(
  data,
  ld.threshold = NULL,
  ldna = FALSE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  # testing
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1
  # keep.gds <- FALSE
  # ldna <- TRUE
  # ld.threshold = 0.2

  timing <- proc.time()
  opt.change <- getOption("width")
  options(width = 70)

  # Check that snprelate is installed
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop('To install SNPRelate:\n
         source("https://bioconductor.org/biocLite.R")
         biocLite("SNPRelate")')
  }

  res <- list()# to store the output

  # dotslist -------------------------------------------------------------------
  radiator.dots <- list(...)
  want <- c("keep.gds")
  unknowned_param <- setdiff(names(radiator.dots), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  keep.gds <- radiator.dots[["keep.gds"]]

  # useful outside this function
  if (is.null(keep.gds)) keep.gds <- FALSE

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")

  # Filename -------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    write.ld <- FALSE
    filename <- stringi::stri_join("radiator_", file.date, ".ld")
  } else {
    write.ld <- TRUE
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, ".ld")
    } else {
      filename <- stringi::stri_join(filename, ".ld")
    }
  }

  filename.gds <- stringi::stri_join(filename, ".gds")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) {
    data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
  }

  markers <- dplyr::distinct(data, MARKERS)

  # Check if data is biallelic -------------------------------------------------
  biallelic <- radiator::detect_biallelic_markers(data = data)
  if (!biallelic) stop("LDna requires biallelic genotypes")

  # Generating SNPRelate data --------------------------------------------------
  message("Generating SNPRelate data")
  res$data.gds <- radiator::write_snprelate(
    data = data,
    biallelic = TRUE,
    filename = filename,
    verbose = FALSE)
  if (keep.gds) {
    message("SNPRelate GDS file generated: ", filename.gds)
    message("To close the connection use SNPRelate::snpgdsClose(filename)")
  }

  # Compute LD -----------------------------------------------------------------
  message("Computing LD matrix...")
  long.distance.ld <- SNPRelate::snpgdsLDMat(
    gdsobj = res$data.gds,
    snp.id = NULL,
    sample.id = NULL,
    slide = -1,
    mat.trim = FALSE,
    method = "r", #composite and corr option are the same with 0, 1, 2 gt coding
    num.thread = parallel.core,
    with.id = TRUE,
    verbose = FALSE) # %$% LD
  # names(long.distance.ld)
  markers$SNPRELATE <- long.distance.ld$snp.id
  long.distance.ld <- long.distance.ld$LD

  # work on the output ---------------------------------------------------------
  # long.distance.ld <- long.distance.ld^2
  colnames(long.distance.ld) <- rownames(long.distance.ld) <- markers$MARKERS

  # Full matrix
  res$ld.full.matrix <- long.distance.ld

  # Lower matrix
  res$ld.lower.matrix <- long.distance.ld
  res$ld.lower.matrix[upper.tri(res$ld.lower.matrix, diag = TRUE)] <- rlang::na_dbl

  if (ldna) {
    res$ldna <- res$ld.lower.matrix^2
  }

  # Upper matrix
  res$ld.upper.matrix <- long.distance.ld
  res$ld.upper.matrix[lower.tri(res$ld.upper.matrix, diag = TRUE)] <- rlang::na_dbl
  long.distance.ld <- NULL

  # LD tibble
  # Use the upper or lower matrix as input
  res$ld.tibble <- ld2df(x = res$ld.upper.matrix)

  # stats ----------------------------------------------------------------------


  # OUTLIERS_LOW = Q25 - (1.5 * IQR)
  # OUTLIERS_HIGH <-  Q75 + (1.5 * IQR)
  # OUTLIERS_LOW_N = length(x[x < OUTLIERS_LOW])
  # OUTLIERS_HIGH_N = length(x[x > OUTLIERS_HIGH])
  # OUTLIERS_TOTAL = OUTLIERS_HIGH_N + OUTLIERS_LOW_N
  # OUTLIERS_PROP = round(OUTLIERS_TOTAL / length(x), 3)

  res$pb.ld <- tibble::tibble(
    x = 1,
    MIN = min(res$ld.tibble$LD, na.rm = TRUE),
    Q25 = stats::quantile(res$ld.tibble$LD, 0.25, na.rm = TRUE),
    MEDIAN = stats::median(res$ld.tibble$LD, na.rm = TRUE),
    Q75 = stats::quantile(res$ld.tibble$LD, 0.75, na.rm = TRUE),
    MAX = max(res$ld.tibble$LD, na.rm = TRUE),
    IQR = stats::IQR(res$ld.tibble$LD, na.rm = TRUE),
    OUTLIERS_LOW = Q25 - (1.5 * IQR),
    OUTLIERS_HIGH =  Q75 + (1.5 * IQR)
  )

  if (res$pb.ld$OUTLIERS_LOW < 0) res$pb.ld$OUTLIERS_LOW <- res$pb.ld$MIN

  element.text <- ggplot2::element_text(size = 10,
                                        family = "Helvetica", face = "bold")

  res$ld.boxplot <- ggplot2::ggplot(data = res$pb.ld, ggplot2::aes(x)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(ymin = MIN, lower = Q25, middle = MEDIAN, upper = Q75, ymax = MAX), stat = "identity") +
    ggplot2::labs(
      y = "Long distance linkage disequilibrium",
      title = "Markers long distance linkage disequilibrium (LD)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = element.text,
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  print(res$ld.boxplot)
  suppressMessages(ggplot2::ggsave(
    filename = "ld.boxplot.pdf",
    # filename = file.path(path.folder.coverage, "plot.coverage.boxplot.pdf"),
    plot = res$ld.boxplot,
    width = 15, height = 20,
    dpi = 300, units = "cm", useDingbats = FALSE))
  # ld.boxplot <- NULL

  res$ld.no.outliers.tibble <- dplyr::filter(res$ld.tibble, LD < res$pb.ld$OUTLIERS_HIGH & LD > res$pb.ld$OUTLIERS_LOW)

  res$pb.ld.no.outliers <- tibble::tibble(
    x = 1,
    MIN = min(res$ld.no.outliers.tibble$LD, na.rm = TRUE),
    Q25 = stats::quantile(res$ld.no.outliers.tibble$LD, 0.25, na.rm = TRUE),
    MEDIAN = stats::median(res$ld.no.outliers.tibble$LD, na.rm = TRUE),
    Q75 = stats::quantile(res$ld.no.outliers.tibble$LD, 0.75, na.rm = TRUE),
    MAX = max(res$ld.no.outliers.tibble$LD, na.rm = TRUE),
    IQR = stats::IQR(res$ld.no.outliers.tibble$LD, na.rm = TRUE),
    OUTLIERS_LOW = Q25 - (1.5 * IQR),
    OUTLIERS_HIGH =  Q75 + (1.5 * IQR)
  )
  if (res$pb.ld.no.outliers$OUTLIERS_LOW < 0) res$pb.ld.no.outliers$OUTLIERS_LOW <- res$pb.ld.no.outliers$MIN

  element.text <- ggplot2::element_text(size = 10,
                                        family = "Helvetica", face = "bold")

  res$ld.boxplot.no.outliers <- ggplot2::ggplot(data = res$pb.ld.no.outliers, ggplot2::aes(x)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(ymin = MIN, lower = Q25, middle = MEDIAN, upper = Q75, ymax = MAX), stat = "identity") +
    ggplot2::labs(
      y = "Long distance linkage disequilibrium",
      title = "Markers long distance linkage disequilibrium (LD)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, family = "Helvetica", face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = element.text,
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  print(res$ld.boxplot.no.outliers)
  suppressMessages(ggplot2::ggsave(
    filename = "ld.no.outliers.boxplot.pdf",
    # filename = file.path(path.folder.coverage, "plot.coverage.boxplot.pdf"),
    plot = res$ld.boxplot.no.outliers,
    width = 15, height = 20,
    dpi = 300, units = "cm", useDingbats = FALSE))

  # Pruning --------------------------------------------------------------------

  if (!is.null(ld.threshold)){
    # ld.threshold <- 0.12
    long.ld.pruning <- SNPRelate::snpgdsLDpruning(
      gdsobj = res$data.gds,
      autosome.only = FALSE,
      remove.monosnp = TRUE,
      maf = NaN, missing.rate = NaN,
      method = "r",
      ld.threshold = ld.threshold,
      num.thread = 1,
      #Warning message:
      # In SNPRelate::snpgdsLDpruning(gdsobj = data, autosome.only = FALSE,  :
      # The current version of 'snpgdsLDpruning()' does not support multi-threading.
      verbose = FALSE)
  }

  # length(long.ld.pruning$chrun)
  # remove SNPRelate GDS object -------------------------------------------------
  if (!keep.gds) {
    res$data.gds <- NULL
    message("Removing SNPRelate GDS file")
    if (file.exists(filename.gds)) file.remove(filename.gds)
  }

  res$markers.ld.info <- markers
  markers <- NULL

  if (unique(stringi::stri_detect_fixed(
    str = sample(x = res$markers.ld.info$MARKERS, size = 10), pattern = "__"))) {
    res$markers.ld.info <- tidyr::separate(
      data = res$markers.ld.info,
      col = MARKERS, into = c("CHROM", "LOCUS", "POS"),
      sep = "__", remove = FALSE)
  }

  if (!is.null(ld.threshold)){
    res$markers.ld.info <- res$markers.ld.info %>%
      dplyr::mutate(
        LONG_DISTANCE_PRUNING = dplyr::if_else(
          SNPRELATE %in% long.ld.pruning$chrun, "whitelist", "blacklist"))
    long.ld.pruning <- NULL

    res$whitelist.snp.ld <- dplyr::filter(res$markers.ld.info,
                                          LONG_DISTANCE_PRUNING == "whitelist") %>%
      dplyr::select(dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS"))) %>%
      readr::write_tsv(x = ., path = "whitelist.snp.long.dist.ld.tsv")

    message("Number of SNPs after pruning for long distance LD: ",
            nrow(res$whitelist.snp.ld))


    res$blacklist.snp.ld <- dplyr::filter(res$markers.ld.info,
                                          LONG_DISTANCE_PRUNING == "blacklist") %>%
      dplyr::select(dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS"))) %>%
      readr::write_tsv(x = ., path = "blacklist.snp.long.dist.ld.tsv")

    message("Number of prunned SNPs based on long distance LD: ",
            nrow(res$blacklist.snp.ld))

    if (nrow(res$blacklist.snp.ld) > 0) {
      res$data <- dplyr::filter(data, MARKERS %in% res$whitelist.snp.ld$MARKERS)
    }
  }
  data <- NULL

  res$manhattan.plot.ld <- res$ld.tibble %>%
    dplyr::mutate(X = "1") %>%
    dplyr::group_by(LD) %>%
    dplyr::tally(.) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(X = "1") %>%
    dplyr::ungroup(.) %>%
    ggplot2::ggplot(data = .,
                    ggplot2::aes(x = X, y = LD, size = n)) +
    ggplot2::geom_jitter(alpha = 0.3, stat = "identity") +
    ggplot2::labs(y = "LD") +
    ggplot2::scale_size_area(name = "Number of marker", max_size = 6) +
    ggplot2::theme_light() +
    ggplot2::theme(
      # legend.position = "none",
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      # panel.grid.major.y = element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica")
    )
  # res$manhattan.plot.ld
  suppressMessages(ggplot2::ggsave(
    filename = "manhattan.plot.ld.pdf",
    plot = res$manhattan.plot.ld,
    width = 15, height = 20,
    dpi = 300, units = "cm", useDingbats = FALSE))

  message("Computation time: ", round((proc.time() - timing)[[3]]), " sec")
  options(width = opt.change)
  return(res)
}#End ld long distance

# Internal nested functions: ---------------------------------------------------

# melt the LD matrice into a data frame --------------------------------------
#' @title ld2df
#' @description melt the LD matrice into a data frame
#' @rdname ld2df
#' @export
#' @keywords internal
ld2df <- function(x) {
  # x <- ld.upper.matrix
  x <- as.matrix(x)
  # diag(x) <- NA
  # x[lower.tri(x)] <- NA
  x <- dplyr::bind_cols(tibble::data_frame(MARKERS_A = rownames(x)),
                        tibble::as_data_frame(x)) %>%
    data.table::as.data.table(.) %>%
    data.table::melt.data.table(
      data = ., id.vars = "MARKERS_A", variable.name = "MARKERS_B", value.name = "LD",
      variable.factor = FALSE) %>%
    tibble::as_data_frame(.) %>%
    dplyr::filter(!is.na(LD)) %>%
    dplyr::arrange(dplyr::desc(LD))
  return(x)
}#End distance2df
