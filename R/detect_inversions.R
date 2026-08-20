#' @title Detect candidate inversions from a VCF using local PCA
#'
#' @description
#' Detect putative chromosomal inversions from SNP data stored in a VCF file,
#' using:
#' \itemize{
#'   \item local PCA on SNP windows with \pkg{lostruct} (\code{eigen_windows}
#'   + \code{pc_dist} + \code{cmdscale}),
#'   \item PCA clustering of individuals (three-class karyotype pattern),
#'   \item heterozygosity profiles across candidate windows,
#'   \item LD structure (mean/median r\eqn{^2}).
#' }
#'
#' This function assumes that \strong{all SNP filtering is done upstream}.
#' Only minimal sanity checks (e.g., removal of monomorphic loci) are
#' performed internally.
#'
#' For each candidate window, three plots are generated and saved to disk:
#' \itemize{
#'   \item PCA of individuals (PC1 vs PC2, coloured by cluster),
#'   \item heterozygosity by cluster,
#'   \item LD heatmap (r\eqn{^2}).
#' }
#' The plots are saved in the folder specified by \code{plot.folder}.
#'
#' @rdname detect_inversions
#'
#' @param vcf (character or \code{vcfR}) Path to a pre-filtered VCF file,
#' or a \code{vcfR} object (from \code{vcfR::read.vcfR()}).
#'
#' @param window.snps (integer) Number of SNPs per window for the local PCA.
#' Default: \code{window.snps = 100L}.
#'
#' @param k.clust (integer) Number of clusters in PCA space. For inversion
#' detection, \code{3} is typical (two homokaryotypes + one heterokaryotype).
#' Default: \code{k.clust = 3L}.
#'
#' @param lostruct.k (integer) Number of axes for the window embedding
#' (MDS space). Default: \code{lostruct.k = 2L}.
#'
#' @param lostruct.outlier.quantile (double) Quantile threshold used to detect
#' outlier windows in lostruct space.
#' Default: \code{lostruct.outlier.quantile = 0.99}.
#'
#' @param ld.method (character) Correlation method used by \code{stats::cor()}
#' to compute LD (r\eqn{^2}).
#' Default: \code{ld.method = "pearson"}.
#'
#' @param parallel.core (integer) Number of cores used by \pkg{lostruct}
#' when computing eigenvectors in windows.
#' Default: \code{parallel.core = 1L}.


#' @inheritParams radiator_common_arguments

#' @return A list of class \code{detect_inversions} containing:
#' \itemize{
#'   \item \code{summary}: tibble of candidate windows,
#'   \item \code{windows}: tibble with window coordinates in lostruct space,
#'   \item \code{results}: list of per-window objects (PCA, clusters,
#'     heterozygosity, LD summaries and LD matrix).
#' }
#'
#' @examples
#' \dontrun{
#' inv <- detect_inversions(
#'   vcf = "filtered_tuna.vcf.gz",
#'   window.snps = 150L,
#'   lostruct.outlier.quantile = 0.995,
#'   parallel.core = 8L
#' )
#'
#' dplyr::filter(inv$summary, support == "strong")
#' }
#'
#' @export
detect_inversions <- function(
    vcf,
    window.snps = 100L,
    k.clust = 3L,
    lostruct.k = 2L,
    lostruct.outlier.quantile = 0.99,
    ld.method = "pearson",
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE
) {


  # no need for this if you look at radiator DESCRIPTION, ggplot2 will be imported...
  # and again... requireNamespace SHOULD NOT be used inside function packages...
  # if (save.plots && !requireNamespace("ggplot2", quietly = TRUE)) {
  #   stop("Package `ggplot2` is required when `save.plots = TRUE`.")
  # }

  #--------------------------------------------------------------------------#
  # 0. Read VCF
  #--------------------------------------------------------------------------#
  if (inherits(vcf, "vcfR")) {
    vcf_r <- vcf
  } else if (is.character(vcf) && file.exists(vcf)) {
    vcf_r <- vcfR::read.vcfR(vcf, verbose = FALSE)
  } else {
    stop("`vcf` must be a vcfR object or a valid VCF path.")
  }

  fix <- tibble::as_tibble(vcf_r@fix)
  chrom <- fix$CHROM
  pos   <- as.numeric(fix$POS)

  #--------------------------------------------------------------------------#
  # 1. Build genotype matrix (0/1/2), individuals × SNPs
  #--------------------------------------------------------------------------#
  gt <- vcfR::extract.gt(vcf_r, element = "GT", as.numeric = FALSE)

  gt_to_numeric <- function(x) {
    x <- gsub("\\|", "/", x)
    x[x %in% c("./.", ".|.", NA)] <- NA
    dplyr::case_when(
      x == "0/0" ~ 0,
      x %in% c("0/1", "1/0") ~ 1,
      x == "1/1" ~ 2,
      TRUE ~ NA_real_
    )
  }

  # variants × samples (this is exactly what eigen_windows wants)
  snps_ls <- t(apply(gt, 1L, gt_to_numeric))
  rownames(snps_ls) <- rownames(gt)
  colnames(snps_ls) <- colnames(gt)

  # individuals × SNPs for PCA/het/LD downstream
  geno_mat <- t(snps_ls)

  indiv  <- rownames(geno_mat)        # sample IDs
  n_ind  <- length(indiv)
  n_snps <- ncol(geno_mat)

  message("VCF: ", n_ind, " individuals, ", n_snps, " SNPs after parsing GT.")

  # Remove monomorphic SNPs (internal safety)
  # geno_mat: individuals (rows) × SNPs (cols)
  is_var <- apply(
    X = geno_mat,
    MARGIN = 2L,
    FUN = function(x) length(unique(x[!is.na(x)])) > 1L
  )

  if (!all(is_var)) {
    message("Dropping ", sum(!is_var), " monomorphic SNPs.")
    geno_mat <- geno_mat[, is_var, drop = FALSE]

    # keep chrom/pos in sync with SNPs
    chrom <- chrom[is_var]
    pos   <- pos[is_var]

    # also subset the SNPs×samples matrix for lostruct
    snps_ls <- snps_ls[is_var, , drop = FALSE]
  }

  n_snps <- ncol(geno_mat)

  #--------------------------------------------------------------------------#
  # 2. Define windows for lostruct (full windows only)
  #--------------------------------------------------------------------------#
  n_windows <- floor(n_snps / window.snps)
  if (n_windows < 1L) {
    stop("Not enough SNPs to form a single full window. ",
         "Increase SNP density or reduce `window.snps`.")
  }

  snp_used_idx <- seq_len(n_windows * window.snps)

  geno_mat_used <- geno_mat[, snp_used_idx, drop = FALSE]
  chrom_used    <- chrom[snp_used_idx]
  pos_used      <- pos[snp_used_idx]

  window_id <- rep(seq_len(n_windows), each = window.snps)

  #--------------------------------------------------------------------------#
  # 3. lostruct local PCA: eigen_windows + pc_dist + cmdscale
  #--------------------------------------------------------------------------#
  message("Running lostruct eigen_windows on ", n_windows,
          " windows of ", window.snps, " SNPs...")

  # lostruct expects matrix with variants in rows, samples in columns
  snps_ls_used <- snps_ls[snp_used_idx, , drop = FALSE]

  pcs <- lostruct::eigen_windows(
    snps_ls_used,
    win = window.snps,
    k   = lostruct.k
  )

  # distance between windows in PC-space
  windist <- lostruct::pc_dist(
    pcs,
    npc = lostruct.k
  )

  # Replace NA distances with a large value = "very far"
  max_d <- max(windist, na.rm = TRUE)
  windist[is.na(windist)] <- max_d * 1.1

  # Identify windows with any NA distance
  na_row <- rowSums(is.na(windist)) > 0
  na_col <- colSums(is.na(windist)) > 0
  bad_windows <- which(na_row | na_col)

  if (length(bad_windows) > 0L) {
    message("Dropping ", length(bad_windows),
            " windows with NA distances before MDS.")
  }

  keep_windows <- setdiff(seq_len(nrow(windist)), bad_windows)

  if (length(keep_windows) < 2L) {
    stop("Too few windows with valid distances after removing NA rows/cols.")
  }

  # MDS embedding of windows
  mds <- stats::cmdscale(
    d   = windist,
    k   = lostruct.k,
    eig = TRUE
  )

  windows_tbl <- tibble::as_tibble(mds$points) |>
    dplyr::rename_with(~ paste0("LS_PC", seq_len(ncol(mds$points)))) |>
    dplyr::mutate(
      window_id = dplyr::row_number(),
      ls_dist   = sqrt(rowSums(dplyr::across(dplyr::starts_with("LS_PC"))^2))
    )

  thr <- stats::quantile(
    x     = windows_tbl$ls_dist,
    probs = lostruct.outlier.quantile,
    na.rm = TRUE
  )

  candidates <- windows_tbl |>
    dplyr::filter(ls_dist >= thr)

  message("Found ", nrow(candidates), " candidate windows (quantile = ",
          lostruct.outlier.quantile, ").")

  # Prepare plot folder
  if (save.plots) {
    if (!dir.exists(plot.folder)) {
      dir.create(plot.folder, recursive = TRUE, showWarnings = FALSE)
    }
  }

  #--------------------------------------------------------------------------#
  # 4. Per-window analyses
  #--------------------------------------------------------------------------#
  analyze_window <- function(wid) {

    snp_idx <- which(window_id == wid)
    if (length(snp_idx) < 10L) return(NULL)

    g <- geno_mat_used[, snp_idx, drop = FALSE]

    # PCA on individuals
    pca <- stats::prcomp(g, center = TRUE, scale. = TRUE)
    scores <- tibble::as_tibble(pca$x[, 1:2, drop = FALSE]) |>
      dplyr::rename(PC1 = 1, PC2 = 2) |>
      dplyr::mutate(INDIVIDUAL = indiv)

    # k-means clustering
    km <- stats::kmeans(scores$PC1, centers = k.clust, nstart = 20L)
    scores <- scores |>
      dplyr::mutate(cluster = factor(km$cluster))

    # Heterozygosity
    het <- rowMeans(g == 1, na.rm = TRUE)
    het_tbl <- tibble::tibble(
      INDIVIDUAL    = indiv,
      window_id     = wid,
      heterozygosity = het,
      cluster       = scores$cluster
    )

    # LD matrix and summaries
    cor_mat <- stats::cor(
      x   = g,
      use = "pairwise.complete.obs",
      method = ld.method
    )
    r2_mat  <- cor_mat^2
    r2_vals <- r2_mat[upper.tri(r2_mat)]
    ld_mean   <- mean(r2_vals, na.rm = TRUE)
    ld_median <- stats::median(r2_vals, na.rm = TRUE)

    # Window coordinates (genomic)
    chr_win <- unique(chrom_used[snp_idx]) |> paste(collapse = ",")
    start   <- min(pos_used[snp_idx], na.rm = TRUE)
    end     <- max(pos_used[snp_idx], na.rm = TRUE)

    # Cluster summary
    clust_summary <- het_tbl |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(
        n_ind   = dplyr::n(),
        mean_het = mean(heterozygosity, na.rm = TRUE),
        .groups  = "drop"
      )

    #--------------------------------------------------------#
    # Plots: PCA, heterozygosity, LD — save to plot.folder
    #--------------------------------------------------------#
    if (save.plots) {

      # PCA plot
      p_pca <- ggplot2::ggplot(
        data = scores,
        ggplot2::aes(x = PC1, y = PC2, colour = cluster)
      ) +
        ggplot2::geom_point(alpha = 0.9, size = 2) +
        ggplot2::labs(
          title = paste0("Window ", wid, " (", chr_win, ":", start, "-", end, ")"),
          subtitle = "PCA of individuals",
          colour = "Cluster"
        ) +
        ggplot2::theme_bw()

      ggplot2::ggsave(
        filename = file.path(plot.folder, paste0("window-", wid, "_pca.png")),
        plot = p_pca,
        width = 6, height = 5, dpi = 300
      )

      # Heterozygosity by cluster
      het_plot_data <- het_tbl |>
        dplyr::left_join(
          scores |> dplyr::select(INDIVIDUAL, cluster),
          by = "INDIVIDUAL"
        )

      p_het <- ggplot2::ggplot(
        data = het_plot_data,
        ggplot2::aes(x = cluster, y = heterozygosity, fill = cluster)
      ) +
        ggplot2::geom_boxplot(alpha = 0.7) +
        ggplot2::labs(
          title = "Heterozygosity by cluster",
          x = "Cluster",
          y = "Heterozygosity"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")

      ggplot2::ggsave(
        filename = file.path(plot.folder, paste0("window-", wid, "_heterozygosity.png")),
        plot = p_het,
        width = 4, height = 4, dpi = 300
      )

      # LD heatmap
      snp_ids <- seq_len(nrow(r2_mat))
      ld_df <- expand.grid(
        SNP1 = snp_ids,
        SNP2 = snp_ids
      )
      ld_df$r2 <- as.vector(r2_mat)

      p_ld <- ggplot2::ggplot(
        data = ld_df,
        ggplot2::aes(x = SNP1, y = SNP2, fill = r2)
      ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(option = "magma", direction = -1) +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title = "LD heatmap (r²)",
          x = "SNP index",
          y = "SNP index",
          fill = "r²"
        ) +
        ggplot2::theme_bw()

      ggplot2::ggsave(
        filename = file.path(plot.folder, paste0("window-", wid, "_ld.png")),
        plot = p_ld,
        width = 5, height = 5, dpi = 300
      )
    }

    list(
      window_id       = wid,
      chrom           = chr_win,
      start_pos       = start,
      end_pos         = end,
      n_snps          = length(snp_idx),
      pca_scores      = scores,
      het             = het_tbl,
      cluster_summary = clust_summary,
      ld_mean         = ld_mean,
      ld_median       = ld_median,
      ld_matrix       = r2_mat
    )
  }

  results <- candidates$window_id |>
    purrr::map(analyze_window) |>
    purrr::compact()

  #--------------------------------------------------------------------------#
  # 5. Summary tibble
  #--------------------------------------------------------------------------#
  summary_tbl <- results |>
    purrr::map_dfr(function(x) {
      tibble::tibble(
        window_id = x$window_id,
        CHROM     = x$chrom,
        START_POS = x$start_pos,
        END_POS   = x$end_pos,
        N_SNPS    = x$n_snps,
        LD_MEAN   = x$ld_mean,
        LD_MEDIAN = x$ld_median
      )
    }) |>
    dplyr::left_join(windows_tbl, by = "window_id") |>
    dplyr::mutate(
      support = dplyr::case_when(
        LD_MEAN > 0.3 & ls_dist >= thr * 1.1 ~ "strong",
        LD_MEAN > 0.2                        ~ "medium",
        TRUE                                 ~ "weak"
      )
    )

  out <- list(
    summary = summary_tbl,
    windows = windows_tbl,
    results = results
  )
  class(out) <- c("detect_inversions", class(out))
  out
}
