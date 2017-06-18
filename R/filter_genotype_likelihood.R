#' @name filter_genotype_likelihood
#' @title Genotype likelihood interactive filter
#' @description Filter genotypes and markers based on genotype likelihood using
#' a tidy data set (long format) and any VCF file both with GL information.
#' The function uses 
#' \code{\link[radiator]{tidy_genomic_data}} and 
#' \code{\link[radiator]{tidy_wide}} to load the file.

# Most arguments are inherited from tidy_genomic_data
#' @inheritParams tidy_genomic_data

#' @param interactive.filter (optional, logical) Do you want the filtering session to 
#' be interactive. With default: \code{interactive.filter == TRUE}, the user is 
#' asked to see figures of distribution before making decisions for filtering 
#' with the genotype likelihood.

#' @param gl.approach Character. By \code{gl.approach = "SNP"} or 
#' by \code{gl.approach = "haplotype"}. 
#' The function will consider the SNP or haplotype GL statistics to filter the marker. 
#' Default: \code{gl.approach = "haplotype"}.

#' @param gl.ind.threshold (Integer, optional if interactive session) 
#' Threshold number of individual's genotype likelihood.
#' @param gl.mean.threshold (Integer, optional if interactive session) 
#' Threshold number of mean genotype likelihood. 
#' @param gl.min.threshold (Integer, optional if interactive session) 
#' Threshold number of min genotype likelihood. 
#' @param gl.diff.threshold (Integer, optional if interactive session)
#' Threshold number of diff genotype likelihood,
#' the difference between the max and min GL over a loci/read/haplotype.

#' @param pop.threshold (optional if interactive session)
#' A threshold number: proportion, percentage
#' or fixed number e.g. 0.50, 50 or 5.
#' 
#' @param percent (logical, optional if interactive session). 
#' Is the pop.threshold argument a percentage? TRUE/FALSE.
#' This argument is necessary to distinguish percentage from integer 
#' for population threshold, (e.g. 5 percent or 5 populations).

#' @param filename (optional) Name of the filtered data set, 
#' written to the working directory.


#' @details There is 4 steps in the interactive version:  
#' 
#' 1. gl_individuals_populations: the user is asked to inspect 
#' the genotype likelihood at the individuals and populations levels"), 
#' 
#' 2. blacklist_genotypes: the option is given to blacklist individual genotypes
#' based on low quality GL, 
#' 
#' 3. gl_markers: the user is asked to inspecting the genotype likelihood 
#' at the marker level and 
#' 
#' 4. blacklist_markers: the option is given to blacklist markers genotypes 
#' based on low quality GL found at the marker level.
#' 
#' Using the haplotype approach: The summary statistics are averaged
#' and nested: SNP -> individuals -> population -> loci. e.g. the mean GL is the average
#' genotype likelihood for all individuals of pop x for loci x.
#' The gl.diff.threshold is the difference between the max and min GL found for 
#' a locus. This argument is given to help your spot big difference in GL for 
#' a loci. Markers with small differences have more stability in the estimate.

#' @return With \code{interactive.filter = FALSE}, the function returns a 
#' filtered tidy vcf data frame inside the global environment. 
#' In the working directory 5 files are written:
#' 1. filters_parameters.tsv is updated or created.
#' The filter parameters, values and returned results 
#' (markers numbers or blacklist genotypes) are inside that file. 
#' 2. blacklist.genotypes.gl.tsv, 3. blacklist.markers.gl.tsv, 
#' 4. whitelist.markers.gl.tsv and 
#' 5. vcf.tidy.paralogs.id.gl.tsv (the same filtered tidy vcf data frame).
#' 
#' With \code{interactive.filter = TRUE}, a list with 20 objects is created. 
#' The information range from summary tables to plots. The objects names are 
#' found by using \code{names(list.name)}. The object can be isolated in separate
#' object outside the list by following the example below.


#' @rdname filter_genotype_likelihood
#' @export

#' @import ggplot2
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub
#' @importFrom dplyr select distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join
#' @importFrom readr write_tsv
#' @importFrom tibble data_frame has_name
#' @importFrom tidyr spread


#' @examples
#' \dontrun{
#' If I have 10 populations, the 3 examples below will give the same output results:
#' 
#' beluga.vcf.tidy.gl <- filter_genotype_likelihood(
#' data = "beluga.vcf", 
#' interactive.filter = TRUE,
#' gl.approach = "haplotype", 
#' gl.ind.threshold = 0,
#' gl.mean.threshold = 10, 
#' gl.min.threshold = 0, 
#' gl.diff.threshold = 200, 
#' pop.threshold = 50, 
#' percent = TRUE)
#' 
#' beluga.vcf.tidy.gl <- filter_genotype_likelihood(
#' data = "beluga.vcf", 
#' interactive.filter = FALSE,
#' gl.approach = "haplotype",
#' gl.ind.threshold = 0,
#' gl.mean.threshold = 10, 
#' gl.min.threshold = 0, 
#' gl.diff.threshold = 200, 
#' pop.threshold = 0.5, 
#' percent = FALSE)
#' 
#' beluga.vcf.tidy.gl <- filter_genotype_likelihood(
#' data = "beluga.vcf", 
#' gl.approach = "haplotype", 
#' gl.ind.threshold = 0,
#' gl.mean.threshold = 10, 
#' gl.min.threshold = 0, 
#' gl.diff.threshold = 100, 
#' pop.threshold = 5, 
#' percent = FALSE)
#' 
#' If interactive.filter = TRUE, a list is created and to view the filtered tidy vcf:
#' tidy.data <- beluga.vcf.tidy.gl$tidy.vcf.filtered.gl
#' 
#' Inside the same list, to isolate the blacklist.genotypes:
#' bg <- beluga.vcf.tidy.gl$blacklist.genotypes
#' }


filter_genotype_likelihood <- function(
  data,
  interactive.filter = TRUE,
  strata = NULL, 
  gl.approach = "haplotype",
  gl.ind.threshold = NULL,
  gl.mean.threshold = NULL,
  gl.min.threshold = NULL,
  gl.diff.threshold = NULL,
  pop.threshold = NULL,
  percent = NULL,
  pop.levels = NULL, 
  pop.labels = NULL, 
  pop.select = NULL,
  blacklist.id = NULL, 
  blacklist.genotype = NULL, 
  whitelist.markers = NULL, 
  monomorphic.out = TRUE, 
  max.marker = NULL,
  snp.ld = NULL, 
  common.markers = FALSE,
  filename = NULL
) {
  
  cat("#######################################################################\n")
  cat("################# radiator: filter_genotype_likelihood ##################\n")
  cat("#######################################################################\n")
  
  
  # manage missing arguments -----------------------------------------------------  
  if (missing(data)) stop("missing input file")
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length (number of groups)")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  if (!is.null(pop.select)) {
    pop.select <- stringi::stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  # Message about steps taken during the process ---------------------------------
  if (interactive.filter) {
    message("With the interactive mode turn on, we will go through 4 steps to visualize and filter the data based on genotype likelihood")
    message("Step 1. gl_individuals_populations: Inspecting the genotype likelihood at the individuals and populations levels")
    message("Step 2. blacklist_genotypes: blacklisting (erasing) or not individual genotypes based on low quality GL")
    message("Step 3. gl_markers: Inspecting the genotype likelihood at the marker level")
    message("Step 4. blacklist_markers: Blacklisting (erasing) or not markers genotypes based on low quality GL found at the marker level")
  }
  # Folder -------------------------------------------------------------------
  # Get date and time to have unique filenaming
  file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
  file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
  file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
  
  path.folder <- stringi::stri_join(getwd(),"/", "filter_gl_", file.date, sep = "")
  dir.create(file.path(path.folder))
  
  message(stringi::stri_join("Folder created: ", path.folder))
  file.date <- NULL #unused object
  
  # Filter parameter file ------------------------------------------------------
  message("Parameters used in this run will be store in a file")
  filters.parameters <- list.files(path = getwd(), pattern = "filters_parameters.tsv", full.names = TRUE)
  if (length(filters.parameters) == 0) {
    filters.parameters <- tibble::data_frame(FILTERS = as.character(), PARAMETERS = as.character(), VALUES = as.integer(), BEFORE = as.character(), AFTER = as.character(), BLACKLIST = as.integer(), UNITS = as.character(), COMMENTS = as.character())
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = FALSE, col_names = TRUE)
    message("Created a file to store filters parameters: filters_parameters.tsv")
  } else {
    message("Using the filters parameters file found in the directory")
  }
  
  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)
  
  # import data ----------------------------------------------------------------
  message("Importing data ...")
  input <- radiator::tidy_genomic_data(
    data = data, 
    vcf.metadata = TRUE,
    blacklist.id = blacklist.id, 
    blacklist.genotype = blacklist.genotype, 
    whitelist.markers = whitelist.markers, 
    monomorphic.out = monomorphic.out, 
    max.marker = max.marker,
    snp.ld = snp.ld, 
    common.markers = common.markers,
    strata = strata, 
    pop.levels = pop.levels, 
    pop.labels = pop.labels, 
    pop.select = pop.select,
    filename = NULL
  )
  
  if (!tibble::has_name(input, "GL") & !tibble::has_name(input, "PL") ) {
    stop("GL or PL information is necessary for filter_genotype_likelihood") 
  }
  
  # create a strata.df
  strata.df <- input %>% 
    dplyr::select(INDIVIDUALS, POP_ID) %>% 
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
  # pop number
  pop.number <- n_distinct(input$POP_ID) # get the number of population
  
  # function--------------------------------------------------------------------
  plot_violinplot_genotype_likelihood_individuals <- function(data) {
    plot <- ggplot(data, aes(x = factor(POP_ID), y = GL, na.rm = TRUE)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.1, fill = "black", outlier.colour = NA) +
      stat_summary(fun.y = "mean", geom = "point", shape = 21, size = 2.5, fill = "white") +
      labs(x = "Sampling sites") +
      labs(y = "Genotype likelihood of individuals") +
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5), 
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold")
      )
    return(plot)
  }
  
  
  ## Step 1: gl_individuals_populations ----------------------------------------
  if (interactive.filter) {
    message("Step 1. gl_individuals_populations: Inspecting the genotype likelihood at the individuals and populations levels")
    # path.folder.step1 <- stringi::stri_join(path.folder, "/01_gl_individuals_populations")
    # dir.create(path.folder.step1)
    summary <- summary_genotype_likelihood(
      data = input, 
      pop.levels = pop.levels, 
      gl.approach = gl.approach, 
      folder = path.folder
    )
  }
  
  
  # plot_1: Violin plot GL individuals and pop
  if (interactive.filter) {
    message("Show the violin plot of individuals genotype likelihood (y/n)): ")
    violinplot <- as.character(readLines(n = 1))
    if (violinplot == "y") {
      message("Generating violin plot may take some time...")
      # plot
      genotype.likelihood.violin.plot.individuals <- plot_violinplot_genotype_likelihood_individuals(data = input)
      print(genotype.likelihood.violin.plot.individuals)
      # save
      ggsave(stringi::stri_join(path.folder, "/genotype.likelihood.violin.plot.individuals.pdf"), width = pop.number, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder, "/genotype.likelihood.violin.plot.individuals.png"), width = pop.number, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.individuals) were saved in this directory: ", path.folder.step1))
    }
  }
  
  
  #manhattan plot
  # manhanttan.ind <- ggplot(data = summary$gl.individuals, aes(x = POP_ID, y = GL_MEAN, colour = POP_ID)) + 
  #   geom_jitter() + 
  #   labs(y = "Individual's Mean Genotype Likelihood") +
  #   labs(x = "Populations") +
  #   labs(colour = "Populations") +
  #   # scale_y_continuous(name = , breaks = )
  #   theme(
  #     legend.position = "none",
  #     axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
  #     axis.text.x = element_text(size = 10, family = "Helvetica"),
  #     axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
  #     axis.text.y = element_text(size = 8, family = "Helvetica")
  #   )
  # manhanttan.ind
  
  manhanttan.makers.pop <- ggplot(data = input, aes(x = POP_ID, y = GL, colour = POP_ID)) + 
    geom_jitter() + 
    labs(y = "Individual's Mean Genotype Likelihood") +
    labs(x = "Populations") +
    labs(colour = "Populations") +
    # scale_y_continuous(name = , breaks = )
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_text(size = 10, family = "Helvetica"),
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 8, family = "Helvetica")
    )
  manhanttan.makers.pop
  
  ## Step 2: blacklist_genotypes -------------------------------------------------
  # Erasing genotypes below threshold
  if (interactive.filter) {
    message("Step 2. blacklist_genotypes: blacklisting (erasing) or not individual genotypes based on low quality GL")
    message("Do you want to erase genotypes below a defined threshold (y/n) ?" )
    erase.genotype <- readLines(n = 1)
    if (erase.genotype == "y") {
      message("Enter the genotype likelihood threshold values, below this threshold (GL < threshold), genotypes are erased:")
      gl.ind.threshold <- as.integer(readLines(n = 1))
    }
  }
  
  if (is.null(gl.ind.threshold)) {
    data.filter.gl.individuals <- input
    blacklist.genotypes <- NULL
  } else {
    
    data.filter.gl.individuals <- input %>% 
      dplyr::mutate(
        BLACKLIST = dplyr::if_else(GL < gl.ind.threshold, "erase", "keep"),
        BLACKLIST = stringi::stri_replace_na(str = BLACKLIST, replacement = "keep")
      )
    
    blacklist.genotypes <- data.filter.gl.individuals %>% 
      dplyr::filter(BLACKLIST == "erase") %>% 
      dplyr::select(CHROM, LOCUS, POS, POP_ID, INDIVIDUALS)
    
    data.filter.gl.individuals <- suppressWarnings(
      data.filter.gl.individuals %>% 
        dplyr::mutate(
          GT = if_else(BLACKLIST == "erase", "000000", GT),
          GL = if_else(BLACKLIST == "erase", as.numeric(NA_character_), GL),
          READ_DEPTH = if_else(BLACKLIST == "erase", as.numeric(NA_character_), READ_DEPTH),
          ALLELE_REF_DEPTH = if_else(BLACKLIST == "erase", as.numeric(NA_character_), ALLELE_REF_DEPTH),
          ALLELE_ALT_DEPTH = if_else(BLACKLIST == "erase", as.numeric(NA_character_), ALLELE_ALT_DEPTH)
        )
    )
    
    if (tibble::has_name(data.filter.gl.individuals, "GT_VCF")) {
      data.filter.gl.individuals <- data.filter.gl.individuals %>% 
        dplyr::mutate(GT_VCF = if_else(BLACKLIST == "erase", "./.", GT_VCF))
    }
    
    if (tibble::has_name(data.filter.gl.individuals, "GT_BIN")) {
      data.filter.gl.individuals <- data.filter.gl.individuals %>% 
        dplyr::mutate(GT_BIN = if_else(BLACKLIST == "erase", NA_character_, GT_BIN))
    }
    
    # interesting stats.
    erased.genotype.number <- length(blacklist.genotypes$INDIVIDUALS)
    total.genotype.number <- length(input$GT[input$GT != "./."])
    percentage <- paste(round(((erased.genotype.number/total.genotype.number)*100), 6), "%", sep = " ")
    cat("######################## ERASING GENOTYPES ############################\n")
    message(stringi::stri_join("Total number of genotypes: ", total.genotype.number))
    message(stringi::stri_join("Blacklisted genotypes bsed on GL: ", erased.genotype.number))
    message(stringi::stri_join("Percentage erased: ", percentage))
    cat("#######################################################################\n")
    
    readr::write_tsv(x = blacklist.genotypes, path = "blacklist.genotypes.gl.tsv", col_names = TRUE)
    message("Writing the blacklist genotypes based on GL information in your working directory\nblacklist.genotypes.gl.tsv")
    
    # Update filters.parameters
    message("Updating the file storing the filters parameters: filters_parameters.tsv")
    filters.parameters <- tibble::data_frame(
      FILTERS = as.character("Genotype likelihood"), 
      PARAMETERS = as.character("gl.ind.threshold"), 
      VALUES = as.integer(gl.ind.threshold), 
      BEFORE = as.character(total.genotype.number), 
      AFTER = as.character(total.genotype.number - erased.genotype.number), 
      BLACKLIST = erased.genotype.number, 
      UNITS = "genotypes", 
      COMMENTS = as.character("NA")
    )
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  }
  
  if (interactive.filter) {
    message("Create summary and show the violin plot of updated individuals genotype likelihood,\ni.e. after erasing lower quality genotypes (y/n)): ")
    summary.blacklist.genotypes <- as.character(readLines(n = 1))
    if (summary.blacklist.genotypes == "y") {
      # updated summary
      summary.blacklist.genotypes <- summary_genotype_likelihood(data = data.filter.gl.individuals, pop.levels = pop.levels, gl.approach = gl.approach, folder = path.folder)
      
      # plot_2: updated plot after erasing genotypes
      message("Generating the updated violin plot may take some time...")
      genotype.likelihood.violin.plot.individuals.updated <- plot_violinplot_genotype_likelihood_individuals(data = data.filter.gl.individuals)
      print(genotype.likelihood.violin.plot.individuals.updated)
      # save
      ggsave(stringi::stri_join(path.folder.step2, "/genotype.likelihood.violin.plot.individuals.updated.pdf"), width = pop.number, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder.step2, "/genotype.likelihood.violin.plot.individuals.updated.png"), width = pop.number, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.individuals.updated) were saved in this directory: ", path.folder.step2))
    }
    # plot_3: before and after figure
    message("Create a combined violin plot of the individuals genotype likelihood with\nbefore and after erasing lower quality genotypes (y/n)): ")
    combined.plot <- as.character(readLines(n = 1))
    if (combined.plot == "y") {
      data.combined <- bind_rows(
        data.select <- input %>% 
          dplyr::select(POP_ID, INDIVIDUALS, GL) %>% 
          dplyr::mutate(GROUP = rep("before", n())),
        data.filter.gl.individuals.select <- data.filter.gl.individuals %>% 
          dplyr::select(POP_ID, INDIVIDUALS, GL) %>% 
          dplyr::mutate(GROUP = rep("after", n()))
      ) %>% 
        dplyr::mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
      
      # combined plot
      message("Generating combined violin plot, may take some time...")
      genotype.likelihood.violin.plot.individuals.combined <- plot_violinplot_genotype_likelihood_individuals(data = data.combined) 
      genotype.likelihood.violin.plot.individuals.combined$facet <- facet_grid(~GROUP)
      print(genotype.likelihood.violin.plot.individuals.combined)
      # save
      ggsave(stringi::stri_join(path.folder.step2, "/genotype.likelihood.violin.plot.individuals.combined.pdf"), width = pop.number*2, height = 10, dpi = 600, units = "cm", useDingbats = F)
      ggsave(stringi::stri_join(path.folder.step2, "/genotype.likelihood.violin.plot.individuals.combined.png"), width = pop.number*2, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.individuals.combined) were saved in this directory: ", path.folder.step2))
    }
  }
  
  ## Step 3: gl_markers -------------------------------------------------------
  # We inspect at the locus level, by pop and generate the figures
  if (interactive.filter) {
    message("Step 3. gl_markers: Inspecting the genotype likelihood at the marker level")
    path.folder.step3 <- stringi::stri_join(path.folder, "/03_gl_markers")
    dir.create(path.folder.step3)
    message("Show the density distribution plot (y/n)): ")
    density.distribution <- as.character(readLines(n = 1))
    
    # plot_4:Density distribution of genotype likelihood summary of loci
    if (density.distribution == "y") {
      genotype.likelihood.density.distribution.figure <- plot_density_distribution_genotype_likelihood(data = summary.blacklist.genotypes$gl.summary.marker.pop, aes.colour = aes(y = ..scaled.., color = GENOTYPE_LIKELIHOOD_GROUP), adjust.bin = 1) + facet_grid(POP_ID~GENOTYPE_LIKELIHOOD_GROUP, scales = "free")
      print(genotype.likelihood.density.distribution.figure)
      ggsave(stringi::stri_join(path.folder.step3, "/genotype.likelihood.density.distribution.figure.pdf"), width = 15, height = pop.number*1.5, dpi = 600, units = "cm", useDingbats = FALSE)
      ggsave(stringi::stri_join(path.folder.step3, "/genotype.likelihood.density.distribution.figure.png"), width = 15, height = pop.number*1.5, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (genotype.likelihood.density.distribution.figure) were saved in this directory: ", path.folder.step3))
    }
    
    # plot_5: violin plot
    message("Show the violin plot (y/n)): ")
    violinplot <- as.character(readLines(n = 1))
    if (violinplot == "y") {
      genotype.likelihood.violin.plot.figure <- plot_boxplot_genotype_likelihood(data = summary$gl.summary.marker.pop) + facet_wrap(~GENOTYPE_LIKELIHOOD_GROUP, nrow = 1, ncol = 5, scales = "free")
      print(genotype.likelihood.violin.plot.figure)
      ggsave(stringi::stri_join(path.folder.step3, "/genotype.likelihood.violin.plot.figure.pdf"), width = pop.number*2, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
      ggsave(stringi::stri_join(path.folder.step3, "/genotype.likelihood.violin.plot.figure.png"), width = pop.number*2, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.figure) were saved in this directory: ", path.folder.step3))
      message("Look for inconsistencies, patterns and trends between your populations")
    }
  }
  
  if (interactive.filter == FALSE) {
    if (gl.approach == "haplotype") {
      message("Approach selected for GL statistics: haplotype")
      data.sum <- data.filter.gl.individuals %>%
        dplyr::group_by(LOCUS, POP_ID) %>% # at the population level
        dplyr::summarise(
          GL_MEAN = mean(GL, na.rm = TRUE),
          GL_MEDIAN = stats::median(GL, na.rm = TRUE),
          GL_MIN = min(GL, na.rm = TRUE),
          GL_MAX = max(GL, na.rm = TRUE),
          GL_DIFF = GL_MAX - GL_MIN
        ) %>% 
        dplyr::group_by(LOCUS, POP_ID)
    } else {
      message("Approach selected for GL statistics: SNP")
      data.sum <- data.filter.gl.individuals %>%
        dplyr::group_by(LOCUS, POS, POP_ID) %>% # at the population level
        dplyr::summarise(
          GL_MEAN = mean(GL, na.rm = TRUE),
          GL_MEDIAN = stats::median(GL, na.rm = TRUE),
          GL_MIN = min(GL, na.rm = TRUE),
          GL_MAX = max(GL, na.rm = TRUE),
          GL_DIFF = GL_MAX - GL_MIN
        ) %>% dplyr::group_by(LOCUS, POS, POP_ID)
    }
  } else {# interactive
    if (gl.approach == "haplotype") {
      data.sum <- summary.blacklist.genotypes$gl.summary.marker.pop %>%
        dplyr::group_by(LOCUS, POP_ID) %>% 
        tidyr::spread(data = ., key = GENOTYPE_LIKELIHOOD_GROUP, value = VALUE)
    } else {
      data.sum <- summary.blacklist.genotypes$gl.summary.marker.pop %>%
        dplyr::group_by(LOCUS, POS, POP_ID) %>% 
        tidyr::spread(data = ., key = GENOTYPE_LIKELIHOOD_GROUP, value = VALUE)
    }
  }
  
  ## Step 4: blacklist_markers ---------------------------------------------------
  # At this point interactive or not, we have the same input data
  if (interactive.filter) {
    message("Step 4. blacklist_markers: Blacklisting (erasing) or not markers genotypes based on low quality GL found at the marker level")
    message("Enter the population threshold (proportion, percentage or fixed) to keep the marker\ni.e. that whatever the locus GL statistics, it will need to pass in x (proportion, percentage or fixed) pop to be kept:")
    pop.threshold <- as.numeric(readLines(n = 1))
    path.folder.step4 <- stringi::stri_join(path.folder, "/04_blacklist_markers")
    dir.create(path.folder.step4)
  }
  
  # percent ?
  if (interactive.filter) {
    message("Is the pop.threshold entered above a percentage (TRUE/FALSE) ?")
    percent <- as.character(readLines(n = 1))
  }
  
  if (stri_detect_fixed(pop.threshold, ".") & pop.threshold < 1) {
    multiplication.number <- 1/pop.number
    message(stringi::stri_join("Using a proportion threshold of , ", pop.threshold))
    threshold.id <- "proportion"
  } else if (stri_detect_fixed(percent, "T")) {
    multiplication.number <- 100/pop.number
    message(stringi::stri_join("Using a percentage threshold of ", pop.threshold))
    threshold.id <- "percent"
  } else {
    multiplication.number <- 1
    message(stringi::stri_join("Using a fixed threshold of ", pop.threshold))
    threshold.id <- "fixed"
  }
  
  # gl.mean.threshold
  if (interactive.filter) {  
    message(stringi::stri_join("Enter the mean genotype likelihood threshold number.\ni.e. markers with gl.mean.threshold >= threshold value in ", pop.threshold, " (", threshold.id, ") ", "pop, will be kept"))
    message("To turn off this filter, enter: off")
    message("gl.mean.threshold:")
    gl.mean.threshold <- readLines(n = 1)
  }
  if (gl.mean.threshold == "off") gl.mean.threshold <- NULL
  
  if (!is.null(gl.mean.threshold)) {
    gl.mean.threshold <- as.integer(gl.mean.threshold)
    data.sum <- data.sum %>%
      dplyr::filter(GL_MIN >= gl.mean.threshold)
    
    # Update filters.parameters
    filters.parameters <- tibble::data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("gl.mean.threshold"), VALUES = as.integer(gl.mean.threshold), BEFORE = as.character("NA"), AFTER = as.character("NA"), BLACKLIST = as.character("NA"), UNITS = "NA", COMMENTS = as.character("NA"))
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  }
  
  
  # gl.min.threshold
  if (interactive.filter) {
    message(stringi::stri_join("Enter the minimum genotype likelihood threshold number.\ni.e. markers with gl.min.threshold >= threshold value in ", pop.threshold, " (", threshold.id, ") ", "pop, will be kept"))
    message("To turn off this filter, enter: off")
    message("gl.min.threshold:")
    gl.min.threshold <- readLines(n = 1)
  }
  
  if (gl.min.threshold == "off") gl.min.threshold <- NULL
  
  if (!is.null(gl.min.threshold)) {
    gl.min.threshold <- as.integer(gl.min.threshold)
    data.sum <- data.sum %>%
      dplyr::filter(GL_MIN >= gl.min.threshold)
    # Update filters.parameters
    filters.parameters <- tibble::data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("gl.min.threshold"), VALUES = as.integer(gl.min.threshold), BEFORE = as.character("NA"), AFTER = as.character("NA"), BLACKLIST = as.character("NA"), UNITS = "NA", COMMENTS = as.character("NA"))
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  }
  
  # gl.diff.threshold
  if (interactive.filter) {
    message(stringi::stri_join("Enter the difference between the max and min genotype likelihood threshold number.\ni.e. markers with gl.diff.threshold <= threshold value in ", pop.threshold, " (", threshold.id, ") ", "pop, will be kept"))
    message("To turn off this filter, enter: off")
    message("gl.diff.threshold:")
    gl.diff.threshold <- readLines(n = 1)
  }
  
  if (gl.diff.threshold == "off") gl.diff.threshold <- NULL
  
  if (!is.null(gl.diff.threshold)) {
    gl.diff.threshold <- as.integer(gl.diff.threshold)
    data.sum <- data.sum %>%
      dplyr::filter(GL_DIFF <= gl.diff.threshold)
    # Update filters.parameters
    filters.parameters <- tibble::data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("gl.diff.threshold"), VALUES = as.integer(gl.diff.threshold), BEFORE = as.character("NA"), AFTER = as.character("NA"), BLACKLIST = as.character("NA"), UNITS = "NA", COMMENTS = as.character("NA"))
    readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  }
  
  # Filter the data set
  if (gl.approach == "haplotype") {
    message("Approach selected: haplotype")
    filter <- data.sum %>%
      dplyr::group_by(LOCUS) %>%
      dplyr::tally(.) %>% # Globally accross loci
      dplyr::filter((n * multiplication.number) >= as.numeric(pop.threshold)) %>%
      dplyr::select(LOCUS) %>% 
      dplyr::left_join(data.filter.gl.individuals, by = "LOCUS") %>%
      dplyr::arrange(LOCUS, POS, POP_ID)
  } else {
    message("Approach selected: SNP")
    filter <- data.sum %>%
      dplyr::group_by(LOCUS, POS) %>%
      dplyr::tally(.) %>% # Globally accross loci
      dplyr::filter((n * multiplication.number) >= as.numeric(pop.threshold)) %>%
      dplyr::select(LOCUS, POS) %>%
      dplyr::left_join(data.filter.gl.individuals, by = c("LOCUS", "POS")) %>%
      dplyr::arrange(LOCUS, POS, POP_ID)
  }
  
  # Update filters.parameters SNP
  filters.parameters <- tibble::data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("pop.threshold"), VALUES = stringi::stri_join(pop.threshold, " (", threshold.id, ")"), BEFORE = as.integer(n_distinct(data.filter.gl.individuals$POS)), AFTER = as.integer(n_distinct(filter$POS)), BLACKLIST = n_distinct(data.filter.gl.individuals$POS) - n_distinct(filter$POS), UNITS = "SNP", COMMENTS = as.character("NA"))
  readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  
  # Update filters.parameters LOCUS
  filters.parameters <- tibble::data_frame(FILTERS = as.character("Genotype likelihood"), PARAMETERS = as.character("pop.threshold"), VALUES = stringi::stri_join(pop.threshold, " (", threshold.id, ")"), BEFORE = as.integer(n_distinct(data.filter.gl.individuals$LOCUS)), AFTER = as.integer(n_distinct(filter$LOCUS)), BLACKLIST = n_distinct(data.filter.gl.individuals$LOCUS) - n_distinct(filter$LOCUS), UNITS = "LOCUS", COMMENTS = as.character("NA"))
  readr::write_tsv(x = filters.parameters, path = "filters_parameters.tsv", append = TRUE, col_names = FALSE)
  
  # saving tidy data 
  if (!is.null(filename)) {
    message("Writing the tidy vcf file filtered with GL information in your working directory...")
    readr::write_tsv(filter, filename, append = FALSE, col_names = TRUE)
  }
  
  # saving whitelist
  whitelist.markers <- dplyr::ungroup(filter) %>%
    dplyr::distinct(CHROM, LOCUS, POS)
  message("Writing the whitelist of markers based on GL information in your working directory\nwhitelist.markers.gl.tsv")
  readr::write_tsv(whitelist.markers, "whitelist.markers.gl.tsv", append = FALSE, col_names = TRUE)
  
  # saving blacklist
  blacklist.markers <- dplyr::ungroup(input) %>%
    dplyr::distinct(CHROM, LOCUS, POS) %>% 
    dplyr::anti_join(whitelist.markers, by = c("CHROM", "LOCUS", "POS"))
  message("Writing the blacklist of markers based on GL information in your working directory\nblacklist.markers.gl.tsv")
  readr::write_tsv(blacklist.markers, "blacklist.markers.gl.tsv", append = FALSE, col_names = TRUE)
  
  if (interactive.filter) {
    message("Summary of individuals genotype likelihood, i.e. after filtering the markers")
    summary.blacklist.markers <- summary_genotype_likelihood(data = filter, pop.levels = pop.levels, gl.approach = gl.approach, folder = path.folder.step4)
    
    
    # before and after figure
    data.combined <- dplyr::bind_rows(
      summary$gl.summary.marker.pop %>% 
        dplyr::mutate(GROUP = rep("before", n())),
      summary.blacklist.markers$gl.summary.marker.pop %>% 
        dplyr::mutate(GROUP = rep("after", n()))
    ) %>% 
      dplyr::mutate(GROUP = factor(GROUP, levels = c("before", "after"), ordered = TRUE))
    
    message("Show the density distribution plot before/after filters (y/n)): ")
    density.distribution <- as.character(readLines(n = 1))
    # plot_6: Density distribution of genotype likelihood summary of loci.
    if (density.distribution == "y") {
      genotype.likelihood.density.distribution.figure.before.after.filters <- plot_density_distribution_genotype_likelihood(data = data.combined, aes.colour = aes(y = ..scaled.., color = GENOTYPE_LIKELIHOOD_GROUP), adjust.bin = 1) + facet_grid(POP_ID+GROUP~GENOTYPE_LIKELIHOOD_GROUP, scales = "free")
      print(genotype.likelihood.density.distribution.figure.before.after.filters)
      ggsave(stringi::stri_join(path.folder.step4, "/genotype.likelihood.density.distribution.figure.before.after.filters.pdf"), width = 20, height = pop.number*3, dpi = 600, units = "cm", useDingbats = FALSE)
      ggsave(stringi::stri_join(path.folder.step4, "/genotype.likelihood.density.distribution.figure.before.after.filters.png"), width = 20, height = pop.number*3, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (genotype.likelihood.density.distribution.figure.before.after.filters) were saved in this directory: ", path.folder.step4))
    }
    
    # plot_7: violin plot 
    message("Show the violin plot before/after filters (y/n)): ")
    violin.plot <- as.character(readLines(n = 1))
    if (violin.plot == "y") {
      # combined plot
      message("Generating combined violin plot, may take some time...")
      genotype.likelihood.violin.plot.figure.before.after.filters <- plot_boxplot_genotype_likelihood(data = data.combined) + facet_grid(GROUP~GENOTYPE_LIKELIHOOD_GROUP, scales = "fixed")
      print(genotype.likelihood.violin.plot.figure.before.after.filters)
      ggsave(stringi::stri_join(path.folder.step4, "/genotype.likelihood.violin.plot.figure.before.after.filters.pdf"), width = pop.number*2, height = 10, dpi = 600, units = "cm", useDingbats = FALSE)
      ggsave(stringi::stri_join(path.folder.step4, "/genotype.likelihood.violin.plot.figure.before.after.filters.png"), width = pop.number*2, height = 10, dpi = 300, units = "cm")
      message(stringi::stri_join("2 versions (pdf and png) of the plot (genotype.likelihood.violin.plot.figure.before.after.filters) were saved in this directory: ", path.folder.step4))
    }
  }
  
  # results --------------------------------------------------------------------
  cat("############################### RESULTS ###############################\n")
  message(stringi::stri_join("The number of markers removed by the GL filter:\nSNP: ", n_distinct(data.filter.gl.individuals$POS) - n_distinct(filter$POS), "\nLOCUS: ", n_distinct(data.filter.gl.individuals$LOCUS) - n_distinct(filter$LOCUS)))
  message("The number of markers before -> after the GL filter")
  message(stringi::stri_join("SNP: ", as.integer(n_distinct(data.filter.gl.individuals$POS)), " -> ", as.integer(n_distinct(filter$POS))))
  message(stringi::stri_join("LOCUS: ", as.integer(n_distinct(data.filter.gl.individuals$LOCUS)), " -> ", as.integer(n_distinct(filter$LOCUS))))
  cat("############################## completed ##############################\n")
  
  if (interactive.filter) {
    res <- list()
    res$gl.summary.individuals <- summary$gl.individuals
    res$gl.summary.marker.pop <- summary$gl.summary.marker.pop
    res$gl.summary.pop <- summary$gl.summary.pop
    res$gl.summary.individuals.blacklist.genotypes <- summary.blacklist.genotypes$gl.individuals
    res$gl.summary.marker.pop.blacklist.genotypes <- summary.blacklist.genotypes$gl.summary.marker.pop
    res$gl.summary.pop.blacklist.genotypes <- summary.blacklist.genotypes$gl.summary.pop
    res$genotype.likelihood.violin.plot.individuals <- genotype.likelihood.violin.plot.individuals
    res$genotype.likelihood.violin.plot.individuals.updated <- genotype.likelihood.violin.plot.individuals.updated
    res$genotype.likelihood.violin.plot.individuals.combined <- genotype.likelihood.violin.plot.individuals.combined
    res$genotype.likelihood.density.distribution.figure <- genotype.likelihood.density.distribution.figure
    res$genotype.likelihood.violin.plot.figure <- genotype.likelihood.violin.plot.figure
    res$gl.summary.individuals.filter.markers <- summary.blacklist.markers$gl.individuals
    res$gl.summary.marker.pop.filter.markers <- summary.blacklist.markers$gl.summary.marker.pop
    res$gl.summary.pop.filter.markers <- summary.blacklist.markers$gl.summary.pop
    res$genotype.likelihood.density.distribution.figure.before.after.filters <- genotype.likelihood.density.distribution.figure.before.after.filters
    res$genotype.likelihood.violin.plot.figure.before.after.filters <- genotype.likelihood.violin.plot.figure.before.after.filters
    res$blacklist.genotypes <- blacklist.genotypes
    res$blacklist.markers <- blacklist.markers
    res$whitelist.markers <- whitelist.markers
    res$tidy.vcf.filtered.gl <- filter
  } else {
    res <- filter
  }
  
  return(res)
}




