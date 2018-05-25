#' @name vcf2dadi
#' @title Create a \code{dadi} SNP input file from a any vcf file.
#' @description This function will create a \code{dadi} SNP input file using a
#' VCF file (Danecek et al. 2011). Missing data can bias demographic inference,
#' \code{vcf2dadi} was created to address this problem, providing a customizable
#' imputation framework specifically designed to work with GBS/RAD data.
#' If your VCF is not filtered, you can supply the function a whitelist of loci and a
#' blacklist of individuals.

#' @inheritParams genomic_converter
#' @inheritParams tidy_genomic_data
#' @inheritParams radiator_imputations_module

#' @param fasta.outgroup (optional) The fasta output file from STACKS. This file is
#' required to use an outgroup. Default: \code{fasta.outgroup = NULL}.

#' @param fasta.ingroup (optional) The fasta output file from STACKS. This file is
#' required to use with an outgroup. Leave empty if no outgroup.
#' Default: \code{fasta.ingroup = NULL}.

#' @param sumstats.outgroup (optional) The sumstats output file from STACKS
#' when running STACKS for the outgroup fasta file. This file is
#' required to use an outgroup. Default: \code{sumstats.outgroup = NULL}.

#' @param sumstats.ingroup (optional) The sumstats output file from STACKS
#' when running STACKS for the ingroup fasta file.This file is
#' required to use with an outgroup. Leave empty if no outgroup.
#' Default: \code{sumstats.ingroup = NULL}.

#' @param dadi.input.filename (optional) Name of the \code{dadi} SNP input file
#' written to the working directory. e.g. \code{dadi.file.tsv}.
#' Default use date and time to make the file. If used, the file extension
#' need to finish with \code{.tsv or .txt}.


#' @details
#' \strong{Input data:}
#'
#' To discriminate the long from the wide format,
#' the function \pkg{radiator} \code{\link[radiator]{tidy_wide}} searches
#' for \code{MARKERS or LOCUS} in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns:
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals),
#' the remaining columns are the markers in separate columns storing genotypes.
#'
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info.
#' (e.g. from a VCF see \pkg{radiator} \code{\link{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID},
#' \code{MARKERS or LOCUS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#'
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 or 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 or 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#'
#'
#' \strong{Imputations details:}
#' The imputations using Random Forest requires more time to compute and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals
#' will require 15 min.

#' @export
#' @rdname vcf2dadi
#' @importFrom parallel detectCores
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join contains
#' @importFrom stringi stri_length stri_pad stri_sub stri_replace_all_fixed stri_replace_na
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom tibble as_data_frame has_name
#' @importFrom readr write_tsv


#' @return The function returns a list with 1 or 2 objects (w/o imputations):
#' `$dadi.no.imputation`
#' `$dadi.imputed`
#' The data frame are accessed form the list with `$`.

#' @examples
#' #See vignette `vignette_vcf2dadi` for more info.
#' \dontrun{
#' dadi.no.imputation <- vcf2dadi(
#' data = "batch_1.vcf",
#' whitelist.markers = "whitelist.loci.txt",
#' strata = "strata.file.tsv",
#' pop.levels = c("PAN", "COS"),
#' common.markers = TRUE,
#' fasta.ingroup = "batch_1.ingroup.fa",
#' fasta.outgroup = "batch_1.outgroup.fa",
#' sumstats.ingroup = "batch_1.sumstats.ingroup.tsv",
#' sumstats.outgroup = "batch_1.sumstats.outgroup.tsv"
#' )
#'
#' With Imputations:
#' dadi.files <- vcf2dadi(
#' data = "batch_1.vcf",
#' whitelist.markers = "whitelist.loci.txt",
#' strata = "strata.file.tsv",
#' pop.levels = c("PAN", "COS"),
#' common.markers = TRUE,
#' fasta.ingroup = "batch_1.ingroup.fa",
#' fasta.outgroup = "batch_1.outgroup.fa",
#' sumstats.ingroup = "batch_1.sumstats.ingroup.tsv",
#' sumstats.outgroup = "batch_1.sumstats.outgroup.tsv",
#' imputation.method = "max",
#' hierarchical.levels = "populations",
#' num.tree = 100,
#' verbose = FALSE,
#' parallel.core = 8
#' )
#' # to get the imputed data frame:
#' dadi.imputed.df <- dadi.files$dadi.imputed
#' }

#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.
#' @references Gutenkunst RN, Hernandez RD, Williamson SH, Bustamante CD (2009)
#' Inferring the Joint Demographic History of Multiple Populations
#' from Multidimensional SNP Frequency Data (G McVean, Ed,).
#' PLoS genetics, 5, e1000695.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841-860.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and
#' Anne-Laure Ferchaud \email{annelaureferchaud@@gmail.com}

vcf2dadi <- function(
  data,
  fasta.ingroup = NULL,
  fasta.outgroup = NULL,
  sumstats.ingroup = NULL,
  sumstats.outgroup = NULL,
  dadi.input.filename = NULL,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  max.marker = NULL,
  strata = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  imputation.method = NULL,
  hierarchical.levels = "populations",
  num.tree = 50,
  pred.mean.matching = 0,
  random.seed = NULL,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1
){

  # dadi unicode character: \u2202
  cat("#######################################################################\n")
  cat("########################## radiator::vcf2dadi ###########################\n")
  cat("#######################################################################\n")

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) stop("Input file missing")

  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...
  # pop.levels, pop.labels, pop.select, strata, etc
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

  # File type detection----------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)


  if (is.null(imputation.method)) {
    message("vcf2dadi: no imputation...")
  } else {
    message("vcf2dadi: with imputations...")
  }

  # Strata argument required for VCF and haplotypes files **********************
  if (data.type == "vcf.file" & is.null(strata)) stop("strata argument is required")
  if (data.type == "haplo.file") stop("This function is for biallelic dataset only")

  # Import input ***************************************************************
  input <- radiator::tidy_genomic_data(
    data = data,
    vcf.metadata = FALSE,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    max.marker = max.marker,
    snp.ld = snp.ld,
    common.markers = common.markers,
    maf.thresholds = maf.thresholds,
    strata = strata,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    pop.select = pop.select,
    filename = NULL
  )

  # create a strata.df
  strata.df <- dplyr::distinct(.data = input, INDIVIDUALS, POP_ID, .keep_all = TRUE)
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels

  # Compute count and Minor Allele Frequency -----------------------------------
  # We split the alleles here to prep for MAF
  # need to compute REF/ALT allele for non VCF file
  if (!tibble::has_name(input, "GT_VCF")) {
    ref.change <- radiator::change_alleles(data = input)$input
    input <- dplyr::left_join(input, ref.change, by = c("MARKERS", "INDIVIDUALS"))
  }

  # MAF = the ALT allele in the VCF
  message("Computing the Allele Frequency Spectrum")

  input.count <- input %>%
    dplyr::group_by(MARKERS, POP_ID, REF, ALT) %>%
    dplyr::summarise(
      N = as.numeric(n()),
      PP = as.numeric(length(GT_VCF[GT_VCF == "0/0"])),
      PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
      QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
    ) %>%
    dplyr::mutate(MAF = ((QQ*2) + PQ)/(2*N)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
      REF = stringi::stri_replace_all_fixed(str = REF, pattern = c("001", "002", "003", "004"), replacement = c("A", "C", "G", "T"), vectorize_all = FALSE),
      ALT = stringi::stri_replace_all_fixed(str = ALT, pattern = c("001", "002", "003", "004"), replacement = c("A", "C", "G", "T"), vectorize_all = FALSE)
    )

  # Function to make dadi input  data format ***********************************
  message("Preparing \u2202a\u2202i input SNP data format")
  write_dadi <- function(input, write.imputation, ...){
    # input <- input.count.imp # testing
    # input <- input.count # testing
    # input.bk <- input # testing

    if (is.null(fasta.ingroup)) {
      dadi.input <- suppressWarnings(
        input %>%
          dplyr::group_by(MARKERS, POP_ID, REF, ALT) %>%
          dplyr::mutate(
            A1 = (2 * PP) + PQ,
            A2 = (2 * QQ) + PQ
          ) %>%
          dplyr::ungroup(.) %>%
          dplyr::select(POP_ID, Allele1 = REF, A1, Allele2 = ALT, A2, MARKERS) %>%
          tidyr::gather(ALLELE_GROUP, COUNT, -c(POP_ID, MARKERS, Allele1, Allele2)) %>%
          tidyr::unite(POP, POP_ID, ALLELE_GROUP, sep = "_") %>%
          dplyr::group_by(MARKERS, Allele1, Allele2) %>%
          tidyr::spread(data = ., key = POP, value = COUNT) %>%
          dplyr::mutate(
            IN_GROUP = rep("---", n()), #version 2
            OUT_GROUP = rep("---", n())
          ) %>%
          dplyr::select(IN_GROUP, OUT_GROUP, Allele1, dplyr::contains("A1"), Allele2, dplyr::contains("A2"), MARKERS) %>%
          dplyr::arrange(MARKERS)
      )
    } # no fasta, no outgroup
    if (!is.null(fasta.ingroup)) {# With outgroup and ingroup fasta file
      message("using outgroup info")
      # Get the list of ref. allele in the vcf of the ingroup
      ref.allele.vcf.ingroup <- input %>%
        dplyr::ungroup(.) %>%
        dplyr::distinct(MARKERS, REF, .keep_all = TRUE) %>%
        dplyr::arrange(MARKERS) %>%
        tidyr::separate(MARKERS, c("CHROM", "LOCUS", "POS"), sep = "__") %>%
        dplyr::distinct(CHROM, LOCUS, POS, REF, .keep_all = TRUE) %>%
        dplyr::mutate(
          CHROM = as.character(stringi::stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1", vectorize_all = FALSE)),
          LOCUS = as.integer(LOCUS),
          POS = as.integer(POS)
        ) %>%
        dplyr::arrange(CHROM, LOCUS, POS)
      # View(ref.allele.vcf.ingroup)

      # keep the list of markers
      markers <- dplyr::ungroup(ref.allele.vcf.ingroup) %>% dplyr::select(-REF)

      # import out- and in- group fasta files ********************************
      fasta.outgroup <- data.table::fread(
        input = fasta.outgroup,
        sep = "\t",
        stringsAsFactors = FALSE,
        showProgress = FALSE,
        verbose = FALSE,
        header = FALSE,
        col.names = "LOCUS"
      ) %>%
        dplyr::mutate(
          ID.FILTER = rep(c("LOCUS", "SEQ"), times = n()/2),
          ANCESTRAL = rep("outgroup", times = n())
        )

      fasta.ingroup <- data.table::fread(
        input = fasta.ingroup,
        sep = "\t",
        stringsAsFactors = FALSE,
        showProgress = FALSE,
        verbose = FALSE,
        header = FALSE,
        col.names = "LOCUS"
      ) %>%
        tibble::as_data_frame() %>%
        dplyr::mutate(
          ID.FILTER = rep(c("LOCUS", "SEQ"), times = n()/2),
          ANCESTRAL = rep("ingroup", times = n())
        )

      fasta.data <- dplyr::bind_rows(fasta.outgroup, fasta.ingroup)

      fasta.ingroup <- NULL
      fasta.outgroup <- NULL

      message("Preparing fasta sequences")
      fasta.prep <- fasta.data %>%
        dplyr::filter(ID.FILTER == "LOCUS") %>%
        dplyr::select(-ID.FILTER) %>%
        dplyr::bind_cols(fasta.data %>%
                           dplyr::filter(ID.FILTER == "SEQ") %>%
                           dplyr::select(-ID.FILTER, -ANCESTRAL, SEQUENCES = LOCUS)
        ) %>%
        tidyr::separate(LOCUS, c("LOCUS", "ALLELE"), sep = "_Sample_", extra = "warn" ) %>%
        tidyr::separate(ALLELE, c("GARBAGE", "ALLELE"), sep = "_Allele_", extra = "warn" ) %>%
        tidyr::separate(ALLELE, c("ALLELE", "INDIVIDUALS"), sep = " ", extra = "warn" ) %>%
        dplyr::mutate(LOCUS = as.integer(stringi::stri_replace_all_fixed(LOCUS, pattern = ">CLocus_", replacement = "", vectorize_all = FALSE))) %>%
        dplyr::select(-GARBAGE, -INDIVIDUALS) %>%
        dplyr::distinct(LOCUS, ALLELE, ANCESTRAL, SEQUENCES, .keep_all = TRUE)

      fasta.data <- NULL

      # import out- and in- group sumstats files*****************************

      # Note: only one sumstats might be necessary to have the info needed, need checking
      sumstats.outgroup <- data.table::fread(
        input = sumstats.outgroup,
        sep = "\t",
        stringsAsFactors = FALSE,
        skip = "Batch ID",
        select = c("Chr", "Locus ID", "BP", "Col"),
        showProgress = TRUE,
        verbose = FALSE
      ) %>%
        tibble::as_data_frame() %>%
        dplyr::select(CHROM = Chr, LOCUS = `Locus ID`, POS = BP, SNP_READ_POS = Col) %>%
        dplyr::mutate(
          CHROM = as.character(stringi::stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1", vectorize_all = FALSE)),
          LOCUS = as.integer(LOCUS),
          POS = as.integer(POS)
        ) %>%
        dplyr::distinct(CHROM, LOCUS, SNP_READ_POS, .keep_all = TRUE) %>%
        dplyr::semi_join(markers, by = c("CHROM", "LOCUS", "POS")) %>%
        dplyr::arrange(CHROM, LOCUS, POS, SNP_READ_POS) %>%
        dplyr::mutate(ANCESTRAL = rep("outgroup", times = n())) %>%
        # the SNP position on the read is +1 for the fasta file
        dplyr::mutate(SNP_READ_POS = SNP_READ_POS + 1)

      sumstats.ingroup <- data.table::fread(
        input = sumstats.ingroup,
        sep = "\t",
        stringsAsFactors = FALSE,
        skip = "Batch ID",
        select = c("Chr", "Locus ID", "BP", "Col"),
        showProgress = TRUE,
        verbose = FALSE
      ) %>%
        tibble::as_data_frame() %>%
        dplyr::select(CHROM = Chr, LOCUS = `Locus ID`, POS = BP, SNP_READ_POS = Col) %>%
        dplyr::mutate(
          CHROM = as.character(stringi::stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1", vectorize_all = FALSE)),
          LOCUS = as.integer(LOCUS),
          POS = as.integer(POS)
        ) %>%
        dplyr::distinct(CHROM, LOCUS, SNP_READ_POS, .keep_all = TRUE) %>%
        dplyr::semi_join(markers, by = c("CHROM", "LOCUS", "POS")) %>%
        dplyr::arrange(CHROM, LOCUS, POS, SNP_READ_POS) %>%
        dplyr::mutate(ANCESTRAL = rep("ingroup", times = n())) %>%
        # the SNP position on the read is +1 for the fasta file
        dplyr::mutate(SNP_READ_POS = SNP_READ_POS + 1)

      markers <- NULL # remove unused object

      # When REF and ALT allele are equal in number, this can be problematic between ANCESTRAL group of alleles.
      # For those loci, we will set them based on both group (out- and in-group), so that there is no differences

      # THINKING: when ingroup REF and ALT equal, set to REF in VCF, and EQUAL in outgroup
      # when outgroup REF and ALT equal, set the REF based on the ingroup VCF

      max.length.read <- stringi::stri_length(str = fasta.prep[1,4])

      message("combining information from the fasta and sumstats file")
      ref.allele.vcf.ingroup.fasta <- ref.allele.vcf.ingroup %>%
        # The sumstats is used to get the position of the markers along the sequence read
        dplyr::left_join(
          sumstats.ingroup %>%
            dplyr::select(-ANCESTRAL)
          , by = c("CHROM", "LOCUS", "POS")
        ) %>%
        dplyr::mutate(SNP_READ_POS = stringi::stri_replace_na(SNP_READ_POS, replacement = "NA")) %>%
        dplyr::filter(SNP_READ_POS != "NA") %>%
        dplyr::mutate(SNP_READ_POS = as.integer(SNP_READ_POS)) %>%
        dplyr::arrange(CHROM, LOCUS, POS) %>%
        dplyr::ungroup(.) %>%
        # get the fasta sequence for those LOCUS...
        dplyr::left_join(
          fasta.prep %>%
            dplyr::filter(ANCESTRAL == "ingroup") %>%
            dplyr::select(-ALLELE, -ANCESTRAL)
          , by = "LOCUS"
        ) %>%
        dplyr::mutate(SNP_READ_POS = stringi::stri_replace_na(SNP_READ_POS, replacement = "NA")) %>%
        dplyr::filter(SNP_READ_POS != "NA") %>%
        # get the flanking bases, left and right, of the SNP
        dplyr::mutate(
          SNP_READ_POS = as.integer(SNP_READ_POS),
          FASTA_REF = stringi::stri_sub(SEQUENCES, from = SNP_READ_POS, to = SNP_READ_POS),
          IN_GROUP = stringi::stri_sub(SEQUENCES, from = SNP_READ_POS - 1, to = SNP_READ_POS + 1)
        ) %>%
        # remove lines with no match between all the alleles in the fasta file and the REF in the VCF
        dplyr::filter(FASTA_REF == REF) %>%
        dplyr::group_by(CHROM, LOCUS, POS, REF) %>%
        dplyr::distinct(CHROM, LOCUS, POS, REF, .keep_all = TRUE) %>%
        dplyr::mutate(
          IN_GROUP = ifelse(SNP_READ_POS == max.length.read, stringi::stri_pad(IN_GROUP, width = 3, side = "right", pad = "-"), IN_GROUP),
          IN_GROUP = ifelse(SNP_READ_POS == 1, stringi::stri_pad(IN_GROUP, width = 3, side = "left", pad = "-"), IN_GROUP)
        )


      ingroup <- dplyr::ungroup(ref.allele.vcf.ingroup.fasta) %>%
        dplyr::select(CHROM, LOCUS, POS, IN_GROUP) %>%
        # dplyr::mutate(
        #   POS = stringi::stri_pad_left(str = POS, width = 8, pad = "0"),
        #   LOCUS = stringi::stri_pad_left(str = LOCUS, width = 8, pad = "0")
        # ) %>%
        dplyr::arrange(CHROM, LOCUS, POS) %>%
        tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__")

      markers.id <- dplyr::ungroup(ref.allele.vcf.ingroup.fasta) %>%
        dplyr::select(CHROM, LOCUS, POS, SNP_READ_POS)

      markers.id.ingroup.nuc <- dplyr::ungroup(ref.allele.vcf.ingroup.fasta) %>%
        dplyr::select(CHROM, LOCUS, POS, IN_GROUP)

      # remove unused objects
      ref.allele.vcf.ingroup <- ref.allele.vcf.ingroup.fasta <- NULL
      sumstats.ingroup <- sumstats.outgroup <- NULL

      # Same thing but with outgroup
      ref.allele.outgroup.fasta <- markers.id %>%
        dplyr::left_join(
          fasta.prep %>%
            dplyr::filter(ANCESTRAL == "outgroup") %>%
            dplyr::select(-ALLELE, -ANCESTRAL)
          , by = "LOCUS"
        ) %>%
        dplyr::mutate(SEQUENCES = stringi::stri_replace_na(SEQUENCES, replacement = "NA")) %>%
        dplyr::filter(SEQUENCES != "NA") %>%
        # get the flanking bases, left and right, of the SNP
        dplyr::mutate(OUT_GROUP = stringi::stri_sub(SEQUENCES, from = SNP_READ_POS - 1, to = SNP_READ_POS + 1)) %>%
        dplyr::select(CHROM, LOCUS, POS, OUT_GROUP, SNP_READ_POS) %>%
        dplyr::group_by(CHROM, LOCUS, POS, OUT_GROUP, SNP_READ_POS) %>%
        dplyr::tally(.) %>%
        dplyr::group_by(CHROM, LOCUS, POS) %>%
        dplyr::filter(n == max(n))

      fasta.prep <- NULL # remove unused object

      # The problem is that some markers in the outgroup will have number of REF and ALT alleles equals...
      # Making the ancestral allele call ambiguous (50/50 chance of differences...)
      ambiguous.ancestral.allele <- ref.allele.outgroup.fasta %>%
        dplyr::ungroup(.) %>%
        dplyr::select(CHROM, LOCUS, POS, SNP_READ_POS) %>%
        dplyr::group_by(CHROM, LOCUS, POS, SNP_READ_POS) %>%
        dplyr::tally(.) %>%
        dplyr::filter(n > 1) %>%
        dplyr::select(CHROM, LOCUS, POS, SNP_READ_POS) %>%
        dplyr::left_join(markers.id.ingroup.nuc, by = c("CHROM", "LOCUS", "POS")) %>%
        dplyr::rename(OUT_GROUP = IN_GROUP) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(CHROM, LOCUS, POS)

      markers.id.ingroup.nuc <- NULL # remove unused object

      outgroup <- ref.allele.outgroup.fasta %>%
        dplyr::select(-n) %>%
        dplyr::anti_join(ambiguous.ancestral.allele, by = c("CHROM", "LOCUS", "POS")) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(CHROM, LOCUS, POS) %>%
        dplyr::bind_rows(ambiguous.ancestral.allele) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(CHROM, LOCUS, POS) %>%
        dplyr::mutate(
          OUT_GROUP = ifelse(SNP_READ_POS == max.length.read, stringi::stri_pad(OUT_GROUP, width = 3, side = "right", pad = "-"), OUT_GROUP),
          OUT_GROUP = ifelse(SNP_READ_POS == 1, stringi::stri_pad(OUT_GROUP, width = 3, side = "left", pad = "-"), OUT_GROUP)
        ) %>%
        dplyr::select(CHROM, LOCUS, POS, OUT_GROUP) %>%
        # dplyr::mutate(
        #   POS = stringi::stri_pad_left(str = POS, width = 8, pad = "0"),
        #   LOCUS = stringi::stri_pad_left(str = LOCUS, width = 8, pad = "0")
        # ) %>%
        dplyr::arrange(CHROM, LOCUS, POS) %>%
        tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__")

      ambiguous.ancestral.allele <- NULL # remove unused object
      ref.allele.outgroup.fasta <- NULL # remove unused object

      # common markers between ingroup and outgroup
      markers.ingroup <- ingroup %>% dplyr::select(MARKERS)
      markers.outgroup <- outgroup %>% dplyr::select(MARKERS)

      markers.ingroup.outgroup.common <- dplyr::bind_rows(markers.ingroup, markers.outgroup) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::tally(.) %>%
        dplyr::filter(n == 2) %>%
        dplyr::arrange(MARKERS) %>%
        dplyr::distinct(MARKERS)

      markers.ingroup <- markers.outgroup <- NULL

      message(stringi::stri_join("Number of markers common between in- and out- group = ", dplyr::n_distinct(markers.ingroup.outgroup.common$MARKERS)))
      dadi.input <- suppressWarnings(
        input %>%
          dplyr::group_by(MARKERS, POP_ID, REF, ALT) %>%
          dplyr::mutate(
            A1 = (2 * PP) + PQ,
            A2 = (2 * QQ) + PQ
          ) %>%
          dplyr::ungroup(.) %>%
          dplyr::select(POP_ID, Allele1 = REF, A1, Allele2 = ALT, A2, MARKERS) %>%
          tidyr::gather(ALLELE_GROUP, COUNT, -c(POP_ID, MARKERS, Allele1, Allele2)) %>%
          tidyr::unite(POP, POP_ID, ALLELE_GROUP, sep = "_") %>%
          dplyr::group_by(MARKERS, Allele1, Allele2) %>%
          tidyr::spread(data = ., key = POP, value = COUNT) %>%
          dplyr::select(Allele1, dplyr::contains("A1"), Allele2, dplyr::contains("A2"), MARKERS) %>%
          dplyr::arrange(MARKERS) %>%
          dplyr::ungroup(.) %>%
          dplyr::semi_join(markers.ingroup.outgroup.common, by = "MARKERS") %>%
          dplyr::inner_join(ingroup, by = "MARKERS") %>%
          dplyr::inner_join(outgroup, by = "MARKERS") %>%
          dplyr::select(IN_GROUP, OUT_GROUP, Allele1, dplyr::contains("A1"), Allele2, dplyr::contains("A2"), MARKERS) %>%
          dplyr::arrange(MARKERS)
      )

      # remove unused objects
      markers.id <- NULL
      markers.ingroup.outgroup.common <- NULL
      ingroup <- NULL
      outgroup <- NULL
    }

    # We need to modify the header to remove "_A1" and "_A2" that were necessary for spreading the info accross columns
    header.dadi <- colnames(dadi.input)
    colnames(dadi.input) <- stringi::stri_replace_all_fixed(
      str = header.dadi, pattern = c("_A1", "_A2"),
      replacement = "", vectorize_all = FALSE
    )

    if (is.null(dadi.input.filename)) {
      # when filename is not provided will save the 'dadi.input' with date and time
      file.date <- format(Sys.time(), "%Y%m%d@%H%M")

      if (write.imputation == FALSE) {
        dadi.input.filename <- stringi::stri_join("dadi_input_", file.date, ".tsv", sep = "")
      }

      if (write.imputation == TRUE) {
        dadi.input.filename.imp <- stringi::stri_join("dadi_input_imputed_", file.date, ".tsv", sep = "")
      }
    } else {
      if (write.imputation == FALSE) {
        dadi.input.filename <- dadi.input.filename
      }

      if (write.imputation == TRUE) {
        dadi.input.filename.imp <- stringi::stri_replace_all_fixed(
          dadi.input.filename,
          pattern = c(".txt", ".tsv"),
          replacement = c("_imputed.txt", "_imputed.tsv"),
          vectorize_all = FALSE
        )
      }
    }
    file.version <- suppressWarnings(stringi::stri_join("#\u2202a\u2202i SNP input file generated with radiator v.", packageVersion("radiator"), sep = ""))
    file.date <- suppressWarnings(stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = ""))
    file.header.line <- suppressWarnings(as.data.frame(stringi::stri_join(file.version, file.date, sep = " ")))

    if (write.imputation == FALSE) {
      readr::write_tsv(x = file.header.line, path = dadi.input.filename, append = FALSE, col_names = FALSE) # write the header line
      readr::write_tsv(x = dadi.input, path = dadi.input.filename, append = TRUE, col_names = TRUE) # write the data frame
      message(stringi::stri_join("\u2202a\u2202i input file name is: ", dadi.input.filename, "\n", "Saved here: ", getwd()))
    }

    if (write.imputation == TRUE) {
      readr::write_tsv(x = file.header.line, path = dadi.input.filename.imp, append = FALSE, col_names = FALSE) # write the header line
      readr::write_tsv(x = dadi.input, path = dadi.input.filename.imp, append = TRUE, col_names = TRUE) # write the data frame
      message(stringi::stri_join("\u2202a\u2202i input file name is: ", dadi.input.filename.imp, "\n", "Saved here: ", getwd()))
    }

    return(dadi.input)
  } # End function write_dadi

  # Create list to store results
  res <- list()

  # without imputations (automatic)
  res$dadi.no.imputation <- write_dadi(input = input.count, write.imputation = FALSE)

  # Imputations-----------------------------------------------------------------

  if (!is.null(imputation.method)) {
    input.imp <- radiator::radiator_imputations_module(
      data = input,
      imputation.method = imputation.method,
      hierarchical.levels = hierarchical.levels,
      num.tree = num.tree,
      pred.mean.matching = pred.mean.matching,
      verbose = verbose,
      parallel.core = parallel.core,
      filename = NULL
    )

    # transform the imputed dataset  -------------------------------------------
    message("Computing the Allele Frequency Spectrum for the imputed data")

    input.count.imp <- dplyr::group_by(.data = input.imp, MARKERS, POP_ID, REF, ALT) %>%
      dplyr::summarise(
        N = as.numeric(n()),
        PP = as.numeric(length(GT_VCF[GT_VCF == "0/0"])),
        PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
        QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
      ) %>%
      dplyr::mutate(MAF = ((QQ*2) + PQ)/(2*N))

    input.imp <- NULL # remove unused object

    # output in dadi
    res$dadi.imputed <- write_dadi(input = input.count.imp, write.imputation = TRUE)
  } # end dadi with imputations
  cat("############################## completed ##############################\n")

  return(res)
} # End vcf2dadi
