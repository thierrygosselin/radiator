# Make genomic input file tidy

#' @name tidy_genomic_data
#' @title Transform common genomic dataset format in a tidy data frame
#' @description Transform genomic data set produced by massive parallel
#' sequencing pipeline (e.g.GBS/RADseq,
#' SNP chip, DArT, etc) into a tidy format. The use of blacklist and whitelist along
#' several filtering options are available to prune the dataset.
#' Several arguments are available to make your data population-wise and easily
#' rename the pop id.
#' Used internally in \href{https://github.com/thierrygosselin/radiator}{radiator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data 12 options: VCF (SNPs or Haplotypes,
#' to make the vcf population ready, see details below),
#' plink, stacks haplotype file, genind (library(adegenet)),
#' genlight (library(adegenet)), gtypes (library(strataG)), genepop, DArT,
#' and a data frame in long/tidy or wide format.
#' \emph{See details} of \code{\link[radiator]{tidy_genomic_data}}.

#' @param vcf.metadata (optional, logical or string) For VCF files only.
#' With \code{logical, TRUE/FALSE}, \code{vcf.metadata = FALSE}, only the genotype
#' information is kept (GT field). With \code{vcf.metadata = TRUE},
#' all the metadata contained in the \code{FORMAT} field will be kept in
#' the tidy data file. radiator is currently keeping/cleaning only these metadata:
#' \code{"DP", "AD", "GL", "PL", "GQ", "GOF", "NR", "NV"}.
#' If you need different one, submit a request.
#' With \code{string}, explicitely ask for a particular FORMAT field you want
#' to keep (except the GT field that is always imported).
#' e.g. \code{vcf.metadata = "PL"} or \code{vcf.metadata = c("DP", "GL")}.
#' Default: \code{vcf.metadata = FALSE}.

#' @param whitelist.markers (optional) A whitelist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by chromosome and/or locus and/or by snp.
#' The whitelist is in the working directory (e.g. "whitelist.txt").
#' de novo CHROM column with 'un' need to be changed to 1.
#' In the VCF, the column ID is the LOCUS identification.
#' Default \code{whitelist.markers = NULL} for no whitelist of markers.

#' @param monomorphic.out (optional) Should the monomorphic
#' markers present in the dataset be filtered out ?
#' Default: \code{monomorphic.out = TRUE}.

#' @param blacklist.genotype (optional) Useful to erase genotype with below
#' average quality, e.g. genotype with more than 2 alleles in diploid likely
#' sequencing errors or genotypes with poor genotype likelihood or coverage.
#' The blacklist has a minimum of 2 column headers (markers and individuals).
#' Markers can be 1 column (CHROM or LOCUS or POS),
#' a combination of 2 (e.g. CHROM and POS or CHROM and LOCUS or LOCUS and POS) or
#' all 3 (CHROM, LOCUS, POS). The markers columns must be designated: CHROM (character
#' or integer) and/or LOCUS (integer) and/or POS (integer). The id column designated
#' INDIVIDUALS (character) columns header. The blacklist must be in the working
#' directory (e.g. "blacklist.genotype.txt"). For de novo VCF, CHROM column
#' with 'un' need to be changed to 1.
#' Default: \code{blacklist.genotype = NULL} for no blacklist of
#' genotypes to erase.

#' @param snp.ld (optional) \strong{For data with locus and SNP info, like VCF and DArT file}.
#' SNP short distance linkage disequilibrium pruning. With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 3 options: \code{"random"} selection, \code{"first"} or
#' \code{"last"} SNP on the same short read/haplotype. For long distance linkage
#' disequilibrium pruning, see details below.
#' Default: \code{snp.ld = NULL}.

#' @param common.markers (optional) Logical. Default: \code{common.markers = TRUE},
#' will only keep markers in common (genotyped) between all the populations.

#' @param maf.thresholds (string, double, optional) String with
#' local/populations and global/overall maf thresholds, respectively.
#' e.g. \code{maf.thresholds = c(0.05, 0.1)} for a local maf threshold
#' of 0.05 and a global threshold of 0.1. Available for VCF, PLINK and data frame
#' files.
#' Default: \code{maf.thresholds = NULL}.

#' @param maf.approach (character, optional).
#' \code{maf.approach = "haplotype"} : looks at the minimum MAF found on the
#' read/haplotype. Using this option will discard all the markers/snp on
#' that read based on the thresholds chosen. This method is only available
#' for VCF and haplotype files, or tidy data frame from those file types.
#' \code{maf.approach = "SNP"} : treats all the SNP on the same
#' haplotype/read as independent. Doesn't work with haplotype file,
#' but does work for all other file type.
#' Default is \code{maf.approach = "SNP"}.

#' @param maf.operator (character, optional) \code{maf.operator = "AND"} or
#' default \code{maf.operator = "OR"}.
#' When filtering over LOCUS or SNP, do you want the local \code{"AND"}
#' global MAF to pass the thresholds, or ... you want the local \code{"OR"}
#' global MAF to pass the thresholds, to keep the marker?

#' @param maf.pop.num.threshold (integer, optional) When maf thresholds are used,
#' this argument is for the number of pop required to pass the maf thresholds
#' to keep the locus. Default: \code{maf.pop.num.threshold = 1}

#' @param max.marker (integer, optional) For large PLINK files,
#' useful to subsample marker number. e.g. if the data set
#' contains 200 000 markers and \code{max.marker = 10000}, 10000 markers are
#' subsampled randomly from the 200000 markers. If you need specific markers,
#' use \code{whitelist.markers} argument.
#' Default: \code{max.marker = NULL}.

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is an object in your
#' global environment or a file in the working directory
#' (e.g. "blacklist.txt").
#' Default: \code{blacklist.id = NULL}.

#' @param pop.levels (optional, string) This refers to the levels in a factor. In this
#' case, the id of the pop.
#' Use this argument to have the pop ordered your way instead of the default
#' alphabetical or numerical order. e.g. \code{pop.levels = c("QUE", "ONT", "ALB")}
#' instead of the default \code{pop.levels = c("ALB", "ONT", "QUE")}.
#' White spaces in population names are replaced by underscore.
#' Default: \code{pop.levels = NULL}.


#' @param pop.labels (optional, string) Use this argument to rename/relabel
#' your pop or combine your pop. e.g. To combine \code{"QUE"} and \code{"ONT"}
#' into a new pop called \code{"NEW"}:
#' (1) First, define the levels for your pop with \code{pop.levels} argument:
#' \code{pop.levels = c("QUE", "ONT", "ALB")}.
#' (2) then, use \code{pop.labels} argument:
#' \code{pop.levels = c("NEW", "NEW", "ALB")}.
#' To rename \code{"QUE"} to \code{"TAS"}:
#' \code{pop.labels = c("TAS", "ONT", "ALB")}.
#' Default: \code{pop.labels = NULL}. If you find this too complicated,
#' there is also the \code{strata} argument that can do the same thing,
#' see below.
#' White spaces in population names are replaced by underscore.


#' @param strata (optional/required) Required for VCF and haplotypes files,
#' optional for the other formats supported.
#'
#' The strata file is a tab delimited file with a minimum of 2 columns headers:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' If a \code{strata} file is specified with all file formats that don't
#' require it, the strata argument will have precedence on the population
#' groupings used internally in those file formats.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.
#' If you have already run
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data,
#' the strata file is similar to a stacks \emph{population map file},
#' make sure you
#' have the required column names (\code{INDIVIDUALS} and \code{STRATA}).
#' Default: \code{strata = NULL}.

#' @param pop.select (string, optional) Selected list of populations for
#' the analysis. e.g. \code{pop.select = c("QUE", "ONT")} to select \code{QUE}
#' and \code{ONT} population samples (out of 20 pops).
#' Default: \code{pop.select = NULL}

#' @param filename (optional) The function uses \code{\link[fst]{write.fst}},
#' to write the tidy data frame in
#' the working directory. The file extension appended to
#' the \code{filename} provided is \code{.rad}.
#' With default: \code{filename = NULL}, the tidy data frame is
#' in the global environment only (i.e. not written in the working directory...).

#' @param parallel.core (optional) The number of core used for parallel
#' execution during vcf import.
#' Default: \code{parallel::detectCores() - 1}.


#' @param verbose (optional, logical) When \code{verbose = TRUE}
#' the function is a little more chatty during execution.
#' Default: \code{verbose = TRUE}.


#' @details
#' \strong{Long distance SNP linkage disequilibrium pruning}
#' If you have markers position on a genome or a linkage map,
#' you can go further in removing linked markers by using \pkg{SNPRelate} or
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK}, \emph{linkage
#' disequilibrium based SNP pruning} option.
#'
#' \strong{Input files:}
#' \enumerate{
#' \item VCF biallelic file (e.g. \code{data = "batch_1.vcf"})
#' To make the VCF population ready, you need the \code{strata} argument.
#'
#' \item VCF haplotypic file (e.g. \code{data = "batch_1.haplotypes.vcf"})
#' To make the VCF population ready, you need the \code{strata} argument.
#'
#' \item haplotype file created in STACKS (e.g. \code{data = "batch_1.haplotypes.tsv"}).
#' To make the haplotype file population ready, you need the \code{strata} argument.
#'
#' \item Data frame
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
#' (e.g. from a VCF see \pkg{radiator} \code{\link[radiator]{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID},
#' \code{MARKERS or LOCUS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#'
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 of 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 of 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#'
#' \emph{How to get a tidy data frame ?}
#' \pkg{radiator} \code{\link[radiator]{tidy_genomic_data}} can transform 11 genomic
#' data formats in a tidy data frame.
#'
#' \item PLINK file in
#' \code{tped/tfam} format (e.g. \code{data =  "data.assignment.tped"}).
#' The \code{tfam} file is found based on the name prefix of the supplied
#' \code{tped} file. The first 2 columns of the \code{tfam} file will be used
#' for the \code{strata} (population identification), unless the \code{strata}
#' argument is provided.
#' Columns 1, 3 and 4 of the \code{tped} are discarded. The remaining columns
#' correspond to the genotype in the format \code{01/04}
#' where \code{A = 01, C = 02, G = 03 and T = 04}. For \code{A/T} format, use
#' PLINK or bash to convert.
#' Use \href{http://vcftools.sourceforge.net/}{VCFTOOLS} with \code{--plink-tped}
#' to convert very large VCF file. For \code{.ped} file conversion to
#' \code{.tped} use \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK}
#' with \code{--recode transpose}
#'
#' \item \code{\link[adegenet]{genind}} object from \code{\link[adegenet]{adegenet}}.
#'
#' \item \code{\link[adegenet]{genlight}} object from \code{\link[adegenet]{adegenet}}.
#'
#' \item \code{\link[strataG]{gtypes}} object from \code{\link[strataG]{strataG}}.
#'
#' \item \href{http://www.diversityarrays.com}{DArT} file.
#'
#' \item genepop data file (e.g. \code{data = "kiwi_data.gen"}). Here, the function can only use
#' alleles encoded with 3 digits.
#' }
#'
#'
#' \strong{GATK VCF files:} Some VCF have an \code{ID} column filled with \code{.},
#' the LOCUS information is all contained along the linkage group in the
#' \code{CHROM} column. To make it work with
#' \href{https://github.com/thierrygosselin/radiator}{radiator},
#' the \code{ID} column is filled with the \code{POS} column info.
#'
#' \strong{platypus VCF files:} Some VCF files don't have an ID filed with values,
#' here the same thing as GATK VCF files above is done.


#' @return The output in your global environment is a tidy data frame.
#' If \code{filename} is provided, the tidy data frame is also
#' written to the working directory with file extension \code{.rad}.

#' @export
#' @rdname tidy_genomic_data
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
#' @importFrom adegenet genind2df
# @importFrom strataG as.data.frame
#' @importFrom stringi stri_join stri_replace_all_fixed stri_extract_all_fixed stri_replace_all_regex stri_sub stri_pad_left stri_count_fixed stri_replace_na
#' @importFrom stats var median quantile
#' @importFrom purrr map flatten keep discard
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom tidyr spread gather unite separate
#' @importFrom utils count.fields
#' @importFrom readr write_tsv read_tsv
#' @importFrom vcfR read.vcfR extract.gt vcf_field_names
#' @importFrom tibble as_data_frame has_name
#' @importFrom parallel detectCores
#' @importFrom pegas VCFloci read.vcf
#' @importFrom fst write.fst

#' @examples
#' \dontrun{
#' tidy.vcf <- tidy_genomic_data(
#' data = "batch_1.vcf",
#' whitelist.markers = "whitelist.vcf.txt",
#' snp.ld = NULL,
#' common.markers = TRUE,
#' blacklist.id = "blacklist.id.treefrog.tsv",
#' strata = "strata.treefrog.tsv",
#' pop.levels = c("PAN", "COS")
#' )
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
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
#' Bender D, et al.
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795
#' @references Jombart T, Devillard S, Balloux F.
#' Discriminant analysis of principal components: a new method for the analysis
#' of genetically structured populations.
#' BMC Genet. 2010:11: 94. doi:10.1186/1471-2156-11-94
#' @references Jombart T, Ahmed I. adegenet 1.3-1: new tools for the analysis
#' of genome-wide SNP data.
#' Bioinformatics. 2011:27: 3070–3071. doi:10.1093/bioinformatics/btr521
#' @references Raymond M. & Rousset F, (1995).
#' GENEPOP (version 1.2): population genetics software for exact tests
#' and ecumenicism.
#' J. Heredity, 86:248-249
#' @references Archer FI, Adams PE, Schneiders BB.
#' strataG: An r package for manipulating, summarizing and analysing population
#' genetic data.
#' Molecular Ecology Resources. 2016. doi:10.1111/1755-0998.12559

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_genomic_data <- function(
  data,
  vcf.metadata = FALSE,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  blacklist.genotype = NULL,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  max.marker = NULL,
  blacklist.id = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  strata = NULL,
  pop.select = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {

  if (verbose) {
    cat("#######################################################################\n")
    cat("##################### radiator::tidy_genomic_data #####################\n")
    cat("#######################################################################\n")
    timing <- proc.time()
  }
  opt.change <- getOption("width")
  options(width = 70)

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

  if (data.type == "haplo.file") {
    if (verbose && !is.null(maf.thresholds)) {
      message("With stacks haplotype file the maf.approach is automatically set to: haplotype")
    }
    maf.approach <- "SNP"
    # confusing, but because the haplotpe file doesn't have snp info, only locus info
    # it's treated as markers/snp info and filtered the same way as the approach by SNP.
    # but it's really by haplotype
  }

  if (maf.approach == "haplotype") {
    if (data.type != "vcf.file" | data.type != "haplo.file") {
      stop("The haplotype approach during MAF filtering is for VCF and
           stacks haplotypes file, only. Use the snp approach for the other file types")
    }
  }


  # Strata argument required for VCF and haplotypes files-----------------------
  if (data.type == "haplo.file" | data.type == "vcf.file") {
    if (is.null(strata)) stop("strata argument is required")
  }

  # Import whitelist of markers-------------------------------------------------
  if (!is.null(whitelist.markers)) {# with Whitelist of markers
    if (is.vector(whitelist.markers)) {
      suppressMessages(whitelist.markers <- readr::read_tsv(whitelist.markers, col_names = TRUE))
    }
    columns.names.whitelist <- colnames(whitelist.markers)
    if ("CHROM" %in% columns.names.whitelist) {
      whitelist.markers$CHROM <- as.character(whitelist.markers$CHROM)
    }
    if ("LOCUS" %in% columns.names.whitelist) {
      whitelist.markers$LOCUS <- as.character(whitelist.markers$LOCUS)
    }
    if ("POS" %in% columns.names.whitelist) {
      whitelist.markers$POS <- as.character(whitelist.markers$POS)
    }
    # haplo.file
    if (data.type == "haplo.file") {
      whitelist.markers <- dplyr::select(.data = whitelist.markers, LOCUS)
      columns.names.whitelist <- colnames(whitelist.markers)
    }
    nrow.before <- nrow(whitelist.markers)
    whitelist.markers <- dplyr::distinct(whitelist.markers)
    nrow.after <- nrow(whitelist.markers)
    duplicate.whitelist.markers <- nrow.before - nrow.after
    if (duplicate.whitelist.markers > 0) {
      message("Whitelist of markers with ", duplicate.whitelist.markers, " duplicated identifiers...")
      message("    Creating unique whitelist")
      message("    Warning: downstream results might be impacted by this, check how you made your VCF file...")
    }
    nrow.before <- nrow.after <- duplicate.whitelist.markers <- NULL
  }
  # Import blacklist id --------------------------------------------------------
  if (!is.null(blacklist.id)) {# With blacklist of ID
    if (is.vector(blacklist.id)) {
      suppressMessages(blacklist.id <- readr::read_tsv(blacklist.id, col_names = TRUE))
    } else {
      if (!tibble::has_name(blacklist.id, "INDIVIDUALS")) {
        stop("Blacklist of individuals should have 1 column named: INDIVIDUALS")
      }
    }
    # not for plink file, where it's done after.
    # see plink section to understand
    blacklist.id$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = blacklist.id$INDIVIDUALS,
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )

    # remove potential duplicate id
    blacklist.id <- dplyr::distinct(.data = blacklist.id, INDIVIDUALS)
    message("Number of individuals in blacklist: ", nrow(blacklist.id))
  }

  # population levels and strata------------------------------------------------
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      # message("strata file: yes")
      number.columns.strata <- max(utils::count.fields(strata, sep = "\t"), na.rm = TRUE)
      col.types <- stringi::stri_join(rep("c", number.columns.strata), collapse = "")
      suppressMessages(strata.df <- readr::read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>%
                         dplyr::rename(POP_ID = STRATA))
    } else {
      # message("strata object: yes")
      colnames(strata) <- stringi::stri_replace_all_fixed(
        str = colnames(strata),
        pattern = "STRATA",
        replacement = "POP_ID",
        vectorize_all = FALSE
      )
      strata.df <- strata
    }

    # filtering the strata if blacklist id available
    if (!is.null(blacklist.id)) {
      strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
    }
    # Remove potential whitespace in pop_id
    strata.df$POP_ID <- stringi::stri_replace_all_fixed(strata.df$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)
    colnames.strata <- colnames(strata.df)
  }

  # Import VCF------------------------------------------------------------------
  if (data.type == "vcf.file") { # VCF
    if (verbose) message("Importing and tidying the VCF...")

    if (is.logical(vcf.metadata) && !vcf.metadata) {
      import.pegas <- TRUE
    } else {
      import.pegas <- FALSE
    }

    # detect stacks
    stacks.vcf <- readr::read_lines(file = data, skip = 2, n_max = 1) %>%
      stringi::stri_detect_fixed(str = ., pattern = "Stacks")

    # import vcf with pegas (fastest, but only GT no metadata)
    if (import.pegas) {
      # change names of columns and CHROM column modif
      input <- pegas::VCFloci(file = data, quiet = verbose) %>%
        dplyr::select(-QUAL, -INFO, -FORMAT) %>%
        dplyr::rename(LOCUS = ID) %>%
        dplyr::mutate(
          CHROM = stringi::stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1")
        ) %>%
        dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = as.character) %>%
        tibble::rownames_to_column(df = ., var = "KEEP") %>%
        dplyr::mutate(KEEP = as.integer(KEEP))
    } else {# import with vcfR
      suppressMessages(vcf.data <- vcfR::read.vcfR(file = data, verbose = FALSE))
      input <- tibble::as_data_frame(vcf.data@fix) %>%
        dplyr::select(-QUAL, -INFO) %>%
        dplyr::rename(LOCUS = ID) %>%
        dplyr::mutate(
          CHROM = stringi::stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1")
        ) %>%
        tibble::rownames_to_column(df = ., var = "KEEP") %>%
        dplyr::mutate(KEEP = as.integer(KEEP))
    }

    # bi- or multi-alllelic VCF
    alt.num <- max(unique(
      stringi::stri_count_fixed(str = unique(input$ALT), pattern = ","))) + 1

    if (alt.num > 1) {
      biallelic <- FALSE
      message("VCF is multi-allelic")
    } else {
      biallelic <- TRUE
      message("VCF is biallelic")
    }


    # Check for duplicate indentifiers
    duplicate.markers <- nrow(input) - nrow(dplyr::distinct(input, CHROM, LOCUS, POS))

    if (duplicate.markers > 0) {
      duplicate.markers <- input %>%
        dplyr::group_by(CHROM, LOCUS, POS) %>%
        dplyr::tally(.) %>%
        dplyr::filter(n > 1) %>%
        dplyr::arrange(CHROM, LOCUS, POS)
      readr::write_tsv(x = duplicate.markers, path = "duplicated.markers.tsv")
      stop("VCF file contains duplicated CHROM, ID and POS\n
             A file named: duplicated.markers.tsv was written in the directory")
    }

    # Scan and filter with FILTER column
    filter.check <- dplyr::distinct(input, FILTER)

    if (nrow(filter.check) > 1) {
      message("Filtering markers based on VCF FILTER column")
      nrow.before <- nrow(input)
      input <- dplyr::filter(input, FILTER %in% "PASS")
      nrow.after <- nrow(input)
      message("    Number of markers before = ", nrow.before)
      message("    Number of markers removed = ", nrow.before - nrow.after)
      message("    Number of markers after = ", nrow.after)
    }
    input <- dplyr::select(input, -FILTER)
    filter.check <- NULL

    # GATK VCF file sometimes have "." in the LOCUS column: replace by CHROM
    # platypus VCF file sometimes have NA in LOCUS column: replace by POS
    weird.locus <- unique(input$LOCUS)
    if (length(weird.locus) <= 1) {
      # if (is.na(weird.locus)) {
      input$LOCUS <- input$POS
      # } else {
      #   input$LOCUS <- input$CHROM
      # }
    }
    weird.locus <- NULL #unused object

    # Unique MARKERS column
    # Since stacks v.1.44 ID as LOCUS + COL (from sumstats) the position of the SNP on the locus.
    # Choose the first 100 markers to scan
    detect.snp.col <- sample(x = unique(input$LOCUS), size = 100, replace = FALSE) %>%
      stringi::stri_detect_fixed(str = ., pattern = "_") %>%
      unique
    # detect.snp.col <- dplyr::n_distinct(input$LOCUS) < length(input$LOCUS)
    if (detect.snp.col && stacks.vcf) {
      if (nrow(input) > 30000) {
        input <- input %>%
          dplyr::mutate(
            SPLIT_VEC = dplyr::ntile(x = 1:nrow(.), n = parallel.core * 3)) %>%
          split(x = ., f = .$SPLIT_VEC) %>%
          .radiator_parallel(
            X = .,
            FUN = split_vcf_id,
            mc.cores = parallel.core
          ) %>%
          dplyr::bind_rows(.) %>%
          dplyr::select(-SPLIT_VEC)
      } else {
        input <- dplyr::rename(input, ID = LOCUS) %>%
          tidyr::separate(data = ., col = ID, into = c("LOCUS", "COL"),
                          sep = "_", extra = "drop", remove = FALSE) %>%
          dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = as.character) %>%
          tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE)
      }
    } else {
      input <- tidyr::unite(
        data = input,
        MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE)
    }

    # Filter with whitelist of markers and FILTER column
    if (!is.null(whitelist.markers)) {

      if (!biallelic) {
        if (ncol(whitelist.markers) >= 3) {
          message("Note: whitelist with CHROM LOCUS POS columns and VCF haplotype:
If the whitelist was not created from this VCF,
the filtering could result in loosing all the markers.
The POS column is different in biallelic and multiallelic file...\n")

          message("Discarding the POS column in the whitelist")
          whitelist.markers <- dplyr::select(whitelist.markers, -POS)
          columns.names.whitelist <- colnames(whitelist.markers)
        }

        if (ncol(whitelist.markers) == 1 && tibble::has_name(whitelist.markers, "MARKERS")) {
          message("Note: whitelist MARKERS column and VCF haplotype:
If the whitelist was not created from this VCF,
the filtering could result in loosing all the markers.
The POS column used in the MARKERS column is different in biallelic and multiallelic file...\n")
        }
      }

      message("Filtering: ", nrow(whitelist.markers), " markers in whitelist")
      input2 <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
test <- input %>% dplyr::filter(LOCUS == "10024")
      if (nrow(input) == 0) stop("No markers left in the dataset, check whitelist...")
    }

    # keep vector
    keep.markers <- dplyr::select(input, KEEP) %>%
      dplyr::arrange(KEEP) %>%
      purrr::flatten_int(.)

    # import genotypes
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "ID", "COL", "REF", "ALT", "INDIVIDUALS", "GT")
    # bk <- input
    if (import.pegas) {
      if (verbose) message("Working on the vcf...")
      input.gt <- pegas::read.vcf(file = data, which.loci = keep.markers, quiet = verbose) %>%
        `colnames<-`(input$MARKERS) %>%
        tibble::as_data_frame(.) %>%
        tibble::rownames_to_column(., var = "INDIVIDUALS") %>%
        data.table::as.data.table(.) %>%
        data.table::melt.data.table(
          data = .,
          id.vars = "INDIVIDUALS",
          variable.name = "MARKERS",
          variable.factor = FALSE,
          value.name = as.character("GT")
        ) %>%
        tibble::as_data_frame(.)

      # input.gt <- pegas::read.vcf(file = data, from = 1, to = nrow(input), quiet = FALSE) %>%
      #   `colnames<-`(input$MARKERS)

      input <- suppressWarnings(
        dplyr::full_join(input.gt, input, by = "MARKERS") %>%
          dplyr::select(dplyr::one_of(want)) %>%
          dplyr::mutate(INDIVIDUALS = stringi::stri_replace_all_fixed(
            str = INDIVIDUALS,
            pattern = c("_", ":"),
            replacement = c("-", "-"),
            vectorize_all = FALSE)
          ) %>%
          dplyr::mutate(GT = stringi::stri_replace_na(
            str = GT, replacement = "./."))
      )
      input.gt <- NULL
    } else {
      # filter the vcf.data
      vcf.data <- vcf.data[keep.markers,]
      keep.markers <- NULL

      input <- suppressWarnings(dplyr::select(
        .data = input,
        dplyr::one_of(want)))#MARKERS, CHROM, LOCUS, POS, REF, ALT)

      filter.check <- NULL
    }


    if (!import.pegas) {
      input <- dplyr::bind_cols(
        input,
        parse_genomic(x = "GT", data = vcf.data, return.alleles = TRUE,
                      verbose = verbose))
      # }

      if (tibble::has_name(input, "ID")) {
        input <- data.table::melt.data.table(
          data = data.table::as.data.table(input),
          id.vars = c("MARKERS", "CHROM", "LOCUS", "POS", "ID", "COL", "REF", "ALT"),
          variable.name = "INDIVIDUALS",
          variable.factor = FALSE,
          value.name = "GT"
        ) %>%
          tibble::as_data_frame()
      } else {
        input <- data.table::melt.data.table(
          data = data.table::as.data.table(input),
          id.vars = c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT"),
          variable.name = "INDIVIDUALS",
          variable.factor = FALSE,
          value.name = "GT"
        ) %>%
          tibble::as_data_frame()
      }
      # metadata
      if (is.logical(vcf.metadata)) {
        overwrite.metadata <- NULL
      } else {
        overwrite.metadata <- vcf.metadata
        vcf.metadata <- TRUE
      }

      if (vcf.metadata) {
        if (verbose) message("Keeping vcf metadata: yes")
        # detect FORMAT fields available
        have <- suppressWarnings(vcfR::vcf_field_names(vcf.data, tag = "FORMAT")$ID)
        # current version doesn't deal well with PL with 3 fields separated with ","
        want <- c("DP", "AD", "GL", "PL", "GQ", "GOF", "NR", "NV")
        if (!is.null(overwrite.metadata)) want <- overwrite.metadata
        parse.format.list <- purrr::keep(.x = have, .p = have %in% want)

        # work on parallelization of this part
        # if (length(parse.format.list) <= 2) {
        input <- dplyr::bind_cols(
          input,
          purrr::map(parse.format.list, parse_genomic, data = vcf.data,
                     gather.data = TRUE, verbose = verbose) %>%
            dplyr::bind_cols(.))

        # some VCF are too big to fit in memory multiple times for parallel processing...
        # } else {
        # if (length(parse.format.list) < parallel.core) {
        # parallel.core.parse <- length(parse.format.list)
        # } else {
        # parallel.core.parse <- parallel.core
        # }
        #   if (verbose) message("Parsing and tidying: ", stringi::stri_join(parse.format.list, collapse = ", "))
        #   parsed.format <- list()
        #   parsed.format <- .radiator_parallel(
        #   # parsed.format <- parallel::mclapply(
        #     X = parse.format.list,
        #     FUN = parse_genomic,
        #     mc.cores = parallel.core.parse,
        #     data = vcf.data,
        #     gather.data = TRUE
        #   ) %>% dplyr::bind_cols(.)
        #   input <- dplyr::bind_cols(input, parsed.format)
        #   parsed.format <- NULL
        # }
      } else {
        if (verbose) message("Keeping vcf metadata: no")
      }
      vcf.data <- NULL
    }#End vcfR GT and metadata import


    # Population levels and strata
    if (verbose) message("Making the vcf population-wise...")
    strata.df$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = strata.df$INDIVIDUALS,
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )

    # Filter blacklisted individuals
    if (!is.null(blacklist.id)) {
      input <- dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS")
      strata.df <- dplyr::anti_join(strata.df, blacklist.id, by = "INDIVIDUALS")
    }

    # check that names match between strata.df and input before going further
    if (!identical(sort(unique(input$INDIVIDUALS)), sort(unique(strata.df$INDIVIDUALS)))) {
      stop("The individuals in the strata file don't match the individuals in the vcf file")
    }

    input <- dplyr::left_join(x = input, y = strata.df, by = "INDIVIDUALS")

    # Using pop.levels and pop.labels info if present
    input <- change_pop_names(
      data = input, pop.levels = pop.levels, pop.labels = pop.labels)

    # Pop select
    if (!is.null(pop.select)) {
      if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
      input$POP_ID <- droplevels(input$POP_ID)
    }

    # Haplotypes or biallelic VCF----------------------------------------------
    # recoding genotype
    if (biallelic) {# biallelic VCF
      if (verbose) message("Recoding bi-allelic VCF...")
      input <- dplyr::rename(input, GT_VCF_NUC = GT)
    } else {#multi-allelic vcf
      if (verbose) message("Recoding VCF haplotype...")
      input <- dplyr::rename(input, GT_HAPLO = GT)
    }

    # split data for 3 rounds of CPU
    # to work, same markers in same split...
    if (verbose) message("Calculating REF/ALT alleles...")
    input <- radiator::change_alleles(
      data = input,
      monomorphic.out = FALSE,
      biallelic = biallelic,
      parallel.core = parallel.core,
      verbose = verbose)$input

    # Re ordering columns
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "ID", "COL", "INDIVIDUALS", "POP_ID",
              "REF", "ALT", "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN",
              "POLYMORPHIC")

    input <- suppressWarnings(
      dplyr::select(input, dplyr::one_of(want), dplyr::everything()))

  } # End import VCF

  # Import PLINK ---------------------------------------------------------------
  if (data.type == "plink.file") { # PLINK
    if (verbose) message("Importing the PLINK files...")
    tfam <- data.table::fread(
      input = stringi::stri_replace_all_fixed(
        str = data,
        pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE
      ),
      sep = " ",
      header = FALSE,
      stringsAsFactors = FALSE,
      verbose = FALSE,
      select = c(1,2),
      colClasses = list(character = c(1,2)),
      col.names = c("POP_ID", "INDIVIDUALS"),
      showProgress = TRUE,
      data.table = FALSE) %>%
      tibble::as_data_frame() %>%
      dplyr::mutate(
        # remove unwanted sep in individual name and replace with "-"
        INDIVIDUALS = stringi::stri_replace_all_fixed(
          str = INDIVIDUALS,
          pattern = c("_", ":"),
          replacement = c("-", "-"),
          vectorize_all = FALSE),
        # remove potential whitespace in tfam pop id column
        POP_ID = stringi::stri_replace_all_fixed(
          POP_ID,
          pattern = " ",
          replacement = "_",
          vectorize_all = FALSE)
      )

    # if no strata tfam = strata.df
    if (is.null(strata)) {
      strata.df <- tfam

      # Check with strata and pop.levels/pop.labels
      if (!is.null(pop.levels)) {
        if (length(levels(factor(strata.df$POP_ID))) != length(pop.levels)) {
          stop("The number of groups in your tfam file must match the number of groups in pop.levels")
        }
      }
    } else {
      # remove unwanted sep in individual name and replace with "-"
      strata.df$INDIVIDUALS <- stringi::stri_replace_all_fixed(
        str = strata.df$INDIVIDUALS,
        pattern = c("_", ":"),
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )
    }

    tped.header.prep <- tfam %>%
      dplyr::select(INDIVIDUALS) %>%
      dplyr::mutate(
        NUMBER = seq(1, n()),
        ALLELE1 = rep("A1", n()), ALLELE2 = rep("A2", n())
      ) %>%
      tidyr::gather(ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, NUMBER)) %>%
      dplyr::arrange(NUMBER) %>%
      dplyr::select(-ALLELES_GROUP) %>%
      tidyr::unite(INDIVIDUALS_ALLELES, c(INDIVIDUALS, ALLELES), sep = "_", remove = FALSE) %>%
      dplyr::arrange(NUMBER) %>%
      dplyr::mutate(NUMBER = seq(from = (1 + 4), to = n() + 4)) %>%
      dplyr::select(-ALLELES)

    tped.header.names <- c("LOCUS", tped.header.prep$INDIVIDUALS_ALLELES)
    tped.header.integer <- c(2, tped.header.prep$NUMBER)

    if (!is.null(blacklist.id)) { # using the blacklist of individuals
      whitelist.id <- tped.header.prep %>%
        dplyr::anti_join(blacklist.id, by = "INDIVIDUALS") %>%
        dplyr::arrange(NUMBER)
      tped.header.names <- c("LOCUS", whitelist.id$INDIVIDUALS_ALLELES)
      tped.header.integer <- c(2, whitelist.id$NUMBER)

      strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
    }

    # import PLINK
    input <- data.table::fread(
      input = data,
      sep = " ",
      header = FALSE,
      stringsAsFactors = FALSE,
      verbose = FALSE,
      select = tped.header.integer,
      col.names = tped.header.names,
      showProgress = TRUE,
      data.table = FALSE) %>%
      tibble::as_data_frame() %>%
      dplyr::mutate(LOCUS = as.character(LOCUS))

    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      if (verbose) message("Filtering with whitelist of markers")
      input <- suppressWarnings(
        dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist)
      )
    }

    # To reduce the size of the dataset we subsample the markers with max.marker
    if (!is.null(max.marker)) {
      if (verbose) message("Using the max.marker to reduce the size of the dataset")
      input <- dplyr::sample_n(tbl = input, size = max(as.numeric(max.marker)), replace = FALSE)

      max.marker.subsample.select <- input %>%
        dplyr::distinct(LOCUS, .keep_all = TRUE) %>%
        dplyr::arrange(LOCUS)

      readr::write_tsv(# save results
        x = max.marker.subsample.select,
        path = "max.marker.subsample.select.tsv")
    }

    # Make tidy
    if (verbose) message("Tidying the PLINK file ...")
    # Filling GT and new separating INDIVIDUALS from ALLELES
    # combining alleles
    input <- data.table::melt.data.table(
      data = data.table::as.data.table(input),
      id.vars = "LOCUS",
      variable.name = "INDIVIDUALS_ALLELES",
      value.name = as.character("GT"),
      variable.factor = FALSE,
      value.factor = FALSE
    ) %>%
      tibble::as_data_frame()

    # detect GT coding
    if (verbose) message("Scanning for PLINK tped genotype coding")
    detect.gt.coding <- unique(sample(x = input$GT, size = 100, replace = FALSE))
    gt.letters <- c("A", "C", "G", "T")

    if (TRUE %in% unique(gt.letters %in% detect.gt.coding)) {
      if (verbose) message("    genotypes coded with letters")
      gt.letters.df <- tibble::data_frame(GT = c("A", "C", "G", "T", "0"), NEW_GT = c("001", "002", "003", "004", "000"))
      input <- dplyr::left_join(
        input,
        gt.letters.df, by = "GT") %>%
        dplyr::select(-GT) %>%
        dplyr::rename(GT = NEW_GT)
      gt.letters.df <- NULL
    } else {
      if (verbose) message("    genotypes coded with integers")
      input <- input %>%
        dplyr::mutate(GT = stringi::stri_pad_left(str = GT, width = 3, pad = "0"))
    }
    detect.gt.coding <- gt.letters <- NULL


    input <- input %>%
      tidyr::separate(
        data = .,
        col = INDIVIDUALS_ALLELES,
        into = c("INDIVIDUALS", "ALLELES"),
        sep = "_") %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula = LOCUS + INDIVIDUALS ~ ALLELES,
        value.var = "GT"
      ) %>%
      tibble::as_data_frame() %>%
      tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>%
      dplyr::select(LOCUS, INDIVIDUALS, GT)

    # population levels and strata
    if (verbose) message("Integrating the tfam/strata file...")

    input <- dplyr::left_join(x = input, y = strata.df, by = "INDIVIDUALS")

    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)

    # Pop select
    if (!is.null(pop.select)) {
      if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
      input$POP_ID <- droplevels(input$POP_ID)
    }

    # removing untyped markers across all-pop
    remove.missing.gt <- input %>%
      dplyr::select(LOCUS, GT) %>%
      dplyr::filter(GT != "000000")

    untyped.markers <- dplyr::n_distinct(input$LOCUS) - dplyr::n_distinct(remove.missing.gt$LOCUS)
    if (untyped.markers > 0) {
      if (verbose) message("Number of marker with 100 % missing genotypes: ", untyped.markers)
      input <- suppressWarnings(
        dplyr::semi_join(input,
                         remove.missing.gt %>%
                           dplyr::distinct(LOCUS, .keep_all = TRUE),
                         by = "LOCUS")
      )
    }

    # Unused objects
    tped.header.prep <- tped.header.integer <- tped.header.names <- remove.missing.gt <- NULL

    # detect if biallelic give vcf style genotypes
    # biallelic <- radiator::detect_biallelic_markers(input)
    input.temp <- radiator::change_alleles(data = input, monomorphic.out = FALSE, verbose = verbose)
    input <- input.temp$input
    biallelic <- input.temp$biallelic
  } # End import PLINK

  # Import genepop--------------------------------------------------------------
  if (data.type == "genepop.file") {
    if (verbose) message("Tidying the genepop file ...")
    data <- radiator::tidy_genepop(data = data, tidy = TRUE)
    data.type <- "tbl_df"
  }

  # Import fst.file ------------------------------------------------------------
  if (data.type == "fst.file") {
    if (verbose) message("Importing the fst.file as a data frame...")
    data <- fst::read.fst(path = data)
    data.type <- "tbl_df"
  }

  # Import DArT ----------------------------------------------------------------
  if (data.type == "dart") {
    if (verbose) message("Tidying DArT data...")
    data <- radiator::tidy_dart(
      data = data,
      strata = strata,
      verbose = FALSE,
      parallel.core = parallel.core)
    data.type <- "tbl_df"
  }

  # Import DF-------------------------------------------------------------------
  if (data.type == "tbl_df") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Importing the data frame ...")
    input <- radiator::tidy_wide(data = data, import.metadata = TRUE)

    # For long tidy format, switch LOCUS to MARKERS column name, if found MARKERS not found
    if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
      input <- dplyr::rename(.data = input, MARKERS = LOCUS)
    }

    # Change individuals names containing special character
    input$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = input$INDIVIDUALS,
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )

    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      if (verbose) message("Filtering with whitelist of markers")
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }

    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      if (verbose) message("Filtering with blacklist of individuals")
      input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }

    # population levels and strata
    if (!is.null(strata)) {
      strata.df$INDIVIDUALS <- stringi::stri_replace_all_fixed(
        str = strata.df$INDIVIDUALS,
        pattern = c("_", ":"),
        replacement = c("-", "-"),
        vectorize_all = FALSE
      )

      input <- input %>%
        dplyr::select(-POP_ID) %>%
        dplyr::left_join(strata.df, by = "INDIVIDUALS")
    }

    # Change potential problematic POP_ID space
    input$POP_ID = stringi::stri_replace_all_fixed(
      input$POP_ID,
      pattern = " ",
      replacement = "_",
      vectorize_all = FALSE
    )

    # Check with strata and pop.levels/pop.labels
    if (!is.null(pop.levels)) {
      if (length(levels(factor(input$POP_ID))) != length(pop.levels)) {
        stop("The number of groups in your POP_ID column file must match the number of groups in pop.levels")
      }
    }

    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)

    # Pop select
    if (!is.null(pop.select)) {
      if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    }

    input.temp <- radiator::change_alleles(data = input, monomorphic.out = FALSE, verbose = verbose)
    input <- input.temp$input
    biallelic <- input.temp$biallelic
  } # End import data frame of genotypes

  # Import haplo---------------------------------------------------------------
  if (data.type == "haplo.file") { # Haplotype file
    if (verbose) message("Importing STACKS haplotype file")
    number.columns <- max(utils::count.fields(data, sep = "\t"))

    input <- data.table::fread(
      input = data,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      colClasses = list(character = 1:number.columns),
      verbose = FALSE,
      showProgress = TRUE,
      data.table = FALSE,
      na.strings = "-"
    ) %>%
      tibble::as_data_frame() %>%
      dplyr::select(-Cnt)

    if (tibble::has_name(input, "# Catalog ID") || tibble::has_name(input, "Catalog ID")) {
      colnames(input) <- stringi::stri_replace_all_fixed(
        str = colnames(input),
        pattern = c("# Catalog ID", "Catalog ID"), replacement = c("LOCUS", "LOCUS"), vectorize_all = FALSE
      )
    }

    if (tibble::has_name(input, "Seg Dist")) {
      input <- dplyr::select(.data = input, -`Seg Dist`)
    }

    n.catalog.locus <- dplyr::n_distinct(input$LOCUS)
    n.individuals <- ncol(input) - 1

    message("\nNumber of loci in catalog: ", n.catalog.locus)
    message("Number of individuals: ", n.individuals)

    # tidy df, individuals in 1 column
    input <- data.table::melt.data.table(
      data = data.table::as.data.table(input),
      id.vars = "LOCUS",
      variable.name = "INDIVIDUALS",
      variable.factor = FALSE,
      value.name = "GT_HAPLO"
    ) %>%
      tibble::as_data_frame()

    number.columns <- NULL

    input$INDIVIDUALS = stringi::stri_replace_all_fixed(
      str = input$INDIVIDUALS,
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )

    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      if (verbose) message("Filtering with whitelist of markers")
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }

    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      if (verbose) message("Filtering with blacklist of individuals")
      input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }

    # remove consensus markers
    if (verbose) message("\nScanning for consensus markers...")
    consensus.markers <- input %>%
      dplyr::filter(GT_HAPLO == "consensus") %>%
      dplyr::distinct(LOCUS, .keep_all = TRUE)

    if (length(consensus.markers$LOCUS) > 0) {
      input <- suppressWarnings(dplyr::anti_join(input, consensus.markers, by = "LOCUS"))
    }
    if (verbose) message("    number of consensus markers removed: ", dplyr::n_distinct(consensus.markers$LOCUS))

    # population levels and strata
    strata.df$INDIVIDUALS = stringi::stri_replace_all_fixed(
      str = strata.df$INDIVIDUALS,
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )

    input <- dplyr::left_join(x = input, y = strata.df, by = "INDIVIDUALS")

    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)

    # Pop select
    if (!is.null(pop.select)) {
      if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    }

    # removing errors and potential paralogs (GT with > 2 alleles)
    if (verbose) message("Scanning for artifactual genotypes...")
    input <- input %>%
      dplyr::mutate(POLYMORPHISM = stringi::stri_count_fixed(GT_HAPLO, "/"))

    blacklist.paralogs <- input %>%
      dplyr::filter(POLYMORPHISM > 1) %>%
      dplyr::select(LOCUS, INDIVIDUALS)

    if (verbose) message("    number of genotypes with more than 2 alleles: ", length(blacklist.paralogs$LOCUS))
    # save.image("testing.haplo.RData")

    if (length(blacklist.paralogs$LOCUS) > 0) {
      input <- input %>%
        dplyr::mutate(GT_HAPLO = replace(GT_HAPLO, which(POLYMORPHISM > 1), NA)) %>%
        dplyr::select(-POLYMORPHISM)

      readr::write_tsv(blacklist.paralogs, "blacklist.genotypes.paralogs.tsv")
    }

    if (verbose) message("Calculating REF/ALT alleles...")

    input <- radiator::change_alleles(
      data = input,
      monomorphic.out = monomorphic.out,
      parallel.core = parallel.core,
      verbose = verbose)$input

    biallelic <- FALSE
    input <- dplyr::rename(input, LOCUS = MARKERS)
  } # End import haplotypes file

  # Import GENIND--------------------------------------------------------------
  if (data.type == "genind") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the genind object ...")
    input <- adegenet::genind2df(data) %>%
      tibble::rownames_to_column("INDIVIDUALS") %>%
      dplyr::rename(POP_ID = pop)

    # scan for the number of character coding the allele
    allele.sep <- input %>% dplyr::select(-INDIVIDUALS, -POP_ID)
    allele.sep <- unique(nchar(allele.sep[!is.na(allele.sep)]))

    if (length(allele.sep) > 1) {
      stop("The number of character/integer string coding the allele is not identical accross markers")
    }

    input <- input %>%
      tidyr::gather(key = LOCUS, value = GT, -c(INDIVIDUALS, POP_ID)) %>%
      tidyr::separate(
        data = ., col = GT, into = c("A1", "A2"),
        sep = allele.sep/2, remove = TRUE, extra = "drop"
      ) %>%
      dplyr::mutate(
        A1 = stringi::stri_pad_left(str = A1, pad = "0", width = 3),
        A2 = stringi::stri_pad_left(str = A2, pad = "0", width = 3)
      ) %>%
      tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>%
      dplyr::mutate(GT = replace(GT, which(GT == "NANA"), "000000"))

    # remove unwanted sep in id and pop.id names
    input <- input %>%
      dplyr::mutate(
        INDIVIDUALS = stringi::stri_replace_all_fixed(
          str = INDIVIDUALS,
          pattern = c("_", ":"),
          replacement = c("-", "-"),
          vectorize_all = FALSE),
        POP_ID = stringi::stri_replace_all_fixed(
          POP_ID,
          pattern = " ",
          replacement = "_",
          vectorize_all = FALSE)
      )

    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      if (verbose) message("Filtering with whitelist of markers")
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }

    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      if (verbose) message("Filtering with blacklist of individuals")
      input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }

    # population levels and strata
    if (!is.null(strata)) {
      input <- input %>%
        dplyr::select(-POP_ID) %>%
        dplyr::mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>%
        dplyr::left_join(strata.df, by = "INDIVIDUALS")
    }

    # Change potential problematic POP_ID space
    input$POP_ID = stringi::stri_replace_all_fixed(
      input$POP_ID,
      pattern = " ",
      replacement = "_",
      vectorize_all = FALSE
    )

    # Check with strata and pop.levels/pop.labels
    if (!is.null(pop.levels)) {
      if (length(levels(factor(input$POP_ID))) != length(pop.levels)) {
        stop("The number of groups in your POP_ID column must match the number of groups in pop.levels")
      }
    }

    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)

    # Pop select
    if (!is.null(pop.select)) {
      if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    }

    # detect if biallelic give vcf style genotypes
    # biallelic <- radiator::detect_biallelic_markers(input)
    input.temp <- radiator::change_alleles(data = input, monomorphic.out = FALSE, verbose = verbose)
    input <- input.temp$input
    biallelic <- input.temp$biallelic

    # Now the genind and genepop are like ordinary data frames
    data.type <- "tbl_df" # for subsequent steps

  } # End tidy genind

  # Import GENLIGHT ------------------------------------------------------------
  if (data.type == "genlight") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the genlight object ...")
    input <- radiator::tidy_genlight(data = data)

    # remove unwanted sep in id and pop.id names
    input <- input %>%
      dplyr::mutate(
        INDIVIDUALS = stringi::stri_replace_all_fixed(
          str = INDIVIDUALS,
          pattern = c("_", ":"),
          replacement = c("-", "-"),
          vectorize_all = FALSE),
        POP_ID = stringi::stri_replace_all_fixed(
          POP_ID,
          pattern = " ",
          replacement = "_",
          vectorize_all = FALSE)
      )

    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      if (verbose) message("Filtering with whitelist of markers")
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }

    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      if (verbose) message("Filtering with blacklist of individuals")
      input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }

    # population levels and strata
    if (!is.null(strata)) {
      input <- input %>%
        dplyr::select(-POP_ID) %>%
        dplyr::mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>%
        dplyr::left_join(strata.df, by = "INDIVIDUALS")
    }

    # Change potential problematic POP_ID space
    input$POP_ID = stringi::stri_replace_all_fixed(
      input$POP_ID,
      pattern = " ",
      replacement = "_",
      vectorize_all = FALSE
    )

    # Check with strata and pop.levels/pop.labels
    if (!is.null(pop.levels)) {
      if (length(levels(factor(input$POP_ID))) != length(pop.levels)) {
        stop("The number of groups in your POP_ID column must match the number of groups in pop.levels")
      }
    }

    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)

    # Pop select
    if (!is.null(pop.select)) {
      if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    }

    # Now the genind and genepop are like ordinary data frames
    data.type <- "tbl_df" # for subsequent steps
    biallelic <- TRUE
  } # End tidy genlight

  # Import STRATAG gtypes ------------------------------------------------------
  if (data.type == "gtypes") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the gtypes object ...")
    input <- tidy_gtypes(data) %>%
      dplyr::mutate(
        INDIVIDUALS = stringi::stri_replace_all_fixed(
          str = INDIVIDUALS,
          pattern = c("_", ":"),
          replacement = c("-", "-"),
          vectorize_all = FALSE),
        POP_ID = stringi::stri_replace_all_fixed(
          POP_ID,
          pattern = " ",
          replacement = "_",
          vectorize_all = FALSE)
      )

    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      if (verbose) message("Filtering with whitelist of markers")
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }

    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      if (verbose) message("Filtering with blacklist of individuals")
      input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }

    # population levels and strata
    if (!is.null(strata)) {
      input <- input %>%
        dplyr::select(-POP_ID) %>%
        dplyr::mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>%
        dplyr::left_join(strata.df, by = "INDIVIDUALS")
    }

    # Change potential problematic POP_ID space
    input$POP_ID = stringi::stri_replace_all_fixed(
      input$POP_ID,
      pattern = " ",
      replacement = "_",
      vectorize_all = FALSE
    )

    # Check with strata and pop.levels/pop.labels
    if (!is.null(pop.levels)) {
      if (length(levels(factor(input$POP_ID))) != length(pop.levels)) {
        stop("The number of groups in your POP_ID column must match the number of groups in pop.levels")
      }
    }

    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)

    # Pop select
    if (!is.null(pop.select)) {
      if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
    }

    # detect if biallelic give vcf style genotypes
    # biallelic <- radiator::detect_biallelic_markers(input)
    input.temp <- radiator::change_alleles(data = input, monomorphic.out = FALSE, verbose = verbose)
    input <- input.temp$input
    biallelic <- input.temp$biallelic

    # Now the gtypes are like ordinary data frames
    data.type <- "tbl_df" # for subsequent steps

  } # End tidy gtypes


  # END IMPORT DATA

  # Arrange the id and create a strata after pop select ------------------------
  input$INDIVIDUALS <- stringi::stri_replace_all_fixed(
    str = input$INDIVIDUALS,
    pattern = c("_", ":"),
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )

  strata.df <- dplyr::ungroup(input) %>%
    dplyr::distinct(POP_ID, INDIVIDUALS)

  # Blacklist genotypes --------------------------------------------------------
  if (is.null(blacklist.genotype)) { # no Whitelist
    if (verbose) message("Erasing genotype: no")
  } else {
    if (verbose) message("Erasing genotype: yes")
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS")
    suppressWarnings(suppressMessages(
      blacklist.genotype <- readr::read_tsv(blacklist.genotype, col_names = TRUE) %>%
        dplyr::mutate(
          INDIVIDUALS = stringi::stri_replace_all_fixed(
            str = INDIVIDUALS,
            pattern = c("_", ":"),
            replacement = c("-", "-"),
            vectorize_all = FALSE
          )
        ) %>%
        dplyr::select(dplyr::one_of(want)) %>%
        dplyr::mutate_all(.tbl = ., .funs = as.character, exclude = NA)))
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)

    if (data.type == "haplo.file") {
      blacklist.genotype <- dplyr::select(.data = blacklist.genotype, INDIVIDUALS, LOCUS)
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }

    # control check to keep only individuals in the strata.df
    blacklist.genotype <- suppressWarnings(
      blacklist.genotype  %>%
        dplyr::filter(INDIVIDUALS %in% strata.df$INDIVIDUALS)
    )

    # control check to keep only whitelisted markers from the blacklist of genotypes
    if (!is.null(whitelist.markers)) {
      blacklist.genotype <- blacklist.genotype
      if (verbose) message("Control check to keep only whitelisted markers present in the blacklist of genotypes to erase.")
      # updating the whitelist of markers to have all columns that id markers
      if (data.type == "vcf.file") {
        whitelist.markers.ind <- input %>% dplyr::distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      } else {
        whitelist.markers.ind <- input %>% dplyr::distinct(LOCUS, INDIVIDUALS)
      }

      # updating the blacklist.genotype
      blacklist.genotype <- suppressWarnings(dplyr::semi_join(whitelist.markers.ind, blacklist.genotype, by = columns.names.blacklist.genotype))
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }

    # Update column names
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)

    input.erase <- dplyr::semi_join(input, blacklist.genotype, by = columns.names.blacklist.genotype) %>%
      dplyr::mutate(GT = rep("000000", n()))
    input <- dplyr::anti_join(input, blacklist.genotype, by = columns.names.blacklist.genotype)
    if (tibble::has_name(input.erase, "GT_VCF")) {
      input.erase <- dplyr::mutate(input.erase, GT_VCF = rep("./.", n()))
    }

    if (tibble::has_name(input.erase, "GT_VCF_NUC")) {
      input.erase <- dplyr::mutate(input.erase, GT_VCF_NUC = rep("./.", n()))
    }

    if (tibble::has_name(input.erase, "GT_BIN")) {
      input.erase <- dplyr::mutate(input.erase, GT_BIN = rep(as.numeric(NA_character_), n()))
    }
    input <- dplyr::bind_rows(input, input.erase)
  } # End erase genotypes

  # dump unused object
  blacklist.id <- whitelist.markers <- whitelist.markers.ind <- NULL
  want <- blacklist.genotype <- NULL
  # SNP LD  --------------------------------------------------------------------
  if (!is.null(snp.ld)) {
    input <- radiator::snp_ld(data = input, snp.ld = snp.ld)
    # if (!tibble::has_name(input, "POS")) {
    #   stop("snp.ld is only available for VCF file, use radiator package for
    #        haplotype file and create a whitelist, for other file type, use
    #        SNPRelate package or PLINK linkage disequilibrium based SNP pruning
    #        option")
    # }
    # if (verbose) message("Minimizing LD...")
    # snp.locus <- input %>% dplyr::distinct(LOCUS, POS)
    #
    # # Random selection
    # if (snp.ld == "random") {
    #   snp.select <- snp.locus %>%
    #     dplyr::group_by(LOCUS) %>%
    #     sample_n(size = 1, replace = FALSE)
    #   message(
    #     "Number of original SNP = ",
    #     dplyr::n_distinct(snp.locus$POS), "\n",
    #     "Number of SNP randomly selected to keep 1 SNP per read/haplotype = ",
    #     dplyr::n_distinct(snp.select$POS), "\n", "Number of SNP removed = ",
    #     dplyr::n_distinct(snp.locus$POS) - dplyr::n_distinct(snp.select$POS))
    # }
    #
    # # Fist SNP on the read
    # if (snp.ld == "first") {
    #   snp.select <- snp.locus %>%
    #     dplyr::group_by(LOCUS) %>%
    #     dplyr::summarise(POS = min(POS))
    #   message(
    #     "Number of original SNP = ",
    #     dplyr::n_distinct(snp.locus$POS), "\n",
    #     "Number of SNP after keeping the first SNP on the read/haplotype = ",
    #     dplyr::n_distinct(snp.select$POS), "\n", "Number of SNP removed = ",
    #     dplyr::n_distinct(snp.locus$POS) - dplyr::n_distinct(snp.select$POS))
    # }
    #
    # # Last SNP on the read
    # if (snp.ld == "last") {
    #   snp.select <- snp.locus %>%
    #     dplyr::group_by(LOCUS) %>%
    #     dplyr::summarise(POS = max(POS))
    #   message(
    #     "Number of original SNP = ", dplyr::n_distinct(snp.locus$POS), "\n",
    #     "Number of SNP after keeping the first SNP on the read/haplotype = ",
    #     dplyr::n_distinct(snp.select$POS), "\n", "Number of SNP removed = ",
    #     dplyr::n_distinct(snp.locus$POS) - dplyr::n_distinct(snp.select$POS))
    # }
    #
    # # filtering the VCF to minimize LD
    # input <- input %>% dplyr::semi_join(snp.select, by = c("LOCUS", "POS"))
    # if (verbose) message("Filtering the tidy VCF to minimize LD by keeping only 1 SNP per short read/haplotype")
  } # End of snp.ld control

  # Unique markers id ----------------------------------------------------------
  # we want to keep LOCUS in the vcf, but not in the other type of input file
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = "LOCUS",
      replacement = "MARKERS",
      vectorize_all = FALSE
    )
  }

  # Removing special characters in markers id ----------------------------------
  input <- input %>%
    dplyr::mutate(
      MARKERS = stringi::stri_replace_all_fixed(
        str = as.character(MARKERS),
        pattern = c("/", ":", "-", "."),
        replacement = "_",
        vectorize_all = FALSE)
    )

  # Markers in common between all populations (optional) -----------------------
  if (common.markers) { # keep only markers present in all pop
    input <- radiator::keep_common_markers(input, verbose = TRUE)$input
  } # End common markers

  # Removing monomorphic markers------------------------------------------------
  if (monomorphic.out) {
    if (verbose) message("Removing monomorphic markers: yes")
    if (tibble::has_name(input, "POLYMORPHIC")) {
      if (verbose) message("Scanning for monomorphic markers...")
      if (verbose) message("    Number of markers before = ", dplyr::n_distinct(input$MARKERS))
      if (tibble::has_name(input, "POS")) {
        mono.markers <- dplyr::filter(input, !POLYMORPHIC) %>%
          dplyr::distinct(MARKERS, CHROM, LOCUS, POS)
      } else {
        mono.markers <- dplyr::filter(input, !POLYMORPHIC) %>%
          dplyr::distinct(MARKERS)
      }

      if (nrow(mono.markers) > 0) {
        if (verbose) message("    Number of monomorphic markers removed = ", nrow(mono.markers))
        input <- dplyr::filter(input, POLYMORPHIC)
        readr::write_tsv(mono.markers, "blacklist.monomorphic.markers.tsv")
        if (verbose) message("    Number of markers after = ", dplyr::n_distinct(input$MARKERS))
      }
      input <- dplyr::select(input, -POLYMORPHIC)
    } else {
      mono.out <- radiator::discard_monomorphic_markers(input, verbose = TRUE)
      mono.markers <- mono.out$blacklist.monomorphic.markers
      if (nrow(mono.markers) > 0) {
        # if (data.type == "haplo.file") {
        #   mono.markers <- dplyr::rename(.data = mono.markers, LOCUS = MARKERS)
        # }
        input <- mono.out$input
        readr::write_tsv(mono.markers, "blacklist.monomorphic.markers.tsv")
      }
      mono.out <- NULL
    }
    mono.markers <- NULL
  } # End monomorphic out

  # Minor Allele Frequency filter ----------------------------------------------
  # maf.thresholds <- c(0.05, 0.1) # test
  if (!is.null(maf.thresholds)) { # with MAF
    maf.info <- radiator_maf_module(
      data = input,
      maf.thresholds = maf.thresholds,
      maf.pop.num.threshold = maf.pop.num.threshold,
      maf.approach = maf.approach,
      maf.operator = maf.operator
    )

    input <- maf.info$input
    # maf.data <- maf.info$maf.data
    maf.info <- NULL
  } # End of MAF filters

  # More VCF cleaning here -----------------------------------------------------
  # check, parse and clean FORMAT columns
  # Note to myself: why not do this at the same time as importing from vcf object?
  # Some software that produce vcf do strange thing and don't follow convention
  # You need the GT field to clean correctly the remaining fields...

  if (data.type == "vcf.file" && vcf.metadata) {

    # for parallel cleaning
    # system.time(split.vec <- dplyr::ntile(x = 1:nrow(input), n = parallel.core * 3))
    n.row <- nrow(input)
    # as.integer is usually twice as light as numeric vector...
    split.vec <- as.integer(floor((parallel.core * 3 * (1:n.row - 1) / n.row) + 1))
    n.row <- NULL

    # chunk <- parallel.core * 3
    # num.row <- nrow(input)
    # split.vec <- rep(1:ceiling(num.row/chunk), each = chunk)[1:num.row]

    # Detect the columns that need a check, parsing and cleaning
    # have <- colnames(input)
    # want <- c("DP", "AD", "GL", "PL", "GQ", "GOF", "NR", "NV")
    # clean.columns <- purrr::keep(.x = have, .p = have %in% want)

    # Fast cleaning columns
    replace_by_na <- function(data, what = ".") {
      replace(data, which(data == what), NA)
    }

    input <-  dplyr::mutate_at(
      .tbl = input, .vars = parse.format.list, .funs =  replace_by_na)

    # Cleaning AD (ALLELES_DEPTH)
    if (tibble::has_name(input, "AD")) {

      if (verbose) message("AD column: splitting coverage info into ALLELE_REF_DEPTH and ALLELE_ALT_DEPTH")

      clean_ad <- function(x) {
        res <- suppressWarnings(
          x %>%
            dplyr::mutate(
              AD = dplyr::if_else(GT_VCF == "./.", NA_character_, AD)) %>%
            tidyr::separate(AD, c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
                            sep = ",", extra = "drop") %>%
            dplyr::mutate(
              ALLELE_REF_DEPTH = as.numeric(
                stringi::stri_replace_all_regex(
                  ALLELE_REF_DEPTH, "^0$", "NA", vectorize_all = TRUE)),
              ALLELE_ALT_DEPTH = as.numeric(
                stringi::stri_replace_all_regex(
                  ALLELE_ALT_DEPTH, "^0$", "NA", vectorize_all = TRUE))
            ) %>%
            dplyr::select(-GT_VCF)
        )
        return(res)
      }#End clean_ad

      input <- dplyr::bind_cols(
        input,
        dplyr::ungroup(input) %>%
          dplyr::select(GT_VCF, AD) %>%
          split(x = ., f = split.vec) %>%
          .radiator_parallel(
            # parallel::mclapply(
            X = ., FUN = clean_ad, mc.cores = parallel.core) %>%
          dplyr::bind_rows(.))
    }#End cleaning AD column

    # Cleaning DP and changing name to READ_DEPTH
    if (tibble::has_name(input, "DP")) {
      if (verbose) message("DP column: cleaning and renaming to READ_DEPTH")
      input <- dplyr::rename(.data = input, READ_DEPTH = DP) %>%
        dplyr::mutate(
          READ_DEPTH = dplyr::if_else(GT_VCF == "./.", as.numeric(NA_character_),
                                      as.numeric(READ_DEPTH))
        )
    }#End cleaning DP column

    # PL for biallelic as 3 values:
    if (tibble::has_name(input, "PL")) {
      if (verbose) message("PL column (normalized, phred-scaled likelihoods for genotypes): separating into PROB_HOM_REF, PROB_HET and PROB_HOM_ALT")
      # Value 1: probability that the site is homozgyous REF
      # Value 2: probability that the sample is heterzygous
      # Value 2: probability that it is homozygous ALT

      clean_pl <- function(x) {
        res <- x %>%
          dplyr::mutate(
            PL = dplyr::if_else(GT_VCF == "./.", NA_character_, PL)) %>%
          tidyr::separate(
            data = ., PL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
            sep = ",", extra = "drop", remove = FALSE) %>%
          dplyr::mutate_at(
            .tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
            .funs = as.numeric) %>%
          dplyr::select(-GT_VCF)
        return(res)
      }#End clean_pl

      input <- dplyr::bind_cols(
        dplyr::select(input, -PL),
        dplyr::ungroup(input) %>%
          dplyr::select(GT_VCF, PL) %>%
          split(x = ., f = split.vec) %>%
          .radiator_parallel(
            # parallel::mclapply(
            X = ., FUN = clean_pl, mc.cores = parallel.core) %>%
          dplyr::bind_rows(.))
    }#End cleaning PL column

    # GL cleaning
    if (tibble::has_name(input, "GL")) {
      if (verbose) message("GL column: cleaning Genotype Likelihood column")
      input <- input %>%
        dplyr::mutate(
          GL = dplyr::if_else(GT_VCF == "./.", NA_character_, GL),
          GL = suppressWarnings(stringi::stri_replace_all_fixed(GL, c(".,.,.", ".,", ",."), c("NA", "", ""), vectorize_all = FALSE))
        )

      # check GL and new stacks version with no GL
      all.missing <- all(is.na(input$GL))

      if (!all.missing) {
        gl.clean <- max(
          unique(stringi::stri_count_fixed(
            str = unique(sample(x = input$GL, size = 100, replace = FALSE)),
            pattern = ",")
          ), na.rm = TRUE
        )

        if (gl.clean == 2) {
          if (verbose) message("GL column: separating into PROB_HOM_REF, PROB_HET and PROB_HOM_ALT")
          # Value 1: probability that the site is homozgyous REF
          # Value 2: probability that the sample is heterzygous
          # Value 2: probability that it is homozygous ALT
          # system.time(input2 <- input %>%
          #   tidyr::separate(data = ., GL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"), sep = ",", extra = "drop", remove = FALSE) %>%
          #   dplyr::mutate_at(.tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"), .funs = as.numeric)
          # )
          clean_gl <- function(x) {
            res <- x %>%
              tidyr::separate(
                data = ., GL, c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
                sep = ",", extra = "drop", remove = FALSE) %>%
              dplyr::mutate_at(
                .tbl = ., .vars = c("PROB_HOM_REF", "PROB_HET", "PROB_HOM_ALT"),
                .funs = as.numeric)
            return(res)
          }
          input <- dplyr::bind_cols(
            dplyr::select(input, -GL),
            dplyr::ungroup(input) %>%
              dplyr::select(GL) %>%
              split(x = ., f = split.vec) %>%
              .radiator_parallel(
                # parallel::mclapply(
                X = ., FUN = clean_gl, mc.cores = parallel.core) %>%
              dplyr::bind_rows(.))

        } else {
          input$GL <- suppressWarnings(as.numeric(input$GL))
        }
        gl.clean <- NULL
      } else {
        input <- dplyr::select(input, -GL)
      }
    }#End cleaning GL column

    # Cleaning GQ: Genotype quality as phred score
    if (tibble::has_name(input, "GQ")) {
      if (verbose) message("GQ column: Genotype Quality")
      input <- dplyr::mutate(
        input,
        GQ = dplyr::if_else(GT_VCF == "./.", as.numeric(NA_character_), as.numeric(GQ))
      )
    }#End cleaning GQ column

    # Cleaning GOF: Goodness of fit value
    if (tibble::has_name(input, "GOF")) {
      if (verbose) message("GOF column: Goodness of fit value")
      input <- dplyr::mutate(
        input,
        GOF = dplyr::if_else(GT_VCF == "./.", as.numeric(NA_character_), as.numeric(GOF))
      )
    }#End cleaning GOF column

    # Cleaning NR: Number of reads covering variant location in this sample
    if (tibble::has_name(input, "NR")) {
      if (verbose) message("NR column: splitting column into the number of variant")
      nr.col <- max(unique(stringi::stri_count_fixed(str = unique(input$NR), pattern = ","))) + 1
      nr.col.names <- stringi::stri_join(rep("NR_", nr.col), seq(1:nr.col))
      nr.col <- NULL

      clean_nr <- function(x, nr.col.names = NULL) {
        res <- tidyr::separate(data = x, col = NR, into = nr.col.names,
                               sep = ",", extra = "drop", remove = FALSE)
        return(res)
      }#End clean_nr

      input <- dplyr::bind_cols(
        dplyr::select(input, -NR),
        dplyr::ungroup(input) %>%
          dplyr::mutate(NR = dplyr::if_else(GT_VCF == "./.", NA_character_, NR)) %>%
          dplyr::select(NR) %>%
          split(x = ., f = split.vec) %>%
          .radiator_parallel_mc(
            # parallel::mclapply(
            X = ., FUN = clean_nr, mc.cores = parallel.core,
            nr.col.names = nr.col.names) %>%
          dplyr::bind_rows(.))
    }#End cleaning NR column

    # Cleaning NV: Number of reads containing variant in this sample
    if (tibble::has_name(input, "NV")) {
      if (verbose) message("NV column: splitting column into the number of variant")
      nv.col <- max(unique(stringi::stri_count_fixed(str = unique(input$NV), pattern = ","))) + 1
      nv.col.names <- stringi::stri_join(rep("NV_", nv.col), seq(1:nv.col))
      nv.col <- NULL

      clean_nv <- function(x, nv.col.names = NULL) {
        res <- tidyr::separate(
          data = x, col = NV, into = nv.col.names,
          sep = ",", extra = "drop", remove = FALSE)
        return(res)
      }

      input <- dplyr::bind_cols(
        dplyr::select(input, -NV),
        dplyr::ungroup(input) %>%
          dplyr::mutate(NV = dplyr::if_else(GT_VCF == "./.", NA_character_, NV)) %>%
          dplyr::select(NV) %>%
          split(x = ., f = split.vec) %>%
          .radiator_parallel_mc(
            # parallel::mclapply(
            X = ., FUN = clean_nv, mc.cores = parallel.core,
            nv.col.names = nv.col.names) %>%
          dplyr::bind_rows(.))
    }#End cleaning NV column

    split.vec <- NULL

    # Re ordering columns
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID",
              "REF", "ALT", "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN")

    input <- suppressWarnings(
      dplyr::select(input, dplyr::one_of(want), dplyr::everything()))

  }# end cleaning columns


  # Write to working directory -------------------------------------------------
  # if (!is.null(filename)) {
  #   if (verbose) message("Writing the tidy data to the working directory: \n", filename)
  #   readr::write_tsv(x = input, path = filename, col_names = TRUE)
  # }

  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".rad")
    message("Writing tidy data set: ", tidy.name)
    # if (!is.null(save.feather)) {
    # feather::write_feather(filter, stri_replace_all_fixed(filename, pattern = ".tsv", replacement = "_feather.tsv", vectorize_all = TRUE))
    # } else {
    fst::write.fst(x = input, path = tidy.name, compress = 85)
    # }
  }

  # Results --------------------------------------------------------------------
  # messages
  n.markers <- dplyr::n_distinct(input$MARKERS)
  if (tibble::has_name(input, "CHROM")) {
    n.chromosome <- dplyr::n_distinct(input$CHROM)
  } else {
    n.chromosome <- "no chromosome info"
  }
  n.individuals <- dplyr::n_distinct(input$INDIVIDUALS)
  n.pop <- dplyr::n_distinct(input$POP_ID)

  if (verbose) {
    cat("############################### RESULTS ###############################\n")

    if (!is.null(filename)) {
      message("Tidy data written in global environment and working directory")
    } else {
      message("Tidy data written in global environment")
    }
    message("Data format: ", data.type)
    if (biallelic) {
      message("Biallelic data")
    } else{
      message("Multiallelic data")
    }
    if (common.markers) {
      message("Number of common markers: ", n.markers)
    } else {
      message("Number of markers: ", n.markers)
    }
    message("Number of chromosome/contig/scaffold: ", n.chromosome)
    message("Number of individuals: ", n.individuals)
    message("Number of populations: ", n.pop)
    timing <- proc.time() - timing
    message("Computation time: ", round(timing[[3]]), " sec")
    cat("################ radiator::tidy_genomic_data completed ################\n")
  }
  res <- input
  options(width = opt.change)
  return(res)
} # tidy genomic data


# Internal nested Function -----------------------------------------------------

#' @title parse_genomic
#' @description function to parse the format field and tidy the results of VCF
#' @rdname parse_genomic
#' @keywords internal
#' @export
parse_genomic <- function(
  x, data = NULL, mask = FALSE, gather.data = FALSE, return.alleles = FALSE,
  verbose = TRUE) {
  format.name <- x

  if (verbose) message("Parsing and tidying: ", format.name)
  x <- tibble::as_data_frame(vcfR::extract.gt(
    x = data, element = x, mask = mask, return.alleles = return.alleles,
    IDtoRowNames = FALSE, convertNA = FALSE))

  if (format.name == "GT") {
    colnames(x) <- stringi::stri_replace_all_fixed(
      str = colnames(x),
      pattern = c("_", ":"),
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )
  }

  if (gather.data) {
    x <- dplyr::mutate(x, ID = seq(1, n()))
    x <- data.table::melt.data.table(
      data = data.table::as.data.table(x),
      id.vars = c("ID"),
      variable.name = "INDIVIDUALS",
      variable.factor = FALSE,
      value.name = format.name
    ) %>%
      tibble::as_data_frame() %>%
      dplyr::select(-ID, -INDIVIDUALS)
  }
  return(x)
}#End parse_genomic


#' @title nuc2integers
#' @description convert nucleotides to integers. Useful while tidying VCF
#' @rdname nuc2integers
#' @keywords internal
#' @export
nuc2integers <- function(x) {

  res <- dplyr::select(x, MARKERS, GT_VCF_NUC) %>%
    dplyr::filter(GT_VCF_NUC != "./.") %>%
    tidyr::separate(col = GT_VCF_NUC, into = c("A1", "A2"), sep = "/") %>%
    tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -MARKERS) %>%
    dplyr::select(-ALLELES_GROUP) %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::tally(.) %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(INTEGERS = seq(0, n() - 1)) %>%
    dplyr::select(-n) %>%
    dplyr::arrange(MARKERS, INTEGERS) %>%
    dplyr::ungroup(.)


  ref.alleles <- res %>%
    dplyr::filter(INTEGERS == 0) %>%
    dplyr::mutate(REF = ALLELES) %>%
    dplyr::distinct(MARKERS, REF)

  alt.alleles <- res %>%
    dplyr::filter(INTEGERS != 0) %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::mutate(ALT = stringi::stri_join(ALLELES, collapse = ",")) %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(MARKERS, ALT)

  res <- dplyr::left_join(
    res,
    dplyr::left_join(ref.alleles, alt.alleles, by = "MARKERS")
    , by = "MARKERS"
  ) %>%
    dplyr::arrange(MARKERS, INTEGERS)
  return(res)
}#End nuc2integers

#' @title nuc2gt
#' @description to apply the conversion of haplo to gt in parallel
#' @rdname nuc2gt
#' @keywords internal
#' @export
nuc2gt <- function(x, conversion.data = NULL, biallelic = TRUE) {
  if (biallelic) {
    res <- dplyr::select(x, MARKERS, GT_VCF_NUC) %>%
      tidyr::separate(data = ., col = GT_VCF_NUC, into = c("A1", "A2"), sep = "/", remove = FALSE) %>%
      dplyr::left_join(dplyr::rename(conversion.data, A1 = ALLELES), by = c("MARKERS", "A1")) %>%
      dplyr::rename(A1_NUC = INTEGERS) %>%
      dplyr::left_join(dplyr::rename(conversion.data, A2 = ALLELES), by = c("MARKERS", "A2")) %>%
      dplyr::rename(A2_NUC = INTEGERS) %>%
      dplyr::mutate(
        GT_VCF = stringi::stri_join(A1_NUC, A2_NUC, sep = "/"),
        GT_VCF = stringi::stri_replace_na(str = GT_VCF, replacement = "./."),
        A1_NUC = as.character(A1_NUC + 1),
        A2_NUC = as.character(A2_NUC + 1),
        A1_NUC = stringi::stri_replace_na(str = A1_NUC, replacement = "0"),
        A2_NUC = stringi::stri_replace_na(str = A2_NUC, replacement = "0"),
        A1_NUC = stringi::stri_pad_left(str = A1_NUC, pad = "0", width = 3),
        A2_NUC = stringi::stri_pad_left(str = A2_NUC, pad = "0", width = 3),
        GT_BIN = as.numeric(stringi::stri_replace_all_fixed(
          str = GT_VCF, pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
          replacement = c("0", "2", "1", "1", NA), vectorize_all = FALSE))
      ) %>%
      tidyr::unite(data = ., col = GT, A1_NUC, A2_NUC, sep = "") %>%
      dplyr::select(-A1, -A2)
  } else {
    res <- dplyr::select(x, MARKERS, GT_VCF_NUC) %>%
      tidyr::separate(data = ., col = GT_VCF_NUC, into = c("A1", "A2"), sep = "/", remove = FALSE) %>%
      dplyr::left_join(dplyr::rename(conversion.data, A1 = ALLELES), by = c("MARKERS", "A1")) %>%
      dplyr::rename(A1_NUC = INTEGERS) %>%
      dplyr::left_join(dplyr::rename(conversion.data, A2 = ALLELES), by = c("MARKERS", "A2")) %>%
      dplyr::rename(A2_NUC = INTEGERS) %>%
      dplyr::mutate(
        GT_VCF = stringi::stri_join(A1_NUC, A2_NUC, sep = "/"),
        A1_NUC = as.character(A1_NUC + 1),
        A2_NUC = as.character(A2_NUC + 1),
        A1_NUC = stringi::stri_replace_na(str = A1_NUC, replacement = "0"),
        A2_NUC = stringi::stri_replace_na(str = A2_NUC, replacement = "0"),
        A1_NUC = stringi::stri_pad_left(str = A1_NUC, pad = "0", width = 3),
        A2_NUC = stringi::stri_pad_left(str = A2_NUC, pad = "0", width = 3)
      ) %>%
      tidyr::unite(data = ., col = GT, A1_NUC, A2_NUC, sep = "") %>%
      dplyr::select(-A1, -A2) %>%
      dplyr::mutate(
        GT_VCF = stringi::stri_replace_na(str = GT_VCF, replacement = "./.")
        # GT_BIN = as.numeric(rep(NA, n()))
      )
  }
  return(res)
}#End nuc2gt


#' @title gt_haplo2gt_vcf_nuc
#' @description to apply the conversion of GT haplotype coding to gt_vcf_nuc in parallel
#' @rdname gt_haplo2gt_vcf_nuc
#' @keywords internal
#' @export

gt_haplo2gt_vcf_nuc <- function(x) {
  res <- x %>%
    dplyr::mutate(
      GT_VCF_NUC = dplyr::if_else(
        stringi::stri_detect_fixed(
          str = GT_HAPLO, pattern = "/"),
        GT_HAPLO,
        stringi::stri_join(GT_HAPLO, GT_HAPLO, sep = "/")),
      GT_VCF_NUC = stringi::stri_replace_na(str = GT_VCF_NUC, replacement = "./.")
    ) %>%
    dplyr::select(-GT_HAPLO)
  return(res)
}#End gt_haplo2gt_vcf_nuc


#' @title split_vcf_id
#' @description split VCF ID in parallel
#' @rdname split_vcf_id
#' @keywords internal
#' @export

split_vcf_id <- function(x) {
  res <- dplyr::rename(x, ID = LOCUS) %>%
    tidyr::separate(data = ., col = ID, into = c("LOCUS", "COL"),
                    sep = "_", extra = "drop", remove = FALSE) %>%
    dplyr::mutate_at(.tbl = ., .vars = c("CHROM", "POS", "LOCUS"), .funs = as.character) %>%
    tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__", remove = FALSE)
  return(res)
}#End split_vcf_id
