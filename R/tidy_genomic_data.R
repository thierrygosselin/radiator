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
#' \code{"DP", "AD", "GL", "PL", "GQ", "HQ", "GOF", "NR", "NV"}.
#' If you need different one, submit a request.
#' With \code{string}, explicitely ask for a particular FORMAT field you want
#' to keep (except the GT field that is always imported).
#' e.g. \code{vcf.metadata = "PL"} or \code{vcf.metadata = c("DP", "GL")}.
#' Default: \code{vcf.metadata = FALSE}.

#' @param whitelist.markers (optional) A whitelist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by chromosome and/or locus and/or by snp.
#' The whitelist is an object in your
#' global environment or a file in the working directory (e.g. "whitelist.txt").
#' Note that \emph{de novo} CHROM column with 'un' need to be changed to 1.
#' In the VCF, the column ID is the LOCUS identification (VCF generated from
#' stacks have the SNP position on the read embedded in the ID,
#' so the ID = no longer represent the LOCUS). Marker names are cleaned of
#' separators that interfere with some packages or codes:
#' \code{/}, \code{:}, \code{-} and \code{.} are changed to an underscore
#' \code{_}.
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
#' INDIVIDUALS (character) columns header.
#' The blacklist is an object in your global environment or
#' a file in the working directory (e.g. "blacklist.genotype.txt").
#' For de novo VCF, CHROM column
#' with 'un' need to be changed to 1.
#' Marker names are cleaned of
#' separators that interfere with some packages or codes:
#' \code{/}, \code{:}, \code{-} and \code{.} are changed to an underscore
#' \code{_}.
#' Ids are also cleaned of separators that interfere with some packages or codes:
#' \code{_} and \code{:} are changed to a dash \code{-}.
#' Default: \code{blacklist.genotype = NULL} for no blacklist of
#' genotypes to erase.

#' @param snp.ld (optional) \strong{For data with locus and SNP info, like VCF and DArT file}.
#' SNP short distance linkage disequilibrium pruning. With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 5 options:
#' \enumerate{
#' \item \code{snp.ld = "random"} for a random selection of 1 SNP on the read,
#' \item \code{snp.ld = "first"} for the first one on the read...,
#' \item \code{snp.ld = "last"} for the last SNP on the read and
#' \item \code{snp.ld = "middle"} for locus with > 2 SNPs/read the option to select at random
#' one SNP between the first and the last SNP on the read. If the locus as <= 2
#' SNPs on the read, the first one is selected. Note that for that last option,
#' the numbers are reported.
#' \item \code{snp.ld = "maf"} will select the SNP on the locus with the maximum global
#' Minor Allele Frequency (MAF).
#' }
#' For long distance linkage disequilibrium pruning, see details below.
#' Default: \code{snp.ld = NULL}, for no pruning.

#' @param common.markers (optional) Logical. Default: \code{common.markers = TRUE},
#' will only keep markers in common (genotyped) between all the populations.

#' @param maf.thresholds (optional) Use for minor allele frequency (maf) or count (mac)
#' filter. Default: \code{maf.thresholds = NULL}, for no filtering.
#'
#' String with 5 values.
#' For example using SNP and frequency approach: \code{maf.thresholds = c("SNP", 0.001, "OR", 0.001, 1)}.
#' For example using locus and count approach: \code{maf.thresholds = c("locus", 3, "OR", 6, 1)}.
#' \enumerate{
#' \item maf approach (character: "SNP"/"locus"):
#' MAF filtering is conducted by SNPs or locus.
#' \code{"SNP"}: will consider all the SNPs on the same locus/read as independent
#' and will be filtered independently of their locus id.
#' \code{"locus"}: looks at the minimum MAF found on the
#' read/locus. Using this option will discard all the markers/snp on
#' that read based on the thresholds chosen. For the locus approach to work, your dataset
#' requires SNP and Locus info (e.g. from a VCF file).
#'
#' \item local threshold (double or integer): For a frequency threshold use
#' a double (e.g. 0.05). For a count threshold, use an integer
#' (e.g. 3 for the number of alternate allele required in a population). This
#' threshold is applied by population.
#' Not sure about the threshold to use, choose the interactive mode argument.
#'
#' \item operator (character: "OR" / "AND"):
#' To consider both the local \code{"AND"} the global thresholds, use \code{"AND"}.
#' To consider the local \code{"OR"} the global thresholds, use \code{"OR"}.
#'
#' \item global threshold (double or integer): For a frequency threshold
#' use a double (e.g. 0.02). For a count threshold, use an integer
#' (e.g. 6 for the number of alternate allele required). This threshold is
#' applied at the dataset level (no population).
#' Not sure about the threshold to use, choose the interactive mode argument.
#'
#'\item maf pop num threshold (integer)
#'The number of pop required to pass the thresholds
#' to keep the marker. Usually, I always use \code{1}, for 1 pop is required to
#' pass thresholds and keep the marker.
#' }

#' @param max.marker (integer, optional) For large PLINK files,
#' useful to subsample marker number. e.g. if the data set
#' contains 200 000 markers and \code{max.marker = 10000}, 10000 markers are
#' subsampled randomly from the 200000 markers. If you need specific markers,
#' use \code{whitelist.markers} argument.
#' Default: \code{max.marker = NULL}.

#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is an object in your
#' global environment or a file in the working directory
#' (e.g. "blacklist.txt"). \code{_} and \code{:} in individual's names
#' are changed to a dash \code{-}.
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
#' \code{pop.labels = c("NEW", "NEW", "ALB")}.
#' To rename \code{"QUE"} to \code{"TAS"}:
#' \code{pop.labels = c("TAS", "ONT", "ALB")}.
#' Default: \code{pop.labels = NULL}. If you find this too complicated,
#' there is also the \code{strata} argument that can do the same thing,
#' see below.
#' White spaces in population names are replaced by underscore.


#' @param strata (optional)
#' The strata file is a tab delimited file with a minimum of 2 columns headers:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' If a \code{strata} file is specified with all file formats that don't
#' require it, the strata argument will have precedence on the population
#' groupings used internally in those file formats. For file formats without
#' population/strata groupings (e.g. vcf, haplotype files) if no strata file is
#' provided, 1 pop/strata grouping will automatically be created.
#' For vcf and haplotypes file, the strata can also be used as a whitelist of id.
#' Samples not in the strata file will be discarded from the data set.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.
#' If you have already run
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data,
#' the strata file is similar to a stacks \emph{population map file},
#' make sure you
#' have the required column names (\code{INDIVIDUALS} and \code{STRATA}).
#' The strata column is cleaned of a white spaces that interfere with some
#' packages or codes: space is changed to an underscore \code{_}.
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
#' execution during import.
#' Default: \code{parallel::detectCores() - 1}.


#' @param verbose (optional, logical) When \code{verbose = TRUE}
#' the function is a little more chatty during execution.
#' Default: \code{verbose = TRUE}.

#' @param ... (optional) To pass further argument for fine-tuning the function.


#' @details
#' \strong{Long distance SNP linkage disequilibrium pruning}
#' If you have markers position on a genome or a linkage map,
#' you can go further in removing linked markers by using
#' \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate} or
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
#' \code{MARKERS or LOCUS} and \code{GT}. The rest are considered metata info.
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
#' written in the working directory with file extension \code{.rad}.
#' The file is written with the
#' \href{https://github.com/fstpackage/fst}{Lightning Fast Serialization of Data Frames for R} package.
#' To read the file back in R use \code{\link[fst]{read.fst}}.

#' @export
#' @rdname tidy_genomic_data
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs
# @importFrom adegenet genind2df
# @importFrom strataG as.data.frame
#' @importFrom stringi stri_join stri_replace_all_fixed stri_extract_all_fixed stri_replace_all_regex stri_sub stri_pad_left stri_count_fixed stri_replace_na
#' @importFrom stats var median quantile
#' @importFrom purrr map flatten keep discard
#' @importFrom data.table fread as.data.table
#' @importFrom tidyr spread gather unite separate
#' @importFrom utils count.fields
#' @importFrom readr write_tsv read_tsv
#' @importFrom tibble as_data_frame has_name
#' @importFrom parallel detectCores
# @importFrom pegas VCFloci read.vcf

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
  max.marker = NULL,
  blacklist.id = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  strata = NULL,
  pop.select = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
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

  # Check pop.levels, pop.labels and pop.select---------------------------------
  check.levels <- check_pop_levels(pop.levels = pop.levels, pop.labels = pop.labels, pop.select = pop.select)
  pop.levels <- check.levels$pop.levels
  pop.labels <- check.levels$pop.labels
  pop.select <- check.levels$pop.select


  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)
  want <- c("keep.allele.names")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }

  radiator.dots <- dotslist[names(dotslist) %in% want]
  keep.allele.names <- radiator.dots[["keep.allele.names"]]

  if (is.null(keep.allele.names)) keep.allele.names <- FALSE


  # File type detection----------------------------------------------------------
  skip.tidy.wide <- FALSE # initiate for data frame below
  data.type <- radiator::detect_genomic_format(data)

  # Import whitelist of markers-------------------------------------------------
  if (!is.null(whitelist.markers)) {# with Whitelist of markers
    if (is.vector(whitelist.markers)) {
      whitelist.markers <- suppressMessages(readr::read_tsv(whitelist.markers, col_names = TRUE) %>%
                                              dplyr::mutate_all(.tbl = ., .funs = as.character))
    }
    columns.names.whitelist <- colnames(whitelist.markers)

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

    whitelist.markers <- dplyr::mutate_all(
      .tbl = whitelist.markers, .funs = clean_markers_names)
  }

  # Import blacklist id --------------------------------------------------------
  if (!is.null(blacklist.id)) {# With blacklist of ID
    if (is.vector(blacklist.id)) {
      suppressMessages(blacklist.id <- readr::read_tsv(
        blacklist.id,
        col_names = TRUE,
        col_types = readr::cols(.default = readr::col_character())))
    } else {
      if (!tibble::has_name(blacklist.id, "INDIVIDUALS")) {
        stop("Blacklist of individuals should have 1 column named: INDIVIDUALS")
      }
      blacklist.id <- dplyr::mutate_all(.tbl = blacklist.id, .funs = as.character)
    }
    # not for plink file, where it's done after.
    # see plink section to understand
    blacklist.id$INDIVIDUALS <- radiator::clean_ind_names(blacklist.id$INDIVIDUALS)

    # remove potential duplicate id
    dup <- dplyr::distinct(.data = blacklist.id, INDIVIDUALS)
    blacklist.id.dup <- nrow(blacklist.id) - nrow(dup)
    if (blacklist.id.dup >1) {
      message("Duplicate id's in blacklist: ", blacklist.id.dup)
      blacklist.id <- dup
    }
    dup <- blacklist.id.dup <- NULL
    message("Number of individuals in blacklist: ", nrow(blacklist.id))
  }

  # Strata and pop levels ------------------------------------------------------
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      strata.df <- readr::read_tsv(
        file = strata, col_names = TRUE,
        col_types = readr::cols(.default = readr::col_character()))
    } else {
      strata.df <- strata
    }

    colnames(strata.df) <- stringi::stri_replace_all_fixed(
      str = colnames(strata.df),
      pattern = "STRATA",
      replacement = "POP_ID",
      vectorize_all = FALSE)


    # Remove potential whitespace in pop_id
    strata.df$POP_ID <- radiator::clean_pop_names(strata.df$POP_ID)
    colnames.strata <- colnames(strata.df)

    # clean ids
    strata.df$INDIVIDUALS <- radiator::clean_ind_names(strata.df$INDIVIDUALS)

    # filtering the strata if blacklist id available
    if (!is.null(blacklist.id)) {
      strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
    }
  } else {
    strata.df <- NULL
  }

  # Import VCF------------------------------------------------------------------
  if (data.type == "vcf.file") { # VCF
    if (verbose) message("Importing and tidying the VCF...")
    input <- radiator::tidy_vcf(
      data = data,
      strata = strata.df,
      vcf.metadata = vcf.metadata,
      parallel.core = parallel.core,
      verbose = verbose,
      whitelist.markers = whitelist.markers,
      blacklist.id = blacklist.id,
      pop.select = pop.select,
      pop.levels = pop.levels,
      pop.labels = pop.labels
    )
    biallelic <- radiator::detect_biallelic_markers(input)
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
      tibble::as_data_frame(.) %>%
      dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                       .funs = clean_ind_names) %>%
      dplyr::mutate_at(.tbl = ., .vars = "POP_ID",
                       .funs = clean_pop_names)

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
      strata.df$INDIVIDUALS <- radiator::clean_ind_names(strata.df$INDIVIDUALS)
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
      tibble::as_data_frame(.) %>%
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
    input <- tidyr::gather(data = input, key = INDIVIDUALS_ALLELES, value = GT, -LOCUS)

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
      dplyr::group_by(LOCUS, INDIVIDUALS) %>%
      tidyr::spread(data = ., key = ALLELES, value = GT) %>%
      dplyr::ungroup(.) %>%
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
      if (is.factor(input$POP_ID)) input$POP_ID <- droplevels(input$POP_ID)
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
    input.temp <- radiator::change_alleles(data = input,
                                           verbose = verbose)
    input <- input.temp$input
    biallelic <- input.temp$biallelic
  } # End import PLINK

  # Import stacks haplotypes----------------------------------------------------
  if (data.type == "haplo.file") { # Haplotype file
    if (verbose) message("Importing STACKS haplotype file")

    strata.df <- strata_haplo(strata = strata.df,
                              data = data, blacklist.id = blacklist.id)

    # import header row
    want <- tibble::data_frame(
      INFO = "CATALOG",
      COL_TYPE = "c") %>%
      dplyr::bind_rows(
        dplyr::select(strata.df, INFO = INDIVIDUALS) %>%
          dplyr::mutate(
            COL_TYPE = rep("c", n()),
            INFO = clean_ind_names(INFO)
          ))

    haplo.col.type <- readr::read_tsv(
      file = data,
      n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character())) %>%
      tidyr::gather(data = .,key = DELETE, value = INFO) %>%
      dplyr::mutate(INFO = clean_ind_names(INFO)) %>%
      dplyr::select(-DELETE) %>%
      dplyr::mutate(INFO = clean_ind_names(INFO)) %>%
      dplyr::left_join(want, by = "INFO") %>%
      dplyr::mutate(COL_TYPE = stringi::stri_replace_na(str = COL_TYPE, replacement = "_")) %>%
      dplyr::select(COL_TYPE)

    haplo.col.type[1,1] <- "c"

    haplo.col.type <- purrr::flatten_chr(haplo.col.type) %>% stringi::stri_join(collapse = "")

    # readr now faster/easier than fread...
    input <- readr::read_tsv(
      file = data, col_names = TRUE, na = "-",
      col_types = haplo.col.type)

    colnames(input) <- stringi::stri_replace_all_fixed(
      str = colnames(input),
      pattern = c("# Catalog ID", "Catalog ID", "# Catalog Locus ID"),
      replacement = c("LOCUS", "LOCUS", "LOCUS"), vectorize_all = FALSE)

    if (tibble::has_name(input, "Seg Dist")) {
      input <- dplyr::select(.data = input, -`Seg Dist`)
    }

    n.catalog.locus <- dplyr::n_distinct(input$LOCUS)
    n.individuals <- ncol(input) - 1

    message("\nNumber of loci in catalog: ", n.catalog.locus)
    message("Number of individuals: ", n.individuals)
    input <- tidyr::gather(
      data = input,
      key = "INDIVIDUALS",
      value = "GT_VCF_NUC", # previously using "GT_HAPLO"
      -LOCUS)

    input$INDIVIDUALS <- radiator::clean_ind_names(input$INDIVIDUALS)

    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      if (verbose) message("Filtering with whitelist of markers")
      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }

    # remove consensus markers
    if (verbose) message("\nScanning for consensus markers...")
    consensus.markers <- dplyr::filter(input, GT_VCF_NUC == "consensus") %>%
      dplyr::distinct(LOCUS)

    if (length(consensus.markers$LOCUS) > 0) {
      input <- suppressWarnings(dplyr::anti_join(input, consensus.markers, by = "LOCUS"))
      readr::write_tsv(consensus.markers, "radiator.tidy.genomic.data.consensus.markers.tsv")
    }
    if (verbose) message("    number of consensus markers removed: ", dplyr::n_distinct(consensus.markers$LOCUS))
    consensus.markers <- NULL

    # population levels and strata
    strata.df$INDIVIDUALS <- radiator::clean_ind_names(strata.df$INDIVIDUALS)

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
      dplyr::mutate(POLYMORPHISM = stringi::stri_count_fixed(GT_VCF_NUC, "/"))

    blacklist.paralogs <- input %>%
      dplyr::filter(POLYMORPHISM > 1) %>%
      dplyr::select(LOCUS, INDIVIDUALS)

    if (verbose) message("    number of genotypes with more than 2 alleles: ", length(blacklist.paralogs$LOCUS))
    if (length(blacklist.paralogs$LOCUS) > 0) {
      input <- input %>%
        dplyr::mutate(GT_VCF_NUC = replace(GT_VCF_NUC, which(POLYMORPHISM > 1), NA)) %>%
        dplyr::select(-POLYMORPHISM)

      readr::write_tsv(blacklist.paralogs, "blacklist.genotypes.paralogs.tsv")
    }
    blacklist.paralogs <- NULL

    if (verbose) message("Calculating REF/ALT alleles...")
    # Prep for REF/ALT alleles and new genotype coding
    # part below could be parallelized if necessary, test with larger dataset for bottleneck...
    input <- input %>%
      dplyr::mutate(
        GT_VCF_NUC = dplyr::if_else(
          POLYMORPHISM == 0,
          stringi::stri_join(GT_VCF_NUC, "/", GT_VCF_NUC), GT_VCF_NUC,
          missing = "./."),
        GT_VCF_NUC = dplyr::if_else(stringi::stri_detect_fixed(GT_VCF_NUC, "N"),
                                    "./.", GT_VCF_NUC)
      ) %>%
      dplyr::select(-POLYMORPHISM)

    input.temp <- radiator::change_alleles(
      data = input,
      biallelic = FALSE,
      parallel.core = parallel.core,
      verbose = verbose)
    input <- input.temp$input
    input.temp <- NULL
    biallelic <- FALSE
    input <- dplyr::rename(input, LOCUS = MARKERS)
  } # End import haplotypes file

  # Import genepop--------------------------------------------------------------
  if (data.type == "genepop.file") {
    if (verbose) message("Tidying the genepop file ...")
    input <- radiator::tidy_genepop(data = data, tidy = TRUE)
    skip.tidy.wide <- TRUE
  }

  # Import DArT ----------------------------------------------------------------
  if (data.type == "dart") {
    if (verbose) message("Tidying DArT data...")
    input <- radiator::tidy_dart(
      data = data,
      strata = strata,
      verbose = FALSE,
      parallel.core = parallel.core)
    skip.tidy.wide <- TRUE

    if (tibble::has_name(strata.df, "NEW_ID")) {
      strata.df <- strata.df %>%
        dplyr::select(-INDIVIDUALS) %>%
        dplyr::rename(INDIVIDUALS = NEW_ID)
    }
  }

  # Import fst.file ------------------------------------------------------------
  if (data.type == "fst.file") {
    if (verbose) message("Importing the fst.file as a data frame...")
    input <- read_rad(data = data)
    skip.tidy.wide <- TRUE
  }

  # Import GENIND--------------------------------------------------------------
  if (data.type == "genind") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the genind object ...")
    input <- radiator::tidy_genind(data = data, keep.allele.names = keep.allele.names)
    data <- NULL
    # remove unwanted sep in id and pop.id names
    input <- input %>%
      dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                       .funs = clean_ind_names) %>%
      dplyr::mutate_at(.tbl = ., .vars = "POP_ID",
                       .funs = clean_pop_names)
    skip.tidy.wide <- TRUE
  } # End tidy genind

  # Import GENLIGHT ------------------------------------------------------------
  if (data.type == "genlight") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the genlight object ...")
    input <- radiator::tidy_genlight(data = data)
    data <- NULL
    # remove unwanted sep in id and pop.id names
    input <- input %>%
      dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                       .funs = clean_ind_names) %>%
      dplyr::mutate_at(.tbl = ., .vars = "POP_ID",
                       .funs = clean_pop_names)
    biallelic <- TRUE
    skip.tidy.wide <- TRUE
  } # End tidy genlight

  # Import STRATAG gtypes ------------------------------------------------------
  if (data.type == "gtypes") { # DATA FRAME OF GENOTYPES
    if (verbose) message("Tidying the gtypes object ...")
    input <- tidy_gtypes(data) %>%
      dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                       .funs = clean_ind_names) %>%
      dplyr::mutate_at(.tbl = ., .vars = "POP_ID",
                       .funs = clean_pop_names)
    data <- NULL
    skip.tidy.wide <- TRUE
  } # End tidy gtypes

  # Import DF-------------------------------------------------------------------
  if (data.type == "tbl_df" || skip.tidy.wide) { # DATA FRAME OF GENOTYPES
    if (verbose) message("Importing the data frame ...")
    if (!skip.tidy.wide) {
      input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
      data <- NULL
    }

    # For long tidy format, switch LOCUS to MARKERS column name, if found MARKERS not found
    if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
      input <- dplyr::rename(.data = input, MARKERS = LOCUS)
    }

    # Change individuals names containing special character
    input$INDIVIDUALS <- radiator::clean_ind_names(input$INDIVIDUALS)

    # initiate ref check
    if (tibble::has_name(input, "REF")) {
      check.ref <- FALSE
    } else {
      check.ref <- TRUE
    }

    # Filter with whitelist of markers
    input$MARKERS <- clean_markers_names(input$MARKERS)

    if (!is.null(whitelist.markers)) {
      if (verbose) message("Filtering with whitelist of markers")

      if (tibble::has_name(input, "MARKERS") && !tibble::has_name(input, "LOCUS")) {
        if (ncol(whitelist.markers) == 1 && colnames(whitelist.markers) == "LOCUS") {
          colnames(whitelist.markers) <- "MARKERS"
          columns.names.whitelist <- "MARKERS"
        }
      }

      input <- suppressWarnings(dplyr::semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }

    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      if (verbose) message("Filtering with blacklist of individuals")
      input <- suppressWarnings(dplyr::anti_join(input, blacklist.id, by = "INDIVIDUALS"))
      check.ref <- TRUE
    }


    # population levels and strata
    if (!is.null(strata)) {
      strata.df$INDIVIDUALS <- radiator::clean_ind_names(strata.df$INDIVIDUALS)

      input <- suppressWarnings(input %>% dplyr::select(-dplyr::one_of(c("POP_ID"))))

      # Check for matching ids before merging info
      check.id <- unique(sort(unique(strata.df$INDIVIDUALS)) %in%
                           sort(unique(input$INDIVIDUALS)))
      if (FALSE %in% check.id) {
        stop("Non-matching INDIVIDUALS between data and strata")
      } else {
        input <- dplyr::left_join(input, strata.df, by = "INDIVIDUALS")
      }
    }

    # Change potential problematic POP_ID space
    input$POP_ID <- radiator::clean_pop_names(input$POP_ID)

    # Check with strata and pop.levels/pop.labels
    if (!is.null(pop.levels)) {
      if (length(levels(factor(input$POP_ID))) != length(pop.levels)) {
        stop("The number of groups in your POP_ID column file must match the number of groups in pop.levels")
      }
    }

    # using pop.levels and pop.labels info if present
    input <- change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)
    if (!identical(pop.levels, pop.labels)) check.ref <- TRUE

    # Pop select
    if (!is.null(pop.select)) {
      if (verbose) message(stringi::stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% dplyr::filter(POP_ID %in% pop.select))
      check.ref <- TRUE
    }

    if (check.ref) {
      input.temp <- radiator::change_alleles(data = input)
      input <- input.temp$input
      biallelic <- input.temp$biallelic
      input.temp <- NULL
    } else {
      biallelic <- radiator::detect_biallelic_markers(data = input)
    }

  } # End import data frame of genotypes

  # END IMPORT DATA-------------------------------------------------------------

  # Arrange the id and create a strata after pop select ------------------------
  input$INDIVIDUALS <- radiator::clean_ind_names(input$INDIVIDUALS)

  strata.df <- dplyr::ungroup(input) %>%
    dplyr::distinct(POP_ID, INDIVIDUALS)

  # Blacklist genotypes --------------------------------------------------------
  if (is.null(blacklist.genotype)) { # no Whitelist
    if (verbose) message("Erasing genotype: no")
  } else {
    if (verbose) message("Erasing genotype: yes")
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS")
    if (is.vector(blacklist.genotype)) {
      suppressWarnings(suppressMessages(
        blacklist.genotype <- readr::read_tsv(blacklist.genotype, col_names = TRUE)))
    }
    suppressWarnings(suppressMessages(
      blacklist.genotype <- blacklist.genotype %>%
        dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
                         .funs = clean_ind_names) %>%
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
      blacklist.genotype <- suppressWarnings(
        dplyr::semi_join(whitelist.markers.ind, blacklist.genotype,
                         by = columns.names.blacklist.genotype))
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }

    # Update column names
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)

    blacklisted.gen.number <- nrow(blacklist.genotype)
    if (blacklisted.gen.number > 0) {
      message("    Number of genotype(s) to erase: ", blacklisted.gen.number)
      input.erase <- dplyr::semi_join(
        input, blacklist.genotype, by = columns.names.blacklist.genotype) %>%
        dplyr::mutate(GT = rep("000000", n()))
      input <- dplyr::anti_join(
        input, blacklist.genotype, by = columns.names.blacklist.genotype)
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
    } else {
      message("There are no genotype left in the blacklist: input file left intact")
    }

    # required because REF/ALT might change after deleting genotypes...
    input <- radiator::change_alleles(data = input)$input
  } # End erase genotypes

  # dump unused object
  blacklist.id <- whitelist.markers <- whitelist.markers.ind <- NULL
  want <- blacklist.genotype <- NULL
  # SNP LD  --------------------------------------------------------------------
  if (!is.null(snp.ld)) {
    input <- radiator::snp_ld(data = input, snp.ld = snp.ld)
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
  input$MARKERS <- radiator::clean_markers_names(input$MARKERS)

  # Markers in common between all populations (optional) -----------------------
  if (common.markers) { # keep only markers present in all pop
    input <- radiator::keep_common_markers(input, verbose = TRUE)$input
  } # End common markers

  # Removing monomorphic markers------------------------------------------------
  if (monomorphic.out) {
    if (tibble::has_name(input, "POLYMORPHIC")) {
      message("Scanning for monomorphic markers...")
      if (verbose) message("    Number of markers before = ", dplyr::n_distinct(input$MARKERS))
      if (tibble::has_name(input, "POS")) {
        mono.markers <- dplyr::filter(input, !POLYMORPHIC) %>%
          dplyr::distinct(MARKERS, CHROM, LOCUS, POS)
      } else {
        mono.markers <- dplyr::filter(input, !POLYMORPHIC) %>%
          dplyr::distinct(MARKERS)
      }

      if (nrow(mono.markers) > 0) {
        message("    Number of monomorphic markers removed = ", nrow(mono.markers))
        input <- dplyr::filter(input, !MARKERS %in% mono.markers$MARKERS)
        readr::write_tsv(mono.markers, "blacklist.monomorphic.markers.tsv")
        message("    Number of markers after = ", dplyr::n_distinct(input$MARKERS))
      }
      input <- dplyr::select(input, -POLYMORPHIC)
    } else {
      mono.out <- radiator::discard_monomorphic_markers(input, verbose = TRUE)
      mono.markers <- mono.out$blacklist.monomorphic.markers
      if (nrow(mono.markers) > 0) {
        input <- mono.out$input
        readr::write_tsv(mono.markers, "blacklist.monomorphic.markers.tsv")
      }
      mono.out <- NULL
    }
    mono.markers <- NULL
  } # End monomorphic out

  # Minor Allele Frequency filter ----------------------------------------------
  if (!is.null(maf.thresholds)) { # with MAF
  input <- radiator::filter_maf(
    data = input,
    interactive.filter = FALSE,
    maf.thresholds = maf.thresholds,
    parallel.core = parallel.core,
    verbose = FALSE)$tidy.filtered.maf
  } # End of MAF filters


  # Write to working directory -------------------------------------------------
  if (!is.null(filename)) {
    tidy.name <- stringi::stri_join(filename, ".rad")
    message("\nWriting tidy data set:\n", tidy.name)
    write_rad(data = input, path = tidy.name)
  }

  # Results --------------------------------------------------------------------
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
  }
  if (verbose) {
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
  }
  message("\nTidy genomic data:")
  if (common.markers) {
    message("    Number of common markers: ", n.markers)
  } else {
    message("    Number of markers: ", n.markers)
  }
  message("    Number of chromosome/contig/scaffold: ", n.chromosome)
  message("    Number of individuals: ", n.individuals)
  message("    Number of populations: ", n.pop)
  if (verbose) timing <- proc.time() - timing
  if (verbose) message("Computation time: ", round(timing[[3]]), " sec")
  if (verbose) {
    cat("################ radiator::tidy_genomic_data completed ################\n")
  }
  res <- input
  options(width = opt.change)
  return(res)
} # tidy genomic data


# Internal nested Function -----------------------------------------------------

#' @title strata_haplo
#' @description Manage strata
#' @rdname strata_haplo
#' @keywords internal
#' @export
strata_haplo <- function(strata = NULL, data = NULL, blacklist.id = NULL) {

  if (is.null(strata)) {
    message("No strata file provided")
    message("    generating a strata with 1 grouping")
    if (is.null(data)) stop("data required to generate strata")
    strata.df <- readr::read_tsv(
      file = data,
      n_max = 1,
      na = "-",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character())) %>%
      tidyr::gather(data = .,key = DELETE, value = INDIVIDUALS) %>%
      dplyr::mutate(INDIVIDUALS = clean_ind_names(INDIVIDUALS)) %>%
      dplyr::select(-DELETE) %>%
      dplyr::filter(!INDIVIDUALS %in% c("Catalog ID", "Cnt")) %>%
      dplyr::distinct(INDIVIDUALS) %>%
      dplyr::mutate(STRATA = rep("pop1", n()))
  } else {
    if (is.vector(strata)) {
      suppressMessages(
        strata.df <- readr::read_tsv(
          file = strata, col_names = TRUE,
          # col_types = col.types
          col_types = readr::cols(.default = readr::col_character())
        ))
    } else {
      strata.df <- strata
    }
  }

  colnames(strata.df) <- stringi::stri_replace_all_fixed(
    str = colnames(strata.df),
    pattern = "STRATA",
    replacement = "POP_ID",
    vectorize_all = FALSE
  )
  # Remove potential whitespace in pop_id
  strata.df$POP_ID <- clean_pop_names(strata.df$POP_ID)
  colnames.strata <- colnames(strata.df)

  # clean ids
  strata.df$INDIVIDUALS <- clean_ind_names(strata.df$INDIVIDUALS)

  # filtering the strata if blacklist id available
  if (!is.null(blacklist.id)) {
    strata.df <- dplyr::anti_join(x = strata.df, y = blacklist.id, by = "INDIVIDUALS")
  }

  strata.df <- dplyr::distinct(strata.df, POP_ID, INDIVIDUALS)
  return(strata.df)
}#End strata_haplo
