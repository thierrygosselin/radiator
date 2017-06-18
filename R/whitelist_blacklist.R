# Create a whitelist from an object with LOCUS column
# whitelist_loci: no visible binding for global variable ‘LOCUS’
# whitelist_loci_snp: no visible binding for global variable ‘LOCUS’
# whitelist_loci_snp: no visible binding for global variable ‘POS’
# whitelist_loci_vcf: no visible binding for global variable ‘POS’
# whitelist_loci_vcf: no visible binding for global variable ‘CHROM’
if(getRversion() >= "2.15.1")  utils::globalVariables(c("LOCUS","POS", "CHROM"))



#' @title Whitelist loci
#' @description This function creates a whitelist of loci, used after applying a filter
#' to keep track of the loci kept by a filter.
#' @param data A tidy vcf or sumstats prep file (using ".tsv") or object in
#' your Environment.
#' @param filename The name of the file written in the directory.
#' @param col.header TRUE and the loci will have a column header 'LOCUS'.
#' @rdname whitelist_loci
#' @export
#' @import dplyr
#' @import readr
#' @importFrom utils write.table
#' @seealso \link{tidy_genomic_data} and  \link{summary_stats_vcf_tidy}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

whitelist_loci <- function(data, filename, col.header) {
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
  } else {
    data <- data
  }
  
  whitelist <- data %>% select(LOCUS) %>% distinct(LOCUS) %>% arrange(LOCUS)
  
    utils::write.table(whitelist, filename, sep = "\t", row.names = F, col.names = col.header,
              quote = F)
invisible(
  cat(
    sprintf(
"Whitelist filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
  whitelist
}




#' @title Whitelist loci and snp
#' @description This function creates a whitelist of loci and snp,
#' useful in the populations module of STACKS.
#' @param data A tidy vcf or sumstats prep file (using ".tsv") or object in
#' your Environment.
#' @param filename The name of the file written in the directory.
#' @param col.header TRUE and the loci will have a column header 'LOCUS'.
#' @rdname whitelist_loci_snp
#' @export
#' @import dplyr
#' @import readr
#' @importFrom utils write.table
#' @seealso \link{tidy_genomic_data} and  \link{summary_stats_vcf_tidy}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

whitelist_loci_snp <- function(data, filename, col.header) {
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
  } else {
    data <- data
  }
  
  whitelist <- data %>%
    select(LOCUS, POS) %>%
    distinct(POS) %>% 
    arrange(LOCUS, POS)
  
  utils::write.table(whitelist, filename, sep = "\t", row.names = F, col.names = col.header,
              quote = F)
  invisible(
  cat(
    sprintf(
"Whitelist filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
  whitelist
}




#' @title Whitelist loci for VCF tools
#' @description This function creates a whitelist of loci for VCF tools.
#' With 2 columns (CHROM and POS).
#' @param data A tidy vcf or sumstats prep file (using ".tsv") or object in
#' your Environment.
#' @param filename The name of the file written in the directory.
#' @rdname whitelist_snp_vcf
#' @export
#' @import dplyr
#' @import readr
#' @importFrom utils write.table
#' @seealso \link{tidy_genomic_data} and  \link{summary_stats_vcf_tidy}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

whitelist_snp_vcf <- function(data, filename) {
  
  if (is.vector(data) == "TRUE") {
    data <- read_tsv(data, col_names = T)
  } else {
    data <- data
  }
  
  whitelist <- data %>% 
    distinct(POS) %>% 
    mutate(CHROM = rep("1", n())) %>% 
    arrange(CHROM, POS) %>% 
    group_by(CHROM) %>%
    select(CHROM, POS)
  
  utils::write.table(whitelist, filename, sep = "\t", row.names = F, col.names = T,
              quote = F)
  invisible(
  cat(
    sprintf(
"VCF whitelist filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
  whitelist
}



#' @title Blacklist loci
#' @description This function creates a blacklist of loci, used after applying a filter
#' to keep track of the loci removed by a filter.
#' @param before.filter.data A tidy vcf or sumstats prep file (using ".tsv") 
#' or object in your Environment. Before the filter you want to test.
#' @param after.filter.data A tidy vcf or sumstats prep file (using ".tsv") 
#' or object in your Environment, after the filter you want to test. 
#' @param filename The name of the file written in the directory.
#' @param col.header TRUE and the loci will have a column header 'LOCUS'.
#' @rdname blacklist_loci
#' @export
#' @import dplyr
#' @import readr
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

blacklist_loci <- function(before.filter.data, after.filter.data, filename, 
                           col.header) {
  
  if (is.vector(before.filter.data) == "TRUE") {
    before.filter.data <- read_tsv(before.filter.data, col_names = T)
  } else {
    before.filter.data <- before.filter.data
  }
  
  if (is.vector(after.filter.data) == "TRUE") {
    after.filter.data <- read_tsv(after.filter.data, col_names = T)
  } else {
    after.filter.data <- after.filter.data
  }
  
  blacklist <- before.filter.data %>%
    distinct(LOCUS) %>% 
    anti_join(after.filter.data %>% 
                distinct(LOCUS),
              by = "LOCUS") %>%
    arrange(LOCUS)
  
  write_tsv(x = blacklist, path = filename, col_names = col.header)
  invisible(
  cat(
    sprintf(
"Blacklist filename:
%s\n
Written in the directory:
%s",
filename, getwd()
)))
  
  blacklist
}


