% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayescan.R
\name{run_bayescan}
\alias{run_bayescan}
\title{Run BayeScan}
\usage{
run_bayescan(
  data,
  n = 5000,
  thin = 10,
  nbp = 20,
  pilot = 5000,
  burn = 50000,
  pr_odds,
  subsample = NULL,
  iteration.subsample = 1,
  parallel.core = parallel::detectCores() - 1,
  bayescan.path = "/usr/local/bin/bayescan"
)
}
\arguments{
\item{data}{(character, path) Read carefully because there's 2 ways.
\enumerate{
\item Path to BayeScan input file.
To get this input rapidly: \code{\link[radiator]{write_bayescan}} if
you already have a tidy data, or use \code{\link[radiator]{genomic_converter}}
\item Path to a tidy data file or object.
This type of input is generated with \code{\link[radiator]{genomic_converter}} or
\code{\link[radiator]{tidy_genomic_data}}.
Use this format if you intend to do subsampling with the
arguments described below \code{subsample} and \code{iteration.subsample}.
\item Remember: you can do both the BayeScan file and tidy data with
\code{\link[radiator]{genomic_converter}}.
}}

\item{n}{(integer) Number of outputted iterations. Default: \code{n = 5000}.}

\item{thin}{(integer) Thinning interval size. Default: \code{thin = 10}}

\item{nbp}{(integer) Number of pilot runs. Default: \code{nbp = 20}.}

\item{pilot}{(integer) Length of pilot runs. Default: \code{pilot = 5000}.}

\item{burn}{(integer) Burn-in length. Default: \code{burn = 50000}.}

\item{pr_odds}{(integer) Prior odds for the neutral model. A \code{pr_odds = 10},
indicates that the neutral model is 10 times more likely than the
model with selection. Larger \code{pr_odds} the more conservative is the results.}

\item{subsample}{(Integer or character)
With \code{subsample = 36}, 36 individuals in each populations are chosen
randomly to represent the dataset. With \code{subsample = "min"}, the
minimum number of individual/population found in the data is used automatically.
Default is no subsampling, \code{subsample = NULL}.}

\item{iteration.subsample}{(Integer) The number of iterations to repeat
subsampling.
With \code{subsample = 20} and \code{iteration.subsample = 10},
20 individuals/populations will be randomly chosen 10 times.
Default: \code{iteration.subsample = 1}.}

\item{parallel.core}{(integer) Number of CPU for parallel computations.
Default: \code{parallel.core = parallel::detectCores() - 1}}

\item{bayescan.path}{(character, path) Provide the FULL path to BayeScan program.
Default: \code{bayescan.path = "/usr/local/bin/bayescan"}. See details.}
}
\value{
For specific \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
output files, see \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
documentation, please read the manual.

radiator::run_bayescan outputs without subsampling:

\enumerate{
\item \code{bayescan}: dataframe with results of BayeScan analysis.
\item \code{selection.summary}: dataframe showing the number of markers in the different group of selections and model choice.
\item \code{whitelist.markers.positive.selection}: Whitelist of markers under diversifying selection and common in all iterations.
\item \code{whitelist.markers.neutral.selection}: Whitelist of neutral markers and common in all iterations.
\item \code{whitelist.markers.neutral.positive.selection}: Whitelist of neutral markers and markers under diversifying selection and common in all iterations.
\item \code{blacklist.markers.balancing.selection}: Blacklist of markers under balancing selection and common in all iterations.
\item \code{markers.dictionary}: BayeScan use integer for MARKERS info. In this dataframe, the corresponding values used inside the function.
\item \code{pop.dictionary}: BayeScan use integer for STRATA info. In this dataframe, the corresponding values used inside the function.
\item \code{bayescan.plot}: plot showing markers Fst and model choice.

Additionnally, if multiple SNPs/locus are detected the object will also have:
\item \code{accurate.locus.summary}: dataframe with the number of accurate locus and the selection types.
\item \code{whitelist.accurate.locus}: whitelist of accurate locus.
\item \code{blacklist.not.accurate.locus}: blacklist of not accurate locus.
\item \code{accuracy.snp.number}: dataframe with the number of SNPs per locus and the count of accurate/not accurate locus.
\item \code{accuracy.snp.number.plot}: the plot showing the proportion of accurate/not accurate locus in relation to SNPs per locus.
\item \code{not.accurate.summary}: dataframe summarizing the number of not accurate locus with selection type found on locus.
}

radiator::run_bayescan outputs WITH subsampling:

\enumerate{
\item \code{subsampling.individuals}: dataframe with indivuals subsample id and random seed number.
\item \code{bayescan.all.subsamples}: long dataframe with combined iterations of bayescan results.
\item \code{selection.accuracy}: dataframe with all markers with selection grouping and number of times observed throughout iterations.
\item \code{accurate.markers}: dataframe with markers attributed the same selection grouping in all iterations.
\item \code{accuracy.summary}: dataframe with a summary of accuracy of selection grouping.
\item \code{bayescan.summary}: dataframe with mean value, averaged accross iterations.
\item \code{bayescan.summary.plot}: plot showing markers Fst and model choice.
\item \code{selection.summary}: dataframe showing the number of markers in the different group of selections and model choice.
\item \code{whitelist.markers.positive.selection}: Whitelist of markers under diversifying selection and common in all iterations.
\item \code{whitelist.markers.neutral.selection}: Whitelist of neutral markers and common in all iterations.
\item \code{blacklist.markers.balancing.selection}: Blacklist of markers under balancing selection and common in all iterations.
\item \code{whitelist.markers.neutral.positive.selection}: Whitelist of neutral markers and markers under diversifying selection and common in all iterations.
\item \code{whitelist.markers.without.balancing.positive}:
Whitelist of all original markers with markers under balancing selection and directional selection removed.
The markers that remains are the ones to use in population structure analysis.
}

Other files are present in the folder and subsampling folder.
}
\description{
\strong{Function highlights:}

\enumerate{
\item \strong{integrated and seamless pipeline: } generate \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
files within radiator and run \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan} inside R!
\item \strong{unbalanced sampling sites impact: } measure and verify genome scan accurary in unbalanced sampling design with subsampling related arguments.
\item \strong{SNP linkage: } detect automatically the presence of multiple SNPs on the same locus and
measure/verify accuracy of genome scan within locus.
\item \strong{summary tables and visualization: } the function generate summary tables and plots of genome scan.
\item \strong{whitelists and blacklists} of markers under different selection identity are automatically generated !
}

This function requires a working
\href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan} program installed
on the computer (\href{http://cmpg.unibe.ch/software/BayeScan/download.html}{install instructions}).
For UNIX machines, please install the 64bits version.
}
\details{
\strong{subsampling:}
During subsampling the function will automatically remove monomorphic
markers that are generated by the removal of some individuals. Also, common markers
between all populations are also automatically detected. Consequently, the number of
markers will change throughout the iterations. The nice thing about the function
is that since everything is automated there is less chance of making an error...

\strong{SNPs data set: }
You should not run BayeScan with SNPs data set that have multiple SNPs on the
same LOCUS. Instead, run radiator::genomic_converter using the \code{snp.ld}
argument to keep only one SNP on the locus. Or run the function by first converting
an haplotype vcf or if your RAD dataset was produced by STACKS, use the
\code{batch_x.haplotypes.tsv} file! If the function detect multiple SNPs on
the same locus, accuracy will be measured automatically.

\strong{UNIX install: } I like to transfer the \emph{BayeScan2.1_linux64bits}
(for Linux) or the \emph{BayeScan2.1_macos64bits} (for MACOs) in \code{/usr/local/bin}
and change it's name to \code{bayescan}. Too complicated ? and you've just
downloaded the last BayeScan version, I would try this :
\code{bayescan.path = "/Users/thierry/Downloads/BayeScan2.1/binaries/BayeScan2.1_macos64bits"}

Make sure to give permission: \code{sudo chmod 777 /usr/local/bin/bayescan}
}
\examples{
\dontrun{
# library(radiator)
# get a tidy data frame and a bayescan file with radiator::genomic_converter:
# to run with a vcf haplotype file
data <- radiator::genomic_converter(
    data = "batch_1.haplotypes.vcf",
    strata = "../../02_project_info/strata.stacks.TL.tsv",
    whitelist.markers = "whitelist.filtered.markers.tsv",
    blacklist.id = "blacklist.id.tsv",
    output = "bayescan",
    filename = "bayescan.haplotypes"
    )
# to run BayeScan:
scan.pops <- radiator::run_bayescan(
    data = "bayescan.haplotypes.txt",
    pr_odds = 1000
    )

# This will use the default values for argument: n, thin, nbp, pilot and burn.
# The number of CPUs will be the number available - 1 (the default).

# To test the impact of unbalance sampling run BayeScan with subsampling,
# for this, you need to feed the function the tidy data frame generated above
# with radiator::genomic_converter:
scan.pops.sub <- radiator::run_bayescan(
    data = data$tidy.data,
    pr_odds = 1000,
    subsample = "min",
    iteration.subsample = 10
    )

# This will run BayeScan 10 times, and for each iteration, the number of individuals
# sampled in each pop will be equal to the minimal number found in the pops
# (e.g. pop1 N = 36, pop2 N = 50 and pop3 N = 15, the subsampling will use 15
# individuals in each pop, taken randomly.
# You can also choose a specific subsample value with the argument.
}
}
\references{
Foll, M and OE Gaggiotti (2008) A genome scan method to identify
selected loci appropriate
for both dominant and codominant markers: A Bayesian perspective.
Genetics 180: 977-993

Foll M, Fischer MC, Heckel G and L Excoffier (2010)
Estimating population structure from
AFLP amplification intensity. Molecular Ecology 19: 4638-4647

Fischer MC, Foll M, Excoffier L and G Heckel (2011) Enhanced AFLP
genome scans detect
local adaptation in high-altitude populations of a small rodent (Microtus arvalis).
Molecular Ecology 20: 1450-1462
}
\seealso{
\href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
}
