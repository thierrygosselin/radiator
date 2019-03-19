[![Travis-CI Build Status](https://travis-ci.org/thierrygosselin/radiator.svg?branch=master)](https://travis-ci.org/thierrygosselin/radiator) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/radiator?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/radiator) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/radiator)](http://cran.r-project.org/package=radiator) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![DOI](https://zenodo.org/badge/14548/thierrygosselin/radiator.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/radiator)

[![packageversion](https://img.shields.io/badge/Package%20version-1.0.0-orange.svg)](commits/master) [![Last-changedate](https://img.shields.io/badge/last%20change-2019--03--18-brightgreen.svg)](/commits/master)

------------------------------------------------------------------------

radiator: an R package for RADseq Data Exploration, Manipulation and Visualization
==================================================================================

This is the development page of the **radiator**, if you want to help, see [contributions section](https://github.com/thierrygosselin/radiator#contributions)

Most genomic analysis look for patterns and trends with various statistics. Bias, noise and outliers can have bounded influence on estimators and interfere with polymorphism discovery. Avoid bad data exploration and control the impact of filters on your downstream genetic analysis. Use radiator to: import, explore, manipulate, visualize, filter, impute and export your GBS/RADseq data.

Installation
------------

To try out the dev version of **radiator**, copy/paste the code below:

``` r
if (!require("devtools")) install.packages("devtools") # to install
devtools::install_github("thierrygosselin/radiator")
library(radiator)
```

Learning radiator
-----------------

See if radiator as the right tools for you:

**1. Prepare a strata file**

-   It's a tab separated file, e.g. `radiator.strata.tsv`.
-   A minimum of 2 columns: `INDIVIDUALS` and `STRATA` is required.
-   The `STRATA` column identifies the individuals stratification, the hierarchical groupings: populations, sampling sites or any grouping you want.
-   It's like *stacks* population map file with header...

To make sure it's going to work properly, try reading it in `R` with:

``` r
strata <- radiator::read_strata("my.strata.tsv")
names(strata)
# more details in with `??radiator::read_strata`
```

**2. Filter your RADseq data: ONE FUNCTION TO RULE THEM ALL**

``` r
data <- radiator::filter_rad(data = "my.vcf", strata = "my.strata.tsv", output = c("genind", "hierfstat"))
```

-   There's a built-in interactive mode that requires users to visualize the data before choosing thresholds.
-   The function is made of modules (see below) that user's can access separately or in combination.
-   Use [magrittr](https://magrittr.tidyverse.org) `%>%` to chain filtering functions together and dig deeper into your data [see vignettes](https://github.com/thierrygosselin/radiator#vignettes-r-notebooks-and-examples)
-   But remember, for 95% of users, `filter_rad` will be enough to start exploring the biology!

Overview
--------

<table style="width:100%;">
<colgroup>
<col width="26%" />
<col width="73%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Caracteristics</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><strong>Import</strong></td>
<td align="left">List of the 11 supported input genomic file formats and their variations:<br> <a href="https://samtools.github.io/hts-specs/">VCF: SNPs and haplotypes</a> (Danecek et al., 2011)<br><a href="http://www.diversityarrays.com">DArT files: genotypes in 1row, alleles counts and coverage in 2 rows and SilicoDArT</a><br><a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr">PLINK: bed/tped/tfam</a> (Purcell et al., 2007)<br><a href="https://github.com/thibautjombart/adegenet">genind</a> (Jombart et al., 2010; Jombart and Ahmed, 2011)<br> <a href="https://github.com/thibautjombart/adegenet">genlight</a> (Jombart et al., 2010; Jombart and Ahmed, 2011)<br><a href="https://github.com/EricArcher/strataG">strataG gtypes</a> (Archer et al., 2016)<br><a href="http://genepop.curtin.edu.au">Genepop</a> (Raymond and Rousset, 1995; Rousset, 2008)<br><a href="http://catchenlab.life.illinois.edu/stacks/">STACKS haplotype file</a> (Catchen et al., 2011, 2013)<br><a href="https://github.com/jgx65/hierfstat">hierfstat</a> (Goudet, 2005)<br><a href="https://github.com/zhengxwen/SeqArray">SeqArray</a> (Zheng et al., 2017)<br><a href="https://github.com/zhengxwen/SNPRelate">SNPRelate</a> (Zheng et al., 2012)<br>Dataframes of genotypes in wide or long/tidy format<br>Reading and tidying is found inside: <code>genomic_converter</code>, <code>tidy_</code> and <code>read_</code> functions</td>
</tr>
<tr class="even">
<td align="left"><strong>Output</strong></td>
<td align="left">25 genomic data formats can be exported out of <strong>radiator</strong>. The module responsible for this is <code>genomic_converter</code>. Separate modules handles the different formats and are also available:<br><code>write_vcf</code>: <a href="https://samtools.github.io/hts-specs/">VCF</a> (Danecek et al., 2011)<br><code>write_plink</code>: <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr">PLINK tped/tfam</a> (Purcell et al., 2007)<br><code>write_genind</code>: <a href="https://github.com/thibautjombart/adegenet">adegenet genind and genlight</a> (Jombart et al., 2010; Jombart and Ahmed, 2011)<br><code>write_genlight</code>: <a href="https://github.com/thibautjombart/adegenet">genlight</a> (Jombart et al., 2010; Jombart and Ahmed, 2011)<br><code>write_gsi_sim</code>: <a href="https://github.com/eriqande/gsi_sim">gsi_sim</a> (Anderson et al. 2008)<br><code>write_gtypes</code>: <a href="https://github.com/EricArcher/strataG">strataG gtypes</a> (Archer et al. 2016)<br><code>write_colony</code>: <a href="https://www.zsl.org/science/software/colony">COLONY</a> (Jones and Wang, 2010; Wang, 2012)<br><code>write_genepop</code>: <a href="http://genepop.curtin.edu.au">Genepop</a> (Raymond and Rousset, 1995; Rousset, 2008)<br><a href="http://catchenlab.life.illinois.edu/stacks/">STACKS haplotype file</a> (Catchen et al., 2011, 2013)<br><code>write_betadiv</code>: <a href="http://adn.biol.umontreal.ca/~numericalecology/Rcode/">betadiv</a> (Lamy, 2015)<br> <code>vcf2dadi</code>: <a href="http://gutengroup.mcb.arizona.edu/software/">δaδi</a> (Gutenkunst et al., 2009)<br> <code>write_structure</code>: <a href="http://pritchardlab.stanford.edu/structure.html">structure</a> (Pritchard et al., 2000)<br> <code>write_arlequin</code>: <a href="http://cmpg.unibe.ch/software/arlequin35/">Arlequin</a> (Excoffier et al. 2005)<br> <code>write_hierfstat</code>: <a href="https://github.com/jgx65/hierfstat">hierfstat</a> (Goudet, 2005)<br> <code>write_snprelate</code>: <a href="https://github.com/zhengxwen/SNPRelate">SNPRelate</a> (Zheng et al. 2012)<br><code>write_seqarray</code>: <a href="https://github.com/zhengxwen/SeqArray">SeqArray</a> (Zheng et al. 2017)<br> <code>write_bayescan</code>: <a href="http://cmpg.unibe.ch/software/BayeScan">BayeScan</a> (Foll and Gaggiotti, 2008)<br><code>write_pcadapt</code>: <a href="https://github.com/bcm-uga/pcadapt">pcadapt</a> (Luu et al. 2017)<br><code>write_hzar</code> (Derryberry et al. 2013) <br><code>write_fineradstructure</code> (Malinsky et al., 2018) <br><code>write_related</code> <a href="https://github.com/timothyfrasier/related">related</a> (Pew et al., 2015)<br><code>write_stockr</code> for stockR package (Foster el al., submitted)<br><code>write_maverick</code> <a href="http://www.bobverity.com/home/maverick/what-is-maverick/">MavericK</a> (Verity &amp; Nichols, 2016)<br><code>write_ldna</code> <a href="https://github.com/petrikemppainen/LDna">LDna</a> (Kemppainen et al. 2015)<br>Dataframes of genotypes in wide or long/tidy format</td>
</tr>
<tr class="odd">
<td align="left"><strong>Conversion function</strong></td>
<td align="left"><code>genomic_converter</code> import/export genomic formats mentioned above. The function is also integrated with usefull filters, blacklist and whitelist.</td>
</tr>
<tr class="even">
<td align="left"><strong>Outliers detection</strong></td>
<td align="left"><code>detect_duplicate_genomes</code>: detect and remove duplicate individuals from your dataset <br><code>detect_mixed_genomes</code>: detect and remove potentially mixed individuals<br><code>stackr::summary_haplotype</code> and <code>filter_snp_number</code>: Discard of outlier markers with <em>de novo</em> assembly artifact (e.g. markers with an extreme number of SNP per haplotype or with irregular number of alleles)</td>
</tr>
<tr class="odd">
<td align="left"><strong>Filters</strong></td>
<td align="left">Targets of filters: alleles, genotypes, markers, individuals and populations and associated metrics and statistics can be filtered and/or selected in several ways inside the main filtering function <code>filter_rad</code> and/or the underlying modules:<br></td>
</tr>
</tbody>
</table>

<br>`filter_rad`: designed for RADseq data, it's the *one function to rule them all*. Best used with unfiltered or very low filtered VCF file as input. The function can handle very large VCF files (e.g. no problem with &gt;2M SNPs, &gt; 30GB files), all within R!! <br>`filter_dart_reproducibility`: blaclist markers under a certain threshold of DArT reproducibility metric. <br>`filter_monomorphic`: blacklist markers with only 1 morph. <br>`filter_common_markers`: keep only markers common between strata. <br>`filter_individuals`: blacklist individuals based on missingness, heterozygosity and/or total coverage. <br>`filter_mac`: blacklist markers based on minor/alternate allele count. <br>`filter_coverage`: blacklist markers based on mean read depth (coverage). <br>`filter_genotype_likelihood`: Discard markers based on genotype likelihood <br>`filter_genotyping`: blacklist markers based on genotyping/call rate. <br>`filter_snp_position_read`: blacklist markers based based on the SNP position on the read/locus. <br>`filter_snp_number`: blacklist locus with too many SNPs. <br>`filter_ld`: blacklist markers based on short and/or long distance linkage disequilibrium. <br>`filter_hwe`: blacklist markers based on Hardy-Weinberg Equilibrium expectations (HWE). <br>`filter_het`: blacklist markers based on the observed heterozygosity (Het obs) <br>`filter_fis`: blacklist markers based on the inbreeding coefficient (Fis) <br>`filter_whitelist`: keep only markers present in a whitelist| |**[ggplot2](http://ggplot2.org)-based plotting**|Visualize distribution of important metric and statistics and create publication-ready figures| |**Parallel**|Codes designed and optimized for fast computations using Genomic Data Structure [GDS](http://zhengxwen.github.io/gdsfmt) file format and data science packages in [tiverse](https://www.tidyverse.org) Works with all OS: Linux, Mac and now PC!|

[More in radiator workflow below](https://github.com/thierrygosselin/radiator#radiator-workflow)

Life cycle
----------

**Missing data**

Visualizing missing data and it's imputations requires special attention that fall outside the scope of *radiator*. Inside my package called [grur](https://github.com/thierrygosselin/grur), users can **visualize patterns of missingness** associated with different variables (lanes, chips, sequencers, populations, sample sites, reads/samples, homozygosity, etc). Several **Map-independent imputations** of missing genotypes are available: **Random Forests** (on-the-fly-imputations or predictive modeling), **Extreme Gradient Tree Boosting**, Strawman imputations (~ max/mean/mode: the most frequently observed, non-missing genotypes is used). Imputations can be conducted **overall samples** or **by populations/strata/grouping**. `radiator::genomic_converter` is integrated with the imputation function of **grur**.

Prerequisite - Suggestions - Troubleshooting
--------------------------------------------

-   **Parallel computing**: follow the steps in this [notebook vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0) to install the packages with OpenMP-enabled compiler and conduct imputations in parallel.
-   [Installation problems.](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0)
-   **Windows users**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
-   The R GUI is unstable with functions using parallel ([more info](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html)), so I recommend using [RStudio](https://www.rstudio.com/products/rstudio/download/) for a better experience.
-   Running codes in chunks inside R Notebook might cause problem, run it outside in the console.

Vignettes, R Notebooks and examples
-----------------------------------

**Vignettes (in development, check periodically for updates):**

-   Vignettes with real data for example in the form of R Notebooks take too much space to be included in package, without CRAN complaining. Consequently, vignettes are gradually being excluded from the package and distributed separately, follow the links below.
-   **installation problems** [notebook vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0)
-   **parallel computing during imputations** [notebook vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0)
-   **vcf2dadi** [Rmd](https://www.dropbox.com/s/bl0mv6kavz97ibz/vignette_vcf2dadi.Rmd?dl=0) or [html](https://www.dropbox.com/s/xbgxk2valwl5o44/vignette_vcf2dadi.html?dl=0)

**R Notebooks:**

-   Missing data visualization and analysis [(html vignette)](https://www.dropbox.com/s/4zf032g6yjatj0a/vignette_missing_data_analysis.nb.html?dl=0) and [(Rmd)](https://www.dropbox.com/s/5fxw2h9w1l1j391/vignette_missing_data_analysis.Rmd?dl=0)

Citation:
---------

To get the citation, inside R:

``` r
citation("radiator")
```

New features
------------

Change log, version, new features and bug history lives in the [NEWS.md file](https://github.com/thierrygosselin/radiator/blob/master/NEWS.md)

Roadmap of future developments:
-------------------------------

-   Updated filters: more efficient, interactive and visualization included: *in progress*.
-   Workflow tutorial that links functions and points to specific vignettes to further explore some problems: *in progress*
-   Integration of several functions with [STACKS](http://catchenlab.life.illinois.edu/stacks/) and [DArT](http://www.diversityarrays.com) database *in progress*.
-   Use Shiny and ggvis (when subplots and/or facets becomes available for ggvis).
-   Until publication **radiator** will change rapidly, stay updated with releases and contribute with bug reports.
-   Suggestions ?

Contributions:
--------------

This package has been developed in the open, and it wouldn’t be nearly as good without your contributions. There are a number of ways you can help me make this package even better:

-   If you don’t understand something, please let me know.
-   Your feedback on what is confusing or hard to understand is valuable.
-   If you spot a typo, feel free to edit the underlying page and send a pull request.

New to pull request on github ? The process is very easy:

-   Click the edit this page on the sidebar.
-   Make the changes using github’s in-page editor and save.
-   Submit a pull request and include a brief description of your changes.
-   “Fixing typos” is perfectly adequate.

radiator workflow
-----------------

The **radiator** package fits currently at the end of the GBS workflow (e.g. after running [STACKS](http://catchenlab.life.illinois.edu/stacks/) inside R with [stackr](https://github.com/thierrygosselin/stackr).

**Table 1: Quality control and filtering RAD/GBS data**

<table style="width:86%;">
<colgroup>
<col width="8%" />
<col width="9%" />
<col width="9%" />
<col width="9%" />
<col width="9%" />
<col width="9%" />
<col width="9%" />
<col width="9%" />
<col width="9%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Parameters</th>
<th align="center">Libraries &amp; Seq.Lanes</th>
<th align="center">Alleles</th>
<th align="center">Genotypes</th>
<th align="center">Individuals</th>
<th align="center">Markers</th>
<th align="center">Sampling sites</th>
<th align="center">Populations</th>
<th align="center">Globally</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Quality</td>
<td align="center">x</td>
<td align="center"></td>
<td align="center"></td>
<td align="center">x</td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
</tr>
<tr class="even">
<td align="left">Assembly and genotyping</td>
<td align="center">x</td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
</tr>
<tr class="odd">
<td align="left">Outliers</td>
<td align="center"></td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
</tr>
<tr class="even">
<td align="left">Missingness</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
</tr>
<tr class="odd">
<td align="left">Coverage</td>
<td align="center"></td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center"></td>
<td align="center">x</td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
</tr>
<tr class="even">
<td align="left">Genotype Likelihood</td>
<td align="center"></td>
<td align="center"></td>
<td align="center">x</td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
</tr>
<tr class="odd">
<td align="left">Prop. Genotyped</td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
</tr>
<tr class="even">
<td align="left">HET &amp; FIS &amp; HWE</td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center"></td>
<td align="center">x</td>
<td align="center"></td>
</tr>
<tr class="odd">
<td align="left">MAF</td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center"></td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
</tr>
<tr class="even">
<td align="left">Missingness</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
<td align="center">x</td>
</tr>
</tbody>
</table>

**Step 1 Quality** Ask yourself these questions:

-   Quality of : DNA, libraries and sequencing lanes ?
-   Please, stop thinking in terms of quantity (e.g. millions of reads returned), and prioritize and think more about the actual quality of your new data.
-   Read about the quality metrics used in available software (e.g. fastqc)

**Step 2 *de novo* assembly and genotyping**

-   This part is conducted outside radiator (e.g. using [STACKS](http://catchenlab.life.illinois.edu/stacks/) inside R with the package [stackr](https://github.com/thierrygosselin/stackr))
-   Software pipelines include: [STACKS](http://catchenlab.life.illinois.edu/stacks/), [pyRAD](http://dereneaton.com/software/), [dDocent](https://ddocent.wordpress.com), [AftrRAD](http://u.osu.edu/sovic.1/downloads/). If you want to develop your own pipeline, there are a multitude of approaches, good luck.
-   At the end of the pipeline, use liberal filter thresholds the go in radiator and do the heavy lifting.

**Step 3 Pattern of missingness**

-   `grur::missing_visualization`: I really like running this function after modifying my RAD data, to make sure bias were not introduce.
-   The trick here is to use the `strata` argument to find patterns associated with different variables of your study (lanes, chips, sequencers, populations, sample sites, reads/samples, etc).
-   Do you see a trend between your missing pattern and reads/samples ? Heterozygosity?
-   Do you need more sequencing? Do you have to re-run some lanes?
-   Usually for this first run I only use the blacklists of ID and markers to start filtering with individuals and markers that won't drag down the polymorphism discovery. I re-introduce individuals at the end of the pipeline and re-run `grur::missing_visualization` to see what the new analysis reveal.

**Step 4 Outliers**

-   Remove replicates (I hope you have some).
-   Remove *de novo* assembly artifact:
-   run `stackr::summary_haplotypes` to automatically generate blacklist of genotypes and whitelist of markers. The function will highlight individuals and locus with more than 2 alleles (outlier individuals and markers).
-   run `filter_snp_number`, function will highlight outlier locus/reads with extreme number of SNP/read or haplotypes
-   Remove potential duplicated samples that went off your radar with `detect_duplicate_genomes`.
-   Remove mixed samples or pooled samples that creates outliers individual's heterozygosity with the function `detect_mixed_individuals`.

**Step 5 Pattern of missingness**

-   Re-run `grur::missing_visualization` with/without your new blacklists (e.g. of genotypes, individuals) and with/without whitelist of markers to examine patterns of missingness in you dataset before more extensive filtering (there is a vignette for this step)
-   The trick here is to use the `strata` argument to find patterns associated with different variables of your study (lanes, chips, sequencers, populations, sample sites, reads/samples, etc).
-   Do you see a trend between your missing pattern and reads/samples ? Heterozygosity?

**Step 6: Metrics and statistics, some thoughts** \* Metrics: what you're observing so far is it *de novo* artefact or a reliable signal of biological polymorphism? \* Statistics: are you going to use haplotype or snp level statistics? Should the statistic you are interested in be consistent throughout the read ? \* Use `snp.ld` argument in several of radiator functions to throughly test the consistensies of SNPs statistics among haplotype.

**Step 7 Coverage and Genotype Likelihood**

-   Coverage is an individual metric. With most software you'll find allele and genotype coverage info.
-   Genotype likelihood is usually a metric based on coverage of the different genotypes found in all of your data. Since it's v.1.45, STACKS no longer output the useful GL metric inside the VCF. It was using only one number to qualify the genotype, while most other pipelline using GL/PL are using 3 numbers for homozygous REF, heterozygous and homozygous ALT genotypes.
-   Good allele coverage is required for reliable genotypes.
-   Reliable genotypes is required for reliable downstream summary statistics.
-   If your data allows it (you have coverage and/or genotype likelihood metrics), explore filtering options in `filter_coverage` and `filter_genotype_likelihood`.

**Step 8 Prop. Genotyped**

-   Use the functions `filter_individual` and `filter_population` to explore if you have enough individuals and enough putative populations for markers filtering.
-   Use blacklist of individuals with different thresholds.
-   Keep different whitelist of markers.
-   Use `common.markers` argument inside most of radiator functions to test the impact of vetting loci based on shared markers. This can be use strategically for Fst calculations.
-   Use imputation methods provided by radiator (inside `tidy_genomic_data` or `genomic_converter`, as a separate module: `radiator_imputations_module`) to assess the impact of lowering or increasing threshold that impact missing data.

**Step 9 HET, Fis, HWE**

-   Overall and/or per populations hwe, heterozygosity and Fis statistics can highlight: *de novo* assembly problems (oversplitting/undermerging), genotyping problems or biological problems.
-   These filters allows to test rapidly if departure from realistic expectations are a problem for downstream analysis ?
-   Choose your threshold wisely and test impact on pipeline.
-   Use `filter_het`, `filter_fis`, `filter_hwe` and look again at the individual's heterozygosity (`detect_mixed_individuals`) for outliers.
-   Hardy-Weinberg Equilibrium: this analysis as [several underlying assumptions](https://en.wikipedia.org/wiki/Hardy–Weinberg_principle). Please do not conduct analysis with sampling sites. Note: most natural populations are violating one or more of the assumptions.

**Step 10 MAF**

-   Remove artifactual and uninformative markers.
-   Use MAF arguments inside several of radiator functions to tailor MAF to your analysis tolerance to minor allelel frequencies.
-   There is also a separate filter in radiator: `filter_maf`
-   I usually run the filter to explore and understand the impact of the different thresholds on the data. Then use different ones inside the different functions of radiator.

**Step 11 Pattern of missingness, yes... again!**

-   Use `grur::missing_visualization` with your latest blacklists (e.g. of genotypes, individuals) and with your latest whitelist of markers to examine patterns of missingness in your dataset after filtering.
-   Hopefully, you will have remove all the bias with the filters.
