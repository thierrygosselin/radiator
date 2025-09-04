# radiator 1.4.0 2025-09-04
* bug fix when using the new DArT data and converting to genind


# radiator 1.3.9 2025-09-02
* Improved DArT file import
* New flavor of strata file import for DArT data


# radiator 1.3.8 2025-07-01
* Improved the use of arrow
* Better reporting


# radiator 1.3.7 2025-06-17

* Improved DArT file import and filtering
* Several functions were not handling generating folder names and numbers inside 
nested function. Presently working on an overhaul. 
VCF and DArT files are impacted by this one. 


# radiator 1.3.6 2025-04-15

* Speed improvements for some functions
* gradually replacing data.table and fst package: always difficult to install to work in parallel
* now using arrow package to store tibble...

# radiator 1.3.5 2024-11-06

* Fix the one column matrices in dplyr::filter warning.
* Fix the error when dplyr::select couldn't find the column STRATA.

# radiator 1.3.4 2024-06-05

* Open the `parallel.core` argument for some internal functions to help windows users with parallel processing problems in R. #188
* Work around for `.DynamicClusterCall` pushed greavess #189


# radiator 1.3.3 2024-05-15

* Fix issue #188 related to coverage in DArT 1row and 2rows format




# radiator 1.3.2 2024-04-30

* works with R 4.3.4
* Fix issue #186 related some particular DArT files
* Fix issue #187 related to sexy_markers and VCF files



# radiator 1.3.1 2024-04-24

* works with R 4.3.3
* Updated DArT code that use COUNT files to check more for problematic markers usually stemming from merged projects


# radiator 1.3.0 2024-02-22

* Bug fix using coverage and DArT data


# radiator 1.2.9 2024-01-25

* Bug fix stemming from genalex files and genind conversion


# radiator 1.2.8 2023-04-03

* Several bug fix while reading VCF from ipyrad



# radiator 1.2.7 2023-03-20

* Additional checks during DArT file reading.


# radiator 1.2.6 2023-03-13

* Bug fix with `pcadapt`, `bayescan`, `genlight` output thanks to @jcaccavo. 
Confusion between `POP_ID` and `STRATA` remained for some less used functions while
the migration towards `STRATA` only inside radiator.


# radiator 1.2.5 2023-01-20

* huge work on `filter_ma` that now incorporate Minor Allele Frequency (MAF), 
Minor Allele Count (MAC) and Minor Allele Depth (MAD).


# radiator 1.2.4 2022-11-17

* bug fix when using some tidy function from GDS, the connection was not closing properly
* work on genind, genlight, genepop, hierfstat, arlequin functions
* updated the vignette


# radiator 1.2.3 2022-10-30

* bug fix with reading VCF, the new update works with SeqArray (>= 1.36.0)
* work on Tidy VCF.


# radiator 1.2.2 2021-12-20

* bug fix when using *arlequin* output file
* bug fix using DArT data

# radiator 1.2.1 2021-10-28

* bug fix when using genepop file and `genomic_converter`
* bug fix in `genomic_converter` when using genlight as input and PLINK as output file
* More GATK VCF problems detections and warnings

# radiator 1.2.0 2021-06-16

* updated to work with R 4.1.0
* bug fix with some format during conversion when necessary genotype format wasn't found
* more safe results when importing VCF files with windows and parallel processing

# radiator 1.1.9 2020-12-14

* bug fix using plink files and future 
* continue to test future backend and carrier
* bug fix when using microsatellites as input and some specific output format


# radiator 1.1.8 2020-10-17

* No longer using Travis CI and AppVeyor to test the package
* R-CMD-check: now using GitHub actions and rhub that test on the 3 OS.
* Taking advantage of future, furrr and carrier packages

# radiator 1.1.7 2020-08-21

* PLINK files: fix a couple of bugs reading the tped
* removed `tidyr::gather` and `tidyr::spread` dependencies (they are deprecated)
* DArT data in 1-row format was not working properly with latest `data.table` melt function. 
Changed to `tidyr::pivot_long`.


# radiator 1.1.6 2020-06-23

* updated radiator so that it work with latest release of SeqArray (v.1.28.1), that
introduced breaking changes.

# radiator 1.1.5 2020-03-07

* `fineRADstructure` rules have changed, I'm no longer doing gymnastic to make the
population work with `genomic_converter` or `write_fineradstructure`, it's not up
to the user to make sure the pop id starts with a letter. 
Some codes might now be broken because of this.


# radiator 1.1.4 2020-02-04

* `read_plink` and `tidy_plink`: new rules to get the best out of PLINK **tped**
and **bed** files.

# radiator 1.1.3 2020-01-21

* `detect_microsatellites`: new radiator function that detect microsatellites using 
GMATA.


# radiator 1.1.0 2019-05-31

* 2 new output formats: `genepopedit` and `rubias`
* 1 new detect function: `detect_paralogs` that copy the method described in
McKinney et al. 2017. This function is the logical step to make available for users
after `detect_mixed_genomes`, `filter_hwe` and `sexy_markers`.
* `sexy_markers`: couple of bug fixed. More testing with the different genomic
format.

# radiator 1.0.0 2019-05-01

* There is still documentation and vignette to fix, but this is 
the release that will be submitted to CRAN
* Major `SeqArray` and `GDS` integration
* Imputation module was removed from `radiator` and now lives exclusively in package `grur`
* `filter_dart` is deprecated. Please use `filter_rad`, the ONE function to rule them all;)
* Worked on travis
* pkgdown website
* vignettes
* better function doc

# radiator 0.0.21 2018-11-14

* `snp.ld`: using missing data now works by chromosome/scaffold.


# radiator 0.0.20 2018-11-12

* `tidy_vcf`, `tidy_genomic_data` and `genomic_converter`: works better with ipyrad vcf's


# radiator 0.0.19 2018-11-07

* `tidy_vcf`, `tidy_genomic_data` and `genomic_converter`: works better with freebayes and stacks vcf


# radiator 0.0.18 2018-10-23

* `tidy_vcf`, `tidy_genomic_data` and `genomic_converter`: work without strata/pop grouping
* `snp_ld`: new argument `ld.threshold` for long.distance linkage disequilibrium.
* `write_vcf`: will now output ID as LOCUS_COL or LOCUS_(POS-1) if COL info is not provided.

# radiator 0.0.17 2018-10-04

* `tidy_vcf`, `tidy_genomic_data` and `genomic_converter`: way faster with huge VCF
* `tidy_vcf` and `write_seqarray`: work better with vcf generated by Stacks
* `write_ldna`: new function that generates a LDna input file from a tidy data frame.

# radiator 0.0.16 2018-09-04

* `tidy_vcf`, `tidy_genomic_data` and `genomic_converter`: way faster with huge VCF
* `write_fineradstructure`: fix bug when data was from DArT


# radiator 0.0.15 2018-08-17

* `genomic_converter`, `tidy_genomic_data`: bug fix when individuals are integers
* `fis_summary`: arguments updated


# radiator 0.0.14 2018-08-04

* `filter_dart` and `filter_rad`: can now opt out of HWE filtering. 
* `filter_hwe`: user can opt to see figures but skip filtering



# radiator 0.0.13 2018-07-09

* working to make radiator work correctly with ggplot2 v.3.0.0
* radiator ready for R 3.5.1 "Feather Spray" released on 2018/07/05


# radiator 0.0.12 2018-06-21

* transferred `write_gsi_sim` from assigner to radiator
* 2 `write` functions: `write_snprelate` and `write_seqarray`


# radiator 0.0.11 2018-04-26

* when individuals in strata file/object and data don't match, an error is generated
* to reduce `radiator` dependencies, several packages were moved in the **Suggests** field.
* `filter_dart`: lots of new stuff. More appropriate filter arguments:
    * `filter.coverage` is deprecated in favour of: filter.markers.coverage
    * `filter.ind.missing.geno` is deprecated in favour of: `filter.markers.missing`
    * `erase.genotypes`: new argument tailored to handle coverage for DArT counts data.
    * `filter.individuals.missing`: allows to blacklist sample early in the filtering
    pipeline.


# radiator 0.0.10 2018-03-06

* New output file: stockr 
  This output file format enables to run the data in the stockR package from Scott Fisher at CSIRO in Hobart.


# radiator 0.0.9 2018-03-01

* New output file: related 
  This output file format enables to run the data in the related R package, which is
  essantially the R version of *COANCESTRY* fortran program developed by Jinliang Wang.


# radiator 0.0.8 2018-02-08

* New output file: fineRADstructure

# radiator 0.0.7 2018-02-06

* Better VCF parsing
* `radiator` is ready for stacks v.2 beta8
* starting to re-introduce `data.table` to increase speed during melting or casting data frames

# radiator 0.0.6 2017-10-02

* `tidy_genomic_data` and `genomic_converter`: bug fix and improvements.
* `write_snprelate` and it's associated output from `genomic_converter` is no longer
available because of difficulties students are having in installing the package,
it's dependencies and/or using the output. The code remains available [here](https://www.dropbox.com/s/7xujizkvpi0ddac/write_snprelate.R?dl=0),
without support.
* `radiator::snp.ld`: the module inside `tidy_genomic_data` and `genomic_converter`
to minimize short linkage disequilibrium as a new option `snp.ld = "middle"`. 
For locus with > 2 SNPs/read the option allows to select at random one SNP between 
the first and the last SNP on the read. If the locus as <= 2 SNPs on the read,
the first one is selected. Note that for that last option, the numbers are reported.
Thanks to Ido Bar for the idea.


# radiator 0.0.5 2017-09-20

* `tidy_dart` and `filter_dart`: 
    * big overhaul.
    * importing 1 and 2 row genotypes (called sometimes binary format) is easier
    and no longer require to prep the DArT file
    (the function parse DArT data with `*` at the beginning of certain lines).
    * filtering is easier and now interactive for those who want.
    * plots are generated automatically. 
* `tidy_genomic_data` and `genomic_converter` can now accept DArT data directly.
* several typos where fix. Thanks to @IdoBar for this.

# radiator 0.0.4 2017-09-07

* new function: `run_bayescan` to run *BayeScan* ... with *radiator*
* new separate `write` functions also included in `genomic_converter`: `write_bayescan`, `write_pcadapt` and `write_hzar`.

# radiator 0.0.3 2017-07-18

* small mod
* writting in specific folder
* writting specific plot

# radiator 0.0.2 2017-07-15

* Work on different functions to prep for official package launch


# radiator 0.0.1 2017-06-18

* First commit
