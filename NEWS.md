# radiator v.0.0.9 2018-03-01

* New output file: related 
  This output file format enables to run the data in the related R package, which is
  essantially the R version of *COANCESTRY* fortran program developed by Jinliang Wang.


# radiator v.0.0.8 2018-02-08

* New output file: fineRADstructure

# radiator v.0.0.7 2018-02-06

* Better VCF parsing
* `radiator` is ready for stacks v.2 beta8
* starting to re-introduce `data.table` to increase speed during melting or casting data frames

# radiator v.0.0.6 2017-10-02

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


# radiator v.0.0.5 2017-09-20

* `tidy_dart` and `filter_dart`: 
    * big overhaul.
    * importing 1 and 2 row genotypes (called sometimes binary format) is easier
    and no longer require to prep the DArT file
    (the function parse DArT data with `*` at the beginning of certain lines).
    * filtering is easier and now interactive for those who want.
    * plots are generated automatically. 
* `tidy_genomic_data` and `genomic_converter` can now accept DArT data directly.
* several typos where fix. Thanks to @IdoBar for this.

# radiator v.0.0.4 2017-09-07

* new function: `run_bayescan` to run *BayeScan* ... with *radiator*
* new separate `write` functions also included in `genomic_converter`: `write_bayescan`, `write_pcadapt` and `write_hzar`.

# radiator v.0.0.3 2017-07-18

* small mod
* writting in specific folder
* writting specific plot

# radiator v.0.0.2 2017-07-15

* Work on different functions to prep for official package launch


# radiator v.0.0.1 2017-06-18

* First commit
