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
