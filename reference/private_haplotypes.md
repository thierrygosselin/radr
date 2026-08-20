# private haplotypes

Minor Allele Frequency filter.

## Usage

``` r
private_haplotypes(data, strata = NULL, verbose = TRUE)
```

## Arguments

- data:

  A tidy genomic dataframe Used internally in
  [radr](https://github.com/thierrygosselin/radr) and might be of
  interest for users. *How to get a tidy data frame ?* Look into
  genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.html).

- strata:

  (optional) A strata file or object in the global environment. The
  strata is a tab-separated data frame with 2 columns: `INDIVIDUALS` and
  `STRATA`. If used, the strata will replace the current STRATA or
  POP_ID in the dataset. Use this argument if you want to find private
  haplotypes on another hierarchical level, other than POP_ID.

## Value

A list with private haplotypes per markers and strata and a summary of
overall number of private haplotypes per strata.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
