# Prune dataset based on LD.

Used internally in [radr](https://github.com/thierrygosselin/radr) Prune
dataset based on LD.

## Usage

``` r
ld_pruning(ld.tibble = NULL, stats = NULL, ld.threshold = 0.8, verbose = TRUE)
```

## Arguments

- ld.tibble:

  (path) The markers LD pairwise data. Default: `ld.tibble = NULL`.

- stats:

  (path) The markers missingness info statistics. Default:
  `stats = NULL`.

- ld.threshold:

  (double) The threshold to prune SNPs in LD. Default:
  `ld.threshold = 0.8`.

- verbose:

  (logical, optional) Default: `verbose = TRUE`.

## Value

A list with blacklisted SNPs. Write the blacklist in the working
directory.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
