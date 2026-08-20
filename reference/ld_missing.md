# Prune dataset based on LD and missingness.

Used internally in [radr](https://github.com/thierrygosselin/radr) Prune
dataset based on LD, uses missingness to keep 1 SNP. This is the
function that does it for 1 chrom. It's then used with purrr::map to run
seriall on all chrom.

## Usage

``` r
ld_missing(
  wl,
  data,
  ld.threshold = seq(0.1, 0.9, by = 0.1),
  denovo = NULL,
  ld.method = "r2",
  ld.figures = TRUE,
  path.folder = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```

## Arguments

- wl:

  The split dataset by chromosome.

- data:

  The GDS object.

- ld.threshold:

  The LD threshold.

- denovo:

  de novo data or not ? (logical).

- ld.method:

  The method to compute LD.

- ld.figures:

  (logical, optional) Default: `ld.figures = TRUE`.

- path.folder:

  (character) The path to the folder.

- parallel.core:

  The number of CPU.

- verbose:

  (logical, optional) Default: `verbose = TRUE`.

## Value

A list with whitelists and blacklists of markers in LD.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
