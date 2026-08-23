# blacklist_hw

blacklist hw

## Usage

``` r
blacklist_hw(
  x,
  unfiltered.data,
  data.temp,
  hw.pop.threshold,
  path.folder = NULL,
  pop.id.levels,
  verbose = TRUE
)
```

## Arguments

- x:

  HWE summary grouped by significance threshold.

- unfiltered.data:

  Original tidy genomic data.

- data.temp:

  Optional data retained outside the HWE calculations.

- hw.pop.threshold:

  Number of strata allowed to depart from HWE.

- path.folder:

  Output directory. Default: `path.folder = NULL`.

- pop.id.levels:

  Ordered strata levels.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
