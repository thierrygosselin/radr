# Compute allele frequencies per markers and populations

Compute allele frequencies per markers and populations. Used internally
in [radr](https://github.com/thierrygosselin/radr) and might be of
interest for users.

## Usage

``` r
allele_frequencies(
  data,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1
)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.html).

- verbose:

  (optional, logical) `verbose = TRUE` to be chatty during execution.
  Default: `verbose = TRUE`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

## Value

A list with allele frequencies in a data frame in long and wide format,
and a matrix. Local (pop) and global minor allele frequency (MAF) is
also computed.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
