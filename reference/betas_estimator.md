# Legacy plural name for population-specific beta

`betas_estimator()` is retained temporarily for compatibility. New code
should use
[`beta_estimator()`](https://thierrygosselin.github.io/radr/reference/beta_estimator.md).

## Usage

``` r
betas_estimator(
  data,
  strata = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```

## Arguments

- data:

  A GDS filepath, an open `SeqVarGDSClass` object, or a tidy data frame
  containing `MARKERS`, `INDIVIDUALS`, `STRATA` (or `POP_ID`), and
  `ALT_DOSAGE`.

- strata:

  Optional strata filepath or data frame used to add or replace
  population assignments. It must contain `INDIVIDUALS` and `STRATA`.
  Default: `strata = NULL`.

- filename:

  Optional output prefix. When supplied, three tab-delimited files are
  written with suffixes `_beta.tsv`, `_within_population.tsv`, and
  `_between_populations.tsv`. Default: `filename = NULL`.

- parallel.core:

  Number of processor cores passed to
  [`genometranslator::read_genome()`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
  when file input must be imported. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress and a population beta summary. Default:
  `verbose = TRUE`.

## Value

The result returned by
[`beta_estimator()`](https://thierrygosselin.github.io/radr/reference/beta_estimator.md).
