# strip_rad

Strip a tidy data set of it's strata and markers meta. Used internally.

## Usage

``` r
strip_rad(
  x,
  m = c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"),
  env.arg = NULL,
  keep.strata = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  The data

- m:

  (character, string) The variables part of the markers metadata.

- env.arg:

  You want to redirect
  [`rlang::current_env()`](https://rlang.r-lib.org/reference/stack.html)
  to this argument.

- keep.strata:

  (logical) Keep the strata in the dataset or remove the info and keep
  only the sample ids.

- verbose:

  (logical) The function will chat more when allowed.
