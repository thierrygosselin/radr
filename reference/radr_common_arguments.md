# Common arguments used by radr

This documentation-only helper centralizes parameters shared by genomic
screening, filtering, and diagnostic functions.

## Usage

``` r
radr_common_arguments(
  interactive.filter = TRUE,
  gds = NULL,
  data = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  random.seed = NULL,
  ...
)
```

## Arguments

- interactive.filter:

  Logical indicating whether an interactive filtering session may
  display diagnostics and ask for thresholds. Default:
  `interactive.filter = TRUE`.

- gds:

  A genome GDS file path or object supported by `genometranslator`.

- data:

  A tidy genomic data frame or another genomic object supported by the
  calling function.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- random.seed:

  Optional integer seed for reproducible operations. Default:
  `random.seed = NULL`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

`NULL`, invisibly. This function exists to share documentation.
