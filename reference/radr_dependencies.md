# Check radr dependencies

Reports the required and optional R packages used by radr and checks
whether the optional `bcftools` and BayeScan executables are available.
The function is diagnostic: it does not install packages or modify a
Conda environment.

## Usage

``` r
radr_dependencies(verbose = TRUE)
```

## Arguments

- verbose:

  Logical. Print guidance for unavailable components. Default:
  `verbose = TRUE`.

## Value

A tibble with the component, source, requirement level, principal use,
and availability.
