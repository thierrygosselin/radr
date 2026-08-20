# Locate and validate BayeScan

Searches an explicit path, the system `PATH`, and optionally a Conda
environment for a BayeScan 2.1 executable.

## Usage

``` r
check_bayescan(
  bayescan.path = NULL,
  conda.env = "genomics",
  conda = NULL,
  verbose = TRUE
)
```

## Arguments

- bayescan.path:

  Explicit path to a BayeScan executable. When `NULL`, the executable is
  discovered automatically. Default: `bayescan.path = NULL`.

- conda.env:

  Conda environment name or prefix containing BayeScan. Default:
  `conda.env = "genomics"`.

- conda:

  Optional path to a Conda, Mamba, or Micromamba executable. Default:
  `conda = NULL`.

- verbose:

  Logical. Display the detected executable. Default: `verbose = TRUE`.

## Value

The normalized BayeScan executable path.
