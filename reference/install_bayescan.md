# Install BayeScan with Conda or Mamba

Creates a dedicated environment and installs the Bioconda BayeScan 2.1
package. Installation is performed only when this function is called
explicitly.

## Usage

``` r
install_bayescan(
  conda.env = "genomics",
  conda = NULL,
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- conda.env:

  Name of the environment to create. Default: `conda.env = "genomics"`.

- conda:

  Optional path to a Conda, Mamba, or Micromamba executable. Default:
  `conda = NULL`.

- force:

  Logical. Re-run the installation even when BayeScan is already
  available in the environment. Default: `force = FALSE`.

- verbose:

  Logical. Display installation progress. Default: `verbose = TRUE`.

## Value

Invisibly returns the installed BayeScan executable path.
