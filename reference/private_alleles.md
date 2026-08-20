# Find private alleles

The function highlight private alleles by strata, using a GDS or tidy
file or object.

## Usage

``` r
private_alleles(data, strata = NULL, verbose = TRUE)
```

## Arguments

- data:

  A tidy genomic data frame or another genomic object supported by the
  calling function.

- strata:

  (path or object) The strata file or object. Additional documentation
  is available in
  [`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.html).
  Use that function to whitelist/blacklist populations/individuals.
  Option to set `pop.levels/pop.labels` is also available.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

## Value

A list with an object highlighting private alleles by markers and strata
and a second object with private alleles summarized by strata.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
corals.private.alleles.by.pop <- radr::private_alleles(data = tidy, strata = strata.pop)
} # }
```
