# Filter data based on DArT reproducibility statistics

This filter removes markers below a certain threshold. Based on the
repoducibility column found in DArT files.

**Filter target:** Markers.

**Statistics**: Reproducibility (established by DArT)

## Usage

``` r
filter_dart_reproducibility(
  data,
  interactive.filter = TRUE,
  filter.reproducibility = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A tidy genomic data frame or another genomic object supported by the
  calling function.

- interactive.filter:

  Logical indicating whether an interactive filtering session may
  display diagnostics and ask for thresholds. Default:
  `interactive.filter = TRUE`.

- filter.reproducibility:

  (double, character) This is best decided after viewing the figures.
  Usually values higher than 0.95 are not uncommon. The value can also
  be character: `filter.reproducibility = "outliers"`. Using this, will
  remove outlier markers using the lower outlier statistics. Default:
  `filter.reproducibility = NULL`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

A list in the global environment with 6 objects:

1.  \$whitelist.markers

2.  \$blacklist.markers

3.  \$filters.parameters

The object can be isolated in separate object outside the list by
following the example below.

## Details

**Interactive version**

There are 2 steps in the interactive version to visualize and filter the
data based on the reproducibility value:

Step 1. Visualization using a box plot

Step 2. Choose the filtering threshold

## Examples

``` r
if (FALSE) { # \dontrun{
spotted.cod <- genometranslator::read_dart(
    data = "Combined_1514and1614_SNP_80Callrate.csv",
    strata = "strata.dart.spotted.cod.tsv"
)
turtle.filtered <- radr::filter_dart_reproducibility(
    data = spotted.cod,
    filter.reproducibility = 0.97
)
} # }
```
