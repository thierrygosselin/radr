# Filter markers mean coverage

This function is designed to remove/blacklist markers based on mean
coverage information.

**Filter target:** Markers.

**Statistics**: mean coverage ( The read depth of individual genotype is
averaged across markers).

## Usage

``` r
filter_coverage(
  data,
  interactive.filter = TRUE,
  filter.coverage = NULL,
  filename = NULL,
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

- filter.coverage:

  (optional, string) 2 options:

  - character string `filter.coverage = "outliers"` will use as
    thresholds the lower and higher outlier values in the box plot.

  - integers string `filter.coverage = c(10, 200)`. For the marker's
    mean coverage lower and upper bound.

  Default: `filter.coverage = NULL`.

- filename:

  (optional, character) Default: `filename = NULL`.

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

With `interactive.filter = FALSE`, a list in the global environment,
with 7 objects:

1.  \$tidy.filtered.mac

2.  \$whitelist.markers

3.  \$blacklist.markers

4.  \$mac.data

5.  \$filters.parameters

With `interactive.filter = TRUE`, a list with 4 additionnal objects are
generated.

1.  \$distribution.mac.global

2.  \$distribution.mac.local

3.  \$mac.global.summary

4.  \$mac.helper.table

## Advance mode

*dots-dots-dots ...* allows to pass several arguments for fine-tuning
the function:

1.  `filter.common.markers` (optional, logical). Default:
    `filter.common.markers = FALSE`, Documented in
    [`filter_common_markers`](https://thierrygosselin.github.io/radr/reference/filter_common_markers.md).

2.  `filter.monomorphic` (logical, optional) Should the monomorphic
    markers present in the dataset be filtered out ? Default:
    `filter.monomorphic = TRUE`. Documented in
    [`filter_monomorphic`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.md).

3.  `path.folder`: to write ouput in a specific path (used internally in
    radr). Default: `path.folder = getwd()`. If the supplied directory
    doesn't exist, it's created.

## Interactive version

To help choose a threshold for the local and global MAF use the
interactive version.

2 steps in the interactive version:

Step 1. Visualization and helper table.

Step 2. Filtering markers based on mean coverage

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# The minumum
} # }
```
