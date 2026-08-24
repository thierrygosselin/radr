# Filter markers mean coverage

This function is designed to remove/blacklist markers based on mean
coverage information.

**Filter target:** Markers.

**Statistics**: Mean marker coverage. Genotype read depth is averaged
across active individuals for each marker.

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

  - numeric vector `filter.coverage = c(10, 200)` for the marker
    mean-coverage lower and upper bounds.

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

The filtered data in the same representation as the input. GDS marker
metadata and active variants are updated in place, so the underlying GDS
file is modified. Coverage tables, figures, marker lists, and filtering
parameters are written to the function output folder when applicable.

## Advanced mode

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

The interactive mode first calculates marker coverage, writes and
displays the coverage distribution and helper plots, and then asks:

1.  `"Choose the min mean coverage threshold (e.g. 7 or 10):"`

2.  `"Choose the max mean coverage threshold (e.g. 100 or 300):"`

Markers outside the inclusive interval are blacklisted. Use
`interactive.filter = FALSE` with an explicit two-value
`filter.coverage` for a reproducible analysis.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
genome <- genometranslator::read_genome("my_genome.gds")

# Inspect marker coverage and choose lower and upper limits interactively.
genome <- radr::filter_coverage(data = genome)

# Alternatively, start from a separate unfiltered GDS for a scripted run.
scripted_genome <- genometranslator::read_genome("my_genome_scripted.gds")
scripted_genome <- radr::filter_coverage(
  data = scripted_genome,
  interactive.filter = FALSE,
  filter.coverage = c(10, 200)
)
} # }
```
