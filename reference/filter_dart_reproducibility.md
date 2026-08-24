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

The filtered data in the same representation as the input. GDS marker
metadata and active variants are updated in place. Diagnostic files,
marker lists, and filtering parameters are written to the output folder.

## Interactive version

The function first displays and writes the DArT reproducibility
distribution and helper outputs. It then asks the following questions:

1.  `"Do you still want to blacklist markers? (y/n):"`

2.  If yes, choose `1` to use the lower boxplot-outlier statistic or `2`
    to enter a threshold.

3.  With option 2, answer
    `"Enter the proportion threshold (0-1), the minimum reproducibility tolerated:"`.

Markers with reproducibility below the selected threshold are
blacklisted. Answering no leaves the data unchanged. Use
`interactive.filter = FALSE` with an explicit `filter.reproducibility`
for a reproducible analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
spotted.cod <- genometranslator::read_dart(
    data = "Combined_1514and1614_SNP_80Callrate.csv",
    strata = "strata.dart.spotted.cod.tsv"
)
turtle.filtered <- radr::filter_dart_reproducibility(
    data = spotted.cod,
    interactive.filter = FALSE,
    filter.reproducibility = 0.97
)
} # }
```
