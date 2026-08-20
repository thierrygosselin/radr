# Filter SNP number per locus/read

This filter removes outlier markers with too many SNP number per
locus/read. The data requires snp and locus information (e.g. from a VCF
file). Having a higher than "normal" SNP number is usually the results
of assembly artifacts or bad assembly parameters. This filter is
population-agnostic, but still requires a strata file if a vcf file is
used as input.

**Filter target:** Markers.

**Statistics**: The number of SNPs per locus.

## Usage

``` r
filter_snp_number(
  data,
  strata = NULL,
  interactive.filter = TRUE,
  filter.snp.number = NULL,
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

- strata:

  Optional strata file or object containing individual and group
  information. Default: `strata = NULL`.

- interactive.filter:

  Logical indicating whether an interactive filtering session may
  display diagnostics and ask for thresholds. Default:
  `interactive.filter = TRUE`.

- filter.snp.number:

  (integer) This is best decided after viewing the figures. If the
  argument is set to 2, locus with 3 and more SNPs will be blacklisted.
  Default: `filter.snp.number = NULL`.

- filename:

  (optional) Name of the filtered tidy data frame file written to the
  working directory (ending with `.tsv`) Default: `filename = NULL`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

A list in the global environment with 6 objects:

1.  \$snp.number.markers

2.  \$number.snp.reads.plot

3.  \$whitelist.markers

4.  \$tidy.filtered.snp.number

5.  \$blacklist.markers

6.  \$filters.parameters

The object can be isolated in separate object outside the list by
following the example below.

## Details

**Interactive version**

There are 2 steps in the interactive version to visualize and filter the
data based on the number of SNP on the read/locus:

Step 1. SNP number per read/locus visualization

Step 2. Choose the filtering thresholds

## Examples

``` r
if (FALSE) { # \dontrun{
turtle.outlier.snp.number <- radr::filter_snp_number(
data = "turtle.vcf",
strata = "turtle.strata.tsv",
max.snp.number = 4,
filename = "tidy.data.turtle.tsv"
)

tidy.data <- turtle.outlier.snp.number$tidy.filtered.snp.number

#Inside the same list, to isolate the markers blacklisted:
blacklist <- turtle.outlier.snp.number$blacklist.markers

} # }
```
