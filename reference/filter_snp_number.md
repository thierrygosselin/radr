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

The filtered data in the same representation as the input. GDS marker
metadata and active variants are updated in place. Diagnostic files,
marker lists, and filtering parameters are written to the output folder.

## Interactive version

The function first displays and writes the distribution of SNPs per
locus and the effect of candidate thresholds. It then asks:

1.  `"Do you still want to blacklist markers? (y/n):"`

2.  If yes, choose `1` to use the upper boxplot-outlier statistic or `2`
    to enter a threshold.

3.  With option 2, answer
    `"Enter the maximum number of SNP per locus allowed:"`.

All SNPs in loci containing more than the selected number are
blacklisted. Answering no leaves the data unchanged. Use
`interactive.filter = FALSE` and provide `filter.snp.number` explicitly
for a reproducible analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
genome <- genometranslator::read_genome(
  data = "turtle.vcf",
  strata = "turtle.strata.tsv"
)

# Inspect the SNP-per-locus distribution interactively.
genome <- radr::filter_snp_number(data = genome)

# Alternatively, use a separate unfiltered GDS for a scripted run that
# retains loci containing at most four SNPs.
scripted_genome <- genometranslator::read_genome("turtle_scripted.gds")
scripted_genome <- radr::filter_snp_number(
  data = scripted_genome,
  interactive.filter = FALSE,
  filter.snp.number = 4
)
} # }
```
