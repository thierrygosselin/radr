# Filter markers based on genotyping / missing rate

Filter markers based on their genotyping (call) rate, i.e. the
proportion of non-missing genotypes per marker. This is a convenient way
to remove SNPs with too much missing data before downstream analyses.

**Filter target:** Markers.

**Statistics**: marker-level missingness (`MISSING_PROP` from
`generate_stats`).

## Usage

``` r
filter_genotyping(
  data,
  interactive.filter = TRUE,
  filter.genotyping = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A Genomic Data Structure (GDS) file or connection (`SeqVarGDSClass` or
  `gds.file`) created by `read_vcf` or other radr import functions.

- interactive.filter:

  (optional, logical). Do you want the filtering session to be
  interactive ? With `interactive.filter = TRUE`, helper tables and
  plots are shown, and the user is prompted for a threshold. Default:
  `interactive.filter = TRUE`.

- filter.genotyping:

  (optional, character or double). Two modes are available:

  - character string: `filter.genotyping = "outliers"` uses the upper
    outlier value from the missingness boxplot as threshold (i.e.
    removes high-missingness outliers);

  - double: `filter.genotyping = 0.2` allows up to `0.2` missing
    genotypes (i.e. keeps markers with `MISSING_PROP <= 0.2`).

  Default: `filter.genotyping = NULL`.

  If `interactive.filter = FALSE` and `filter.genotyping = NULL`, the
  function returns `data` unchanged.

- filename:

  (optional, character). Basename for some output files (plots, tables).
  If `NULL`, a name is generated automatically from the function and
  date. Default: `filename = NULL`.

- parallel.core:

  (optional, integer). Number of CPU cores to use in helper statistics.
  Default: `parallel.core = parallel::detectCores() - 1`.

- verbose:

  (optional, logical). Show messages and progress. Default:
  `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

The filtered GDS connection with updated `markers.meta`. The underlying
GDS file is updated on disk.

Side-effects:

- helper tables and plots are written under
  `filter_genotyping_YYYYMMDD@HHMM`;

- `markers.meta` inside the GDS is updated and synchronised;

- the `filters.parameters` file is updated to record the filter.

## Details

Internally, the function uses `generate_stats` to compute marker-level
missingness (`MISSING_PROP`). Markers with missing proportion greater
than the chosen threshold are blacklisted and have their `FILTERS`
column updated to `"filter.genotyping"` in `markers.meta`.

With `interactive.filter = TRUE`, a helper table and plot are generated
to assist in choosing a threshold:

- `genotyping.helper.table.tsv`: number of markers kept/removed for
  thresholds from 0 to 1 by 0.1;

- if strata are present, `markers.pop.missing.helper.table.tsv`
  summarises missingness per population;

- a PDF/PNG helper plot summarising these patterns.

This function operates on the GDS representation of the data and is
meant for genotyping/missingness-based pruning *after* VCF import (e.g.
after `read_vcf`), in contrast to VCF-level slimming using tools like
[`filter_monomorphic_vcf`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic_vcf.md).

## Interactive version

The interactive mode first displays and writes the marker-missingness
distribution and helper tables. It then asks
`"Choose the maximum missing proportion allowed:"`, with a value from 0
to 1. Markers with `MISSING_PROP` greater than the selected value are
blacklisted. For example, 0.20 permits at most 20 percent missing
genotypes per marker. Use `interactive.filter = FALSE` and provide
`filter.genotyping` explicitly for a reproducible analysis.

## See also

[`filter_common_markers`](https://thierrygosselin.github.io/radr/reference/filter_common_markers.md),
[`filter_monomorphic`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.md),
[`filter_ma`](https://thierrygosselin.github.io/radr/reference/filter_ma.md),
`read_vcf`,
[`explore_genomes`](https://thierrygosselin.github.io/radr/reference/explore_genomes.md).

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
genome <- genometranslator::read_genome("populations.snps.vcf")

# Inspect marker missingness and choose the threshold interactively.
genome <- radr::filter_genotyping(data = genome)

# Alternatively, start from a separate unfiltered GDS for a scripted run.
scripted_genome <- genometranslator::read_genome("populations_scripted.gds")
scripted_genome <- radr::filter_genotyping(
  data = scripted_genome,
  interactive.filter = FALSE,
  filter.genotyping = 0.20
)
} # }
```
