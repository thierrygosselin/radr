# Summarise genomic context around genome-scan signals

A regional peak is easier to interpret when the analyst can ask whether
it overlaps low recombination, unusual marker density, increased
missingness, a sequencing or mapping batch, a candidate inversion, or a
signal supported by a method based on a different genomic signature.

Use this function after basic quality control but before interpreting
genome scans. When candidate inversion regions exist, compare at least
three views: the complete genome, a collinear sensitivity analysis
excluding candidate regions, and an inversion-specific analysis.
Candidate regions are annotated but never silently excluded.

## Usage

``` r
genome_scan_context(
  data,
  window.snps = 250L,
  step.snps = window.snps,
  window.bp = NULL,
  step.bp = window.bp,
  inversion.regions = NULL,
  scan.statistics = NULL,
  ld.max.snps = 250L,
  filename = "genome_scan_context",
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A GDS filepath or open `SeqVarGDSClass` object.

- window.snps:

  Number of SNPs per fixed-SNP window. Default: `window.snps = 250`.

- step.snps:

  Number of SNPs between window starts. Default:
  `step.snps = window.snps`.

- window.bp:

  Optional physical window size in base pairs. Supplying it replaces
  fixed-SNP windows. Default: `window.bp = NULL`.

- step.bp:

  Physical distance between window starts. Default:
  `step.bp = window.bp`.

- inversion.regions:

  Optional table, or the result from
  [`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md).
  A table requires `chromosome`, `start`, and `end`; `candidate_id` and
  `candidate_class` are retained when present. Default:
  `inversion.regions = NULL`.

- scan.statistics:

  Optional marker-level table with chromosome, position, and one or more
  scan statistics. Column names are matched without regard to case.
  Numeric statistics are summarised by window mean, median, and maximum.
  Default: `scan.statistics = NULL`.

- ld.max.snps:

  Maximum evenly spaced SNPs used for LD in each window. Default:
  `ld.max.snps = 250`.

- filename:

  Output filename stem. Default: `filename = "genome_scan_context"`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

- ...:

  Standard `radr` workflow arguments, including `path.folder`.

## Value

A list containing the window context table, normalised region table,
plots, data-source information, and output paths.

## Details

Build a window-level context table for interpreting regional genome-scan
peaks. The table places marker density, call rate, heterozygosity, minor
allele frequency, depth when available, local LD, candidate inversion or
structural-region annotations, and user-supplied scan statistics beside
one another. It does not decide whether a region is under selection.

The context table is descriptive. Overlap with an inversion-like
haploblock, centromere, low-recombination region, assembly gap, or
technical anomaly changes the interpretation of a peak but does not by
itself validate or invalidate selection. Repeat analyses after excluding
putative heterokaryotypes when arrangement calls are available, and seek
support from methods based on different summaries.

## References

Booker TR, Yeaman S, Whitlock MC (2020). Variation in recombination rate
affects detection of outliers in genome scans under neutrality.
Molecular Ecology, 29, 4274-4279.

Faria R, Johannesson K, Butlin RK, Westram AM (2019). Evolving
inversions. Trends in Ecology & Evolution, 34, 239-248.
