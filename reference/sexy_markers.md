# Screen a GDS for candidate sex-linked markers

`sexy_markers()` compares marker presence, heterozygosity, and, when
available, normalized read depth between known females and males. It is
a read-only screen: active GDS sample and variant selections are
respected and restored, and no filter is applied or written to the GDS.

## Usage

``` r
sexy_markers(
  data,
  strata = NULL,
  sex.column = "STRATA",
  presence.threshold = 0.9,
  absence.threshold = 0.1,
  min.heterozygosity.difference = 0.2,
  coverage.ratio.threshold = 1.5,
  coverage.threshold = 1,
  min.samples.per.sex = 5L,
  fdr.threshold = 0.05,
  require.significance = FALSE,
  chunk.size = 2000L,
  save.plots = TRUE,
  plot.formats = c("png", "pdf"),
  folder.name = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A GDS filepath or an open `SeqVarGDSClass` object.

- strata:

  A data frame or tab-delimited file containing `INDIVIDUALS` and
  `sex.column`. If `NULL`, individual metadata are read from the GDS.

- sex.column:

  Metadata column containing sex. `F`, `female`, `M`, and `male` are
  recognized without regard to case; other values are unknown.

- presence.threshold:

  Minimum presence in the expected sex for Y/W.

- absence.threshold:

  Maximum presence in the other sex for Y/W.

- min.heterozygosity.difference:

  Minimum absolute female-minus-male heterozygosity difference for X/Z
  candidates.

- coverage.ratio.threshold:

  Minimum larger-to-smaller normalized depth ratio for X/Z coverage
  candidates.

- coverage.threshold:

  Minimum read depth considered present. Without depth data, a
  non-missing genotype is considered present.

- min.samples.per.sex:

  Minimum number of known females and males.

- fdr.threshold:

  Maximum FDR when significance is required.

- require.significance:

  Require effect size and method-specific FDR.

- chunk.size:

  Number of variants read from the GDS at a time.

- save.plots:

  Save diagnostic plots in the result folder.

- plot.formats:

  One or both of `"png"` and `"pdf"`.

- folder.name:

  Optional result-folder stem. Default: `sexy_markers`.

- verbose:

  Display progress and summary messages.

- ...:

  Common arguments: `path.folder` and `internal`.

## Value

A `sexy_markers` object containing marker statistics, candidates, sample
summaries, a metadata audit, plots, and the result-folder path.

## Details

Y-like markers are present mainly in males, W-like markers mainly in
females, X-like markers have greater heterozygosity or normalized depth
in females, and Z-like markers show the reverse. These labels are
candidates, not proof of chromosomal location. Sex-specific dropout,
population structure, plate or library effects, paralogy, and mapping
bias can produce similar patterns.

## Statistical tests

Presence and heterozygosity use two-sample score tests for proportions.
Normalized depth uses Welch tests on `log2(1 + depth)`. P-values are
adjusted separately by method with Benjamini-Hochberg FDR. Candidate
selection always uses explicit effect-size thresholds;
`require.significance = TRUE` also requires the method-specific FDR to
pass `fdr.threshold`.

## Read-depth normalization

When depth is available, each sample is divided by its mean positive
depth across active markers. This reduces sample-wide sequencing-depth
differences but cannot remove marker-by-batch interactions or
confounding between sex and plate, lane, library, population, or
sampling group.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- radr::sexy_markers(
  data = "study.gds",
  strata = "sample_metadata.tsv",
  sex.column = "SEX"
)
result$candidates
} # }
```
