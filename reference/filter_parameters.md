# Track radr filtering parameters

A package-facing wrapper around
[`genometranslator::genome_parameters()`](https://thierrygosselin.github.io/genometranslator/reference/genome_parameters.html)
that maintains the common genomic operation history during radr
filtering. For GDS input, every completed operation also writes an
identifier-level audit beside the usual parameter file: complete marker
and individual tables for records newly removed and records still
active, plus `filter_audit_manifest.tsv`. Consequently this audit also
covers state-changing functions whose names do not start with `filter_`,
including
[`detect_duplicate_genomes()`](https://thierrygosselin.github.io/radr/reference/detect_duplicate_genomes.md)
and
[`detect_mixed_genomes()`](https://thierrygosselin.github.io/radr/reference/detect_mixed_genomes.md).

## Usage

``` r
filter_parameters(...)
```

## Arguments

- ...:

  Arguments passed to
  [`genometranslator::genome_parameters()`](https://thierrygosselin.github.io/genometranslator/reference/genome_parameters.html).
