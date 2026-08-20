# Extract the reference genome filename/path used for a dataset

Retrieve the reference genome filename or path associated with a
dataset, when available.

This helper is intended for internal or downstream use. It does *not*
attempt to determine whether the dataset is reference-guided (see
[`detect_ref_genome`](https://thierrygosselin.github.io/radr/reference/detect_ref_genome.md)
for that). Instead, it simply attempts to recover the reference
information from:

**1. VCF files** Parsed via `check_header_source_vcf()`, reading the
`##reference=` header line when present.

**2. SeqArray / radr GDS files** Extracted using
`SeqArray::seqSummary(gds, "$reference")`. If this header field exists,
it is returned directly.

If no reference path is detected, the function returns `NULL` quietly.

## Usage

``` r
extract_ref_genome(data = NULL, verbose = FALSE)
```

## Arguments

- data:

  A VCF file path, a radr GDS object, a SeqArray GDS object, or a path
  to a GDS file. Default: `data = NULL`.

- verbose:

  (logical) When `TRUE`, the function prints optional diagnostics
  explaining why no reference path could be found. Default:
  `verbose = FALSE`.

## Value

A character scalar containing the reference genome filename or path, or
`NULL` if no reference could be detected.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# From a VCF:
ref.path <- radr::extract_ref_genome("variants.vcf.gz")

# From a radr GDS:
gds <- genometranslator::read_genome("my_data.gds")
ref.path <- radr::extract_ref_genome(gds)
} # }
```
