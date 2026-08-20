# Detect whether a dataset is reference-guided or de novo assembled

Determine if genomic data were produced using a **reference genome**
(reference-assisted assembly) or using a **de novo assembly**.

The function integrates multiple sources of evidence depending on the
input provided:

**1. Radr GDS (preferred)** If a radr-generated GDS is supplied, the
function checks, in order:

- the `/radr/reference.genome` node (logical value);

- the `/radr/markers.meta/CHROM` field;

- the `SeqArray` `chromosome` node.

**2. VCF imported with SeqArray or read by radr** When a VCF has been
converted to GDS via SeqArray, or when a VCF file is provided directly,
the function may detect reference-guided assembly from:

- the VCF header’s `##reference=` field;

- the presence of `##contig=<ID=...>` definitions;

- the structure of chromosome labels.

**3. Data-source-specific heuristics** Additional rules are applied for
specific pipelines:

- **Stacks**: optionally, detection of `.` or `+` in `annotation/id`
  (see internal comments in the function body);

- **GATK**: chromosome names containing `"contig"`;

- **FreeBayes**: inspection of reference and contig metadata in the VCF
  header.

## Usage

``` r
detect_ref_genome(data = NULL, verbose = TRUE)
```

## Arguments

- data:

  A radr GDS object, a SeqArray GDS object, or a path to a file. If
  provided, the function extracts information from radr metadata,
  SeqArray header fields, `markers.meta`, and other pipeline-specific
  indicators. Default: `data = NULL`.

- verbose:

  (logical, optional) When `TRUE`, the function prints messages
  describing the decision process and final classification. Default:
  `verbose = TRUE`.

## Value

A single logical value:

- `TRUE` – reference-assisted assembly;

- `FALSE` – de novo assembly.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Using a radr GDS file
gds <- genometranslator::read_genome("my_data.gds")
ref.genome <- radr::detect_ref_genome(data = gds)

# Using a VCF file directly
ref.genome <- radr::detect_ref_genome(data = "variants.vcf.gz")
} # }
```
