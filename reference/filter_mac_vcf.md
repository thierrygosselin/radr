# Filter low-MAC variants in a VCF using bcftools (AC-based)

Wrapper around `bcftools view` and `bcftools query` to:

- keep only variants with at least one ALT allele having a minor-allele
  count (MAC) above a user-defined threshold, and

- generate a blacklist of low-MAC variants.

The filter is implemented using the `INFO/AC` field written by the
variant caller. For multi-ALT sites, the bcftools expression
`INFO/AC >= mac.threshold` evaluates to `TRUE` if any ALT allele has
`AC >= mac.threshold`. Variants that do not satisfy this condition are
considered low-MAC and are blacklisted.

The function assumes that `INFO/AC` (and `INFO/AN` if present) come from
a recent calling step and are still valid (i.e. no sample subsetting has
occurred since calling). It does not recalculate AC/AN; it only uses the
INFO fields already in the VCF.

This is a lightweight, VCF-level pre-filter that is particularly useful
to slim down large, noisy VCFs (e.g. per-scaffold outputs) before
merging or importing into a GDS. For GDS-level MAC filtering based on
genotypes, see
[`filter_ma`](https://thierrygosselin.github.io/radr/reference/filter_ma.md).

## Usage

``` r
filter_mac_vcf(
  vcf,
  mac.threshold = 4L,
  out.vcf = NULL,
  blacklist.file = NULL,
  bcftools.path = "bcftools",
  compress = TRUE,
  index = TRUE,
  verbose = TRUE
)
```

## Arguments

- vcf:

  (character) Path to the input VCF (`.vcf` or `.vcf.gz`).

- mac.threshold:

  (integer) Minimum ALT allele count required to keep a variant.
  Variants where all ALT alleles have `AC < mac.threshold` are removed.
  Default: `mac.threshold = 4L`.

- out.vcf:

  (optional, character) Output VCF filename for MAC-filtered variants.
  If `NULL`, the function generates a filename from the input VCF root
  by appending `".mac.vcf"` and `".gz"` when `compress = TRUE`. Default:
  `out.vcf = NULL`.

- blacklist.file:

  (optional, character) Path to a TSV file storing low-MAC variants
  (`CHROM`, `POS`, `ID`). If `NULL`, the function generates a filename
  from the VCF root by appending `".mac_blacklist.tsv"`. Default:
  `blacklist.file = NULL`.

- bcftools.path:

  (character) Path or name of the `bcftools` executable (e.g.
  `"bcftools"` or a full path inside a conda environment). Default:
  `bcftools.path = "bcftools"`.

- compress:

  (logical) Compress the output VCF using bgzip (`-Oz`). Default:
  `compress = TRUE`.

- index:

  (logical) Index the output VCF with tabix when `compress = TRUE`.
  Default: `index = TRUE`.

- verbose:

  (logical) Show progress messages and bcftools command lines. Default:
  `verbose = TRUE`.

## Value

Invisibly returns a list with:

- `mac.vcf` – path to the MAC-filtered VCF;

- `blacklist` – path to the low-MAC variant list;

- `log` – path to the bcftools log file.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
filter_mac_vcf(
  vcf            = "chr1_freebayes_poly.vcf.gz",
  mac.threshold  = 4L,
  out.vcf        = NULL,   # will append ".mac.vcf.gz" by default
  blacklist.file = NULL,
  bcftools.path  = "bcftools"
)
} # }
```
