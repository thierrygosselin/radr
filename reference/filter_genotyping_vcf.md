# Filter SNPs in a VCF based on genotyping / missing rate (bcftools)

Wrapper around `bcftools +fill-tags` and `bcftools view/query` to:

- keep only markers with missing proportion below or equal to a
  user-defined threshold, and

- generate a blacklist of markers above that threshold.

The missing proportion per site is taken from the `F_MISSING` INFO tag,
computed by `bcftools +fill-tags`. A variant is kept if:


    F_MISSING <= genotyping.threshold

This makes the function well suited as a lightweight, VCF-level pruning
step (e.g. before merging scaffold-level VCFs), without opening a GDS.

**Filter target:** VCF variants.

## Usage

``` r
filter_genotyping_vcf(
  vcf,
  genotyping.threshold,
  out.vcf = NULL,
  blacklist.file = NULL,
  bcftools.path = "bcftools",
  compress = TRUE,
  index = TRUE,
  keep.filled.vcf = FALSE,
  verbose = TRUE
)
```

## Arguments

- vcf:

  (character) Path to the input VCF (`.vcf` or `.vcf.gz`).

- genotyping.threshold:

  (double) Maximum allowed missing proportion per site (i.e.
  `F_MISSING`). For example, `genotyping.threshold = 0.5` keeps markers
  with at most 50 percent missing genotypes. Must be between 0 and 1.

- out.vcf:

  (optional, character) Output VCF filename for filtered markers. If
  `NULL`, the function generates a filename from the input VCF root by
  appending `".genotyping.vcf"` and `".gz"` when `compress = TRUE`.
  Default: `out.vcf = NULL`.

- blacklist.file:

  (optional, character) Path to a TSV file storing blacklisted markers
  (`CHROM`, `POS`, `ID`). If `NULL`, the function generates a filename
  from the VCF root by appending `".genotyping_blacklist.tsv"`. Default:
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

- keep.filled.vcf:

  (logical) Keep the intermediate VCF with `F_MISSING` filled-in (the
  output of `bcftools +fill-tags`)? If `FALSE`, the function removes
  this file and its index at the end. Default:
  `keep.filled.vcf = FALSE`.

- verbose:

  (logical) Show progress messages and bcftools command lines. Default:
  `verbose = TRUE`.

## Value

Invisibly returns a list with:

- `genotyping.vcf` – path to the filtered VCF;

- `blacklist` – path to the blacklist of high-missing markers;

- `log` – path to the bcftools log file;

- `filled.vcf` – path to the intermediate VCF with `F_MISSING` (may have
  been deleted if `keep.filled.vcf = FALSE`).

## Details

The function operates entirely at the VCF level:

1.  `bcftools +fill-tags` is used to compute `F_MISSING` for each site
    (if not already present);

2.  `bcftools view -i 'F_MISSING<=X'` keeps well-genotyped markers;

3.  `bcftools query -i 'F_MISSING>X'` writes a blacklist of removed
    markers.

This is conceptually similar to
[`filter_genotyping`](https://thierrygosselin.github.io/radr/reference/filter_genotyping.md),
but:

- `filter_genotyping_vcf` works on the raw VCF using `F_MISSING`
  (per-site missing proportion),

- [`filter_genotyping`](https://thierrygosselin.github.io/radr/reference/filter_genotyping.md)
  works on a GDS and uses the full radr statistics and plotting
  machinery.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Keep markers with at most 50% missing genotypes
filter_genotyping_vcf(
  vcf               = "chr1_freebayes.vcf.gz",
  genotyping.threshold = 0.5,
  out.vcf           = "chr1_freebayes_geno.vcf.gz",
  blacklist.file    = "chr1_freebayes_geno_blacklist.tsv",
  bcftools.path     = "bcftools"
)
} # }
```
