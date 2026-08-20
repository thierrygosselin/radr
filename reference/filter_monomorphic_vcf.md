# Filter monomorphic SNPs in a VCF using bcftools (AC/AN-based)

Wrapper around `bcftools view` and `bcftools query` to:

- keep only polymorphic markers in a VCF, and

- generate a blacklist of monomorphic sites.

A site is considered polymorphic in the sample if at least one ALT
allele has an allele count strictly between 0 and the total allele count
(`INFO/AN`). In bcftools expression syntax:


    INFO/AC>0 && INFO/AC<INFO/AN

This works for multi-ALT sites: as soon as one ALT allele is segregating
(some REF and some ALT\\\_i\\), the variant is kept. Sites fixed for REF
or fixed for a single ALT allele are dropped as monomorphic.

The function is intended for use on freshly created VCFs where `INFO/AC`
and `INFO/AN` are trustworthy (e.g. direct output from FreeBayes or
bcftools). It does not inspect genotypes or GDS, only the VCF INFO
fields.

## Usage

``` r
filter_monomorphic_vcf(
  vcf,
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

- out.vcf:

  (optional, character) Output VCF filename for polymorphic markers. If
  `NULL`, the function generates a filename from the input VCF root by
  appending `".polymorphic.vcf"` and `".gz"` when `compress = TRUE`.
  Default: `out.vcf = NULL`.

- blacklist.file:

  (optional, character) Path to a TSV file storing monomorphic markers
  (`CHROM`, `POS`, `ID`). If `NULL`, the function generates a filename
  from the VCF root by appending `".monomorphic_blacklist.tsv"`.
  Default: `blacklist.file = NULL`.

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

- `polymorphic.vcf` – path to the filtered VCF;

- `blacklist` – path to the monomorphic marker list;

- `log` – path to the bcftools log file.

## Details

**Important distinction — allele-level vs genotype-level polymorphism**

`filter_monomorphic_vcf` identifies polymorphic sites solely from VCF
INFO fields (`AC` and `AN`), using the bcftools expression:


    INFO/AC>0 && INFO/AC<INFO/AN

A variant is kept if **at least one ALT allele** has a non-zero allele
count and is not fixed in the sample. This is an **allele-level
definition of polymorphism**:

If some chromosomes carry REF and some carry ALT\\\_i\\, the variant is
considered polymorphic—even if all individuals share the same genotype
(e.g., all REF/ALT heterozygotes).

This pre-filtering step is intentionally lightweight and designed for
**VCF-level slimming before merging for example scaffolds**. It does not
inspect genotypes and does not attempt to detect genotype-level
monomorphism.

As a result, the set of monomorphic markers detected here may differ
from those detected later inside a GDS by
[`filter_monomorphic`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.md).

## Note

**Why results differ between `filter_monomorphic_vcf` and
[`filter_monomorphic`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.md)**

These functions intentionally implement two different biological
definitions:

- `filter_monomorphic_vcf` removes sites that are
  **allele-monomorphic**, based on INFO/AC and INFO/AN.

- [`filter_monomorphic`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.md)
  removes sites that are **genotype-monomorphic**, based on the
  distribution of genotype phenotypes in the GDS.

A variant can contain both REF and ALT alleles (allele-level
polymorphism) but still have only one genotype phenotype across all
individuals. Therefore, it is expected and correct that
[`filter_monomorphic`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.md)
may remove additional markers.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
filter_monomorphic_vcf(
  vcf            = "chr1_freebayes.vcf.gz",
  out.vcf        = "chr1_freebayes_poly.vcf.gz",
  blacklist.file = "chr1_freebayes_mono.tsv",
  bcftools.path  = "bcftools"
)
} # }
```
