# Filter monomorphic markers

Filter monomorphic markers. This filter will remove from the dataset
markers with just *one genotype phenotype*:

- genotypes are ALL homozygotes REF/REF (pp)

- genotypes are ALL heterozygotes REF/ALT, ALT/REF (pq or qp)

- genotypes are ALL homozygotes ALT/ALT (qq)

**Filter target:** Markers.

**Statistics**: the number of genotype phenotypes

Used internally in [radr](https://github.com/thierrygosselin/radr) and
might be of interest for users who wants to keep only polymorphic
markers in their dataset.

## Usage

``` r
filter_monomorphic(
  data,
  filter.monomorphic = TRUE,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A tidy genomic data frame or another genomic object supported by the
  calling function.

- filter.monomorphic:

  (optional, logical) Default: `filter.monomorphic = TRUE`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress messages. Default: `verbose = FALSE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

A list with the filtered input, whitelist and blacklist of markers..

## Details

**Important distinction — genotype-level monomorphism**

`filter_monomorphic` evaluates monomorphism based on the actual
genotypes stored inside a GDS file. A marker is considered monomorphic
when all non-missing individuals display the **same genotype
phenotype**.

Internally, this is assessed using the variant-level alternate allele
dosage (`$dosage_alt`):


    length(unique(g[!is.na(g)])) == 1

This captures cases such as:

- all REF/REF (dosage = 0)

- all REF/ALT (dosage = 1)

- all ALT/ALT (dosage = 2)

Even if the allele frequency in the population is neither 0 nor 1 (e.g.,
all individuals are REF/ALT heterozygotes), the variant is still
considered genotype-monomorphic.

This is a **genotype-level** definition of polymorphism, which is more
conservative and more appropriate for downstream population genomics.

Consequently, `filter_monomorphic` will often identify additional
monomorphic markers that were not removed earlier by
[`filter_monomorphic_vcf`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic_vcf.md),
which uses allele-level logic.

**This discrepancy is by design**.

## Note

**Why results differ between
[`filter_monomorphic_vcf`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic_vcf.md)
and `filter_monomorphic`**

These functions intentionally implement two different biological
definitions:

- [`filter_monomorphic_vcf`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic_vcf.md)
  removes sites that are **allele-monomorphic**, based on INFO/AC and
  INFO/AN.

- `filter_monomorphic` removes sites that are **genotype-monomorphic**,
  based on the distribution of genotype phenotypes in the GDS.

A variant can contain both REF and ALT alleles (allele-level
polymorphism) but still have only one genotype phenotype across all
individuals. Therefore, it is expected and correct that
`filter_monomorphic` may remove additional markers.

## See also

[`explore_genomes`](https://thierrygosselin.github.io/radr/reference/explore_genomes.md),
[`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html),
[`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.html).

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
require(SeqArray) # when using gds
mono <- radr::filter_monomorphic(data = "my.radr.gds.rad", verbose = TRUE)
} # }
```
