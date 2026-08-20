# Estimate population-specific beta

Estimates population-specific \\\beta_i\\ following Weir and Goudet
(2017). The implementation uses diploid, biallelic alternate-allele
dosage (`ALT_DOSAGE`: 0, 1, 2, or `NA`) and works directly from a GDS
file, an open `SeqVarGDSClass` object, or a tidy genotype table.

`beta_estimator()` does not filter individuals or markers. Missing
genotypes are omitted within each marker-population combination. A
marker contributes only when it has observed genotypes in at least two
populations and its between-population diversity is finite.

## Usage

``` r
beta_estimator(
  data,
  strata = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```

## Arguments

- data:

  A GDS filepath, an open `SeqVarGDSClass` object, or a tidy data frame
  containing `MARKERS`, `INDIVIDUALS`, `STRATA` (or `POP_ID`), and
  `ALT_DOSAGE`.

- strata:

  Optional strata filepath or data frame used to add or replace
  population assignments. It must contain `INDIVIDUALS` and `STRATA`.
  Default: `strata = NULL`.

- filename:

  Optional output prefix. When supplied, three tab-delimited files are
  written with suffixes `_beta.tsv`, `_within_population.tsv`, and
  `_between_populations.tsv`. Default: `filename = NULL`.

- parallel.core:

  Number of processor cores passed to
  [`genometranslator::read_genome()`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
  when file input must be imported. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress and a population beta summary. Default:
  `verbose = TRUE`.

## Value

A named list containing:

- beta:

  One row per population with `BETA`, the number of contributing
  markers, and the summed within- and between-population diversities.

- within_population:

  Marker-population allele counts, frequencies, and `HW`.

- between_populations:

  Marker-level number of represented populations and `HB`.

## Details

For marker \\l\\ and population \\i\\, within-population gene diversity
is estimated as \$\$H\_{W,li} = \frac{n\_{li}}{n\_{li}-1} \left(1 -
p\_{li}^2 - (1-p\_{li})^2\right),\$\$ where \\n\_{li}\\ is the number of
observed gene copies and \\p\_{li}\\ is the alternate-allele frequency.

Between-population diversity is calculated for every marker from the
populations with observed genotypes. Population-specific beta is then
\$\$\beta_i = 1 - \frac{\sum_l H\_{W,li}}{\sum_l H\_{B,l}}.\$\$
Consequently, populations can be compared using different numbers of
loci when their missingness patterns differ. Review `N_MARKERS` in the
returned summary before interpreting differences among populations.

## References

Weir, B. S. and Goudet, J. (2017). A unified characterization of
population structure and relatedness. *Genetics*, 206, 2085–2103.
[doi:10.1534/genetics.116.198424](https://doi.org/10.1534/genetics.116.198424)

Goudet, J., Kay, T. and Weir, B. S. (2018). How to estimate kinship.
*Molecular Ecology*, 27, 4121–4135.
[doi:10.1111/mec.14833](https://doi.org/10.1111/mec.14833)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
genome <- genometranslator::read_genome(
  data = "individuals.vcf.gz",
  strata = "strata.tsv"
)

beta <- beta_estimator(genome)
beta$beta
beta$within_population
beta$between_populations
} # }
```
