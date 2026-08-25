# radr

`radr` explores, diagnoses, visualizes, and filters individual genomic
data. It works primarily with GDS files and objects created by
[`genometranslator`](https://thierrygosselin.github.io/genometranslator/).

The two packages have deliberately different responsibilities:

- `genometranslator` reads, standardizes, and writes genomic formats;
- `radr` investigates data quality and applies explicit filters.

[`explore_genomes()`](https://thierrygosselin.github.io/radr/reference/explore_genomes.md)
offers a guided first exploration. It is not a universal filtering
recipe: established analyses should use selected `detect_*()` and
`filter_*()` functions in an order justified for the dataset.

## Installation

Starting from a basic R installation, install the CRAN installer and
required Bioconductor foundation first:

``` r

install.packages(c("BiocManager", "remotes"))

BiocManager::install(c(
  "gdsfmt",
  "Rsamtools",
  "SeqArray"
))

remotes::install_github("thierrygosselin/tgbase")
remotes::install_github("thierrygosselin/genometranslator")
remotes::install_github("thierrygosselin/radr")
```

The `Remotes` field in radr’s `DESCRIPTION` records its GitHub
dependencies, but the explicit sequence above makes a clean installation
easier to diagnose.

Check the installation without changing it:

``` r

radr::radr_dependencies()
```

The returned table distinguishes required components from optional
components and states which workflow uses each optional dependency.

### Optional R packages

Install only what is needed for the planned analysis:

``` r

# LD, linkage pruning, and IBS calculations on GDS
BiocManager::install("SNPRelate")

# Tidy-data distances, sex markers, and fast IBM PNG rendering
install.packages(c("amap", "quantreg", "ragg"))
```

Function documentation identifies its additional dependencies. An
optional package is not required merely to install or load radr.

### Optional command-line tools with Conda

Some VCF-level filters use `bcftools`, while
[`run_bayescan()`](https://thierrygosselin.github.io/radr/reference/run_bayescan.md)
uses the BayeScan executable. These are programs, not R packages. A
shared Conda environment can provide both:

``` bash
conda create --name genomics --channel conda-forge --channel bioconda bcftools bayescan=2.1
conda activate genomics
bcftools --version
bayescan --help
```

For an existing environment:

``` bash
conda activate genomics
conda install --channel conda-forge --channel bioconda bcftools bayescan=2.1
```

Start R or RStudio from the activated environment, then verify
visibility:

``` r

Sys.which(c("bcftools", "bayescan"))
radr::radr_dependencies()
radr::check_bayescan()
```

## A minimal workflow

Import and standardize a genomic file with `genometranslator`, then
diagnose and filter the resulting GDS with radr:

``` r

genome <- genometranslator::read_genome(
  data = "individuals.vcf.gz",
  strata = "strata.tsv"
)

# Preserve the original sample and marker order for the first missingness view
ibm <- radr::detect_ibm(
  data = genome,
  filename = "initial_missingness.png"
)

# Guided exploration for a new dataset
screened <- radr::explore_genomes(data = genome)
```

Filtering order should follow what is known about the project rather
than a fixed recipe. Filtering individuals first changes marker
statistics, while filtering markers first changes individual statistics.

The [getting-started
vignette](https://thierrygosselin.github.io/radr/articles/using_radr.html#build-a-tailored-workflow)
develops two contrasting examples: a marker-first workflow for a noisy
callset and a sample-first workflow for known sequencing failures. It
also explains how to return to guided exploration after correcting a
known problem and how to compare alternative filtering orders
reproducibly.

## Citation

``` r

citation("radr")
packageVersion("radr")
```

Until a dedicated publication or DOI is available, cite the version and,
for a development build, record the Git commit and access date:

> Gosselin, T. (2026). *radr: Explore, diagnose and filter genomic
> data*. R package version 0.0.0.9000.
> <https://github.com/thierrygosselin/radr>. Accessed 2026-08-26.

## Website and support

Documentation and articles are available at
<https://thierrygosselin.github.io/radr/>. Report problems or request
features through the [GitHub issue
tracker](https://github.com/thierrygosselin/radr/issues).
