# Detect heterozygotes outliers and estimate miscall rate

Explore departure from H-W equilibrium in bi-allelic RADseq data.
Highlight excess of homozygotes present in numeros RADseq studies. The
function estimate the genotyping error rate and heterozygote miscall
rate. The model focus on heterozygotes being incorrectly called as
homozygotes. See details below for more info.

## Usage

``` r
detect_het_outliers(
  data,
  nreps = 2000,
  burn.in = NULL,
  strata = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A GDS filename or open `SeqVarGDSClass` object.

- nreps:

  (integer, optional) The number of MCMC sweeps to do. Default:
  `nreps = 2000`.

- burn.in:

  (integer, optional) The number of MCMC burn-in reps. With default,
  during execution, you will be asked to enter the nuber of burn-in. For
  this, a plot showing the heterozygote miscall rate for all the MCMC
  sweeps will be printed. This plot will help pinpoint the number of
  burn-in. The remaining MCMC sweeps will be used to average the
  heterozygote miscall rate. e.g. of common value `burn.in = 500`.
  Default: `burn.in = NULL`.

- strata:

  Optional strata file or object containing individual and group
  information. Default: `strata = NULL`.

- filename:

  Optional prefix used for generated output files. Default:
  `filename = NULL`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

A folder generated automatically with date and time, the file
`het.summary.tsv` contains the summary statistics. The file
`markers.genotypes.boundaries.pdf` is the plot with boundaries. The
overall genotyping and heterozygotes miscall rate is writen in the file
`overall_error_rate.tsv`. The function also returns a list inside the
global environment with 8 objects:

1.  input the input data, cleaned if filters were used during import.

2.  outlier.summary a list with a tibble and plot of genotypes
    frequencies and boundaries (also written in the folder).

3.  summary.alt.allele a tibble summarizing the number of markers with:

    - no homozygote for the alternate allele (NO_HOM_ALT)

    - no heterozygote genotype (NO_HET)

    - one homozygote for the alternate allele(ONE_HOM_ALT)

    - one heterozygote genotype (ONE_HET)

    - one homozygote for the alternate allele only (ONE_HOM_ALT_ONLY)

    - one heterozygote genotype only (ONE_HET_ONLY)

    - one homozygote for the alternate allele and one heterozygote
      genotype only (ONE_HOM_ALT_ONE_HET_ONLY)

4.  m.nreps A tibble with the heterozygote miscall rate for each MCMC
    replicate

5.  overall.genotyping.error.rate The overall genotyping error rate

6.  overall.m The overall heterozygote miscall rate

7.  simmed_genos The simulated genotypes

The statistics are summarized per population and overall, the grouping
is found in the last column called `POP_ID`.

## Details

**Before using the function:**

1.  Don't use raw RADseq data, this function will work best with
    filtered data

2.  Remove duplicate
    [`detect_duplicate_genomes`](https://thierrygosselin.github.io/radr/reference/detect_duplicate_genomes.md).

3.  Remove mixed samples
    [`detect_mixed_genomes`](https://thierrygosselin.github.io/radr/reference/detect_mixed_genomes.md).

4.  Look at other filters in radr package...

**During import:**

By default the function will keep only polymorphic markers and markers
common between all populations. If you supply a tidy data frame or a
`.rad` file, the function skip all the filters, pop selection, etc. It
will however scan and remove monomorphic markers automatically.

**Keep track of the data:**

Use the argument filename to write the imported (and maybe further
filtered) tidy genomic data set inside the folder. The filename will be
automatically appended `.rad` to it. This file can be used again
directly inside this function and other radr functions. See
`read_genome`.

## Author

Eric Anderson <eric.anderson@noaa.gov> and Thierry Gosselin
<thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
het.prob <- radr::detect_het_outliers(
data = "tuna.vcf", strata = "tuna.strata.tsv", nreps = 2000)
} # }
```
