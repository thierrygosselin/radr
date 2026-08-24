# MAC, MAF and MAD filter

Remove or blacklist markers using the global minor allele count (MAC),
minor allele frequency (MAF), or minor allele depth (MAD). These
statistics can help remove weakly supported variants, sequencing or
genotyping noise, and markers with very low polymorphism.

Filtering rare alleles also changes the allele-frequency spectrum. This
can bias analyses that depend on that spectrum, including population
structure, demographic inference, differentiation, and assignment.
Select a threshold that matches the intended analysis, and report it
explicitly.

**Filter target:** Markers.

**Statistics**: Minor allele count, frequency, or read depth per marker.
The statistic is computed globally across the active individuals. Strata
or population labels are not used.

## Usage

``` r
filter_ma(
  data,
  interactive.filter = TRUE,
  ma.stats = "mac",
  filter.ma = NULL,
  calibrate.alleles = "depth",
  keep.biallelic = TRUE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A tidy genomic data frame or another genomic object supported by the
  calling function.

- interactive.filter:

  Logical indicating whether an interactive filtering session may
  display diagnostics and ask for thresholds. Default:
  `interactive.filter = TRUE`.

- ma.stats:

  (optional, character) The statistic of the alternate (minor) allele to
  keep a SNP (see details). Options are: `"mac", "maf", "mad"`. Default:
  `ma.stats = "mac"`.

- filter.ma:

  (optional, integer or double) The threshold to blacklist a marker (see
  Details). For example, `ma.stats = "mac"` and `filter.ma = 3`
  blacklist markers with a minor allele count strictly below 3. The
  threshold is applied globally and independently to every SNP. Default:
  `filter.ma = NULL`.

- calibrate.alleles:

  (character) When ancestral allele information is not available, REF /
  ALT alleles can be calibrated using alleles count or depth (sequencing
  coverage for REF and ALT alleles). Removing individuals will impact
  REF/ALT alleles. Ideally, it should be based on observed read depth,
  when this information is available for both alleles. By selecting
  `calibrate.alleles = "depth"`, radr will check if the information is
  available, if not, it will use `calibrate.alleles = "count"`. Using
  `calibrate.alleles = "ancestral"` will turn off calibration and use
  the information preset available in the data. Default:
  `calibrate.alleles = "depth"`.

- keep.biallelic:

  (logical) Keep only biallelic variants (REF + 1 ALT) when computing
  allele counts from GDS. Default: `keep.biallelic = TRUE`.

- filename:

  (optional, character) Output filename prefix when a tidy tibble is
  supplied. The filtered table is written as an Arrow Parquet file. GDS
  input is updated through its marker metadata and is not duplicated by
  this argument. Default: `filename = NULL`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

The filtered data in the same representation as the input. For GDS
input, the function updates the `FILTERS` field in the marker metadata,
synchronizes the active variants, and returns the GDS connection. The
GDS file on disk is therefore modified. For tidy input, it returns a
filtered tibble. Diagnostic tables, figures, blacklists, whitelists, and
filtering parameters are written to the function output folder when
applicable.

## Note

Thanks to Charles Perrier and Jeremy Gaudin for very useful comments on
previous version of this function.

## Why is the statistic computed globally?

Applying a separate threshold within each stratum makes marker retention
depend on stratum sample size, missingness, and allele-frequency
differences. It can preferentially remove variants that are rare in one
population but common in another, erase locally informative alleles, or
enrich the retained data for particular patterns of differentiation.
These effects can bias downstream comparisons among populations.

For this reason, `filter_ma()` applies one criterion to each marker
using all active individuals and without consulting population labels.
This does not make every threshold unbiased. Global filtering still
changes the allele-frequency spectrum, and a globally retained allele
may occur mostly in one stratum. The appropriate threshold therefore
depends on the data and the intended analysis.

## MAC, MAF or MAD?

**MAC** is the preferred general-purpose choice in radr. A fixed count
threshold asks for the same minimum number of observed minor alleles at
every marker. In contrast, a fixed MAF threshold corresponds to
different allele counts when the number of successfully genotyped
individuals varies among markers.

This is especially relevant for RADseq data. Restriction-site
polymorphism, variation in read depth, allele dropout, library effects,
and genotype calling can produce marker-specific missingness. Because
MAF uses the number of observed chromosomes as its denominator, two
markers with the same minor allele count can fall on opposite sides of a
MAF threshold solely because their missing-data rates differ.

Consider 36 diploid individuals and three SNPs with the same minor
allele count but different numbers of genotyped individuals:

- `SNP: genotyped individuals: major/minor allele counts: MAF`

- SNP1: 36: 69/3: 0.0417

- SNP2: 30: 57/3: 0.0500

- SNP3: 24: 45/3: 0.0625

Each marker has three copies of the minor allele. With
`filter.ma = 0.05`, SNP1 is removed because its MAF is below 0.05,
whereas SNP2 and SNP3 are retained. With a MAC threshold of 3, all three
are treated identically and retained. MAC does not preserve the
allele-frequency spectrum, but it avoids changing the effective count
criterion from one marker to another as missingness changes.

**MAF** remains useful when a frequency threshold is scientifically
required or when marker call rates are sufficiently uniform. Its
dependence on the observed denominator should be considered and
reported.

**MAD** uses the total read support for the less abundant allele. It can
identify variants that pass a count threshold but whose minor allele is
supported by very little sequence coverage. For example, two markers may
both have MAC = 4. If the minor allele has total depths of 4 and 32
reads, respectively, a MAD threshold of 5 removes the first marker and
retains the second. MAD is therefore a read-support filter, not a
replacement for MAC or MAF. It is available only when allele-depth
information is present, and its interpretation depends on library depth
and genotype-calling properties.

## Advanced mode

*dots-dots-dots ...* allows to pass several arguments for fine-tuning
the function:

1.  `filter.common.markers` (optional, logical). Default:
    `filter.common.markers = FALSE`, Documented in
    [`filter_common_markers`](https://thierrygosselin.github.io/radr/reference/filter_common_markers.md).

2.  `filter.monomorphic` (logical, optional) Should the monomorphic
    markers present in the dataset be filtered out ? Default:
    `filter.monomorphic = TRUE`. Documented in
    [`filter_monomorphic`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.md).

3.  `path.folder`: to write ouput in a specific path (used internally in
    radr). Default: `path.folder = getwd()`. If the supplied directory
    doesn't exist, it's created.

## Interactive version

The interactive mode exposes the decisions made by the filter before it
modifies marker status. It has three steps:

1.  Calculate the available global minor-allele statistics, write the
    helper tables and figures, and display the distributions and
    consequences of candidate thresholds.

2.  Ask `"Choose the statistic:"`. The available answers are drawn from
    `"mac"`, `"maf"`, and `"mad"`; MAD is offered only when allele-depth
    information is available.

3.  Ask `"Choose the <statistic> threshold:"`, where `<statistic>` is
    the choice made in step 2. Markers with a value strictly below the
    selected threshold are blacklisted.

Use `interactive.filter = FALSE` for a reproducible scripted analysis,
and supply both `ma.stats` and `filter.ma` explicitly.

## References

Linck, E., & Battey, C. J. (2019). Minor allele frequency thresholds
strongly affect population structure inference with genomic data sets.
*Molecular Ecology Resources*, 19(3), 639-647.
[doi:10.1111/1755-0998.12995](https://doi.org/10.1111/1755-0998.12995)

Roesti, M., Salzburger, W., & Berner, D. (2012). Uninformative
polymorphisms bias genome scans for signatures of selection. *BMC
Evolutionary Biology*, 12, 94.
[doi:10.1186/1471-2148-12-94](https://doi.org/10.1186/1471-2148-12-94)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Open a GDS created by genometranslator.
genome <- genometranslator::read_genome("my_genome.gds")

# Explore the distributions and choose the statistic and threshold
# interactively. The questions and available diagnostics are described in
# the Interactive version section.
genome <- radr::filter_ma(data = genome)

# Reproducible non-interactive MAC filtering. A threshold of 4 retains a
# marker only when at least four copies of its minor allele are observed.
# Four copies could be carried by four heterozygous individuals, two
# minor-allele homozygotes, or one minor-allele homozygote plus two
# heterozygotes. Consequently, at least two diploid individuals must carry
# the minor allele.

genome <- radr::filter_ma(
  data = genome,
  interactive.filter = FALSE,
  ma.stats = "mac",
  filter.ma = 4
)

# Alternative starting from a separate, unfiltered copy of the source GDS:
# use MAF when a proportional threshold is appropriate for the analysis.
# Here, markers with global MAF below 0.01 are blacklisted.
maf_genome <- genometranslator::read_genome("my_genome_for_maf.gds")
maf_genome <- radr::filter_ma(
  data = maf_genome,
  interactive.filter = FALSE,
  ma.stats = "maf",
  filter.ma = 0.01
)

# MAD requires allele-depth information. This example retains markers whose
# minor allele is supported by at least five reads across active samples.
mad_genome <- genometranslator::read_genome("my_genome_for_mad.gds")
mad_genome <- radr::filter_ma(
  data = mad_genome,
  interactive.filter = FALSE,
  ma.stats = "mad",
  filter.ma = 5
)
} # }
```
