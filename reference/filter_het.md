# Heterozygosity filter

Observed Heterozygosity based filtering. The filter arguments of
`filter_het` allows you to test rapidly if departure from realistic
expectations of heterozygosity statistics are a problem in downstream
analysis.

1.  Highlight outliers individual's observed heterozygosity for a quick
    diagnostic of mixed samples or poor polymorphism discovery. The
    statistic is also contrasted with missing data to help differentiate
    problems. (also computed in
    [detect_mixed_genomes](https://thierrygosselin.github.io/radr/reference/detect_mixed_genomes.md))

2.  The observed heterozygosity in the dataset: an assembly artefact, a
    genotyping problem, a problem of population groupings or a reliable
    signal of biological polymorphism? Detect assembly artifact or
    genotyping problem (e.g. under/over-splitting loci) by looking at
    marker's observed heterozygosity statistics by population or
    overall.

3.  Use haplotype or snp level statistics. When the haplotype approach
    is selected, consistencies of heterozygosity statistics along the
    read are highlighted.

4.  Interactive approach help by visualizing the data before making a
    decision on thresholds.

## Usage

``` r
filter_het(
  interactive.filter = TRUE,
  data,
  strata = NULL,
  ind.heterozygosity.threshold = NULL,
  het.approach = c("SNP", "overall"),
  het.threshold = 1,
  het.dif.threshold = 1,
  outlier.pop.threshold = 1,
  helper.tables = FALSE,
  coverage.info = FALSE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- interactive.filter:

  (optional, logical) Do you want the filtering session to be
  interactive. With the default, the user is asked to see figures of
  distribution before making decisions for filtering with heterozygosity
  statistics. Default: `interactive.filter = TRUE`.

- data:

  A GDS filename or open `SeqVarGDSClass` object.

- strata:

  Optional strata file or object containing individual and group
  information. Default: `strata = NULL`.

- ind.heterozygosity.threshold:

  (string, double, optional) Blacklist individuals based on observed
  heterozygosity (averaged across markers).

  The string contains 2 thresholds values (min and max). The values are
  proportions (0 to 1), where 0 turns off the min threshold and 1 turns
  off the max threshold. Individuals with mean observed heterozygosity
  higher (\>) or lower (\<) than the thresholds will be blacklisted.

  A value of `NULL` turns off filtering while retaining the plots and
  heterozygosity tables. Default: `ind.heterozygosity.threshold = NULL`.

- het.approach:

  (Character string). First value, `"SNP"` or `"haplotype"`. The
  haplotype approach considers the statistic consistency on the read
  (locus/haplotype). The major difference: the haplotype approach
  results in blacklisting the entire locus/haplotype with all the SNPs
  on the read. With the SNP approach, SNPs are independently analyzed
  and blacklisted. The second value, will use the statistics by
  population `"pop"` or will consider the data overall, as 1 large group
  `"overall"`. Default: `het.approach = c("SNP", "overall")`.

- het.threshold:

  Number Biallelic markers usually max 0.5. But departure from this
  value is common. The higher the proportion threshold, the more relaxed
  the filter is. With default there is no filtering. Default:
  `het.threshold = 1`.

- het.dif.threshold:

  Number (0 - 1). For `het.approach = "haplotype"` only. You can set a
  threshold for the difference in het along your read. Set the number
  your willing to tolerate on the same read/haplotype. e.g. if you have
  2 SNP on a read/haplotype and on as a het of 0.9 and the other 0.1 and
  you set `het.dif.threshold = 0.3`, this markers will be blacklisted.
  You should strive to have similar statistics along short read like
  RAD. The higher the proportion threshold, the more relaxed the filter
  is. With default, there is no filtering and all the range of
  differences are allowed. Default: `het.dif.threshold = 1`.

- outlier.pop.threshold:

  (integer, optional) Useful to incorporate problematic populations
  dragging down polymorphism discovery, but still wanted for analysis.
  Use this threshold to allow variance in the number of populations
  passing the thresholds described above. e.g. with
  `outlier.pop.threshold = 2`, you tolerate a maximum of 2 populations
  failing the `het.threshold` and/or `het.dif.threshold`. Manage outlier
  markers, individuals and populations downstream with blacklists and
  whitelists produced by the function. See details for more information.
  Default: `outlier.pop.threshold = 1`.

- helper.tables:

  (logical) Output tables that show the number of markers blacklisted or
  whitelisted based on a series of automatic thresholds to guide
  decisions. When `interactive.filter == TRUE`, helper tables are
  written to the directory. Default: `helper.tables = FALSE`.

- coverage.info:

  (optional, logical) Use `coverage.info = TRUE`, if you want to
  visualize the relationshio between locus coverage and locus observed
  heterozygosity statistics. Coverage information is required (e.g. in a
  vcf file...). Default: `coverage.info = FALSE`.

- filename:

  (optional) The function uses
  [`write.fst`](http://www.fstpackage.org/reference/write_fst.md), to
  write the tidy data frame in the folder created in the working
  directory. The file extension appended to the `filename` provided is
  `.rad`. Default: `filename = NULL`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

The function returns inside the global environment a list with 15
objects, the objects names are found by using `names(your.list.name)`:

1.  filtered tidy data frame: `$tidy.filtered.het`

2.  whitelist of markers:`$whitelist.markers`

3.  the strata:`$strata`

4.  the filters parameters used:`$filters.parameters`

5.  the individual's heterozigosity:`$individual.heterozigosity`

6.  the heterozygosity statistics per populations and
    overall:`$heterozygosity.statistics`

7.  the blacklisted individuals based on the individual's
    heterozigosity:`$blacklist.ind.het`

8.  a list containing the helper tables:`$helper.table.het`

9.  the boxplot of individual observed
    heterozygosity:`$individual.heterozygosity.boxplot`

10. the manhattan plot of individual
    heterozygosity:`$individual.heterozygosity.manhattan.plot`

11. the boxplot of observed heterozygosity averaged across markers and
    pop:`$markers.pop.heterozygosity.boxplot`

12. the density plot of observed heterozygosity averaged across markers
    and pop:`$markers.pop.heterozygosity.density.plot`

13. the manhattan plot of observed heterozygosity averaged across
    markers and pop:`$markers.pop.heterozygosity.manhattan.plot`

14. the relationship between the number of SNPs on a locus and locus
    observed heterozygosity statistics:`$snp.per.locus.het.plot`

15. the relationship between locus coverage (read depth) and locus
    observed heterozygosity statistics:`$het.cov.plot`

16. whitelist of markers:`$whitelist.markers`

17. whitelist of markers:`$whitelist.markers`

In the working directory, output is conditional to `interactive.filter`
argument

## Details

**Interactive version**

There are 4 steps in the interactive version to visualize and filter the
data based heterozygosity statistics:

Step 1. Individual's observed heterozygosity: outliers that might
represent mixed samples Step 2. Blacklist outliers based on a proportion
threshold of mean observed heterozygosity Step 3. Observed
heterozygosity statistics per populations and overall Step 4: Blacklist
markers based on observed heterozygosity

**outlier.pop.threshold**

If your a regular radr user, you've seen the `pop.num.threshold`.
`outlier.pop.threshold`, is different and requires more thinking,
because the number of populations genotyped potentially vary across
markers, which makes the use of `pop.num.threshold` less optimal. e.g.
If only 4 populations out of 10 are genotyped, for a marker you want to
keep, using `pop.num.threshold = 0.6` will lead to unwanted results and
inconsistensis. It's easier to tolerate outliers with this new approach:
`outlier.pop.threshold = 2`.

**Individual observed heterozygosity (averaged across markers):** To
help discard an individual based on his observed heterozygosity
(averaged across markers), use the manhanttan plot to:

1.  contrast the individual with population and overall samples.

2.  visualize the impact of missigness information (based on population
    or overall number of markers) and the individual observed
    heterozygosity. The larger the point, the more missing genotypes.

**Outlier above average:**

- potentially represent two samples mixed together (action: blacklist),
  or...

- a sample with more sequecing effort (point size small): did you merge
  your replicates fq files ? (action : keep and monitor)

- a sample with poor sequencing effort (point size large) where the
  genotyped markers are all heterozygotes, verify this with missingness
  (action: discard)

In all cases, if there is no bias in the population sequencing effort,
the size of the point will usually be "average" based on the population
or overall number of markers.

You can visualize individual observed heterozygosity, choose thresholds
and then visualize, choose thresholds and filter markers based on
observed heterozygosity in one run with: radr `filter_het`.

**Outlier below average:**

- A point with a size larger than the population or overall average (=
  lots of missing): the poor polymorphism discovery of the sample is
  probably the result of bad DNA quality, a bias in sequencing effort,
  etc. (action: blacklist)

- A point with a size that looks average (not much missing): this sample
  requires more attention (action: blacklist) and conduct more tests.
  e.g. for biallelic data, look for coverage imbalance between ALT/REF
  allele. At this point you need to distinguish between an artifact of
  poor polymorphism discovery or a biological reason (highly inbred
  individual, etc.).

## See also

[plot_density_distribution_het](https://thierrygosselin.github.io/radr/reference/plot_density_distribution_het.md)
