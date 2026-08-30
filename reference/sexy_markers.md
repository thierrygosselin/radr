# Identify sex-linked markers and reassign genetic sex

Identifies markers putatively located on heterogametic or homogametic
sex chromosomes and can reassign genetic sex using candidate Y- or
W-linked markers. The function works best in: DArT silico (counts) \>
DArT counts or RADseq with allele read depth \> DArT silico (genotypes)
\> RADseq (genotypes) and DArT (1-row, 2-rows genotypes).

## Usage

``` r
sexy_markers(
  data,
  silicodata = NULL,
  strata = NULL,
  boost.analysis = FALSE,
  coverage.thresholds = 1,
  filters = TRUE,
  interactive.filter = TRUE,
  folder.name = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A genomic input supported by
  [`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
  that resolves to a GDS object. An open GDS object is preferred for
  repeated analyses. VCF, DArT, PLINK, and other supported inputs can be
  translated automatically. For format-specific import control, call the
  corresponding `genometranslator::read_*()` function first.

- silicodata:

  (optional, file) A silico (dominant marker) DArT file
  `.csv or .tsv`.  
  This can be count or genotyped data. Note that both `data` and
  `silicodata` can be used at the same time.  
  Default: `silicodata = NULL`.

- strata:

  (file) A tab delimited file with a minimum of 2 columns
  (`INDIVIDUALS, STRATA`) for VCF files and 3 columns for DArT files
  (`TARGET_ID, INDIVIDUALS, STRATA`).

  - `TARGET_ID:` it's the column header of the DArT file.

  - `STRATA:` is the grouping column, here the sex info. 3 values: `M`
    for male, `F` for female and `U` for unknown. Anything else is
    converted to `U`.

  - `INDIVIDUALS:` Is how you want your samples to be named.

  Default: `strata = NULL`, the function will look for sex info in the
  tidy data or GDS data (individuals.meta section). You can easily build
  the strata file by starting with the output of these functions:
  [`extract_dart_target_id`](https://thierrygosselin.github.io/genometranslator/reference/extract_dart_target_id.html)
  and `extract_individuals_vcf`

- boost.analysis:

  (optional, logical) This method uses machine learning approaches to
  find sex markers and re-assign samples in sex group.  
  The approach is currently under construction. Default:
  `boost.analysis = FALSE`.

- coverage.thresholds:

  (optional, integer) The minimum coverage required to call a marker
  absent. For silico genotype data this must be \< 1.  
  Default: `coverage.thresholds = 1`.

- filters:

  (optional, logical) When `filters = TRUE`, the data will be filtered
  for monomorphic loci, missingness of individuals and heterozygosity of
  individuals.  
  CAUTION: we advice to use these filter, since not filtering or
  filtering too stringently will results in false positive or false
  negative detections.  
  Default: `filters = TRUE`.

- interactive.filter:

  (optional, logical) When `interactive.filter = TRUE` the function will
  ask for your input to define thresholds. If
  `interactive.filter = FALSE` the function expects additional arguments
  (see Advance mode).  
  Default: `interactive.filter = TRUE`.

- folder.name:

  (optional,character) Name of the folder to store the results. Default:
  `folder.name = NULL`. The name sexy_markers_datetime will be
  generated.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

The created object contains:

1.  A list with (1) the summarised SNP data per sex and (2) the
    summarised silico data per sex. This should allow you to re-create
    the various plots.

2.  A vector with the names of the sex-linked marker. One vector for
    each sex method.

3.  A dataframe with a summary of the sex-linked markers and their
    sequence (if available).

## Details

This function takes DArT and RAD-type data to find markers that have a
specific pattern that is linked to sex.  
The function hypothesizes the presence of sex-chromosomes in your
species/population. The tests are designed to identify markers that are
located on putative heterogametic (Y or W) or homogametic (X or Z)
chromosomes.  
*Note:* Violating Assumptions or Prerequisites (see below) can lead to
false positive or the absence of detection of sex-linked markers.

## Assumptions

1.  **Genetic Sex Determination System** over a e.g.
    environmental-sex-determination system.

2.  **Genome coverage:** restriction sites randomly spread throughout
    the whole genome.

3.  **Mutations:** Processes such as sex-specific mutations in the
    restriction sites could lead to false positive results.

4.  **Deletions/duplications:** Processes such as sex-specific deletions
    or duplications could lead to false positive results.

5.  **Homology:** The existence of homologous sequences between the
    homogametic and heterogametic chromosomes could lead to false
    negative results.

6.  **Absence of population signal**

## Prerequisites

1.  **Sample size:** Ideally, the data must have enough individuals (n
    \> 100).

2.  **Batch effect:** Sex should be randomized on lanes/chips during
    sequencing.

3.  **Sex ratio:** Dataset with equal ratio work best.

4.  **Genotyping rate**: for DArT data, if the minimum call rate is
    greater than 0.5, ask DArT to lower its filtering threshold. RADseq
    data, lower markers missingness thresholds during filtering (e.g.
    stacks `r` and `p`).

5.  **Identity-by-Missingness:** Absence of artifactual pattern of
    missingness ([see missing
    visualization](https://github.com/thierrygosselin/grur))

6.  **Low genotyping error rate:** see
    [`detect_het_outliers`](https://thierrygosselin.github.io/radr/reference/detect_het_outliers.md)
    and [whoa](https://github.com/eriqande/whoa).

7.  **Low heterozygosity miscall rate:** see
    [`detect_het_outliers`](https://thierrygosselin.github.io/radr/reference/detect_het_outliers.md)
    and [whoa](https://github.com/eriqande/whoa).

8.  **Absence of pattern of heterozygosity driven by missingness:** see
    [`detect_mixed_genomes`](https://thierrygosselin.github.io/radr/reference/detect_mixed_genomes.md).

9.  **Absence of paralogous sequences in the data**.

## Sex methods

**Heterogametic sex-markers:**

- **Presence/Absence method**: To identify markers on Y or W
  chromosomes, we look at the presence or absence of a marker between
  females and males. More specifically, if a marker is always present in
  males but never in females, they are putatively located in the
  Y-chromosome; and vice versa for the W-linked markers.

**Homogametic sex-markers:** We have two different methods to identify
markers on X or Z chromosomes:

- **Heterozygosity method:** By looking at the heterozygosity of a
  marker between sexes, we can identify markers that are always
  homozygous in one sex (e.g. males for an XY system), while exhibiting
  an intermediate range of heterozygosity in the other sex (0.1 - 0.5).

- **Coverage method:** If the data includes count information, this
  function will look for markers that have double the number of counts
  for either of the sexes. For example if an XY/XX system is present,
  females are expected to have double the number of counts for markers
  on the X chromosome.

## Advance mode

*dots-dots-dots ...* allows to pass several arguments for fine-tuning
the function:

- `species`: To give your figures some meanings. Default
  `species = NULL`.

- `population`: To give your figures some meanings. Default
  `species = NULL`.

- `tau`: The quantile used in regression to distinguish homogametic
  markers with the **heterozygosity method**. See
  [`rq`](https://rdrr.io/pkg/quantreg/man/rq.html).  
  Default `tau = 0.03`.

- `mis.threshold.data`: Threshold to filter the SNP data on missingness.
  Only if `interactive.filter = FALSE`.

- `mis.threshold.silicodata`: Threshold to filter the silico data on
  missingness. No default. Only if `interactive.filter = FALSE`.

- `threshold.y.markers`: Threshold to select heterogametic sex-linked
  markers from the SNP data with the **presence/absence method**.No
  default. Only if `interactive.filter = FALSE`.

- `threshold.x.markers.qr`: Threshold to select homogametic sex-linked
  markers from the SNP data with the **heterozygosity method**. No
  default. Only if `interactive.filter = FALSE`.

- `zoom.data`: Threshold to subset the F/M ratio of mean SNP coverage.
  Used to improve the histogram resolution to select a better
  `threshold.x.markers.RD` threshold. No default. Only if
  `interactive.filter = FALSE`.

- `threshold.x.markers.RD`: Threshold to select homogametic sex-linked
  markers from the SNP data with the **coverage method**.No default.
  Only if `interactive.filter = FALSE`.

- `zoom.silicodata`: Threshold to subset the F/M ratio of mean silico
  coverage. Used to improve the histogram resolution to select a better
  `threshold.x.markers.RD.silico` threshold. No default. Only if
  `interactive.filter = FALSE`.

- `threshold.x.markers.RD.silico`: Threshold to select homogametic
  sex-linked markers from the silico data with the **coverage method**.
  No default. Only if `interactive.filter = FALSE`.

- `sex.id.input`: (integer) `sex.id.input = c(1, 2 or 3)` to recalculate
  the sex based on (1) 'visual', (2) 'genetic SNP' or (3) 'genetic
  SILICO' sexID. No default. Only if `interactive.filter = FALSE`.

- `het.qr.input`: (integer) `het.qr.input = c(1 or 2)` to plot the
  heterozygosity residual plot for (1) X-linked markers (heterozygous
  for females), or (2) Z-linked markers (heterozygous for males). No
  default. Only if `interactive.filter = FALSE`.

## Life cycle

Machine Learning approaches (Random Forest and Extreme Gradient Boosting
Trees) are currently been tested. They usually show a lower discovery
rate but tend to perform better with new samples.

## See also

Eric Anderson's [whoa](https://github.com/eriqande/whoa) package.

## Author

Floriaan Devloo-Delva <Floriaan.Devloo-Delva@csiro.au> and with help
from Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# The minimum
sex.markers <- radr::sexy_markers(
    data = "shark.dart.data.csv",
    strata = "shark.strata.tsv")
# This will use the default: interactive version and a list is created and to view the sex markers
} # }
```
