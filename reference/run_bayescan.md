# Run a BayeScan genome scan

Run BayeScan, import its results, classify loci using a chosen
false-discovery-rate threshold, and generate diagnostic tables, figures,
whitelists, and blacklists.

**Function highlights:**

1.  **integrated and seamless pipeline:** generate
    [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) files within
    radr and run [BayeScan](http://cmpg.unibe.ch/software/BayeScan/)
    inside R!

2.  **unbalanced sampling sites impact:** measure and verify genome scan
    accurary in unbalanced sampling design with subsampling related
    arguments.

3.  **SNP linkage:** detect automatically the presence of multiple SNPs
    on the same locus and measure/verify accuracy of genome scan within
    locus.

4.  **summary tables and visualization:** the function generate summary
    tables and plots of genome scan.

5.  **whitelists and blacklists** of markers under different selection
    identity are automatically generated !

This function requires a working
[BayeScan](http://cmpg.unibe.ch/software/BayeScan/) program installed on
the computer ([install
instructions](http://cmpg.unibe.ch/software/BayeScan/download.md)). For
UNIX machines, please install the 64bits version.

## Usage

``` r
run_bayescan(
  data,
  n = 5000,
  thin = 10,
  nbp = 20,
  pilot = 5000,
  burn = 50000,
  pr_odds,
  fdr = 0.05,
  subsample = NULL,
  iteration.subsample = 1,
  parallel.core = parallel::detectCores() - 1,
  bayescan.path = NULL,
  conda.env = "genomics",
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  BayeScan input file, GDS file or open GDS object, or tidy genomic
  data.

  1.  Path to BayeScan input file. Generate this with
      [`write_bayescan`](https://thierrygosselin.github.io/genometranslator/reference/write_bayescan.html).

  2.  GDS file or open GDS object. Population assignments are obtained
      from the individual metadata stored in the GDS.

  3.  Tidy genomic data file or object. Tidy data and GDS input can be
      used with `subsample` and `iteration.subsample`.

- n:

  Integer. Number of output iterations. Default: `n = 5000`.

- thin:

  Integer. Thinning interval. Default: `thin = 10`.

- nbp:

  Integer. Number of pilot runs. Default: `nbp = 20`.

- pilot:

  Integer. Length of each pilot run. Default: `pilot = 5000`.

- burn:

  Integer. Burn-in length. Default: `burn = 50000`.

- pr_odds:

  Numeric. Prior odds for the neutral model. A `pr_odds = 10` indicates
  that the neutral model is 10 times more likely than the model with
  selection. Larger values make the analysis more conservative. This
  argument is required and has no default.

- fdr:

  Numeric false-discovery-rate threshold used to classify loci from
  BayeScan q-values. A locus is classified as diversifying or balancing
  only when its q-value is at or below this threshold; the sign of
  `ALPHA` determines the direction. BayeScan already calculates these
  Bayesian q-values; do not pass its posterior probabilities or
  posterior odds to the qvalue package as if they were p-values. A value
  of 0.05 is a common candidate-discovery threshold, whereas 0.01 is
  more conservative. The choice should reflect the cost of false
  discoveries and be examined together with prior-odds sensitivity.
  Default: `fdr = 0.05`.

- subsample:

  Integer, proportion, or `"min"`. With `subsample = 36`, 36 individuals
  in each populations are chosen randomly to represent the dataset. With
  `subsample = "min"`, the minimum number of individual/population found
  in the data is used automatically. Default: `subsample = NULL`.

- iteration.subsample:

  Integer. Number of repeated subsampling iterations. subsampling. With
  `subsample = 20` and `iteration.subsample = 10`, 20
  individuals/populations will be randomly chosen 10 times. Default:
  `iteration.subsample = 1`.

- parallel.core:

  Integer. Number of CPU cores available to BayeScan. Default:
  `parallel.core = parallel::detectCores() - 1`.

- bayescan.path:

  Character. Full path to the BayeScan executable. When `NULL`,
  [`check_bayescan()`](https://thierrygosselin.github.io/radr/reference/check_bayescan.md)
  searches the system `PATH` and the environment specified by
  `conda.env`. Default: `bayescan.path = NULL`. See details.

- conda.env:

  Conda environment name or prefix used when discovering BayeScan.
  Default: `conda.env = "genomics"`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

- ...:

  Additional arguments, including `path.folder` for the parent results
  directory.

## Value

For specific [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) output
files, see [BayeScan](http://cmpg.unibe.ch/software/BayeScan/)
documentation, please read the manual.

radr::run_bayescan outputs without subsampling:

1.  `bayescan`: dataframe with results of BayeScan analysis.

2.  `selection.summary`: dataframe showing the number of markers in the
    different group of selections and model choice.

3.  `whitelist.markers.positive.selection`: Whitelist of markers under
    diversifying selection and common in all iterations.

4.  `whitelist.markers.neutral.selection`: Whitelist of neutral markers
    and common in all iterations.

5.  `whitelist.markers.neutral.positive.selection`: Whitelist of neutral
    markers and markers under diversifying selection and common in all
    iterations.

6.  `blacklist.markers.balancing.selection`: Blacklist of markers under
    balancing selection and common in all iterations.

7.  `markers.dictionary`: BayeScan use integer for MARKERS info. In this
    dataframe, the corresponding values used inside the function.

8.  `pop.dictionary`: BayeScan use integer for STRATA info. In this
    dataframe, the corresponding values used inside the function.

9.  `bayescan.plot`: plot showing markers Fst and model choice.

    Additionnally, if multiple SNPs/locus are detected the object will
    also have:

10. `accurate.locus.summary`: dataframe with the number of accurate
    locus and the selection types.

11. `whitelist.accurate.locus`: whitelist of accurate locus.

12. `blacklist.not.accurate.locus`: blacklist of not accurate locus.

13. `accuracy.snp.number`: dataframe with the number of SNPs per locus
    and the count of accurate/not accurate locus.

14. `accuracy.snp.number.plot`: the plot showing the proportion of
    accurate/not accurate locus in relation to SNPs per locus.

15. `not.accurate.summary`: dataframe summarizing the number of not
    accurate locus with selection type found on locus.

radr::run_bayescan outputs WITH subsampling:

1.  `subsampling.individuals`: dataframe with indivuals subsample id and
    random seed number.

2.  `bayescan.all.subsamples`: long dataframe with combined iterations
    of bayescan results.

3.  `selection.accuracy`: dataframe with all markers with selection
    grouping and number of times observed throughout iterations.

4.  `accurate.markers`: dataframe with markers attributed the same
    selection grouping in all iterations.

5.  `accuracy.summary`: dataframe with a summary of accuracy of
    selection grouping.

6.  `bayescan.summary`: dataframe with mean value, averaged accross
    iterations.

7.  `bayescan.summary.plot`: plot showing markers Fst and model choice.

8.  `selection.summary`: dataframe showing the number of markers in the
    different group of selections and model choice.

9.  `whitelist.markers.positive.selection`: Whitelist of markers under
    diversifying selection and common in all iterations.

10. `whitelist.markers.neutral.selection`: Whitelist of neutral markers
    and common in all iterations.

11. `blacklist.markers.balancing.selection`: Blacklist of markers under
    balancing selection and common in all iterations.

12. `whitelist.markers.neutral.positive.selection`: Whitelist of neutral
    markers and markers under diversifying selection and common in all
    iterations.

13. `whitelist.markers.without.balancing.positive`: Whitelist of all
    original markers with markers under balancing selection and
    directional selection removed. The markers that remains are the ones
    to use in population structure analysis.

Other files are present in the folder and subsampling folder.

## Details

**Selection classification:** A q-value at or below `fdr` is required
before a locus is classified as a selection candidate. Among supported
candidates, positive `ALPHA` indicates diversifying selection and
negative `ALPHA` indicates balancing or purifying selection. Loci above
the FDR threshold are neutral regardless of the sign of `ALPHA`.

BayeScan q-values estimate the minimum Bayesian false discovery rate
incurred when a locus is included among the candidates. They are not
ordinary p-values. Consequently, applying `qvalue::qvalue()` a second
time to BayeScan posterior probabilities, posterior odds, or q-values is
not recommended. The qvalue package remains appropriate for methods that
produce valid p-values, including a typical pcadapt workflow.

**Limitations and complementary genome scans:** BayeScan is an
FST-outlier method. Its results depend on how well the model represents
population history and sampling. Hierarchical population structure,
isolation by distance, range expansion, bottlenecks, unequal effective
population sizes, admixture, linked markers, low-information variants,
and unbalanced sampling can alter power or increase false discoveries.
The neutral-model prior odds also influence posterior support,
especially for weakly informative loci.

Foll and Gaggiotti (2008) demonstrated this directly with a spatial
human expansion model. Including isolated populations that had undergone
severe bottlenecks increased false positives, particularly for
directional selection. Excluding those populations substantially reduced
the problem. This result does not justify removing inconvenient
populations after seeing the scan. Instead, it motivates analyses and
simulations based on plausible demographic histories, with population
inclusion rules defined from independent biological and historical
information.

Statistical power is also asymmetric. In the simulations of Foll and
Gaggiotti (2008), detecting balancing selection with biallelic AFLP or
SNP markers was nearly impossible when neutral FST was at or below 0.05.
Power depended strongly on genetic differentiation, the number of
populations, and the sample size. They found that approximately 30
individuals per population were generally sufficient when at least six
populations were analysed under their simulated conditions. This is a
study-specific result, not a universal sampling rule.

Finally, outlier behaviour is not unique evidence of selection.
Differences in mutation rate among loci can also produce outliers, a
concern emphasized by Foll and Gaggiotti (2008) for microsatellites.
They recommended separate analyses for marker classes with different
mutation processes, such as di-, tri-, and tetranucleotide
microsatellites. For SNP data, marker quality, ascertainment, linkage,
and demographic history remain important alternative explanations for an
apparent selection signal.

Candidate inversions and other low-recombination haploblocks need
explicit treatment. A broad BayeScan peak can reflect linked selection,
arrangement- frequency differences, recombination suppression, a
centromere, another structural variant, or a technical regional effect.
Screen for candidates before LD pruning with
[`detect_inversions`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md),
then compare the complete genome, a collinear sensitivity dataset, and
candidate-region or homokaryotype-only analyses. Use
[`genome_scan_context`](https://thierrygosselin.github.io/radr/reference/genome_scan_context.md)
to place peaks beside marker density, missingness, LD, and regional
annotations. No candidate region is removed automatically.

Treat significant loci as candidates rather than confirmed targets of
selection. Examine sensitivity to `pr_odds`, `fdr`, population grouping,
sample-size balance, linkage pruning, and marker filtering. When
possible, compare methods based on different summaries and assumptions.
For example, pcadapt detects markers excessively associated with
principal components and does not require predefined populations;
prepare its input with
[`write_pcadapt`](https://thierrygosselin.github.io/genometranslator/reference/write_pcadapt.html).
Environmental association or haplotype-aware methods may provide
additional evidence when suitable covariates or genomic information are
available.

Comparisons among genome-scan methods are conditional on the scenarios
used to evaluate them. For example, the elevated false-discovery rate
and reduced BayeScan power in the presence of admixed individuals
reported by Luu et al. (2017) came from simulations and should not be
generalized to every dataset. Conversely, analyses of whole-genome 1000
Genomes data by Meisner et al. (2021) showed that pcadapt can itself
produce inflated statistics under discrete population structure or when
principal components reflect sequencing and genotype-calling artefacts.
Empirical whole-genome scans can recover known adaptive regions, but
they cannot by themselves provide a complete false-discovery benchmark
because the true selected loci are not fully known.

Agreement among methods can strengthen a candidate's priority, but
disagreement is also informative because methods target different forms
of selection and respond differently to demographic history. Do not
require a simple intersection of every candidate list. Report
method-specific results, inspect genomic clustering and biological
context, and validate important candidates with independent data or
simulations tailored to the study design.

**subsampling:** During subsampling the function will automatically
remove monomorphic markers that are generated by the removal of some
individuals. Also, common markers between all populations are also
automatically detected. Consequently, the number of markers will change
throughout the iterations. The nice thing about the function is that
since everything is automated there is less chance of making an error...

**SNPs data set:** You should not run BayeScan with SNPs data set that
have multiple SNPs on the same LOCUS. Instead, run
genometranslator::genome_translator using the `snp.ld` argument to keep
only one SNP on the locus. Or run the function by first converting an
haplotype vcf or if your RAD dataset was produced by STACKS, use the
`batch_x.haplotypes.tsv` file! If the function detect multiple SNPs on
the same locus, accuracy will be measured automatically.

**UNIX install:** I like to transfer the *BayeScan2.1_linux64bits* (for
Linux) or the *BayeScan2.1_macos64bits* (for MACOs) in `/usr/local/bin`
and change it's name to `bayescan`. Too complicated ? and you've just
downloaded the last BayeScan version, I would try this :
`bayescan.path = "/Users/thierry/Downloads/BayeScan2.1/binaries/BayeScan2.1_macos64bits"`

Make sure to give permission: `sudo chmod 777 /usr/local/bin/bayescan`

## References

Foll, M., & Gaggiotti, O. E. (2008). A genome-scan method to identify
selected loci appropriate for both dominant and codominant markers: a
Bayesian perspective. *Genetics*, 180, 977-993.
[doi:10.1534/genetics.108.092221](https://doi.org/10.1534/genetics.108.092221)

Foll M, Fischer MC, Heckel G and L Excoffier (2010) Estimating
population structure from AFLP amplification intensity. Molecular
Ecology 19: 4638-4647

Fischer MC, Foll M, Excoffier L and G Heckel (2011) Enhanced AFLP genome
scans detect local adaptation in high-altitude populations of a small
rodent (Microtus arvalis). Molecular Ecology 20: 1450-1462

Excoffier, L., Hofer, T., & Foll, M. (2009). Detecting loci under
selection in a hierarchically structured population. *Heredity*, 103,
285-298. [doi:10.1038/hdy.2009.74](https://doi.org/10.1038/hdy.2009.74)

de Villemereuil, P., Frichot, E., Bazin, E., Francois, O., & Gaggiotti,
O. E. (2014). Genome scan methods against more complex models: when and
how much should we trust them? *Molecular Ecology*, 23, 2006-2019.
[doi:10.1111/mec.12705](https://doi.org/10.1111/mec.12705)

Lotterhos, K. E., & Whitlock, M. C. (2014). Evaluation of demographic
history and neutral parameterization on the performance of FST outlier
tests. *Molecular Ecology*, 23, 2178-2192.
[doi:10.1111/mec.12725](https://doi.org/10.1111/mec.12725)

Luu, K., Bazin, E., & Blum, M. G. B. (2017). pcadapt: an R package to
perform genome scans for selection based on principal component
analysis. *Molecular Ecology Resources*, 17, 67-77.
[doi:10.1111/1755-0998.12592](https://doi.org/10.1111/1755-0998.12592)

Meisner, J., Albrechtsen, A., & Hanghoj, K. (2021). Detecting selection
in low-coverage high-throughput sequencing data using principal
component analysis. *BMC Bioinformatics*, 22, 470.
[doi:10.1186/s12859-021-04375-2](https://doi.org/10.1186/s12859-021-04375-2)

## See also

[BayeScan](http://cmpg.unibe.ch/software/BayeScan/),
[`write_pcadapt`](https://thierrygosselin.github.io/genometranslator/reference/write_pcadapt.html),
and [pcadapt](https://bcm-uga.github.io/pcadapt/)

## Examples

``` r
if (FALSE) { # \dontrun{
# library(radr)
# get a tidy data frame and a bayescan file with genometranslator::genome_translator:
# to run with a vcf haplotype file
data <- genometranslator::genome_translator(
    data = "batch_1.haplotypes.vcf",
    strata = "../../02_project_info/strata.stacks.TL.tsv",
    whitelist.markers = "whitelist.filtered.markers.tsv",
    blacklist.id = "blacklist.id.tsv",
    output = "bayescan",
    filename = "bayescan.haplotypes"
    )
# to run BayeScan:
scan.pops <- radr::run_bayescan(
    data = "bayescan.haplotypes.txt",
    pr_odds = 1000
    )

# This will use the default values for argument: n, thin, nbp, pilot and burn.
# The number of CPUs will be the number available - 1 (the default).

# To test the impact of unbalance sampling run BayeScan with subsampling,
# for this, you need to feed the function the tidy data frame generated above
# with genometranslator::genome_translator:
scan.pops.sub <- radr::run_bayescan(
    data = data$tidy.data,
    pr_odds = 1000,
    subsample = "min",
    iteration.subsample = 10
    )

# This will run BayeScan 10 times, and for each iteration, the number of individuals
# sampled in each pop will be equal to the minimal number found in the pops
# (e.g. pop1 N = 36, pop2 N = 50 and pop3 N = 15, the subsampling will use 15
# individuals in each pop, taken randomly.
# You can also choose a specific subsample value with the argument.
} # }
```
