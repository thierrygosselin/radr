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
  `ALPHA` determines the direction. Default: `fdr = 0.05`.

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

Foll, M and OE Gaggiotti (2008) A genome scan method to identify
selected loci appropriate for both dominant and codominant markers: A
Bayesian perspective. Genetics 180: 977-993

Foll M, Fischer MC, Heckel G and L Excoffier (2010) Estimating
population structure from AFLP amplification intensity. Molecular
Ecology 19: 4638-4647

Fischer MC, Foll M, Excoffier L and G Heckel (2011) Enhanced AFLP genome
scans detect local adaptation in high-altitude populations of a small
rodent (Microtus arvalis). Molecular Ecology 20: 1450-1462

## See also

[BayeScan](http://cmpg.unibe.ch/software/BayeScan/)

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
