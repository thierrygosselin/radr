# Visual diagnostic of coverage imbalance

GBS data and STACKS pipeline sometimes output REF and ALT alleles that
have coverage imbalance, i.e. the coverage is not equal and skewed
towards the REF or ALT alleles. Thw density distribution figure of
coverage imbalance between REF and ALT alleles will highlight the
problem in you data.

## Usage

``` r
plot_coverage_imbalance_diagnostic(
  tidy.vcf.file,
  pop.levels,
  read.depth.threshold,
  aes.colour,
  adjust.bin
)
```

## Arguments

- tidy.vcf.file:

  The tidy VCF file created with tidy_genome.

- pop.levels:

  Character string defining your ordered populations.

- read.depth.threshold:

  Define the threshold you wish to analyse.

- aes.colour:

  GGPLOT2 aesthetics, e.g. aes(y = ..count..).

- adjust.bin:

  Adjust GGPLOT2 bin size (0 to 1).

## Value

4-plots highlighting the different combination of under or over the
coverage threshold and mean genotype likelihood. Y- axis show the
distribution of genotypes. X- axis show the coverage imbalance the
Negative ratio (left x axis) : REF \> ALT. Positive ratio (right x axis)
: ALT \> REF.

## Details

The figures shows the results of the of coverage threshold selected and
mean genotype likelihood. You can test different threshold to inspect
your data. Ideally the lower left pannel of the 4-plot should be empty.
If it is, this shows that setting the threshold of the genotype
likelihood filter to the mean or close to it take care of the allelic
coverage imbalance. \#' e.g. fig \<- plot_coverage_imbalance_diagnostic(
tidy.vcf.file, pop.levels, read.depth.threshold, aes.colour, adjust.bin)
Use ( fig + facet_grid(GROUP_GL ~ GROUP_COVERAGE)). The ratio is
calculated : (read depth ALT allele - read depth REF allele)/(read depth
ALT allele + read depth REF allele).
