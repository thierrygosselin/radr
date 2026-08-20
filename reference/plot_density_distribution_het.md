# Figure density distribution of the observed heterozygosity summary statistics.

Create density distribution of the observed heterozygosity summary
statistics.

## Usage

``` r
plot_density_distribution_het(
  data,
  pop.levels,
  het.group,
  aes.colour,
  adjust.bin,
  x.title
)
```

## Arguments

- data:

  sumstats or tidy vcf files.

- pop.levels:

  Character string defining your ordered populations.

- het.group:

  = aes(x = HET_MAX, na.rm = F)

- aes.colour:

  GGPLOT2 aesthetics colour, e.g. aes(y = ..scaled.., color = GROUP).

- adjust.bin:

  Adjust GGPLOT2 bin size (0 to 1).

- x.title:

  Title of the x-axis.
