# Figure density distribution of minor allele frequency (MAF) summary statistics.

Create density distribution of MAF summary statistics.

## Usage

``` r
plot_density_distribution_maf(
  data,
  maf.group,
  aes.colour = ggplot2::aes(y = ..scaled.., color = POP_ID),
  adjust.bin = 1,
  x.title
)
```

## Arguments

- data:

  sumstats or tidy vcf summarised files.

- maf.group:

  The GGPLOT2 aes (e.g. aes(x = FREQ_ALT, na.rm = F)).

- aes.colour:

  GGPLOT2 aesthetics colour, e.g. aes(y = ..scaled.., color = POP_ID) or
  aes(y = ..count.., color = POP_ID)

- adjust.bin:

  Adjust GGPLOT2 bin size (0 to 1).

- x.title:

  Title of the x-axis. e.g. "MAF distribution"

## Details

turn off the legend using `fig + theme(legend.position = "none")` or
zoom in a section with 'fig + coord_cartesian(xlim = c(0, 0.1), ylim =
c(0, 1))'. Save the figure with :
`ggsave("figure name.pdf", width = 40, height = 20, dpi = 600, units = "cm", useDingbats = F)`.

## See also

[tidy_genome](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.html)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
