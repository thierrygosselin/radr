# Figure of distribution (histogram)

Create density distribution of coverage summary statistics. Use the
coverage summary file created with coverage_summary function.

Create density distribution of coverage summary statistics. Use the
coverage summary file created with coverage_summary function.

## Usage

``` r
plot_histogram(data, variable, x.axis.title, y.axis.title)

plot_density_distribution_coverage(data, aes.colour, adjust.bin)
```

## Arguments

- data:

  Coverage summary file.

- aes.colour:

  GGPLOT2 aesthetics colour, e.g. aes(y = ..scaled.., color =
  COVERAGE_GROUP).

- adjust.bin:

  Adjust GGPLOT2 bin size (0 to 1).
