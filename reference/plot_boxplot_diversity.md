# Box plot of the diversity (Gene and Haplotypes)

GGPLOT2 Box plot of the diversity (Gene and Haplotypes).

GGPLOT2 radr helper plot used in MAC, MAF and MAD filter. Lines and
points, to help decide filter threshold.

## Usage

``` r
plot_boxplot_diversity(data, aes.x.y, y.title)

radr_helper_plot(
  data,
  strata = FALSE,
  stats,
  x.axis.title = NULL,
  x.breaks = NULL,
  text.orientation = "horizontal",
  plot.filename = NULL
)
```

## Arguments

- data:

  containing the stats

- aes.x.y:

  The GGPLOT2 aesthetics, e.g. aes.x.y = aes(x = factor(POP_ID), y =
  GENE_DIVERSITY, na.rm = T).

- y.title:

  Title of the y-axis.

- strata:

  Containing strata info for the visual. Default: `strata = FALSE`.

- stats:

  e.g. MAC, MAF or MAD.

- x.axis.title:

  Title for the x axis

- text.orientation:

  Orientation of the x axis text. Default:
  `text.orientation = "horizontal"horizontal`.

- plot.filename:

  To save the figure provide a filename. Default: `plot.filename = NULL`
