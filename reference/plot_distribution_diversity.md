# Density distribution of diversity (Gene and Haplotypes)

GGPLOT2 Density distribution of diversity (Gene and Haplotypes).

## Usage

``` r
plot_distribution_diversity(data, aes.x, aes.colour, x.title, y.title)
```

## Arguments

- data:

  The hapstats summary file or object.

- aes.x:

  GGPLOT2 aesthetics, e.g. aes.x = aes(x = GENE_DIVERSITY, na.rm = T).

- aes.colour:

  GGPLOT2 aesthetics colour, e.g. aes.colour = aes(y = ..scaled..,
  colour = POP_ID).

- x.title:

  Title of the x-axis.

- y.title:

  Title of the y-axis.
