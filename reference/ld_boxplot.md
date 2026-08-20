# ld_boxplot

Generate a boxplot of LD using subsampled markers. Return the matrix of
LD. Write to directory the stats and boxplot.

## Usage

``` r
ld_boxplot(
  gds,
  subsample.markers.stats = 0.2,
  ld.method = "r2",
  path.folder = NULL,
  parallel.core = parallel::detectCores() - 2,
  verbose = TRUE
)
```
