# minor_allele_stats

Generate Minor Allele necessary statistics

## Usage

``` r
minor_allele_stats(
  data,
  calibrate.alleles = c("count", "depth", "ancestral"),
  keep.biallelic = TRUE,
  blacklist.markers = NULL,
  markers.meta = NULL,
  path.folder = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```
