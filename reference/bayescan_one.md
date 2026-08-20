# bayescan one iteration

bayescan_one

## Usage

``` r
bayescan_one(
  x = NULL,
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
  path.folder,
  file.date,
  bayescan.path = "/usr/local/bin/bayescan"
)
```
