# heterozygotes miscall rate

estimate the heterozygotes miscall rate from a simple model

## Usage

``` r
estimate_m(
  data,
  nreps = 200,
  m_init = stats::runif(1),
  a0 = 0.5,
  a1 = 0.5,
  sm = 0.005
)
```

## Arguments

- data:

  A tidy data frame (biallelic only).

- nreps:

  number of MCMC sweeps to do. Default: `nreps = 200`.

- m_init:

  initial starting value for m must be between 0 and 1

- a0:

  beta parameter for reference alleles

- a1:

  beta parameter for alternate alleles

- sm:

  standard devation of proposal distribution for m

## Author

Eric Anderson <eric.anderson@noaa.gov>
