# simulate_genos_from_posterior

simulate values for the genotypes given the observed genotype, then
estimate allele frequencies, and the genotyping error rate.

This is a helper function for the estimate_m function

## Usage

``` r
simulate_genos_from_posterior(D, p, m)
```

## Arguments

- D:

  an 012,-1 matrix of observed genotypes

- p:

  the estimated allele freqs

- m:

  the genotyping error rate (must be a scalar)

## Author

Eric Anderson <eric.anderson@noaa.gov>
