# haplotype_reconstruction

Reconstruct haplotypes

## Usage

``` r
haplotype_reconstruction(data, parallel.core = parallel::detectCores() - 1)
```

## Arguments

- data:

  A tidy genomic data frame or another genomic object supported by the
  calling function.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.
