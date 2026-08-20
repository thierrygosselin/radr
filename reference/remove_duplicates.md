# Read tidy genomic data file ending .rad

Used internally in [radr](https://github.com/thierrygosselin/radr) To
remove duplicate individuals based on threshold established from the
visualization figures.

## Usage

``` r
remove_duplicates(
  data = "individuals.pairwise.dist.tsv",
  stats = "genotyped.statistics.tsv",
  dup.threshold = 0.25,
  diff.pop.remove = TRUE,
  path.folder = NULL
)
```

## Arguments

- data:

  (path) The individual's pairwise data. Default:
  `data = "individuals.pairwise.dist.tsv"`.

- stats:

  (path) The genotype statistics Default:
  `stats = "genotyped.statistics.tsv"`.

- dup.threshold:

  (double) The threshold to filter out duplicates Default:
  `dup.threshold = 0.25`.

- diff.pop.remove:

  Remove all individuals in pairs from different pop. Both samples are
  potentially problems. With defautl, the function will not keep one
  sample in the duplicate pair. Default: `diff.pop.remove = TRUE`.

## Value

A list with blacklisted duplicates. Write the blacklist in the working
directory.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
