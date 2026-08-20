# Fis filter

TODO

## Usage

``` r
filter_fis(
  data,
  approach = "haplotype",
  fis.min.threshold,
  fis.max.threshold,
  fis.diff.threshold,
  pop.threshold,
  percent,
  filename,
  verbose = TRUE
)
```

## Arguments

- data:

  A GDS filename or open `SeqVarGDSClass` object.

- approach:

  Character. By `"SNP"` or by `"haplotype"`. The function will consider
  the SNP or haplotype statistics to filter the marker. Default:
  `approach = "haplotype"`.

- fis.min.threshold:

  Number.

- fis.max.threshold:

  Number.

- fis.diff.threshold:

  Number (0 - 1)

- pop.threshold:

  Fixed number of pop required to keep the locus.

- percent:

  Is the threshold a percentage ? TRUE or FALSE.

- filename:

  (optional) The function uses
  [`write.fst`](http://www.fstpackage.org/reference/write_fst.md), to
  write the tidy data frame in the folder created in the working
  directory. The file extension appended to the `filename` provided is
  `.rad`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

## Details

The Fis calculated uses the ratio of averages (1-mean(Ho)/mean(Hs)) and
NOT THE AVERAGE OF RATIOS (Nei 1987).

## References

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University
Press.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
