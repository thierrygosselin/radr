# MA diagnostic

Minor Allele diagnostic, help choose a filter threshold.

## Usage

``` r
diagnostic_ma(data, group.rank, filename = NULL)
```

## Arguments

- data:

  A file in the working directory or object in the global environment in
  wide or long (tidy) formats. To import, the function uses
  [radr](https://github.com/thierrygosselin/radr) `read_genome`. *See
  details of this function for more info*.

- group.rank:

  (Number) The number of group to class the MAF.

- filename:

  (optional) Name of the file written to the working directory.

## Details

Highly recommended to look at the distribution of MAF

## See also

[filter_ma](https://thierrygosselin.github.io/radr/reference/filter_ma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
problem <- radr::diagnostic_ma(
data = tidy.salmon.data, group.rank = 10)
} # }
```
