# read blacklist of genotypes

Read a blacklist object or file.

Used internally in [radr](https://github.com/thierrygosselin/radr) and
might be of interest for users.

## Usage

``` r
read_blacklist_genotypes(blacklist.genotypes = NULL, verbose = FALSE, ...)
```

## Arguments

- blacklist.genotypes:

  (path or object) The blacklist is an object in your global environment
  or a file in the working directory (e.g. "blacklist.geno.tsv"). The
  dataframe contains at least these 2 columns: `MARKERS, INDIVIDUALS`.
  Additional columns are allowed: `CHROM, LOCUS, POS`.

  Useful to erase genotypes with bad QC, e.g. genotype with more than 2
  alleles in diploid likely sequencing errors or genotypes with poor
  genotype likelihood or coverage.

  Columns are cleaned of separators that interfere with some packages or
  codes, detailed in `clean_markers_names` and `clean_ind_names`
  Default: `blacklist.genotypes = NULL`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Life cycle

This function arguments will be subject to changes. Currently the
function uses erase.genotypes, but using the `dots-dots-dots ...`
arguments allows to pass `erase.genotypes and masked.genotypes`. These
arguments do exactly the same thing and only one can be used.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
bl <- radr::read_blacklist_genotypes(data = data,
    blacklist.genotypes = "blacklist.geno.iguana.tsv")
} # }
```
