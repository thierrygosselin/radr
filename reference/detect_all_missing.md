# Detect markers with all missing genotypes

Detect if markers in tidy dataset have no genotypes at all. Used
internally in [radr](https://github.com/thierrygosselin/radr) and might
be of interest for users.

## Usage

``` r
detect_all_missing(data, verbose = FALSE)
```

## Arguments

- data:

  A tidy data frame object in the global environment.

- verbose:

  Logical. Display progress messages. Default: `verbose = FALSE`.

## Value

The filtered dataset if problematic markers were found. Otherwise, the
untouch dataset.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
