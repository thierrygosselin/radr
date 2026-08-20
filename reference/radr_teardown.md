# Common teardown helper for radr functions

Companion to
[`radr_startup()`](https://thierrygosselin.github.io/radr/reference/radr_startup.md).
Restores working directory and options, stops the timer, and prints a
closing header.

Intended to be called from
[`on.exit()`](https://rdrr.io/r/base/on.exit.html) in the calling
function.

## Usage

``` r
radr_teardown(start.obj)
```

## Arguments

- start.obj:

  (list) The object returned by
  [`radr_startup()`](https://thierrygosselin.github.io/radr/reference/radr_startup.md).

## Value

Invisibly returns `NULL`. Called for its side effects.
