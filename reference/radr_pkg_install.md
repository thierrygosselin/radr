# Legacy radr dependency helper

`radr_pkg_install()` no longer installs or updates packages. Automatic
installation made library changes difficult to anticipate and obscured
the distinction between required and optional components. Use
[`radr_dependencies()`](https://thierrygosselin.github.io/radr/reference/radr_dependencies.md)
for a read-only diagnostic and follow the installation instructions in
the radr README.

## Usage

``` r
radr_pkg_install(check = TRUE, minimal.install = FALSE)
```

## Arguments

- check:

  Retained for compatibility and ignored. Default: `check = TRUE`.

- minimal.install:

  Retained for compatibility and ignored. Default:
  `minimal.install = FALSE`.

## Value

The dependency table returned by
[`radr_dependencies()`](https://thierrygosselin.github.io/radr/reference/radr_dependencies.md).
