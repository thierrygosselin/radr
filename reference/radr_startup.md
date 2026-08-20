# Common startup helper for radr functions

Perform standard radr initialisation steps inside a function:

- print a header with
  [`tgbase::function_header()`](https://rdrr.io/pkg/tgbase/man/function_header.html);

- record execution date/time;

- record and temporarily modify global options (e.g. `width`,
  `future.globals.maxSize`);

- start a timing object with
  [`tgbase::tic()`](https://rdrr.io/pkg/tgbase/man/tic.html).

The function returns a small list with everything needed for a matching
teardown helper. You should typically pair this with
[`radr_teardown()`](https://thierrygosselin.github.io/radr/reference/radr_teardown.md)
inside an [`on.exit()`](https://rdrr.io/r/base/on.exit.html) call in the
calling function.

## Usage

``` r
radr_startup(f.name, verbose = TRUE, width = 70L)
```

## Arguments

- f.name:

  (character) Name of the calling function (e.g. `"read_vcf"`). Used
  only for logging in
  [`tgbase::function_header()`](https://rdrr.io/pkg/tgbase/man/function_header.html).

- verbose:

  (logical) When `TRUE`, the helper prints the execution date/time and
  passes verbosity to the teardown helper. Default: `verbose = TRUE`.

- width:

  (integer) Temporary value for `getOption("width")`. Default:
  `width = 70`.

## Value

A named list containing:

- `file.date` – character timestamp `"YYYYMMDD@HHMM"`;

- `old.dir` – working directory to restore;

- `opt.width` – original `options("width")`;

- `opt.future` – original `options("future.globals.maxSize")`;

- `timing` – timer object from
  [`tgbase::tic()`](https://rdrr.io/pkg/tgbase/man/tic.html);

- `f.name` – function name;

- `verbose` – verbosity flag.
