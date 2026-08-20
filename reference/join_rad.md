# join_rad

Join back the parts stripped in strip_rad. Used internally.

## Usage

``` r
join_rad(x, s, m, g, env.arg = NULL)
```

## Arguments

- x:

  The data

- s:

  The strata metadata.

- m:

  The markers metadata.

- g:

  The genotypes metadata.

- env.arg:

  You want to redirect
  [`rlang::current_env()`](https://rlang.r-lib.org/reference/stack.html)
  to this argument.
