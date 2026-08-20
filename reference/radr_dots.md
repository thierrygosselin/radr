# radr_dots

Parse and manage `...` for radr functions:

- assign supported arguments ("keepers") into a target environment

- fill in internal defaults for missing keepers

- detect deprecated and unknown arguments

- return a compact summary tibble for logging/debugging.

## Usage

``` r
radr_dots(
  func.name = as.list(sys.call())[[1]],
  fd = NULL,
  args.list = NULL,
  dotslist = NULL,
  keepers = c("subsample.markers.stats", "force.stats", "id.stats", "subsample",
    "filter.reproducibility", "filter.individuals.missing",
    "filter.individuals.heterozygosity", "filter.individuals.coverage.total",
    "filter.individuals.coverage.median", "filter.individuals.coverage.iqr",
    "filter.common.markers", "filter.monomorphic", "filter.ma", "ma.stats",
    "filter.coverage", "dp", "filter.genotyping", "filter.snp.position.read",
    "filter.snp.number", "filter.short.ld", "filter.long.ld", "long.ld.missing",
    "ld.method", 
     "ld.figures", "detect.mixed.genomes",
    "ind.heterozygosity.threshold", "detect.duplicate.genomes", "filter.hwe",
    "filter.strands", "random.seed", "path.folder", "filename", "parameters",
    "blacklist.genotypes", "erase.genotypes", "gt", "alt.dosage", "gt.vcf", "gt.vcf.nuc",
    "pop.levels", "pop.labels", "pop.select", "blacklist.id", "markers.info",
    "keep.allele.names", "keep.gds", "calibrate.alleles", "vcf.metadata", "vcf.stats",
    "wide", "whitelist.markers", "write.tidy", "missing.memory", "dart.sequence", 
    
    "internal", "heatmap.fst", "tidy.check", "tidy.vcf", "tidy.dart", "species",
    "population", "tau", "threshold.y.markers", "threshold.y.silico.markers",
    "sex.id.input", "threshold.x.markers.qr", "threshold.x.markers.RD",
    "threshold.x.markers.RD.silico", "mis.threshold.data", "mis.threshold.silicodata",
    "zoom.data", "zoom.silicodata", "sex.id.input", "het.qr.input"),
  deprecated = c("maf.thresholds", "common.markers", "max.marker", "monomorphic.out",
    "snp.ld", "filter.call.rate", "filter.markers.coverage", "filter.markers.missing",
    "number.snp.reads", "mixed.genomes.analysis", "duplicate.genomes.analysis",
    "ref.calibration"),
  env = parent.frame(),
  assign = TRUE,
  verbose = TRUE
)
```

## Arguments

- func.name:

  Name of the calling function (used only in messages). Default:
  `as.list(sys.call())[[1]]`.

- fd:

  (optional) The formal argument names of the calling function. Default:
  [`rlang::fn_fmls_names()`](https://rlang.r-lib.org/reference/fn_fmls.html).

- args.list:

  (optional) Named list of current arguments from the calling
  environment. Default: `as.list(environment())`.

- dotslist:

  Captured `...` arguments, usually from
  `rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)`.

- keepers:

  Character vector of argument names that are valid inside `...` and
  should be assigned into the calling environment.

- deprecated:

  Character vector of deprecated argument names.

- env:

  Environment where keeper arguments and defaults should be assigned.
  Default: [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html).

- assign:

  (logical) If `TRUE`, assign keeper and default values into `env`. If
  `FALSE`, no assignment is performed and the function only returns the
  summary tibble. Default: `assign = TRUE`.

- verbose:

  (logical) When `TRUE`, messages about function call arguments, `...`
  resolution, deprecated and unknown arguments are printed. Default:
  `verbose = TRUE`.

## Value

A tibble with columns:

- `ARGUMENTS`: argument name

- `VALUES`: character representation of the value or class

- `GROUPS`: one of `"fct.call.args"`, `"fct.call..."`, `"default..."`,
  `"deprecated..."`, `"unknowned..."`.

The side-effect (when `assign = TRUE`) is to assign resolved keeper
arguments and their defaults into `env`.
