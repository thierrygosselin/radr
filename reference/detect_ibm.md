# Detect identity-by-missingness structure

Computes a binary identity-by-missingness matrix directly from a GDS
object and returns a heatmap together with marker- and individual-level
missingness summaries. The heatmap reveals whether missing genotypes are
concentrated in particular individuals, markers, or groups, helping
users decide whether individual- or marker-level filtering should be
investigated first.

## Usage

``` r
detect_ibm(
  data,
  strata = NULL,
  strata.select = "STRATA",
  sort.individuals = "input",
  sort.markers = "input",
  sample.max = NULL,
  marker.max = 50000,
  title = "Identity-by-missingness",
  filename = NULL,
  image.width = 1800L,
  image.height = 2400L,
  image.res = 150,
  facet = TRUE,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A GDS filepath or an open `SeqVarGDSClass` object.

- strata:

  (optional) Path to a strata file with a minimum of 2 columns:
  `INDIVIDUALS` and the grouping column selected in `strata.select`.
  Default: `strata = NULL`.

- strata.select:

  (character) Column used to facet individuals in the heatmap. Default:
  `strata.select = "STRATA"`.

- sort.individuals:

  (character) Sorting method for individuals. Choices are: `"input"`,
  `"missingness"`, `"strata"`. Keeping the input order is recommended
  for the initial diagnostic plot because it can reveal sequencing,
  plate, library-preparation, or sample- processing structure. Default:
  `sort.individuals = "input"`.

- sort.markers:

  (character) Sorting method for markers. Choices are: `"input"`,
  `"missingness"`, `"position"`. Keeping the input order is recommended
  for the initial diagnostic plot. Default: `sort.markers = "input"`.

- sample.max:

  (optional, integer) Maximum number of individuals plotted. Default:
  `sample.max = NULL`.

- marker.max:

  (optional, integer) Maximum number of markers plotted. Default:
  `marker.max = 50000`.

- title:

  Optional heatmap title. Use `NULL` to omit it. Default:
  `title = "Identity-by-missingness"`.

- filename:

  Optional name for writing the heatmap directly to a PNG file inside
  the function results folder. When no `.png` extension is supplied, it
  is added automatically. The plot is not displayed by this argument; it
  remains available in the returned `heatmap` component. Default:
  `filename = NULL`.

- image.width:

  Width of the PNG in pixels. Default: `image.width = 1800`.

- image.height:

  Height of the PNG in pixels. Default: `image.height = 2400`.

- image.res:

  PNG resolution in pixels per inch. Default: `image.res = 150`.

- facet:

  (logical) Should multiple groups be labelled above the heatmap using
  `strata.select`? A redundant single-group label is omitted. Default:
  `facet = TRUE`.

- parallel.core:

  (optional, integer) Number of cores. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- ...:

  (optional) Further arguments passed to advanced mode. Use
  `path.folder` to select the parent results directory.

## Value

A list with:

- `heatmap`: a `ggplot2` object containing the raster

- `ibm.matrix`: integer matrix, 0 = missing, 1 = genotyped

- `raster`: raster representation used by the heatmap

- `image.file`: path to the PNG, or `NULL` when no file was written

- `individuals.missingness`: tibble with per-individual missingness

- `markers.missingness`: tibble with per-marker missingness

- `plot.data`: `NULL`; retained for compatibility because the raster
  renderer does not generate a long plotting table

## Details

Missing genotypes are shown in black and observed genotypes in grey.
Broad vertical bands indicate individuals with elevated missingness,
whereas broad horizontal bands indicate markers with elevated
missingness. Structure restricted to a stratum can indicate a group- or
batch-specific issue. These patterns can guide whether individual or
marker filtering should be examined first; they are diagnostic evidence
rather than an automatic filtering rule. The heatmap is intended to
reveal dataset-level patterns and guide filtering strategy. It is not
intended to identify specific markers or individuals; consequently,
individual and marker labels are deliberately omitted.

For the first inspection, the defaults preserve the current GDS order
for both individuals and markers. This is intentional: ordering
inherited from sequencing, plates, libraries, or sample processing can
expose technical structure that would be obscured by sorting on
missingness. The sorted modes are best used as complementary views after
examining the input-order plot.

## Further exploration of missing data

`detect_ibm()` is intentionally a fast, global diagnostic. It helps
reveal broad missingness structure and decide whether marker- or
individual-level filtering should be investigated first.

For a more detailed investigation, see the the [grur
package](https://thierrygosselin.github.io/grur/), particularly
[`grur::missing_visualization()`](https://thierrygosselin.github.io/grur/reference/missing_visualization.html).
It explores missingness at the individual, marker, and population
levels; examines relationships with biological and technical metadata;
and uses identity-by-missingness ordinations such as PCoA/MDS and RDA to
reveal additional patterns. The accompanying [missing-data analysis
vignette](https://thierrygosselin.github.io/grur/articles/vignette_missing_data_analysis.html)
provides a complete workflow.
