# Summarise genomic data

Calculate summary statistics from tidy genomic data or a GDS file.
Statistics are calculated by population and marker and include REF and
ALT allele frequencies, observed and expected heterozygosity, and the
inbreeding coefficient (FIS).

## Usage

``` r
summarise_genomic_data(data, path.folder = NULL, digits = 4, verbose = TRUE)
```

## Arguments

- data:

  Tidy genomic data or a Genomic Data Structure (GDS) file or object:

  - tidy data

  - Genomic Data Structure (GDS)

  *How to get GDS and tidy data ?* Use
  [`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
  to import supported formats and
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.html)
  when a tidy table is needed.

- path.folder:

  (path, optional) By default will print results in the working
  directory. Default: `path.folder = NULL`.

- digits:

  (integer, optional). Default: `digits = 4`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.
