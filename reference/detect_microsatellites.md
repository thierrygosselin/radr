# Detect microsatellites

Detect Simple Sequence Reapeats (SSR) commonly known as
microsatellites... **radr** is not re-inventing the wheel here, it uses
the software [GMATA: Genome-wide Microsatellite Analyzing Toward
Application](https://github.com/XuewenWangUGA/GMATA).

## Usage

``` r
detect_microsatellites(data, gmata.dir = NULL, verbose = TRUE, ...)
```

## Arguments

- data:

  (path or object) Object in your global environment or a file in the
  working directory. The tibble must contain 2 columns named: `MARKERS`
  and `SEQUENCE`. When RADseq data from DArT is used,
  [`explore_genomes`](https://thierrygosselin.github.io/radr/reference/explore_genomes.md)
  generates automatically this file under the name
  `whitelist.markers.tsv`.

- gmata.dir:

  (path) For the function to work, the path to the directory with GMATA
  software needs to be given. If not found or `NULL`, the function
  download GMATA from github in the working directory. Default:
  `gmata.dir = NULL`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

6 files are returned in the folder: detect_microsatellites:

1.  ".fa.fms": the fasta file of sequences

2.  ".fa.fms.sat1": the summary of sequences analysed (not important)

3.  ".fa.ssr": The microsatellites found per markers (see GMATA doc)

4.  ".fa.ssr.sat2": Extensive summary (see GMATA doc).

5.  "blacklist.microsatellites.tsv": The list of markers with
    microsatellites.

6.  "whitelist.microsatellites.tsv": The whitelist of markers with NO
    microsatellites.

In the global environment, the object is a list with the blacklist and
the whitelist.

## Note

Thanks to Peter Grewe for the idea of including this type of filter
inside radr.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# The simplest way to run the function when the raw data was DArT:
mic <- radr::detect_microsatellites(data = "my_whitelist.tsv")

# With stacks pipeline, the populations module need to be run with --fasta-loci
# You could prepare the file this way (uncomment the function):
#
# prep_stacks_fasta <- function(fasta.file) {
#   fasta <- suppressWarnings(
#     vroom::vroom(
#      file = fasta.file,
#      delim = "\t",
#      col_names = "DATA",
#      col_types = "c",
#      comment = "#"
#    ) %>%
#      dplyr::mutate(MARKERS = stringi::stri_sub(str = DATA, from = 3, to = 7)) %>%
#      tidyr::separate(data = ., col = DATA, into = c("SEQUENCE", "LOCUS"), sep = "_")
#  )
#
#  fasta <- dplyr::bind_cols(
#   dplyr::filter(fasta, MARKERS == "Locus") %>%
#   dplyr::select(LOCUS),
#   dplyr::filter(fasta, MARKERS != "Locus") %>%
#   dplyr::select(SEQUENCE)
#  ) %>%
#  dplyr::mutate(LOCUS = as.numeric(LOCUS))
#   return(fasta)
#  } #prep_stacks_fasta

} # }
```
