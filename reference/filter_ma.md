# MAC, MAF and MAD filter

The Allele Frequency is the relative frequency of an allele for a
markers/SNP. Calculations usually involve classifying the alleles into
Major/Minor, Reference (REF)/Alternate (ALT) and/or Ancestral or derived
allele. radr provides a filter based of Minor Allele Count, frequency
and depth function. Remove/blacklist markers based on Minor/Alternate
Allele Count (MAC), Frequency (MAF) or Depth of coverage (MAD). Use it
to remove noise, sequencing errors or low marker polymorphism. Some
analysis performs better with the full spectrum of allele frequency,
consequently it's important to understand the limite of using the
thresholds to avoid generating biases.

**Filter targets**: Marker's alternate allele(s)

**Statistics**: count, frequency or depth of the allele per markers.
This filter as NO concept of strata/populations. It is computed
globally.

## Usage

``` r
filter_ma(
  data,
  interactive.filter = TRUE,
  ma.stats = "mac",
  filter.ma = NULL,
  calibrate.alleles = "depth",
  keep.biallelic = TRUE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A tidy genomic data frame or another genomic object supported by the
  calling function.

- interactive.filter:

  Logical indicating whether an interactive filtering session may
  display diagnostics and ask for thresholds. Default:
  `interactive.filter = TRUE`.

- ma.stats:

  (optional, character) The statistic of the alternate (minor) allele to
  keep a SNP (see details). Options are: `"mac", "maf", "mad"`. Default:
  `ma.stats = "mac"`.

- filter.ma:

  (optional, integer or double) The threshold to blacklist Minor Allele
  (see details). e.g. with `ma.stats = "mac"` and `filter.ma = 3` will
  blacklist markers based on minor allele count less or equal to 3. This
  argument as no concept of locus, local or strata or pop. It's applied
  globally, by SNPs. Default: `filter.ma = NULL`.

- calibrate.alleles:

  (character) When ancestral allele information is not available, REF /
  ALT alleles can be calibrated using alleles count or depth (sequencing
  coverage for REF and ALT alleles). Removing individuals will impact
  REF/ALT alleles. Ideally, it should be based on observed read depth,
  when this information is available for both alleles. By selecting
  `calibrate.alleles = "depth"`, radr will check if the information is
  available, if not, it will use `calibrate.alleles = "count"`. Using
  `calibrate.alleles = "ancestral"` will turn off calibration and use
  the information preset available in the data. Default:
  `calibrate.alleles = "depth"`.

- keep.biallelic:

  (logical) Keep only biallelic variants (REF + 1 ALT) when computing
  allele counts from GDS. Default: `keep.biallelic = TRUE`.

- filename:

  (optional, character) Write to folder the MAF filtered tidy dataset.
  The name will be appended `.rad`. With default, the filtered data is
  only in the global environment. Default: `filename = NULL`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level screening or filtering
  functions.

## Value

With `interactive.filter = FALSE`, a list in the global environment,
with 7 objects:

1.  \$tidy.filtered.mac

2.  \$whitelist.markers

3.  \$blacklist.markers

4.  \$mac.data

5.  \$filters.parameters

With `interactive.filter = TRUE`, a list with 4 additionnal objects are
generated.

1.  \$distribution.mac.global

2.  \$distribution.mac.local

3.  \$mac.global.summary

4.  \$mac.helper.table

**mac.helper.table:**

First and second variables (in columns) represents POP_ID and sample
size (n).

The last 2 rows are the local MAF (suggested based on the lowest pop
value) and the TOTAL/GLOBAL observations.

Columns starting with ALT are the variable corresponding to the number
of alternative (ALT) allele (ranging from 1 to 20). The observations in
the ALT allele variable columns are the local (for the pop) and global
(last row) MAF of your dataset. e.g. ALT_3 can potentially represent 3
heterozygote individuals with the ALT allele or 1 homozygote individuals
for the ALT allele and 1 heterozygote individual. And so on...

## Note

Thanks to Charles Perrier and Jeremy Gaudin for very useful comments on
previous version of this function.

## MAC, MAF or MAD ?

Using count or frequency to remove a SNPs ? The preferred choice in radr
as changed from frequency to count, because we think the filtering
should not alter the spectrum and this is only achieved if the same
criteria is applied for each SNP.

Even small differences in missing data between RADseq markers generates
differences in MAF frequency thresholds applied.

Example with a datset consisting of N = 36 individuals and 3 SNPs with
varying level of missing genotypes:

- `SNP number : number samples genotypes : REF/ALT counts`

- SNP1 : 36 : 69/3

- SNP2 : 30 : 65/3

- SNP3 : 24 : 45/3

Each SNPs have the same alternate allele count, corresponding to 2
individuals with the polymorphism: 1 homozygote + 1 heterozygote.
Applying a MAF threshold of 0.05 would mean that SNP3 would be
blacklisted (`24 * 2 * 0.05 = 2.4 alt alleles required to pass`).

**Using count instead of frequency allows each RADseq markers, with
varying missing data, to be treated equally.**

## Advance mode

*dots-dots-dots ...* allows to pass several arguments for fine-tuning
the function:

1.  `filter.common.markers` (optional, logical). Default:
    `filter.common.markers = FALSE`, Documented in
    [`filter_common_markers`](https://thierrygosselin.github.io/radr/reference/filter_common_markers.md).

2.  `filter.monomorphic` (logical, optional) Should the monomorphic
    markers present in the dataset be filtered out ? Default:
    `filter.monomorphic = TRUE`. Documented in
    [`filter_monomorphic`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.md).

3.  `path.folder`: to write ouput in a specific path (used internally in
    radr). Default: `path.folder = getwd()`. If the supplied directory
    doesn't exist, it's created.

## Interactive version

To help choose a threshold for the local and global MAF use the
interactive version.

2 steps in the interactive version:

Step 1. Visualization and helper table.

Step 2. Minor Allele statistic and threshold selection

Step 3. Filtering markers based on MAC, MAF or MAD

## References

Good reading: Linck, E., Battey, C. (2019). Minor allele frequency
thresholds strongly affect population structure inference with genomic
data sets. Molecular Ecology Resources 19(3), 639-647.
https://dx.doi.org/10.1111/1755-0998.12995

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# The minumum
mac <- radr::filter_ma(data = turtle.tidy.data)
# This will use the default: interactive version,
# a list is created and to view the filtered tidy data:
mac.tidy.data <- mac$tidy.filtered.mac

# No user interaction

# Using filter.ma = 4, is a good practice to remove mostly sequencing errors
# and assembly artifacts, because it requires the markers to be genotyped in
# 4 heterozygote individuals or 2 homozygote individuals for the alternate
# allele or 2 heterozygote individuals and 1 homozygote individual for the
# alternate allele. Overall, never less than 2 indiduals are required.

mac <- radr::filter_ma(
        data = turtle.gds, # using gds object
        ma.stats = "mac",
        filter.ma = 4,
        filename = "turtle.mac")

# This will remove monomorphic markers (by default filter.monomorphic = TRUE)
# The filtered data will be written in the directory under the name:
# turtle.maf.rad
} # }
```
