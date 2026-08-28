# Detect candidate inversion-associated genomic regions

Scan a diploid, biallelic GDS for genomic windows whose local population
structure differs from the genomic background. The implementation
follows the local-PCA principle: each window is represented by a
low-rank covariance matrix among individuals, distances are calculated
between those matrices, and classical multidimensional scaling (MDS) is
used to identify unusual windows. Candidate regions are then summarised
using regional PCA, heterozygosity, and linkage disequilibrium (LD).

## Usage

``` r
detect_inversions(
  data,
  chromosome = NULL,
  window.snps = 100L,
  step.snps = window.snps,
  window.bp = NULL,
  step.bp = window.bp,
  window.method = c("snps", "bp", "ld"),
  ld.window.threshold = 0.1,
  ld.window.min.snps = 50L,
  ld.window.max.snps = 500L,
  sensitivity.window.snps = c(100L, 250L, 500L, 1000L),
  n.pcs = 2L,
  mds.axes = 2L,
  outlier.quantile = 0.99,
  min.window.snps = max(10L, n.pcs + 2L),
  min.call.rate = 0.8,
  min.candidate.windows = 1L,
  cluster.k = 3L,
  arrangement.labels = c("AA", "AB", "BB"),
  known.regions = NULL,
  ld.max.snps = 500L,
  return.ld = FALSE,
  save.plots = TRUE,
  plot.formats = c("png", "pdf"),
  random.seed = 42L,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A GDS filepath or an open `SeqVarGDSClass` object.

- chromosome:

  Optional chromosome or scaffold names to scan. By default, all
  chromosomes represented by at least one complete window are scanned.

- window.snps:

  Number of SNPs per window.

- step.snps:

  Number of SNPs between consecutive window starts. Defaults to
  `window.snps` (non-overlapping windows).

- window.bp:

  Optional fixed physical window size in base pairs. When supplied,
  physical windows are used instead of fixed-SNP windows.

- step.bp:

  Distance in base pairs between physical window starts. Defaults to
  `window.bp`.

- window.method:

  Window construction method: `"snps"` for a fixed number of SNPs,
  `"bp"` for fixed physical windows, or `"ld"` for experimental
  LD-scaled windows. Supplying `window.bp` selects `"bp"`. Default:
  `window.method = "snps"`.

- ld.window.threshold:

  Adjacent-marker r-squared threshold used to end an experimental
  LD-scaled window after the minimum number of markers. Default:
  `ld.window.threshold = 0.1`.

- ld.window.min.snps:

  Minimum markers in an LD-scaled window. Default:
  `ld.window.min.snps = 50`.

- ld.window.max.snps:

  Maximum markers examined in an LD-scaled window. This bounds memory
  and computation. Default: `ld.window.max.snps = 500`.

- sensitivity.window.snps:

  Additional fixed-SNP window sizes used for a sensitivity analysis.
  These runs summarise PC1 variance and LD along the genome without
  replacing the primary candidate scan. Use `NULL` to skip this
  additional work.

- n.pcs:

  Number of local covariance axes retained per window.

- mds.axes:

  Number of MDS axes used to score unusual windows.

- outlier.quantile:

  Quantile of the robust window score used as the candidate threshold.

- min.window.snps:

  Minimum number of usable polymorphic SNPs required in a window after
  missing-data and variance checks.

- min.call.rate:

  Minimum genotype call rate required for a SNP within a window.

- min.candidate.windows:

  Minimum number of consecutive candidate windows required to form a
  candidate region.

- cluster.k:

  Number of regional PCA clusters. The biological expectation for a
  common polymorphic inversion is often three, representing the two
  homokaryotypes and their heterokaryotype, but this is diagnostic
  rather than proof.

- arrangement.labels:

  Three labels, ordered from the lowest to highest regional PC1 cluster,
  used when `cluster.k = 3`. Default:
  `arrangement.labels = c("AA", "AB", "BB")`.

- known.regions:

  Optional data frame describing centromeres, regions of low
  recombination, assembly gaps, or other annotations. It must contain
  `chromosome`, `start`, `end`, and `type` columns. Overlapping
  annotation types are reported for each candidate but are not used to
  select or score candidate windows.

- ld.max.snps:

  Maximum number of evenly spaced SNPs used for each regional LD matrix.
  This bounds memory use without changing the GDS input.

- return.ld:

  Logical indicating whether sampled regional LD matrices are retained
  in the result.

- save.plots:

  Logical indicating whether standard PDF and PNG diagnostic figures are
  written to the results folder.

- plot.formats:

  One or more of `"png"` and `"pdf"`.

- random.seed:

  Integer seed used for reproducible k-means clustering.

- verbose:

  Logical indicating whether progress messages are printed.

- ...:

  Further arguments for the standard radr workflow. Use `path.folder` to
  choose the parent results directory.

## Value

An object of class `detect_inversions` containing:

- `windows`: coordinates, quality statistics, MDS coordinates, robust
  scores, and candidate flags for all analysed windows;

- `candidates`: one row per contiguous candidate region;

- `diagnostics`: regional PCA scores, cluster assignments,
  heterozygosity summaries, and optional LD matrices;

- `arrangement.genotypes`: one row per individual and candidate with
  putative arrangement genotype and relative assignment confidence;

- `homokaryotype.whitelist`: candidate-specific `AA` and `BB`
  individuals, plus `homokaryotype.all.candidates` for the intersection;

- `sensitivity`: optional summaries for additional fixed-SNP window
  sizes;

- `path.folder` and `output.files`: locations of written results;

- `settings`: the effective analysis settings.

## Details

Run this screen before LD pruning. LD pruning can remove the extended
correlation pattern that makes an inversion-associated haploblock
detectable. After candidates have been reviewed, repeat downstream
analyses with the complete genome, with candidate regions excluded, and
within each candidate region or inferred arrangement.

This is a screening method. A candidate region is not proof of a
physical inversion, and the returned coordinates describe an
inversion-associated haploblock rather than validated breakpoints.
Long-read, read-pair, split-read, cytogenetic, or genetic-map evidence
is needed for breakpoint confirmation.

## Chromosome-specific windows

Windows are constructed independently within each chromosome, linkage
group, or scaffold. A window can therefore never contain markers from
two linkage groups. Local PCA is not performed once on an entire linkage
group: it is performed separately for every SNP window within that
linkage group. The resulting window-level covariance summaries are
compared across all valid windows in the requested scan. Use
`chromosome` to restrict that comparison to one or more linkage groups.
Standard genomic-position figures are faceted by chromosome or linkage
group in the results folder.

## Missing genotypes and RADseq data

SNPs below `min.call.rate` are excluded separately within each window.
For the remaining SNPs, each missing genotype is replaced temporarily by
that SNP's mean dosage before covariance PCA is calculated. Mean
imputation is equivalent to giving a missing sample the observed
allele-frequency expectation at that SNP: after centring, it contributes
zero rather than an invented homozygous or heterozygous deviation. The
GDS is not modified, and observed calls are retained for heterozygosity
calculations.

This simple imputation keeps all samples in a common PCA space, but it
does not correct non-random missingness. In RADseq and other
reduced-representation data, allele dropout at restriction-site
polymorphisms, uneven depth, library quality, lanes, plates,
populations, marker panels, or alignment quality may correlate with
biological groups. Mean imputation can then shrink affected samples
toward the window centre, alter covariance, weaken a real signal, or
create a batch-associated local signal. Before interpreting candidates,
plot call rate and depth against chromosome position, compare PCA
clusters with batch and plate metadata, repeat the scan at stricter
`min.call.rate` values, and confirm that candidate windows remain after
removing problematic samples or markers. Imputation makes the matrix
computable; it does not make missing data unbiased.

LD uses observed genotypes with pairwise-complete correlations; missing
LD genotypes are not mean-imputed.

## Candidate evidence summary

Regional k-means clustering is treated as a hypothesis, not as evidence
by itself. `three_cluster_evidence` requires three groups with at least
three samples each, a smallest-cluster frequency of at least 0.05, and a
minimum adjacent-centre separation of one pooled within-cluster standard
deviation. The candidate table also reports cluster compactness, PC1
variance, heterozygosity excess in the middle cluster, LD within
inferred arrangement groups, LD in flanking windows, boundary contrasts,
and the largest internal score transition.

`evidence_score` is a transparent screening heuristic from zero to five.
One point is assigned for quantitative three-cluster support, positive
middle- cluster heterozygosity excess, a positive candidate-to-flank
score contrast, regional LD above flanking LD, and continuity across at
least two windows. Scores of 0–2 are labelled `weak`, 3–4 `moderate`,
and 5 `strong`. These labels prioritise review; they do not convert a
candidate into a structurally confirmed inversion. Known-region overlaps
are reported separately and do not increase or decrease the evidence
score.

Every candidate is described conservatively. A local-PCA signal can
reflect a putative inversion-associated haploblock, but it can also
arise near a centromere, in a region of low recombination, from assembly
or mapping problems, introgression, population-specific missingness, or
another form of structural variation. The `candidate_class` and
`alternative_explanations` columns make these alternatives explicit. The
function never reports a structurally confirmed inversion.

## Arrangement genotypes and sensitivity datasets

For a three-cluster regional PCA, clusters are ordered along PC1 and
labelled `AA`, `AB`, and `BB` by default. These are putative arrangement
genotypes, not sequence-level breakpoint genotypes. The individual table
includes the arrangement call, its numeric dosage (0, 1, or 2), and a
relative assignment-confidence score.

For every candidate, the function writes whitelists for each arrangement
and a homokaryotype-only whitelist containing `AA` and `BB`. It also
returns a combined whitelist for individuals classified as an outer
arrangement in every candidate region. These datasets support analyses
that exclude putative heterokaryotypes. This is useful before haplotype
scans: `selscan`, for example, does not accept missing genotype or
haplotype data.

## Output and plotting

Following other `radr` detection functions, each call creates a dated
`detect_inversions` results folder in the working directory (or below
the parent supplied with `path.folder`). It records the function
arguments, window and candidate tables, individual PCA scores, cluster
summaries, LD summaries, and standard diagnostic plots. PNG and PDF are
written by default.

## References

Li H, Ralph P (2019). Local PCA shows how the effect of population
structure differs along the genome. Genetics, 211, 289-304.
[doi:10.1534/genetics.118.301747](https://doi.org/10.1534/genetics.118.301747)
.

Faria R, Johannesson K, Butlin RK, Westram AM (2019). Evolving
inversions. Trends in Ecology & Evolution, 34, 239-248.
[doi:10.1016/j.tree.2018.12.005](https://doi.org/10.1016/j.tree.2018.12.005)
.

Wellenreuther M, Bernatchez L (2018). Eco-evolutionary genomics of
chromosomal inversions. Trends in Ecology & Evolution, 33, 427-440.
[doi:10.1016/j.tree.2018.04.002](https://doi.org/10.1016/j.tree.2018.04.002)
.

## Author

Thierry Gosselin <Thierry.Gosselin@csiro.au>
