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
  strata = NULL,
  chromosome = NULL,
  chromosome.lengths = NULL,
  reference.genome = NULL,
  window.snps = 100L,
  step.snps = window.snps,
  window.bp = NULL,
  step.bp = window.bp,
  window.method = c("snps", "bp", "ld"),
  ld.window.threshold = 0.1,
  ld.window.min.snps = 50L,
  ld.window.max.snps = 500L,
  sensitivity.window.snps = NULL,
  n.pcs = 2L,
  mds.axes = 2L,
  outlier.quantile = 0.99,
  min.window.snps = max(10L, n.pcs + 2L),
  min.call.rate = 0.8,
  min.candidate.windows = 1L,
  cluster.k = 3L,
  stability.replicates = 0L,
  stability.fraction = 0.8,
  parallel.core = 1L,
  chromosome.pca = TRUE,
  chromosome.pca.max.snps = 2000L,
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

- strata:

  Optional sample metadata supplied as a data frame or a tabular file.
  It must contain `INDIVIDUALS`; all other columns are retained. Rows
  act as a sample whitelist, and metadata such as `STRATA`, sequencing
  batch, library, lane, plate, or caller can be compared with regional
  PCA and putative arrangement assignments. The strata file remains
  sample metadata; it is not copied into or used to modify the GDS.

- chromosome:

  Optional chromosome or scaffold names to scan. By default, all
  chromosomes represented by at least one complete window are scanned.

- chromosome.lengths:

  Optional chromosome-length information supplied as a named numeric
  vector, a data frame with `CHROM` and `LENGTH` columns, or a tabular
  filepath containing those columns. Explicit values override lengths
  stored in the GDS.

- reference.genome:

  Optional reference FASTA or FASTA-index (`.fai`) filepath used to
  recover chromosome lengths. When a FASTA is supplied, its accompanying
  `FASTA.fai` file must already exist. This is unnecessary when the GDS
  retains a VCF contig dictionary with declared sequence lengths.

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
  replacing the primary candidate scan. This additional work is opt-in;
  the default `NULL` skips it. A useful focused set is
  `c(50, 100, 250, 500)`.

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

- stability.replicates:

  Number of SNP-resampling replicates used to test
  arrangement-assignment stability within every candidate. Use `0` to
  skip. Each replicate samples regional SNPs with replacement, repeats
  PCA and clustering, and compares assignments after matching cluster
  order on PC1.

- stability.fraction:

  Fraction of regional SNPs sampled in each stability replicate.
  Default: `stability.fraction = 0.8`.

- parallel.core:

  Number of independent R workers used for the local-PCA window scan.
  The default `1` is memory-conscious. Values greater than one open a
  separate read-only GDS connection in each worker, avoiding concurrent
  use of one connection. Each worker holds its own sample covariance
  matrix, so increase this value gradually on large datasets.

- chromosome.pca:

  Logical indicating whether an independent PCA is calculated for every
  chromosome or linkage group and coloured afterward with each candidate
  region's arrangement assignments. This diagnoses whether the same
  groups separate locally or throughout the genome.

- chromosome.pca.max.snps:

  Maximum number of evenly distributed SNPs used for each chromosome
  PCA. This bounds computation while retaining chromosome-wide coverage.
  Default: `2000`.

- arrangement.labels:

  Three labels, ordered from the lowest to highest regional PC1 cluster,
  used when `cluster.k = 3` and quantitative three-cluster support is
  present. Default: `arrangement.labels = c("AA", "AB", "BB")`.

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

- `candidate.summaries`: concise narrative summaries of the evidence,
  arrangement groups, regional LD, and boundary caution for each
  candidate;

- `diagnostics`: regional PCA scores, cluster assignments,
  heterozygosity, coverage, allele balance, cluster-number comparison,
  assignment stability, sample-metadata audit, marker loadings,
  arrangement differentiation, and optional LD matrices;

- `arrangement.genotypes`: one row per individual and candidate with
  putative arrangement genotype and relative assignment confidence;

- `homokaryotype.whitelist`: candidate-specific `AA` and `BB`
  individuals, plus `homokaryotype.all.candidates` for the intersection;

- `sensitivity`: optional summaries for additional fixed-SNP window
  sizes;

- `chromosome.pca`: independent chromosome or linkage-group PCA scores
  and summaries. Candidate-region arrangement labels are joined only for
  plotting and never influence these chromosome-wide PCAs;

- `chromosome.lengths`: declared or estimated chromosome lengths,
  observed marker maxima, and the provenance of each length;

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

## Chromosome length and candidate extent

Candidate tables report physical span, chromosome length, chromosome
fraction and percentage, and the left and right flanking lengths.
Declared chromosome lengths are preferred. The function first uses
explicit `chromosome.lengths`, then `reference.genome`, then the VCF
contig dictionary retained in the GDS. Chromosome names such as `14` and
`chr14` are matched. If no declared length is available, the largest
observed marker position is used as a clearly labelled underestimate of
chromosome length; the resulting candidate percentage may therefore be
overestimated.

## Recommended staged workflow

Begin with the default genome-wide scan on a carefully filtered dataset,
before LD pruning. Use the default result to locate unusual chromosomes
and inspect chromosome-wide PCA, regional PCA, missingness, depth,
heterozygosity, LD, cluster support, and sample metadata.

When a candidate is found, rerun its chromosome separately. Compare
several primary `window.snps` and `step.snps` combinations, and
optionally physical or LD-scaled windows. Record the candidate start,
end, span, chromosome percentage, and interval overlap for every run.
`sensitivity.window.snps` provides complementary PC1-variance and LD
summaries at additional scales, but it does not formally recall
candidate boundaries. Boundary stability therefore requires separate
focused calls with different primary window settings. Repeat promising
candidates under stricter call-rate, sample, relatedness, and batch
filters before seeking linkage, recombination, long-read, read-pair,
assembly, or breakpoint confirmation.

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

## Sample metadata and technical confounders

`strata` may contain `STRATA` and any additional sample-level variables
such as sequencing batch, library, lane, plate, extraction method,
genotype caller, or sampling year. Its rows first act as a whitelist.
For every candidate, the function then measures the association of each
metadata variable with regional PC1 and with arrangement assignments. A
strong association does not automatically reject a candidate, but it
identifies a biological or technical alternative that must be checked.

Candidate outputs also compare call rate, mean read depth, and
heterozygote allele balance among putative arrangements when those
quantities are stored in the GDS. Standard `DP` and biallelic `AD` nodes
are recognized, as are the genometranslator genotype-metadata nodes
`READ_DEPTH`, `ALLELE_REF_DEPTH`, and `ALLELE_ALT_DEPTH`. Consequently,
these diagnostics can be available for a VCF, DArT two-row count file,
or another source that retained read counts. Depth or allele-balance
shifts can indicate mapping artifacts, paralogy, copy-number variation,
or other structural variation. They should not be interpreted
automatically as inversion support.

## Cluster number and assignment stability

The requested `cluster.k` controls the arrangement calls, but the
function also compares one-, two-, and three-cluster descriptions of
regional PC1 using an approximate Gaussian BIC. This prevents three
clusters from being accepted merely because three centres were
requested.

With `stability.replicates > 0`, regional SNPs are resampled and
arrangement calls are repeated. The individual and overall agreement
values quantify how dependent the calls are on the exact SNP set. PC
direction is arbitrary, so reversed cluster labels are aligned before
agreement is measured.

## Chromosome-wide PCA context

When `chromosome.pca = TRUE`, the function calculates an independent PCA
for every chromosome or linkage group, using at most
`chromosome.pca.max.snps` evenly distributed SNPs. The arrangement
labels inferred from each candidate interval are then used only to
colour these chromosome-wide PCA panels. They do not influence the PCA
or force three groups on any chromosome.

Separation concentrated on the candidate linkage group supports a
localized haploblock interpretation. Similar separation across many
linkage groups instead suggests genome-wide population structure,
relatedness, admixture, or technical confounding. Absence of separation
in a whole-linkage-group PCA does not reject a shorter candidate:
unrelated SNPs elsewhere on the linkage group can dilute a strong
regional signal. Use this overview with the local window scan and the
candidate-specific technical diagnostics.

## Arrangement differentiation and marker loadings

Regional PC1 loadings are returned for every SNP and ranked by absolute
magnitude. Pairwise allele-frequency differences and Hudson's FST are
also reported between inferred arrangement classes. These summaries
describe the regional genetic contrast and help select diagnostic
markers. They do not make the arrangement classes populations.

Dxy is deliberately not calculated between `AA`, `AB`, and `BB`. In
particular, `AB` is an inferred heterokaryotype rather than an
independently sampled lineage. Population-level Dxy may be appropriate
later when genuine populations or evolutionary lineages are compared.

## Arrangement genotypes and sensitivity datasets

For a three-cluster regional PCA, clusters are ordered along PC1 and
labelled `AA`, `AB`, and `BB` by default when quantitative three-cluster
support is present. `AA` and `BB` identify the two outer PC1 clusters;
neither label establishes the reference, ancestral, derived, or
physically inverted arrangement. `AB` is the intermediate cluster and is
interpreted as a putative heterokaryotype. Candidates without
quantitative three-cluster support retain their numeric cluster IDs but
use neutral `Group 1`, `Group 2`, and `Group 3` labels and do not
receive arrangement dosages.

These are putative arrangement genotypes, not sequence-level breakpoint
genotypes. For supported candidates, the individual table includes the
arrangement call, its numeric dosage (0, 1, or 2), the original numeric
cluster ID, and a relative assignment-confidence score.

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

Hudson RR, Slatkin M, Maddison WP (1992). Estimation of levels of gene
flow from DNA sequence data. Genetics, 132, 583-589.

Bhatia G, Patterson N, Sankararaman S, Price AL (2013). Estimating and
interpreting FST: the impact of rare variants. Genome Research, 23,
1514-1521.
[doi:10.1101/gr.154831.113](https://doi.org/10.1101/gr.154831.113) .

## Author

Thierry Gosselin <Thierry.Gosselin@csiro.au>
