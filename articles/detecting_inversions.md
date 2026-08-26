# Detecting candidate chromosomal inversions with radr

## Scope

[`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md)
is a screening tool for finding genomic regions whose local population
structure is unusual relative to the rest of a GDS dataset. It is
designed for diploid, biallelic SNP data, including carefully filtered
RADseq and similar reduced-representation datasets.

The function detects **candidate inversion-associated haploblocks**. SNP
data can show the expected population-genomic footprint of an inversion,
but they do not directly demonstrate reversed chromosome orientation or
locate physical breakpoints. Use long reads, split reads, discordant
read pairs, linkage maps, comparative assemblies, or cytogenetics for
structural confirmation.

## How `radr` builds on previous inversion work

[`radr::detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md)
builds on analytical principles explored by Li and Ralph (2019), who
formalised local PCA for detecting changes in population structure along
a genome; Huang et al. (2020), who demonstrated that candidate
inversions can be recovered from reduced-representation SNP data; and
Pearse et al. (2019) and Akopyan et al. (2025), who showed the value of
combining regional genotype groups, linkage disequilibrium,
recombination, diversity, and structural evidence when interpreting
inversion-associated haploblocks.

The function is an independent, implementation rather than a wrapper
around another inversion package. It was designed to make these
transferable ideas practical for the quality-controlled RADseq and
similar datasets already used in `radr`. In particular,
[`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md):

- reads genotypes directly from GDS;
- constructs every window within a chromosome, linkage group, or
  scaffold;
- combines local covariance PCA with contiguous-window candidate
  detection;
- evaluates regional PCA groups quantitatively instead of treating a
  requested three-group k-means solution as evidence by itself;
- summarises heterozygosity, overall and within-group LD, candidate
  boundaries, flanking windows, internal transitions, and window-size
  sensitivity;
- accepts centromere, low-recombination, repeat, and assembly-gap
  annotations without allowing them to determine the candidate calls;
- explicitly records the consequences of missing-data filtering and
  temporary mean imputation for RADseq interpretation; and
- writes reproducible tables and publication-ready `ggplot2` diagnostics
  to a standard results folder while returning all results for further
  exploration.

This integrated workflow is the main advantage of using `radr`: it moves
from a broad local-structure scan to an auditable candidate evidence
table without claiming more resolution than the marker data provide. The
output is intended to guide the next experiment, not replace it. Linkage
mapping can test recombination suppression; haplotagging or other
linked-read approaches can phase the alternative haplotypes and refine
structural hypotheses; and ONT or PacBio long reads, breakpoint PCR,
comparative assemblies, or cytogenetics can provide physical
confirmation. Haplotagging has recovered long population haplotypes and
inversion-associated barcode-sharing patterns at scale (Meier et
al. 2021), but exact or repetitive breakpoints may still require
continuous long reads or a junction-specific assay.

## What is a chromosomal inversion?

An inversion occurs when a chromosome segment is reversed. Individuals
may carry two copies of one arrangement, one copy of each arrangement,
or two copies of the alternative arrangement. These are often described
as two homokaryotypes and one heterokaryotype.

Recombination is commonly reduced between alternative arrangements in
heterokaryotypes. Consequently, alleles across a large interval can
remain associated as a haplotype. This can protect combinations of
locally adaptive alleles from being broken apart by recombination,
especially when populations exchange migrants. Inversions can therefore
maintain genomic differentiation even when much of the genome shows weak
population structure.

Not every inversion is adaptive, and not every long differentiated
haploblock is an inversion. Selection, centromeres, low-recombination
regions, recent admixture, family structure, paralogs, assembly errors,
marker-density changes, and technical batches can generate partially
similar patterns.

## Evidence expected from SNP data

A convincing candidate normally combines several lines of evidence:

1.  **Local structure:** contiguous windows show a population structure
    that differs from the genomic background.
2.  **Three genotype groups:** regional PCA may show two homokaryotype
    groups and an intermediate heterokaryotype group.
3.  **Heterozygosity:** the intermediate group is often more
    heterozygous within the region than the two outer groups.
4.  **Extended LD:** SNPs across the region show elevated linkage
    disequilibrium.
5.  **Genomic continuity:** the signal spans adjacent markers and does
    not merely follow an isolated outlier SNP.
6.  **Technical independence:** the signal does not track plates, lanes,
    libraries, missingness, coverage, or other processing variables.
7.  **Population differentiation:** when biologically defined
    populations differ in arrangement frequency, a broad windowed
    $`F_{ST}`$ peak overlaps the local structure signal and helps
    prioritise the region for follow-up.

None of these criteria alone proves an inversion. For example, LD may be
uneven inside an old inversion because gene conversion and double
crossovers permit some gene flux between arrangements.

## Lessons from empirical fish examples

Published fish studies illustrate both the value and the limits of an
inversion screen.
[`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md)
does not reproduce their complete analyses: several relied on
whole-genome sequence, linkage maps, pedigrees, long reads, or phenotype
data that are not present in a typical RADseq GDS. Instead, the function
uses the parts of their reasoning that transfer to marker data and
reports candidates for further investigation.

### Rainbow trout: a structurally validated double inversion

Pearse et al. (2019) characterized a roughly 55-Mb double-inversion
supergene on rainbow trout chromosome Omy05. Linkage mapping showed
almost complete recombination suppression in heterokaryotypic parents,
while comparative genomics and long-read assembly helped resolve the
structure. The two arrangements were associated with migratory tendency,
with effects that depended on sex and dominance.

This example motivates several choices in
[`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md):

- evidence should extend across a contiguous haploblock rather than a
  single strongly associated marker;
- regional genotype groups, LD, heterozygosity, and boundary changes
  should be interpreted together;
- one broad signal may contain adjacent or nested rearrangements and
  should not automatically be described as one simple inversion; and
- association with sex, migration, or another phenotype can support
  biological relevance, but it does not establish chromosome
  orientation.

The outer candidate coordinates returned by
[`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md)
therefore describe a marker-supported interval. Internal changes in
window scores, LD, or clustering may justify examining more than one
rearrangement inside that interval.

### Chinook salmon: a major-effect locus is not necessarily an inversion

Thompson et al. (2020) showed that variation near `GREB1L` is strongly
associated with adult migration timing in Chinook salmon. This is an
important counterexample for interpreting a genomic peak: a locus with a
large phenotypic effect, strong differentiation, or extended haplotypes
is not by itself evidence of a chromosomal inversion.

Consequently,
[`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md)
does not use phenotype association or an isolated differentiation peak
as sufficient evidence. A narrow selected region may be biologically
important without producing the extended local-PCA, three-group,
heterozygosity, and LD pattern expected from an inversion haploblock.

### Atlantic silverside: distinguish inversions from centromeres

Akopyan et al. (2025) compared 168 Atlantic silverside genomes from four
populations using whole-genome variation and recombination maps. Large,
abruptly bounded differentiation haploblocks coincided with known
inversions and showed the characteristic three tight regional-PCA
clusters. Narrower differentiation peaks frequently coincided with
putative centromeres. Regional PCA in centromeric regions could also
reflect reduced recombination, but the individuals were more dispersed
and retained more haplotype variation than in the inversion regions.

The study also demonstrates why relative differentiation alone can
mislead. Elevated $`F_{ST}`$ near centromeres resulted partly from low
within-population diversity, whereas absolute sequence divergence
($`d_{XY}`$) was elevated in the large inversions but reduced near
centromeres. For `radr`, the transferable diagnostic lessons are to:

1.  compare broad, abrupt and contiguous haploblocks with narrow or
    gradual peaks;
2.  examine the compactness and separation of the three inferred
    genotype groups, not merely the existence of a regional PCA pattern;
3.  compare candidates with known or putative centromeres, recombination
    maps, marker density, repeats, and assembly gaps; and
4.  avoid calling a low-recombination region an inversion without
    independent structural or linkage evidence.

Sparse, ascertained RADseq markers generally do not support the same
robust windowed $`d_{XY}`$ analysis as whole-genome sequence.
[`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md)
therefore does not manufacture an absolute-divergence statistic from
insufficient data. Where dense sequence data and defensible population
groups are available, $`F_{ST}`$, diversity, and $`d_{XY}`$ provide
valuable downstream validation alongside the `radr` candidate scan.

### Cross-reference candidates with chromosome-resolved $`F_{ST}`$

Population differentiation is most informative when it is plotted along
each linkage group rather than reduced to one genome-wide value. The
complementary `assigner::fst_WC84()` function calculates Weir and
Cockerham’s (1984) $`F_{ST}`$ from independently defined populations
and, when `CHROM` and numeric `POS` metadata are available, returns:

- `fst.linkage.groups`, with one summary per chromosome or linkage
  group;
- `fst.windows`, with WC84 variance components combined in physical
  windows; and
- `fst.genome.plot`, a chromosome-concatenated Manhattan-style figure.

The physical window size should reflect marker density and the expected
scale of the candidate. For sparse RADseq data, a 25-kb window may
contain too few markers; larger windows and a sensitivity analysis
across several window sizes are usually more defensible. Use the same
genome assembly, chromosome labels, position metadata, sample filters,
and marker filters in both analyses.

``` r

library(ggplot2)

strata <- readr::read_tsv("strata.tsv", show_col_types = FALSE)

fst <- assigner::fst_WC84(
  data = genome,
  strata = strata,
  linkage.group = TRUE,
  window.size = 1e6,
  window.method = "fixed",
  filename = "inversion_fst_crosscheck"
)

# Whole-genome overview supplied by assigner
fst$fst.genome.plot

# Linkage-group view with radr candidate intervals shaded
candidate_intervals <- inv$candidates |>
  dplyr::transmute(
    CHROM = as.character(chromosome),
    xmin = start,
    xmax = end
  )

ggplot(fst$fst.windows, aes(WINDOW_MID, FST_WC84)) +
  geom_rect(
    data = candidate_intervals,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#E69F00",
    alpha = 0.18
  ) +
  geom_line(colour = "grey45") +
  geom_point(aes(colour = OUTLIER), size = 1.4) +
  facet_wrap(~ CHROM, scales = "free_x") +
  scale_x_continuous(
    labels = scales::label_number(scale = 1e-6, suffix = " Mb")
  ) +
  labs(
    x = "Genomic position",
    y = expression("Windowed WC84 " * F[ST]),
    colour = "Exploratory outlier"
  ) +
  theme_bw()
```

An overlapping, spatially coherent $`F_{ST}`$ peak strengthens the case
that a candidate haploblock contributes to differentiation among the
sampled populations. Its absence does not reject an inversion: both
arrangements may occur at similar frequencies among populations.
Conversely, an $`F_{ST}`$ peak can arise from local selection, a
centromere, low recombination, reduced diversity, marker ascertainment,
or technical structure without an inversion.

Population strata must be specified independently of the inversion scan.
Calculating $`F_{ST}`$ between the regional PCA clusters inferred from
the same candidate genotypes is circular and should not be presented as
independent support. Where whole-genome sequence is available, compare
relative differentiation with $`d_{XY}`$, within-population diversity,
recombination, and structural evidence, as illustrated by Akopyan et
al. (2025).

### Rainbow trout across populations: inversion signals are context-dependent

Campbell, Anderson, Garza, and Pearse (2021) compared independently
landlocked rainbow trout populations with an anadromous source
population. All landlocked populations showed an increased frequency of
the large Omy05 inversion, while three of four showed an increased
frequency of the Omy20 inversion. Outside these major rearrangements,
however, the genomic response to similar selection was much less
parallel.

This work, including the population-genomic contribution of Eric C.
Anderson, reinforces two principles used in
[`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md).
First, a long inversion-associated haplotype can represent standing
adaptive variation that changes frequency repeatedly across populations.
Second, neither the strength of a genomic signal nor its association
with migration should be assumed to be constant across the species
range. Candidate clusters should therefore be cross-tabulated against
population and environment, and biological associations should be tested
rather than assigned from the candidate region alone.

## How the local analysis is organised

Windows are created **within** chromosome, linkage-group, or scaffold
boundaries. If LG1 ends with 40 unused SNPs and LG2 follows in the GDS,
those 40 SNPs are never combined with the first 60 SNPs of LG2. This
prevents an artificial window spanning unrelated genomic regions.

The analysis is not one PCA per linkage group. Instead:

1.  each linkage group is divided into SNP windows;
2.  a separate local covariance PCA is calculated for every window;
3.  each window is represented by its leading covariance axes;
4.  distances between window representations are calculated;
5.  multidimensional scaling (MDS) places similar windows together and
    unusual windows farther from the genomic background;
6.  contiguous unusual windows are joined into candidate regions;
7.  a regional PCA, clustering, heterozygosity summary, and LD summary
    are calculated for each candidate.

When several linkage groups are scanned, their window summaries are
compared in the same MDS analysis. The returned table retains the
linkage-group label, so chromosome-position plots should normally be
faceted by linkage group. A regional PCA belongs to one candidate on one
linkage group and therefore does not itself require an LG facet.

## Run a scan from GDS

The examples are not evaluated because they refer to user files.

``` r

library(radr)

genome <- genometranslator::read_genome("filtered_dataset.gds")

inv <- detect_inversions(
  data = genome,
  window.snps = 100,
  step.snps = 100,
  sensitivity.window.snps = c(100, 250, 500, 1000),
  min.call.rate = 0.90,
  min.candidate.windows = 2,
  return.ld = TRUE
)

inv
inv$candidates
inv$path.folder
```

Known centromeres, assembly gaps, or other low-recombination annotations
can be supplied without allowing those annotations to determine which
windows are selected:

``` r

known_regions <- data.frame(
  chromosome = c("LG5", "LG14", "LG14"),
  start = c(22000000, 18000000, 41000000),
  end = c(26000000, 23000000, 42500000),
  type = c("putative_centromere", "low_recombination", "assembly_gap")
)

inv <- detect_inversions(
  data = genome,
  window.snps = 100,
  step.snps = 50,
  known.regions = known_regions
)
```

The `known_region_overlap` column lists overlapping annotation types. An
overlap is a warning for interpretation, not evidence for or against an
inversion, and it does not change candidate selection or the evidence
score.

To investigate a known chromosome-14 signal without letting other
linkage groups define the background:

``` r

chr14 <- detect_inversions(
  data = genome,
  chromosome = "14",
  window.snps = 100,
  step.snps = 50,
  min.call.rate = 0.90,
  min.candidate.windows = 2
)
```

Overlapping windows (`step.snps < window.snps`) can help localise
transitions, but adjacent results are then strongly dependent.
Coordinates should not be reported as precise breakpoints merely because
a score changes between two overlapping windows.

For physical windows instead of equal-SNP windows, supply `window.bp`.
This is particularly useful for checking whether a candidate is robust
to variation in RAD marker density:

``` r

chr14_mb <- detect_inversions(
  data = genome,
  chromosome = "14",
  window.bp = 1e6,
  step.bp = 5e5,
  min.window.snps = 20
)
```

## Plot window evidence by linkage group

[`detect_inversions()`](https://thierrygosselin.github.io/radr/reference/detect_inversions.md)
automatically saves the standard window-score, MDS, LD, call-rate,
regional-PCA, heterozygosity, and candidate-LD figures. It also returns
the `ggplot2` objects in `inv$output.files$plots`, so they can be
changed without rerunning the genomic calculations.

``` r

library(ggplot2)

ggplot(
  inv$windows,
  aes(x = (start + end) / 2, y = robust_score,
      colour = candidate_window)
) +
  geom_point() +
  geom_hline(
    yintercept = inv$settings$score.threshold,
    linetype = 2
  ) +
  facet_wrap(~ chromosome, scales = "free_x") +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6,
                                                    suffix = " Mb")) +
  labs(
    x = "Window midpoint",
    y = "Robust local-structure score",
    colour = "Candidate"
  ) +
  theme_bw()
```

This is the appropriate LG-faceted view: window score versus genomic
position. The MDS coordinates can also be inspected, but faceting them
by linkage group answers a different question because MDS axes represent
similarity among windows rather than physical position.

``` r

ggplot(inv$windows, aes(MDS1, MDS2, colour = chromosome,
                        shape = candidate_window)) +
  geom_point(size = 2) +
  theme_bw()
```

## Inspect candidate genotype groups

Each candidate has a matching entry in `inv$diagnostics`. A regional PCA
plot shows individuals, not windows:

``` r

candidate1 <- inv$diagnostics[[1]]

ggplot(candidate1$scores, aes(PC1, PC2, colour = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(colour = "Putative arrangement genotype") +
  theme_bw()

candidate1$cluster_summary
```

Three clusters are biologically suggestive only when the central PCA
cluster also has the expected heterozygosity and the pattern is not
explained by population, family, or batch. Two clusters may occur when
one arrangement is rare or absent, while continuous PCA scores may
indicate ordinary population structure rather than a polymorphic
inversion.

The LD diagnostic uses the upper triangle for all individuals and the
lower triangle for the more common outer PCA cluster, treated
provisionally as a homokaryotype. This comparison is useful because
mixing arrangements can create strong LD even when LD is lower within
one arrangement. It is not independent validation: the subset was
inferred from the same candidate-region genotypes. The group-specific
values are available in `candidate1$ld_summary`.

## Mean imputation: what it does

Within each window, SNPs below `min.call.rate` are excluded. At each
remaining SNP, a missing genotype is temporarily replaced by the
observed mean allele dosage for that SNP. The value is used only to make
covariance PCA matrices complete; the GDS is not changed. After marker
centring, an imputed value contributes zero deviation at that SNP. It
therefore supplies no evidence that the individual carries either
arrangement.

Only covariance PCA uses these imputed values. LD is calculated from
observed genotypes using pairwise-complete correlations; missing LD
genotypes are not mean-imputed.

This is intentionally simple and transparent. It prevents PCA from
discarding every individual with one missing call, but it is not a
correction for biased missingness.

## Why RADseq missingness needs special attention

RADseq missingness is often structured rather than random. Relevant
causes include:

- polymorphism at restriction sites and consequent allele dropout;
- uneven sequencing depth among samples or libraries;
- lane, plate, library-preparation, or genotyping batches;
- population differences in reference divergence or mapping quality;
- separate marker discovery or genotype calling among sample groups;
- paralogous loci and collapsed repetitive regions;
- marker density that varies with restriction sites, assembly quality,
  or filtering.

If one population or batch is missing preferentially at a set of linked
SNPs, mean imputation pulls those samples toward the local PCA centre.
Depending on the pattern, this can weaken a biological cluster,
manufacture an intermediate group, or make a batch-specific genomic
region look unusual.

For every candidate, perform sensitivity analyses:

1.  plot marker and individual call rate, depth, and heterozygosity;
2.  cross-tabulate inferred clusters against population, plate, lane,
    library, family, and sampling date;
3.  repeat the scan with stricter `min.call.rate` values;
4.  repeat after excluding the most incomplete samples;
5.  compare fixed-SNP and fixed-base-pair windows when marker density is
    uneven;
6.  inspect whether the signal remains when close relatives are removed;
7.  verify that marker density and assembly gaps do not define the
    candidate boundaries.

A robust biological signal should not disappear under every reasonable
QC choice, and its clusters should not be synonymous with one technical
batch.

## Results folder and customised plots

Each call creates a dated `detect_inversions` folder following the
standard `radr` workflow. It contains the recorded call, TSV tables, and
PNG and PDF figures. Timestamped folders prevent a normal run from
overwriting an earlier analysis. The returned object also contains the
window table, candidate table, PCA scores, cluster summaries, LD
summaries, and plot objects. A standard plot can therefore be adjusted
and saved under an additional name:

``` r

p <- inv$output.files$plots$window_scores +
  labs(title = "Candidate inversion scan after batch QC")

ggsave(file.path(inv$path.folder, "window_scores_annotated.png"), p,
       width = 10, height = 7, dpi = 300)
```

Use `save.plots = FALSE` when only tables and returned objects are
wanted. The results folder and reproducibility tables are still created.

## Interpreting coordinates

The candidate `start` and `end` values are the outer marker coordinates
of the selected windows. They depend on window size, step size, marker
density, missingness filters, and the reference assembly. Report them as
an approximate inversion-associated interval. With RADseq data, there
may be a substantial gap between the last marker outside the signal and
the first marker inside it.

A defensible report might state:

> A candidate inversion-associated haploblock was detected from
> approximately X to Y Mb on LG14, supported by concordant local-PCA,
> regional genotype-group, heterozygosity, and LD patterns. These
> coordinates describe the marker-defined interval and not validated
> structural breakpoints.

## Interpreting the candidate evidence table

The candidate table combines continuous diagnostics with a deliberately
simple screening grade. Important columns include:

- `cluster_separation`: the smallest gap between adjacent regional-PC1
  cluster centres, divided by the pooled within-cluster standard
  deviation;
- `cluster_compactness`: the proportion of regional PC1 variation
  explained by the inferred clusters;
- `smallest_cluster_n` and `smallest_cluster_frequency`: protection
  against treating a few outlying individuals as an inversion
  arrangement;
- `middle_heterozygosity_excess`: middle-cluster heterozygosity minus
  the mean of the two outer clusters;
- `regional_mean_ld_r2`, `homokaryotype_mean_ld_r2`, and
  `heterokaryotype_mean_ld_r2`: LD across all individuals and within
  inferred genotype groups;
- `flanking_mean_ld_r2` and `boundary_contrast`: comparison with the
  immediate noncandidate windows;
- `internal_transition_max`: the largest score change between candidate
  windows, useful for finding complex, adjacent, or nested signals; and
- `known_region_overlap`: any user-supplied centromere,
  low-recombination, assembly-gap, or other annotation intersecting the
  candidate.

`three_cluster_evidence` is not merely a record that k-means was run
with `cluster.k = 3`. It additionally requires all three inferred groups
to contain at least three samples, the smallest group to represent at
least 5% of samples, and adequate separation relative to within-group
spread.

`evidence_score` assigns one point for each of five observations:
quantitative three-cluster support, positive middle-cluster
heterozygosity excess, positive boundary contrast, regional LD above
flanking LD, and continuity across at least two windows. Scores of 0–2,
3–4, and 5 are labelled `weak`, `moderate`, and `strong`, respectively.
This grade prioritises candidates for inspection; it is not a posterior
probability and does not prove an inversion.

## Recommended workflow

Use the scan as one part of a staged analysis:

1.  complete sample, marker, batch, depth, and missingness QC;
2.  scan chromosome-specific windows for unusual local structure;
3.  inspect contiguous candidates rather than isolated outlier windows;
4.  evaluate regional PCA groups, heterozygosity, and LD jointly;
5.  test sensitivity to filtering, window size, step size, and
    relatedness;
6.  compare candidates with centromeres, assembly gaps, repeats, and
    gene annotations;
7.  use linkage maps or haplotagging to test recombination and phase the
    alternative arrangements; and
8.  seek breakpoint-spanning long reads, junction PCR, comparative
    assembly, or cytogenetic evidence before calling a structurally
    confirmed inversion.

## References

Akopyan M, Tigano A, Jacobs A, Wilder AP, Therkildsen NO (2025) Genetic
differentiation is constrained to chromosomal inversions and putative
centromeres in locally adapted populations with higher gene flow.
*Molecular Biology and Evolution*, 42, msaf092.

Campbell MA, Anderson EC, Garza JC, Pearse DE (2021) Polygenic basis and
the role of genome duplication in adaptation to similar selective
environments. *Journal of Heredity*, 112, 614-625.

Faria R, Johannesson K, Butlin RK, Westram AM (2019) Evolving
inversions. *Trends in Ecology & Evolution*, 34, 239-248.

Gosselin T, Anderson EC, Bradbury I (2020) assigner: assignment analysis
with GBS/RAD data using R. R package.
<https://doi.org/10.5281/zenodo.592677>.

Hoffmann AA, Rieseberg LH (2008) Revisiting the impact of inversions in
evolution: from population genetic markers to drivers of adaptive shifts
and speciation? *Annual Review of Ecology, Evolution, and Systematics*,
39, 21-42.

Huang K, Andrew RL, Owens GL, Ostevik KL, Rieseberg LH (2020) Multiple
chromosomal inversions contribute to adaptive divergence of a dune
sunflower ecotype. *Molecular Ecology*, 29, 2535-2549.

Kess T, Bentzen P, Lehnert SJ, Sylvester EVA, Lien S, Kent MP,
Sinclair-Waters M, Morris CJ, Regular P, Fairweather R, Bradbury IR
(2019) A migration-associated supergene reveals loss of biocomplexity in
Atlantic cod. *Science Advances*, 5, eaav2461.

Li H, Ralph P (2019) Local PCA shows how the effect of population
structure differs along the genome. *Genetics*, 211, 289-304.

Meier JI, Salazar PA, Kučka M, Davies RW, Dréau A, Aldás I, Power OB,
Nadeau NJ, Bridle JR, Rolian C, Barton NH, McMillan WO, Jiggins CD, Chan
YF (2021) Haplotype tagging reveals parallel formation of hybrid races
in two butterfly species. *Proceedings of the National Academy of
Sciences of the United States of America*, 118, e2015005118.

Mérot C (2020) Making the most of population genomic data to understand
the importance of chromosomal inversions for adaptation and speciation.
*Molecular Ecology*, 29, 2513-2516.

Pearse DE, Barson NJ, Nome T, Gao G, Campbell MA, Abadía-Cardoso A,
Anderson EC, Rundio DE, Williams TH, Naish KA, Moen T, Liu S, Kent M,
Moser M, Minkley DR, Rondeau EB, Brieuc MSO, Sandve SR, Miller MR,
Cedillo L, Baruch K, Hernandez AG, Ben-Zvi G, Shem-Tov D, Barad O,
Kuzishchin K, Garza JC, Lindley ST, Koop BF, Thorgaard GH, Palti Y, Lien
S (2019) Sex-dependent dominance maintains migration supergene in
rainbow trout. *Nature Ecology & Evolution*, 3, 1731-1742.

Thompson NF, Anderson EC, Clemento AJ, Campbell MA, Pearse DE, Hearsey
JW, Kinziger AP, Garza JC (2020) A complex phenotype in salmon
controlled by a simple change in migratory timing. *Science*, 370,
609-613.

Weir BS, Cockerham CC (1984) Estimating F-statistics for the analysis of
population structure. *Evolution*, 38, 1358-1370.

Wellenreuther M, Bernatchez L (2018) Eco-evolutionary genomics of
chromosomal inversions. *Trends in Ecology & Evolution*, 33, 427-440.
