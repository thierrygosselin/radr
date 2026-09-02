# Detecting and validating candidate sex-linked markers with radr

## Scope

Sex-linked markers can help identify the sex-determination system of a
species, describe sex-chromosome differentiation, prevent sex-linked
variation from confounding population-genetic analyses, and assign
genetic sex when external sex characters are absent or unreliable.

[`radr::sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md)
screens an already prepared GDS dataset for markers whose presence,
heterozygosity, or normalized read depth differs between recorded
females and males. It supports diploid, biallelic SNP datasets and
dominant SilicoDArT presence/absence markers stored with their marker
type in GDS. These include RADseq, DArTseq, and related
reduced-representation data. The function does not alter the GDS. It
respects the active sample and marker selections, restores them before
returning, and writes an auditable results folder.

**A candidate is not proof of a sex-chromosome location.** Sex-specific
dropout, plate or library effects, population structure, paralogy,
mapping bias, and errors in recorded sex can reproduce the expected
patterns. Discovery should be followed by technical, genomic, and
biological validation.

## Biological background

In an XX/XY system, females are the homogametic sex and males are the
heterogametic sex. In a ZZ/ZW system, males are homogametic and females
are heterogametic. Sex chromosomes can range from nearly homomorphic
chromosomes with a small non-recombining region to strongly
differentiated chromosomes with extensive sequence loss and divergence.
Consequently, no single marker pattern is expected in every species.

Reduced-representation data can reveal three complementary signals.

### Sex-specific presence

A sequence confined to the heterogametic chromosome should be present
mainly in one sex:

- a Y-like marker is expected mainly in males;
- a W-like marker is expected mainly in females.

These markers are especially useful for genetic sex assignment. However,
a similar presence/absence pattern can be produced by low read depth,
sex-confounded batches, allelic dropout, or a polymorphic restriction
site.

### Sex-biased heterozygosity

Markers in homologous or gametologous regions can exhibit different
genotype patterns between sexes. Under a simple model, an X-like marker
may show greater heterozygosity in females than males, whereas a Z-like
marker may show the reverse. The direction and magnitude depend on the
divergence and mapping behaviour of the chromosome copies; a biological
sex-linked marker does not have to satisfy this simple expectation.

### Sex-biased read depth

Copy-number differences can also produce depth differences. Under an
idealized model, an X-linked region can have greater normalized depth in
XX females than XY males, while a Z-linked region can have greater
normalized depth in ZZ males than ZW females. Library size, capture
efficiency, paralogy, and mapping bias can all affect this signal, so
[`sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md)
normalizes depth within samples and still treats the result as candidate
evidence.

## Sex-linked markers can create false population structure

Sex-linked markers matter even when sex determination is not the study
target. Benestan et al. (2017) demonstrated this with RADseq and GBS
data from American Lobster (*Homarus americanus*) and Arctic Char
(*Salvelinus alpinus*). Their multivariate analyses initially separated
individuals by sex rather than by the expected geographic or ecological
groups. When the sex ratio differed among sampling groups, a small set
of sex-linked markers created false differentiation in the panmictic
lobster comparison and inflated differentiation in Arctic Char.

The study used a complementary allele-frequency approach:

1.  DAPC was used to diagnose the unexpected male-female clustering.
2.  BayeScan treated females and males as the two groups and identified
    markers with unusually high differentiation between sexes.
3.  Sex ratios were deliberately varied among artificial sampling groups
    to measure their effect on FST.
4.  Candidate markers were removed in descending order of sex
    differentiation and FST was recalculated with
    `assigner::fst_WC84()`.
5.  LD, BLAST, transcriptome annotation, comparative genome placement,
    and a linkage map were used to investigate the candidate regions.

Only 12 of 1,717 SNPs in American Lobster (0.7%) and 94 of 6,147 SNPs in
Arctic Char (1.5%) were identified as sex outliers, yet they materially
changed the population-genetic conclusions. Removing 11 of 12 candidates
eliminated the false lobster differentiation in the most extreme
sampling scenario. In Arctic Char, FST reached a plateau after
approximately 80 candidates were removed and was more than threefold
smaller than before their removal. In the paper’s literature review,
only 5 of 52 marine and diadromous population-genomic studies reported
individual sex.

This result has direct consequences for a `radr` workflow:

- record and report sex whenever it is known;
- cross-tabulate sex by population, site, year, plate, and sequencing
  batch;
- colour PCA, DAPC, or other ordinations by sex before interpreting
  clusters;
- repeat population-genetic summaries after removing or stratifying
  candidate sex-linked markers; and
- do not assume that thousands of autosomal markers will dilute a small
  number of highly differentiated sex-linked markers in a weakly
  structured species.

The Benestan et al. (2017) strategy is valuable but tests a different
signal from
[`sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md).
BayeScan detects unusually large allele-frequency differences between
sex groups. It does not directly model sex-specific marker presence or
normalized read depth, and a sex outlier is not automatically on a sex
chromosome. Population structure, batch effects, family structure, or
unequal missingness can also differentiate the recorded sexes.

When sex is unknown, Benestan et al. (2017) suggested first identifying
genetic clusters unrelated to geography, treating those clusters as
groups in BayeScan, and inspecting the heterozygosity of the outliers.
This is an exploratory diagnostic, not a general replacement for
recorded sex. The authors noted that the heterozygosity pattern is
informative only when sex-linked regions are sufficiently
differentiated. It can miss a small or weakly differentiated
sex-determining region and can mistake another source of structure for
sex.

The study also illustrates good wet-laboratory design: Arctic Char
samples from each location were distributed across at least six
sequencing multiplexes to reduce confounding between sampling location
and library. Computational adjustment is much less convincing when sex,
population, and technical batch were never separated experimentally.

## Origin of the method

The original `sexy_markers` method was developed by Floriaan
Devloo-Delva, Thierry Gosselin, Robin B. Thomson, and collaborators and
described by Devloo-Delva et al. (2024). The study introduced a
genotyped-data workflow based on the three signals above and evaluated
it with 558 White Sharks (*Carcharodon carcharias*) and 23,393 SNPs. It
recovered nine Y-linked and 406 X-linked candidates, supported sex
assignment, and led to a PCR assay for independent validation.

The study also showed why study design matters. In that White Shark
dataset, subsamples with fewer than 100 individuals, fewer than 10,000
markers or half of the available markers, and unbalanced sex ratios
produced more false positives. These values are empirical results from
one species and dataset, not universal thresholds. They are nevertheless
a useful warning that a small, unbalanced discovery set should be
interpreted conservatively.

The current `radr` implementation is a modern reimplementation rather
than an exact copy of the published `radiator` function. It uses GDS
input only, reads markers in chunks, applies explicit effect-size rules,
calculates method-specific tests and false-discovery rates, audits
metadata for confounding with recorded sex, and produces a separate Y/W
panel for downstream assignment. Cite Devloo-Delva et al. (2024) for the
biological method and original validation, and record the `radr` version
used for every new analysis.

## Prepare the discovery dataset

### What data are best for discovery?

The best discovery dataset is not defined by one sequencing format. It
combines a study design in which recorded sex is reliable and separable
from population and technical batch with data that preserve the three
signals tested by
[`sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md):
genotype, marker presence, and read depth. Ideally, it also includes
marker sequences, chromosome-scale coordinates, and an independent set
of individuals reserved for validation.

Individual-level whole-genome sequencing provides the most complete
discovery data when sufficient depth and a suitable reference assembly
are available. It samples the genome more evenly than restriction-site
methods and is better suited to locating candidates, measuring copy
number, and resolving a sex-linked region. Its cost and computational
requirements can nevertheless limit sample size, and a small or
sex-confounded whole-genome dataset is not automatically better than a
well-balanced reduced-representation dataset.

RADseq, ddRAD, GBS, and DArTseq are cost-effective discovery data and
closely match the setting in which the original method was developed.
They can provide many individuals, genotypes, read depth, and short
marker sequences. Their limitations are sparse genomic sampling,
restriction-site ascertainment, allelic dropout, and the possibility
that no marker samples the sex-determining region. Preserve allele
counts or retained depth and marker sequences during conversion whenever
they are available.

SilicoDArT can be particularly informative for strongly differentiated
Y- or W-linked sequence. It records whether a restriction fragment is
observed and can therefore retain a sex-specific sequence that has no
shared SNP genotype, is absent from the reference assembly, or maps
poorly. It also works without a reference genome. Its limitation is
equally important: absence can reflect a true missing sequence,
insufficient depth, restriction-site polymorphism, library failure, or
allelic dropout. Dominant SilicoDArT observations cannot distinguish
heterozygotes from presence homozygotes and must not be interpreted with
a diploid heterozygosity model.

A DArT file should be imported together with its sample metadata so that
names, population, sex, and technical variables are embedded in the GDS
from the start. Keep sex in a dedicated column rather than replacing
population in `STRATA`:

``` r

silicodart_gds <- genometranslator::read_genome(
  data = "study.silicodart.csv",
  strata = "study.strata.tsv"
)

silicodart_result <- radr::sexy_markers(
  data = silicodart_gds,
  strata = NULL,
  sex.column = "SEX"
)
```

[`genometranslator::read_dart()`](https://thierrygosselin.github.io/genometranslator/reference/read_dart.html)
detects SilicoDArT automatically. It stores binary absence/presence in
the GDS, labels the markers with `MARKER_TYPE = SILICODART`, and
preserves original counts as read depth when a count-form SilicoDArT
report is supplied.
[`sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md)
consequently uses these markers for Y/W-like presence tests, disables
their heterozygosity tests, and does not assign X/Z-like genotype
classes from their synthetic GDS dosage. VCF SNP data remain richer for
codominant genotype, heterozygosity, LD, mapping, and regional analyses.
When both products are available, analyse them as complementary
discovery datasets and prioritize candidates that replicate across data
types or independent assays.

Other data have more specialized roles:

- a SNP array or targeted panel is usually better for validation than
  discovery because its markers were selected in advance and may exclude
  sex-specific or poorly mapping variation;
- RNA-seq can support transcript- or gametolog-based methods such as
  SDpop, but expression differences among tissues, stages, and sexes are
  additional sources of signal and confounding; and
- pooled sequencing can reveal sex-biased allele frequencies or
  coverage, but it does not retain the individual-level observations and
  metadata required by the current
  [`radr::sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md)
  implementation.

For SNP analysis in `radr`, the ideal input is a diploid, biallelic
SeqArray GDS with individual genotypes, retained per-sample depth or
allele counts, marker sequences, reference coordinates when available,
and complete sample metadata. Dominant SilicoDArT is the supported
exception and is identified explicitly in the GDS rather than presented
as a biological diploid genotype. GDS is the analysis format, not a
restriction on the original data source. `genometranslator` can
standardize VCF, DArT, PLINK BED/TPED, GENEPOP, FSTAT, tidy or wide
genotype tables, Parquet, and existing `genind`, `genlight`, and
`gtypes` objects before a GDS is written. Use a format-specific
`genometranslator::read_*()` function when non-default import choices
matter, then write the standardized data with
[`genometranslator::write_gds()`](https://thierrygosselin.github.io/genometranslator/reference/write_gds.html).

Do not apply ordinary marker filtering mechanically. Remove unequivocal
assay failures and paralogous or technically invalid markers, but avoid
aggressive call-rate filtering before examining sex-specific
missingness: genuine Y- or W-like presence is expressed as missingness
in the other sex. Similarly, do not LD-prune the initial discovery scan,
because a linked block of candidates is valuable genomic evidence. A
useful sensitivity analysis compares a minimally filtered screen with
the conventionally filtered dataset and reports which candidates
persist.

Complete the ordinary sample and marker quality-control workflow before
the final screen, but do not automatically discard every marker with
sex-biased missingness: sex-specific absence may be the biological
signal of interest. Inspect missingness before and after filtering and
retain the original audit files.

At minimum:

1.  verify individual identifiers and recorded sex;
2.  remove failed samples, obvious duplicates, and unresolved mixtures;
3.  inspect genotype call rate, depth, allele balance, and paralogy;
4.  record plate, library, lane, sequencing run, population, collection
    site, project, and year when available;
5.  check whether recorded sex is confounded with any technical or
    biological group; and
6.  retain enough known females and males for discovery and independent
    validation.

The metadata must contain `INDIVIDUALS` and the column named by
`sex.column`. The values `F`, `female`, `M`, and `male` are recognized
without regard to case. Other values are treated as unknown and are not
used to discover candidates.

``` r

sample_metadata <- readr::read_tsv("sample_metadata.tsv")

table(sample_metadata$SEX, useNA = "ifany")
table(sample_metadata$SEX, sample_metadata$PLATE, useNA = "ifany")
table(sample_metadata$SEX, sample_metadata$POPULATION, useNA = "ifany")
```

If sex and population or batch are inseparable, the data cannot
distinguish a sex effect from that other effect. More modelling cannot
recover information that the sampling design does not contain. Add
balanced samples or validate the markers in independent plates,
populations, or projects.

## Run the screen

``` r

sex_markers <- radr::sexy_markers(
  data = "study.filtered.gds",
  strata = sample_metadata,
  sex.column = "SEX",
  presence.threshold = 0.90,
  absence.threshold = 0.10,
  min.heterozygosity.difference = 0.20,
  coverage.ratio.threshold = 1.50,
  coverage.threshold = 1,
  min.samples.per.sex = 5,
  fdr.threshold = 0.05,
  require.significance = TRUE,
  path.folder = getwd()
)

sex_markers
sex_markers$candidates
```

The permissive default `require.significance = FALSE` selects markers
from the effect-size rules alone. This is useful for exploration and for
very small datasets in which formal tests have low power. For a final
discovery panel, `require.significance = TRUE` is usually preferable: a
marker must then satisfy both the effect-size rule and the FDR threshold
for its method. Significance is not a substitute for a meaningful effect
size or external validation.

### Thresholds

`presence.threshold` and `absence.threshold` define Y-like and W-like
presence/absence patterns. For example, `0.90` and `0.10` require the
marker to be detected in at least 90% of the expected sex and no more
than 10% of the other sex.

`min.heterozygosity.difference` is the minimum absolute difference
between female and male heterozygosity. `coverage.ratio.threshold` is
the minimum ratio between the larger and smaller mean normalized depths.
Thresholds should be specified before examining the candidate identities
when possible, then tested in sensitivity analyses.

`coverage.threshold` defines read presence when depth is available. If
retained depth is unavailable, the function uses a non-missing genotype
as presence and does not run the normalized-depth method.

### Statistical tests

Presence and heterozygosity are tested with two-sample score tests for
proportions. Normalized depth is tested with Welch tests on
`log2(1 + normalized depth)`. Benjamini-Hochberg FDR correction is
applied separately to each method because the methods test different
observations and hypotheses.

Within each sample, read depth is divided by its mean positive depth
over the active markers. This reduces sample-wide sequencing-depth
differences. It does not remove marker-by-batch interactions or
confounding between sex and plate, lane, library, population, or
collection group.

## Understand the result

The returned `sexy_markers` object contains:

``` r

names(sex_markers)

sex_markers$statistics
sex_markers$candidates
sex_markers$assignment_panel
sex_markers$sample_summary
sex_markers$metadata_audit
sex_markers$plots
```

The direction of every statistic is female minus male:

| Quantity | Positive value | Negative value |
|----|----|----|
| `presence_difference` | more frequent in females | more frequent in males |
| `heterozygosity_difference` | more heterozygous in females | more heterozygous in males |
| `log2(coverage_ratio_female_male)` | greater normalized depth in females | greater normalized depth in males |

Therefore, a strong negative presence difference is Y-like and a strong
positive difference is W-like. Positive heterozygosity or depth evidence
is labelled X-like, and negative evidence is labelled Z-like. These are
operational labels for expected patterns, not direct chromosome
annotations.

A marker can be selected by more than one method. Agreement among
methods is useful evidence, but disagreement is not automatically an
error: young, homomorphic, diverged, repetitive, or incompletely
assembled sex chromosomes can generate different combinations of
signals.

## Files written

Each run creates a dated results folder containing:

| File | Purpose |
|----|----|
| `radr_sexy_markers_args_*.tsv` | complete argument record |
| `sex_marker_statistics.tsv` | statistics and candidate flags for every tested marker |
| `candidate_sex_markers.tsv` | markers passing at least one candidate rule |
| `sex_assignment_panel.tsv` | assignment-ready Y-like and W-like markers |
| `sex_sample_summary.tsv` | numbers and mean positive depth by recorded sex |
| `sex_metadata_audit.tsv` | association of recorded sex with metadata columns |
| `candidate_sex_markers.fasta` | candidate sequences, when sequence metadata exist |
| `sex_markers_presence.*` | presence effect versus FDR plot |
| `sex_markers_heterozygosity.*` | heterozygosity effect versus FDR plot |
| `sex_markers_coverage.*` | normalized-depth effect versus FDR plot, when available |

The metadata audit reports Cramer’s V between recorded sex and each
available categorical metadata variable. A strong association is a
warning about possible confounding; it is not proof that the candidate
is technical, and it is not an adjusted association test. Review the
underlying contingency table and the experimental design.

## Validate candidates

Discovery and validation answer different questions. Reproducing the
recorded sex of the same individuals used to discover a panel measures
resubstitution, not out-of-sample accuracy.

A strong validation plan includes several levels:

1.  **Technical validation:** inspect per-sample genotypes, depth,
    allele balance, missingness, and sequence uniqueness. Repeat
    libraries when possible.
2.  **Design validation:** test the panel in different plates, lanes,
    projects, sampling years, and laboratories.
3.  **Population validation:** test independent populations because sex
    chromosomes and restriction-site polymorphisms can vary
    geographically.
4.  **Genomic validation:** map sequences to a chromosome-scale assembly
    and test whether candidates co-localize in a plausible sex-linked
    region.
5.  **Biological validation:** genotype independently sexed animals and
    report sensitivity, specificity, ambiguous calls, and confidence
    intervals.
6.  **Assay validation:** develop and test an independent PCR, amplicon,
    or SNP assay for the final operational marker panel.

Use a holdout set or cross-validation during panel development, but
retain a fully independent final validation set whenever the intended
application is management, monitoring, or clinical-quality sex
assignment.

## Pass the panel to assigner

`radr` discovers and audits candidate markers. Genetic sex assignment
belongs in `assigner`, which keeps prediction separate from discovery
and does not overwrite recorded sex.

``` r

assignment <- assigner::assign_genetic_sex(
  data = "study.filtered.gds",
  panel = sex_markers$assignment_panel,
  strata = sample_metadata,
  sex.column = "SEX",
  path.folder = getwd()
)

assignment$assignments
assignment$comparison
```

The assignment panel intentionally contains Y-like and W-like presence
markers, which provide direct directional evidence for the heterogametic
sex. X-like and Z-like candidates remain valuable for biological
interpretation and custom models, but they are not silently converted
into presence-based assignment rules.

## Alternative methods

No method is best for every input type or sex-chromosome system.
Applying more than one method can be informative when their assumptions
and input data are understood.

| Method | Suitable input and main signal | When it is useful |
|----|----|----|
| **RADSex** | Demultiplexed RADseq reads; sex-biased marker presence and depth distributions | A fast, purpose-built workflow starting from reads, without requiring an existing SNP GDS |
| **dartR / dartR.sexlinked** | Reduced-representation SNPs in `genlight`; sex-linked genotype and call-rate patterns | An R workflow for identifying, filtering, and using multiple classes of sex-linked loci |
| **SDpop** | Genotypes or RNA-seq from natural populations; hierarchical models of sex-linkage and gametologs | Model-based inference of XY or ZW systems and posterior evidence for sex-linked genes, especially with mapped genomic data |
| **BayeScan by sex** | Genotyped SNPs with females and males treated as groups; locus-specific FST outliers | Detecting strongly differentiated sex-associated markers and testing whether they bias population-genetic analyses |
| **OutFLANK** | Genotyped SNPs partitioned by sex; allele-frequency differentiation | A generic outlier scan that can reveal differentiated homologous regions but is not a sex-marker-specific model |
| **pcadapt** | Genotyped SNPs; association with principal components | A complementary structure-based outlier scan, but results can reflect population or technical structure unrelated to sex |

RADSex begins closer to raw RADseq reads and is especially attractive
when a SNP dataset has not already been produced. `dartR.sexlinked`
provides a broader genlight-based workflow that can identify sex-linked
classes and infer sex. SDpop is a probabilistic population-genomic
method with different data and model requirements. BayeScan, OutFLANK,
and pcadapt are generic outlier methods: they can detect
allele-frequency differences between recorded sexes but do not by
themselves establish sex linkage. Benestan et al. (2017) provide a
particularly useful demonstration of BayeScan as both a discovery method
and a sensitivity analysis for sex-ratio bias in estimates of population
differentiation.

### Why the current radr workflow differs from dartR.sexlinked

The comparison below is based on `dartR.sexlinked` 1.2.2, source commit
`4bb887a` dated 23 March 2026. Software changes, so the package version
should always accompany a methodological comparison. Floriaan
Devloo-Delva is listed as an author of `dartR.sexlinked`; the present
comparison concerns the current code and inferential workflow, not
contributor credit.

`dartR.sexlinked` is useful when data already live in the dartR
`genlight` ecosystem, the investigator is prepared to declare an XY or
ZW system in advance, and identifying several genotype-based classes for
filtering or assignment is the main objective. Its published companion
workflow was tested in birds and a mammal and is an important
alternative to `radr`.

For a new `radr` analysis, however,
[`radr::sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md)
has several deliberate advantages:

- **The hypotheses are not forced by a declared system.** `radr` reports
  Y-, W-, X-, and Z-like patterns from their observed directions. The
  current dartR functions require `system = "xy"` or `system = "zw"`
  before screening.
- **Effect sizes are explicit and configurable.** `radr` requires
  user-recorded presence, absence, heterozygosity-difference, and
  depth-ratio thresholds. In the reviewed dartR code, Y/W detection uses
  a fixed call rate of 0.10 in the expected absent sex and FDR of 0.01,
  without also requiring a high call rate in the expected present sex.
  X/Z and gametolog labels require FDR at most 0.01 and the expected
  direction, but no minimum heterozygosity difference.
- **Observed zero counts remain zero.** The reviewed dartR
  implementation changes every zero cell to one before its Fisher or
  chi-square tests. This prevents numerical failures but changes the
  observed contingency table, effect estimate, and p-value. `radr` does
  not add this pseudocount.
- **Read depth is an independent evidence stream.** `radr` can use
  retained GDS depth, normalized within samples, to detect X/Z-like
  copy-number differences. The reviewed dartR screen uses genotype call
  rate and heterozygosity but does not implement the read-depth scenario
  from Devloo-Delva et al. (2024).
- **Technical confounding is reported.** `radr` audits the association
  between recorded sex and plate, lane, library, population, project,
  and other supplied metadata. The reviewed dartR workflow does not
  perform an equivalent audit or adjust its marker tests for these
  variables.
- **Large datasets do not require one full R genotype matrix.** `radr`
  reads the active GDS selection in chunks and restores it unchanged.
  The reviewed dartR implementation converts the complete `genlight`
  object to a conventional R matrix before analysis, which can become
  memory intensive.
- **The GDS requirement does not limit the upstream data format.**
  [`radr::sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md)
  deliberately analyses only a standardized SeqArray GDS, while
  `genometranslator` can create that GDS from VCF, DArT, PLINK, GENEPOP,
  FSTAT, tabular or Parquet data, and `genind`, `genlight`, or `gtypes`
  objects. By contrast, the reviewed `dartR.sexlinked` functions accept
  a dartR `genlight` object. The accompanying `dartR.base` package
  exposes dedicated readers for DArT SNP, SilicoDArT, VCF, PLINK, CSV,
  and FASTA inputs before analysis. Thus, neither sex-marker function
  directly accepts all of these source formats; the difference is the
  breadth of the conversion layer leading to its required analysis
  object.
- **Discovery is separated from assignment.** `radr` writes a versioned,
  marker-level Y/W panel. `assigner::assign_genetic_sex()` applies
  explicit agreement and margin thresholds, retains recorded sex for
  comparison, and returns `U` when evidence is tied or insufficient. The
  reviewed `dartR.sexlinked::gl.infer.sex()` combines presence calls
  with two K-means classifications; its result depends on a random seed
  unless one is supplied, and disagreement is returned as a starred
  majority assignment for manual review rather than an explicit
  unresolved class. In the reviewed code, equal support defaults to the
  first ordered label, and three missing preliminary assignments can
  consequently produce `*F` instead of an unresolved result.
- **The complete analysis is auditable.** `radr` writes the arguments,
  all marker statistics, candidate table, assignment panel, sample
  summary, metadata audit, sequences, and diagnostics to a dated folder.
  It also has automated tests for GDS selection restoration,
  classification direction, metadata handling, and ambiguous assignment.
  The reviewed `dartR.sexlinked` repository does not contain an
  automated test suite.

These are not cosmetic differences. Replacing zero counts, omitting
minimum effect sizes, or ignoring a sex-by-batch association can change
which markers are called. Conversely, `radr` is not automatically
correct because it is more auditable: its thresholds still require
biological justification, and every panel still requires independent
validation.

Devloo-Delva et al. (2024) compared the original `sexy_markers` workflow
with `dartR::gl.report.sexlinked()`, OutFLANK, and pcadapt in the White
Shark data. The methods recovered markedly different candidate sets,
which is expected because they target different signals. Candidate count
alone is therefore not a measure of accuracy. Compare genomic
localization, replication, assignment performance, and independent assay
validation.

## Reporting checklist

Report enough detail for the screen to be reproduced and challenged:

- sample counts by recorded sex, population, and technical batch;
- how sex was recorded and how uncertain records were treated;
- the starting and active numbers of samples, loci, and SNPs;
- all upstream filters and the complete marker audit;
- whether genotype calls or read depth defined marker presence;
- all effect-size, FDR, depth, and minimum-sample thresholds;
- the `radr` version and saved argument file;
- candidate counts by evidence class and overlap among classes;
- metadata associations and known design confounding;
- validation data that were not used for discovery; and
- sensitivity, specificity, ambiguous-call rate, and population
  transferability for any operational sex-assignment panel.

## References

Benestan L, Moore J-S, Sutherland BJG, Le Luyer J, Maaroufi H, Rougeux
C, Normandeau E, Rycroft N, Atema J, Harris LN, Tallman RF, Greenwood
SJ, Clark FK, Bernatchez L (2017) Sex matters in massive parallel
sequencing: evidence for biases in genetic parameter estimation and
investigation of sex determination systems. *Molecular Ecology*, 26,
6767-6783. <https://doi.org/10.1111/mec.14217>

Devloo-Delva F, Gosselin T, Butcher PA, Grewe PM, Huveneers C, Thomson
RB, Werry JM, Feutry P (2024) An R-based tool for identifying sex-linked
markers from restriction site-associated DNA sequencing with
applications to elasmobranch conservation. *Conservation Genetics
Resources*, 16, 11-16. <https://doi.org/10.1007/s12686-023-01331-5>

Feron R, Pan Q, Wen M, Imarazene B, Jouanno E, Anderson J, Herpin A,
Journot L, Parrinello H, Klopp C, Kottler VA, Roco AS, Du K, Kneitz S,
Adolfi M, Wilson CA, McCluskey B, Amores A, Desvignes T, Goetz FW,
Takanashi A, Kawaguchi M, Detrich HW, Oliveira MA, Nóbrega RH, Sakamoto
T, Nakamoto M, Wargelius A, Karlsen Ø, Wang Z, Stöck M, Waterhouse RM,
Braasch I, Postlethwait JH, Schartl M, Guiguen Y (2021) RADSex: a
computational workflow to study sex determination using restriction
site-associated DNA sequencing data. *Molecular Ecology Resources*, 21,
1715-1731. <https://doi.org/10.1111/1755-0998.13360>

Foll M, Gaggiotti OE (2008) A genome-scan method to identify selected
loci appropriate for both dominant and codominant markers: a Bayesian
perspective. *Genetics*, 180, 977-993.
<https://doi.org/10.1534/genetics.108.092221>

Gruber B, Unmack PJ, Berry OF, Georges A (2018) dartr: an R package to
facilitate analysis of SNP data generated from reduced representation
genome sequencing. *Molecular Ecology Resources*, 18, 691-699.
<https://doi.org/10.1111/1755-0998.12745>

Käfer J, Lartillot N, Marais GAB, Picard F (2021) Detecting sex-linked
genes using genotyped individuals sampled in natural populations.
*Genetics*, 218, iyab053. <https://doi.org/10.1093/genetics/iyab053>

Luu K, Bazin E, Blum MGB (2017) pcadapt: an R package to perform genome
scans for selection based on principal component analysis. *Molecular
Ecology Resources*, 17, 67-77. <https://doi.org/10.1111/1755-0998.12592>

Mijangos JL, Gruber B, Berry O, Pacioni C, Georges A (2022) dartR v2: an
accessible genetic analysis platform for conservation, ecology and
agriculture. *Methods in Ecology and Evolution*, 13, 2150-2158.
<https://doi.org/10.1111/2041-210X.13918>

Robledo-Ruiz DA, Austin L, Amos JN, Castrejón-Figueroa J, Harley DKP,
Magrath MJL, Sunnucks P, Pavlova A (2025) Easy-to-use R functions to
separate reduced-representation genomic datasets into sex-linked and
autosomal loci, and conduct sex assignment. *Molecular Ecology
Resources*, 25, e13844. <https://doi.org/10.1111/1755-0998.13844>

Whitlock MC, Lotterhos KE (2015) Reliable detection of loci responsible
for local adaptation: inference of a null model through trimming the
distribution of FST. *The American Naturalist*, 186, S24-S36.
<https://doi.org/10.1086/682949>
