# Exploring and filtering genomic data with radr

## Responsibilities of the packages

`radr` begins after genomic data have been imported and standardized.
[`genometranslator`](https://thierrygosselin.github.io/genometranslator/)
handles file formats and creates the GDS representation used throughout
this article. `radr` diagnoses the data and performs deliberate
filtering. Writers in `genometranslator` do not silently perform
quality-control filtering.

The examples are not evaluated because they refer to user files.

## Start from sample metadata and an unfiltered GDS

Maintain one authoritative sample metadata table and generate the
compact strata file from it. The strata file should contain stable
`INDIVIDUALS` names and the grouping variable required for the current
analysis. Consult the [genometranslator
vignette](https://thierrygosselin.github.io/genometranslator/articles/using_genometranslator.html)
for identifier, `NEW_ID`, and DArT `TARGET_ID` recommendations.

``` r

library(radr)

genome <- genometranslator::read_genome(
  data = "individuals.vcf.gz",
  strata = "strata.tsv"
)

genometranslator::summary_gds(genome)
```

Keep the original genomic source file untouched. A GDS is stateful:
filtering changes the active samples or markers and records those
changes in the same backing file. An open GDS object is a connection to
that file, not an in-memory copy of all its contents. Consequently,
assigning the connection to another R name does not create an
independent dataset:

``` r

genome_copy <- genome
```

Here, `genome` and `genome_copy` refer to the same physical GDS.
Filtering through either name updates the marker or individual `FILTERS`
metadata and the active selection seen through both names. Similarly,
assigning the output of a filter to a new object does not protect the
input GDS from modification.

Before continuing an analysis, inspect both the active dimensions and
the filters recorded in the metadata:

``` r

genometranslator::summary_gds(genome, check.sync = TRUE)
genometranslator::list_filters(genome, history = TRUE)
```

This is particularly important after an interrupted interactive
analysis, after reopening a previously used GDS, or before assuming that
an object is an unfiltered starting point.

Each filtering stage creates a dated results folder. Depending on the
function, it contains the function-call arguments, diagnostic figures
and tables, blacklists, and summaries. Filtering history is stored
inside the GDS and each dated `filters_parameters_*.tsv` is a cumulative
snapshot through that stage. It records operation order, parameter
values, before-and-after dimensions, and blacklist information. This
works identically for piped and sequential calls: the shared backing
GDS, rather than the R expression, connects the history.

## Inspect missingness before choosing what to remove

A common first question is whether to filter markers or individuals
first. There is no universal answer because each choice changes
statistics calculated for the other dimension.

[`detect_ibm()`](https://thierrygosselin.github.io/radr/reference/detect_ibm.md)
displays the raw variant-by-sample presence/absence pattern. The goal is
to reveal broad dataset structure, not to identify one labelled marker
or sample. For the first view, preserve input order: sequencing plates,
libraries, batches, and processing order can create patterns that
sorting would hide.

``` r

ibm <- detect_ibm(
  data = genome,
  sort.individuals = "input",
  sort.markers = "input",
  filename = "initial_missingness.png"
)
```

Interpret broad patterns cautiously:

- vertical bands suggest samples with unusual missingness;
- horizontal bands suggest marker sets with unusual missingness;
- rectangular blocks can indicate batches, strata, libraries, or marker
  panels;
- isolated patterns can reflect technical problems, biology, or both.

Use marker- or individual-missingness sorting only as a complementary
second view. A pattern is evidence for investigation, not an automatic
blacklist.

## Guided exploration

[`explore_genomes()`](https://thierrygosselin.github.io/radr/reference/explore_genomes.md)
provides an interactive overview for a new dataset and for users
learning how candidate thresholds affect it:

``` r

screened <- explore_genomes(
  data = genome,
  interactive.filter = TRUE
)
```

The function creates a dated folder containing arguments, figures,
helper tables, filters, and summaries. It is opinionated and convenient,
but it is not “one function to rule them all.” Marker type, sequencing
design, organism, population structure, and analysis objective determine
which filters are valid and in which order they should run.

## Build a tailored workflow

After learning the dataset, call individual functions with recorded
thresholds. The filtering order should express a hypothesis about how
the data became problematic. The two examples below deliberately use
opposite orders.

### Computational strategy: remove clear failures before expensive comparisons

Runtime and memory depend primarily on the number of active markers,
number of individuals, availability of depth fields, requested figures,
and number of cores. The order therefore matters computationally as well
as biologically. A filter applied early to 400,000 markers can make
every subsequent genotype, coverage, LD, distance, and local-PCA
calculation substantially cheaper.

The following timing classes are intentionally approximate for
contemporary desktop hardware. Tens of thousands of markers and hundreds
of individuals may finish near the lower end; several hundred thousand
markers and 1,000-2,000 individuals can move a step into the next class.

| Operation | Typical scale | Main cost and strategy |
|----|---:|----|
| Metadata validation, existing VCF filters, reproducibility fields, and monomorphic checks | seconds to a few minutes | Usually inexpensive; run first when applicable. |
| [`filter_genotyping()`](https://thierrygosselin.github.io/radr/reference/filter_genotyping.md) | seconds to a few minutes | Very fast even for large datasets and often removes many markers. It is usually the first marker-targeting filter to run. |
| [`filter_ma()`](https://thierrygosselin.github.io/radr/reference/filter_ma.md) | seconds to several minutes | Low memory cost with substantial potential to remove low-information markers. It is usually the second marker-targeting filter. Prefer MAC when preservation of the allele-frequency spectrum matters. |
| [`filter_snp_number()`](https://thierrygosselin.github.io/radr/reference/filter_snp_number.md) | seconds to a few minutes | Fast, but usually removes fewer markers. Use only when `LOCUS` represents a genuine shared read or RAD locus; it is generally not informative for one-record-per-locus reference VCFs. |
| [`filter_coverage()`](https://thierrygosselin.github.io/radr/reference/filter_coverage.md) | minutes to tens of minutes | Can remove a substantial number of markers, but reading and summarizing DP/AD is computationally and memory intensive. Run it last among the routine marker-targeting filters, after the dataset has already been reduced. |
| Individual missingness, heterozygosity, and coverage | minutes to tens of minutes | Recalculate after consequential marker filtering; remove known failed libraries earlier when project QC already identifies them. |
| Short-distance LD | seconds to minutes | Cheap only when genuine shared RAD/read `LOCUS` identifiers exist; it is skipped for one-record-per-locus reference VCFs. |
| Long-distance LD | minutes to hours | One of the expensive filters; cost grows with active marker count, chromosome size, missingness-aware selection, and requested figures. Run after basic marker and sample QC. |
| Pairwise genome comparisons, duplicate/mixed-genome diagnostics, HWE by groups, local PCA, and inversion diagnostics | minutes to hours | Potentially expensive in markers, samples, groups, windows, or pairs. Reserve for the reduced dataset unless the scientific question requires unfiltered data. |

A useful default computational order for the routine marker filters is:

1.  Validate samples, caller provenance, metadata, and already assigned
    VCF filters. Remove samples known from laboratory or sequencing QC
    to have failed completely.
2.  Run
    [`filter_genotyping()`](https://thierrygosselin.github.io/radr/reference/filter_genotyping.md).
    It is quick, works well on very large datasets, and frequently
    provides the largest early reduction in marker number.
3.  Run
    [`filter_ma()`](https://thierrygosselin.github.io/radr/reference/filter_ma.md)
    with a defensible global MAC threshold. It has a relatively small
    memory footprint and can remove many low-information markers.
    Low-count alleles contribute little information to many analyses and
    can make later calculations expensive or unstable. Preserve the
    unfiltered master when the downstream analysis needs the rare-allele
    spectrum.
4.  Run
    [`filter_snp_number()`](https://thierrygosselin.github.io/radr/reference/filter_snp_number.md)
    when loci genuinely group several SNP records. This filter is fast,
    although it commonly removes fewer markers than the two preceding
    filters. Skip it when every VCF record has a unique `LOCUS`, as is
    common for reference-based FreeBayes data.
5.  Run
    [`filter_coverage()`](https://thierrygosselin.github.io/radr/reference/filter_coverage.md)
    last among the routine marker-targeting filters. Coverage and
    allele-depth information can identify many technical failures, but
    DP/AD extraction and summarization are substantially more expensive
    than the preceding filters.
6.  Recalculate individual statistics and remove supported sample
    failures.
7.  Run short-distance pruning only when loci genuinely group several
    SNP records. Run long-distance LD and the more expensive pairwise or
    windowed diagnostics after the dataset has been reduced.

In compact form, the recommended marker-filtering sequence is therefore:

``` r

filtered <- genome |>
  radr::filter_genotyping() |>
  radr::filter_ma() |>
  radr::filter_snp_number() |> # omit when LOCUS does not group genuine reads
  radr::filter_coverage()
```

This is not a universal biological filter order. For example, severely
failed samples can inflate marker missingness and should be removed
first, whereas widespread bad markers can make good samples appear poor.
Use the initial diagnostics to decide which problem dominates, record
before-and-after dimensions and elapsed time, and recalculate statistics
affected by each major decision. A fast filter is not automatically a
valid filter, and computational convenience should not erase variation
required by the intended analysis.

### Chain filters with a pipe

Most radr filters take `data` as their first argument, update the active
GDS, and return its connection. They can therefore be chained directly
with the native R pipe:

``` r

filtered <- genome |>
  radr::filter_genotyping() |>
  radr::filter_ma()
```

The result returned by
[`filter_genotyping()`](https://thierrygosselin.github.io/radr/reference/filter_genotyping.md)
becomes the first argument of
[`filter_ma()`](https://thierrygosselin.github.io/radr/reference/filter_ma.md).
No placeholder is required. This works in interactive mode too: the
first filter completes its questions and updates the GDS before the
second filter begins.

For a reproducible scripted analysis, provide the thresholds explicitly:

``` r

filtered <- genome |>
  radr::filter_genotyping(
    interactive.filter = FALSE,
    filter.genotyping = 0.20
  ) |>
  radr::filter_ma(
    interactive.filter = FALSE,
    filter.ma = 2
  )
```

If an explicit placeholder is useful, the native `|>` pipe uses `_`, not
`.`:

``` r

filtered <- radr::filter_genotyping(data = genome) |>
  radr::filter_ma(data = _)
```

The `.` placeholder belongs to the magrittr pipe `%>%`:

``` r

filtered <- radr::filter_genotyping(data = genome) %>%
  radr::filter_ma(data = .)
```

Prefer the first form when `data` is the first argument. Remember that
piping does not create a new physical GDS: each filter updates the same
backing file, then passes its connection to the next function. In the
example above, `filtered` and `genome` remain connections to that same
modified file; the new object name is not a snapshot of the result.

To compare alternative thresholds or filter orders, create separate
physical GDS files from an unfiltered master *before* running either
workflow. Close its open connection before copying the file:

``` r

# Create and preserve this master before applying any radr filter.
master_path <- "genome_unfiltered.gds"

SeqArray::seqClose(genome)

stopifnot(file.copy(master_path, "genome_workflow_A.gds"))
stopifnot(file.copy(master_path, "genome_workflow_B.gds"))

genome_A <- genometranslator::read_genome("genome_workflow_A.gds")
genome_B <- genometranslator::read_genome("genome_workflow_B.gds")

genometranslator::list_filters(genome_A)
genometranslator::list_filters(genome_B)
```

Copying a GDS after filters have already been applied also copies their
recorded state. If no clean master exists, recreate one from the
original VCF, DArT, or other source file.
[`genometranslator::reset_filters()`](https://thierrygosselin.github.io/genometranslator/reference/reset_filters.html)
can deliberately reactivate selected or all markers and individuals, but
it modifies the same GDS and is not a substitute for independent files
when comparing workflows.

### What “SNPs per locus” means for STACKS and FreeBayes

The interpretation of
[`filter_snp_number()`](https://thierrygosselin.github.io/radr/reference/filter_snp_number.md)
depends on whether the input preserves a genuine RAD/read locus
identifier. STACKS reconstructs RAD loci and associates multiple SNPs
with the same locus. With that representation, asking whether a roughly
70-100-bp sequence contains three or more SNPs can be informative. An
unusually dense locus may reflect paralogy, poor locus assembly,
misalignment, unexpectedly high local diversity, or unsuitable assembly
parameters. The appropriate threshold still depends on read length,
organism, sequencing design, and expected diversity.

FreeBayes uses the word *haplotype* differently. It is haplotype-aware
while evaluating nearby alleles and can emit a SNP, multi-nucleotide
polymorphism, indel, or complex allele as a single VCF record. A
FreeBayes record is therefore not automatically a reconstructed RAD
locus. In particular, the length of its REF or ALT allele is not the
number of SNPs observed on one read.

For a reference-based VCF,
[`filter_snp_number()`](https://thierrygosselin.github.io/radr/reference/filter_snp_number.md)
counts VCF records that share the stored `LOCUS`. If the upstream
workflow did not retain a common RAD/read locus identifier, each VCF
record has a unique locus. The original read boundaries cannot be
recovered unambiguously from genomic coordinates alone, so radr reports
one SNP per locus and skips this filter. Short-distance pruning in
[`filter_ld()`](https://thierrygosselin.github.io/radr/reference/filter_ld.md)
is skipped for the same reason, while genomic long-distance LD pruning
remains available.

For FreeBayes data without retained RAD/read grouping, place greater
emphasis on genotype quality, allele balance, depth and coverage, excess
heterozygosity, mapping diagnostics, and genomic LD. Variant density can
also be summarized in windows such as 70, 100, or 150 bp, but such
windows are genomic diagnostics: they should not be described
automatically as reads or RAD loci. If SNP-per-read filtering is
required, preserve a shared locus or fragment identifier upstream and
carry it into the GDS metadata.

### Marker noise is the dominant problem

Imagine that genotyping parameters were not tailored to the project, or
that files genotyped separately were combined after variant calling.
Differences in coverage, callers, parameter settings, reference
versions, or marker discovery can introduce a large number of
inconsistently genotyped SNPs. These markers can make otherwise
acceptable samples appear to have high missingness.

Remove clear marker failures first, then recalculate individual
statistics:

``` r

markers_first <- genome |>
  filter_genotyping(
    interactive.filter = FALSE,
    filter.genotyping = 0.20
  ) |>
  filter_ma(
    interactive.filter = FALSE,
    filter.ma = 2
  ) |>
  filter_individuals(
    interactive.filter = FALSE,
    filter.individuals.missing = 0.30
  )
```

This is a downstream repair strategy. When source data are available, a
consistent joint-genotyping workflow is preferable. At minimum, verify
that combined files used compatible references, marker definitions,
callers, filters, and genotype fields. Marker filtering cannot prove
that batch effects have been removed.

### Poor samples are the dominant problem

Now imagine that the project records already identify failed or weak
libraries: some samples received few reads, had poor base quality or
adapter contamination, mapped badly, or had abnormally low coverage.
Those individuals can introduce missing calls at many markers that
perform well in the rest of the project.

Remove the failed samples first, then return to
[`explore_genomes()`](https://thierrygosselin.github.io/radr/reference/explore_genomes.md)
with the individuals that remain. This resumes the normal guided
discovery workflow without assuming which marker filter should come
next:

``` r

samples_first <- genome |>
  filter_individuals(
    interactive.filter = FALSE,
    filter.individuals.missing = 0.50,
    filter.individuals.coverage.total = 1e6
  ) |>
  explore_genomes()
```

Ideally, these samples are identified before genotyping through FASTQ
and library QC: read yield, per-base quality, adapter content,
duplication, contamination, mapping rate, and depth should be reviewed
with suitable tools and project-specific expectations. Removing them in
radr is a downstream confirmation and safeguard, not a replacement for
upstream QC.

The numeric thresholds in both examples are illustrations. Derive
thresholds from diagnostic figures, project design, sequencing
expectations, and the requirements of the intended analysis. When
comparing alternative orders, create separate physical copies of the
original GDS file or re-import the untouched source data. Assigning the
same GDS connection to a different R variable does not create an
independent copy.

At each stage, ask what changed:

- Did filtering individuals alter marker missingness or allele counts?
- Did marker filtering change which individuals appear unusual?
- Is a signal restricted to one stratum, sequencing batch, or marker
  class?
- Does the chosen statistic match the intended downstream analysis?

Re-run the relevant diagnostic after a consequential filter. Avoid
choosing a threshold solely because it is conventional in another study.

## Review active filters

The parameter file documents how the dataset reached its current state.
To see which individual and marker filters are currently active in the
GDS metadata, use:

``` r

genometranslator::list_filters(genome)

# Include the chronological history and export a cumulative snapshot
filters <- genometranslator::list_filters(
  genome,
  history = TRUE,
  filename = "filters_parameters_complete.tsv"
)
```

The active-filter list answers *what is applied now*. `filters$history`
answers *what was done, in which order, and with which thresholds*,
including an operation that removed no markers or a filter that was
later reset.

Older GDS files do not yet contain this embedded history. Their existing
parameter files can be imported once, in chronological order, using
either the files themselves or their result folders:

``` r

genometranslator::import_filter_history(
  genome,
  paths = c(
    "06_filter_genotyping",
    "07_filter_ma"
  )
)
```

### Identifier-level filtering audit

For GDS workflows, every radr operation using the common filtering
mechanism also writes a complete identifier-level audit in its results
folder. This is not limited to functions named `filter_*`:
state-changing branches of
[`detect_duplicate_genomes()`](https://thierrygosselin.github.io/radr/reference/detect_duplicate_genomes.md)
and
[`detect_mixed_genomes()`](https://thierrygosselin.github.io/radr/reference/detect_mixed_genomes.md)
are audited in the same way.

The `filter_audit_manifest.tsv` file gives, for each operation, the
number of active markers and individuals before filtering, newly
removed, and kept. It points to four corresponding tables:

- `audit_*_markers_removed.tsv` and `audit_*_markers_kept.tsv`;
- `audit_*_individuals_removed.tsv` and `audit_*_individuals_kept.tsv`.

Empty removal tables are intentional: they prove that a completed
operation did not newly remove that entity type. The compact
chronological history stays inside the GDS, while these potentially
large tables remain in the operation folder. Diagnostic functions that
only *suggest* a blacklist without changing the GDS are not recorded as
applied filters. Likewise, VCF-to-VCF utilities audit their external
output files rather than the GDS, and genotype-masking functions alter
calls rather than the marker or individual filter state.

## Detection does not necessarily mean filtering

The `detect_*()` functions investigate particular structures or
problems:

``` r

duplicates <- detect_duplicate_genomes(data = filtered)
mixed <- detect_mixed_genomes(data = filtered)
paralogs <- detect_paralogs(data = filtered)
```

Review returned tables and figures together with sample metadata. A
duplicate, sex-linked marker, family group, or biologically divergent
sample may be important information rather than an error.

### Screen for candidate sex-linked markers

[`sexy_markers()`](https://thierrygosselin.github.io/radr/reference/sexy_markers.md)
accepts only a GDS filepath or an open GDS object. It respects the
active marker and sample selections, reads the GDS in chunks, and
restores the selections when it finishes. It does not run filters or
alter the GDS. See
[`vignette("detecting_sex_linked_markers", package = "radr")`](https://thierrygosselin.github.io/radr/articles/detecting_sex_linked_markers.md)
for the biological background, statistical interpretation, alternative
methods, and validation workflow.

``` r

sex.markers <- sexy_markers(
  data = filtered,
  strata = sample.metadata,
  sex.column = "SEX",
  require.significance = TRUE
)

sex.markers$candidates
```

The screen compares three complementary signals between known females
and males:

- marker presence identifies Y-like and W-like candidates;
- heterozygosity identifies X-like and Z-like candidates;
- normalized read depth provides a second X-like or Z-like signal when
  depth is available in the GDS.

The complete `sex_marker_statistics.tsv` file reports effect sizes,
p-values, and FDR values for every active marker. The smaller
`candidate_sex_markers.tsv` file applies the thresholds recorded in the
argument file. Candidate labels describe the direction of a signal; they
do not prove that a marker is physically located on a sex chromosome.

For classification, `sex_assignment_panel.tsv` contains only
assignment-ready Y-like and W-like presence/absence candidates. Validate
this panel in samples that were not used for marker discovery, then pass
it to `assigner::assign_genetic_sex()`. The assignment step belongs in
`assigner` because discovery accuracy, individual classification, and
validation are separate inferential tasks. Recorded sex is retained for
comparison and is never overwritten.

Include plate, lane, library, population, and other relevant columns in
`strata`. The `sex_metadata_audit.tsv` file measures their association
with recorded sex. Strong confounding can make technical dropout look
sex-linked. Within-sample depth normalization reduces global
sequencing-depth differences, but it cannot remove marker-specific batch
effects. Confirm important markers in independent samples and, where
possible, with genome placement or targeted validation.

## Export only after filtering decisions

Use `genometranslator` to write the final GDS to a downstream format.
The writer does not silently improve or filter the data:

``` r

genometranslator::write_genome(
  data = filtered,
  output = "vcf",
  filename = "analysis_ready"
)
```

Confirm the destination method’s assumptions before export. For example,
BayeScan, structure-like programs, relatedness software, and LD analyses
can require different marker, linkage, ploidy, and population
treatments.

## Dependencies and external tools

Inspect requirements at any time with:

``` r

radr_dependencies()
```

The diagnostic separates the required package foundation from optional
components. `SNPRelate` supports LD and IBS workflows, and `ragg`
accelerates PNG output from
[`detect_ibm()`](https://thierrygosselin.github.io/radr/reference/detect_ibm.md).

VCF-level filtering functions may call the external `bcftools`
executable, and
[`run_bayescan()`](https://thierrygosselin.github.io/radr/reference/run_bayescan.md)
requires BayeScan. Install these through the shared Conda `genomics`
environment described in the README, then start R or RStudio from that
activated environment.

## Reproducibility checklist

For a defensible filtering workflow, retain:

- the original genomic and sample metadata files;
- the standardized strata file;
- radr argument and filter-parameter logs;
- generated blacklists and whitelists;
- diagnostic figures before and after important decisions;
- the random seed where sampling is used;
- R, radr, genometranslator, tgbase, and executable versions;
- the final GDS and the exact export command.

During active development, record the Git commit in addition to the
package version.
