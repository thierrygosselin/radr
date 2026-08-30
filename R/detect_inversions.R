#' Detect candidate inversion-associated genomic regions
#'
#' Scan a diploid, biallelic GDS for genomic windows whose local population
#' structure differs from the genomic background. The implementation follows
#' the local-PCA principle: each window is represented by a low-rank covariance
#' matrix among individuals, distances are calculated between those matrices,
#' and classical multidimensional scaling (MDS) is used to identify unusual
#' windows. Candidate regions are then summarised using regional PCA,
#' heterozygosity, and linkage disequilibrium (LD).
#'
#' Run this screen before LD pruning. LD pruning can remove the extended
#' correlation pattern that makes an inversion-associated haploblock
#' detectable. After candidates have been reviewed, repeat downstream analyses
#' with the complete genome, with candidate regions excluded, and within each
#' candidate region or inferred arrangement.
#'
#' This is a screening method. A candidate region is not proof of a physical
#' inversion, and the returned coordinates describe an inversion-associated
#' haploblock rather than validated breakpoints. Long-read, read-pair, split-read,
#' cytogenetic, or genetic-map evidence is needed for breakpoint confirmation.
#'
#' @section Chromosome-specific windows:
#' Windows are constructed independently within each chromosome, linkage group,
#' or scaffold. A window can therefore never contain markers from two linkage
#' groups. Local PCA is not performed once on an entire linkage group: it is
#' performed separately for every SNP window within that linkage group. The
#' resulting window-level covariance summaries are compared across all valid
#' windows in the requested scan. Use `chromosome` to restrict that comparison
#' to one or more linkage groups. Standard genomic-position figures are faceted
#' by chromosome or linkage group in the results folder.
#'
#' @section Chromosome length and candidate extent:
#' Candidate tables report physical span, chromosome length, chromosome
#' fraction and percentage, and the left and right flanking lengths. Declared
#' chromosome lengths are preferred. The function first uses explicit
#' `chromosome.lengths`, then `reference.genome`, then the VCF contig dictionary
#' retained in the GDS. Chromosome names such as `14` and `chr14` are matched.
#' If no declared length is available, the largest observed marker position is
#' used as a clearly labelled underestimate of chromosome length; the resulting
#' candidate percentage may therefore be overestimated.
#'
#' @section Recommended staged workflow:
#' Begin with the default genome-wide scan on a carefully filtered dataset,
#' before LD pruning. Use the default result to locate unusual chromosomes and
#' inspect chromosome-wide PCA, regional PCA, missingness, depth,
#' heterozygosity, LD, cluster support, and sample metadata.
#'
#' When a candidate is found, rerun its chromosome separately. Compare several
#' primary `window.snps` and `step.snps` combinations, and optionally physical
#' or LD-scaled windows. Record the candidate start, end, span, chromosome
#' percentage, and interval overlap for every run. `sensitivity.window.snps`
#' provides complementary PC1-variance and LD summaries at additional scales,
#' but it does not formally recall candidate boundaries. Boundary stability
#' therefore requires separate focused calls with different primary window
#' settings. Repeat promising candidates under stricter call-rate, sample,
#' relatedness, and batch filters before seeking linkage, recombination,
#' long-read, read-pair, assembly, or breakpoint confirmation.
#'
#' @section Missing genotypes and RADseq data:
#' SNPs below `min.call.rate` are excluded separately within each window. For the
#' remaining SNPs, each missing genotype is replaced temporarily by that SNP's
#' mean dosage before covariance PCA is calculated. Mean imputation is
#' equivalent to giving a missing sample the observed allele-frequency
#' expectation at that SNP: after centring, it contributes zero rather than an
#' invented homozygous or heterozygous deviation. The GDS is not modified, and
#' observed calls are retained for heterozygosity calculations.
#'
#' This simple imputation keeps all samples in a common PCA space, but it does not
#' correct non-random missingness. In RADseq and other reduced-representation
#' data, allele dropout at restriction-site polymorphisms, uneven depth, library
#' quality, lanes, plates, populations, marker panels, or alignment quality may
#' correlate with biological groups. Mean imputation can then shrink affected
#' samples toward the window centre, alter covariance, weaken a real signal, or
#' create a batch-associated local signal. Before interpreting candidates, plot
#' call rate and depth against chromosome position, compare PCA clusters with
#' batch and plate metadata, repeat the scan at stricter `min.call.rate` values,
#' and confirm that candidate windows remain after removing problematic samples
#' or markers. Imputation makes the matrix computable; it does not make missing
#' data unbiased.
#'
#' LD uses observed genotypes with pairwise-complete correlations; missing LD
#' genotypes are not mean-imputed.
#'
#' @section Candidate evidence summary:
#' Regional k-means clustering is treated as a hypothesis, not as evidence by
#' itself. `three_cluster_evidence` requires three groups with at least three
#' samples each, a smallest-cluster frequency of at least 0.05, and a minimum
#' adjacent-centre separation of one pooled within-cluster standard deviation.
#' The candidate table also reports cluster compactness, PC1 variance,
#' heterozygosity excess in the middle cluster, LD within inferred arrangement
#' groups, LD in flanking windows, boundary contrasts, and the largest internal
#' score transition.
#'
#' `evidence_score` is a transparent screening heuristic from zero to five. One
#' point is assigned for quantitative three-cluster support, positive middle-
#' cluster heterozygosity excess, a positive candidate-to-flank score contrast,
#' regional LD above flanking LD, and continuity across at least two windows.
#' Scores of 0--2 are labelled `weak`, 3--4 `moderate`, and 5 `strong`. These
#' labels prioritise review; they do not convert a candidate into a structurally
#' confirmed inversion. Known-region overlaps are reported separately and do
#' not increase or decrease the evidence score.
#'
#' Every candidate is described conservatively. A local-PCA signal can reflect
#' a putative inversion-associated haploblock, but it can also arise near a
#' centromere, in a region of low recombination, from assembly or mapping
#' problems, introgression, population-specific missingness, or another form of
#' structural variation. The `candidate_class` and `alternative_explanations`
#' columns make these alternatives explicit. The function never reports a
#' structurally confirmed inversion.
#'
#' @section Sample metadata and technical confounders:
#' `strata` may contain `STRATA` and any additional sample-level variables such
#' as sequencing batch, library, lane, plate, extraction method, genotype caller,
#' or sampling year. Its rows first act as a whitelist. For every candidate, the
#' function then measures the association of each metadata variable with
#' regional PC1 and with arrangement assignments. A strong association does not
#' automatically reject a candidate,
#' but it identifies a biological or technical alternative that must be checked.
#'
#' Candidate outputs also compare call rate, mean read depth, and heterozygote
#' allele balance among putative arrangements when those quantities are stored
#' in the GDS. Standard `DP` and biallelic `AD` nodes are recognized, as are the
#' genometranslator genotype-metadata nodes `READ_DEPTH`, `ALLELE_REF_DEPTH`, and
#' `ALLELE_ALT_DEPTH`. Consequently, these diagnostics can be available for a
#' VCF, DArT two-row count file, or another source that retained read counts.
#' Depth or allele-balance shifts can indicate
#' mapping artifacts, paralogy, copy-number variation, or other structural
#' variation. They should not be interpreted automatically as inversion support.
#'
#' @section Cluster number and assignment stability:
#' The requested `cluster.k` controls the arrangement calls, but the function
#' also compares one-, two-, and three-cluster descriptions of regional PC1
#' using an approximate Gaussian BIC. This prevents three clusters from being
#' accepted merely because three centres were requested.
#'
#' With `stability.replicates > 0`, regional SNPs are resampled and arrangement
#' calls are repeated. The individual and overall agreement values quantify how
#' dependent the calls are on the exact SNP set. PC direction is arbitrary, so
#' reversed cluster labels are aligned before agreement is measured.
#'
#' @section Chromosome-wide PCA context:
#' When `chromosome.pca = TRUE`, the function calculates an independent PCA for
#' every chromosome or linkage group, using at most
#' `chromosome.pca.max.snps` evenly distributed SNPs. The arrangement labels
#' inferred from each candidate interval are then used only to colour these
#' chromosome-wide PCA panels. They do not influence the PCA or force three
#' groups on any chromosome.
#'
#' Separation concentrated on the candidate linkage group supports a localized
#' haploblock interpretation. Similar separation across many linkage groups
#' instead suggests genome-wide population structure, relatedness, admixture,
#' or technical confounding. Absence of separation in a whole-linkage-group PCA
#' does not reject a shorter candidate: unrelated SNPs elsewhere on the linkage
#' group can dilute a strong regional signal. Use this overview with the local
#' window scan and the candidate-specific technical diagnostics.
#'
#' @section Arrangement differentiation and marker loadings:
#' Regional PC1 loadings are returned for every SNP and ranked by absolute
#' magnitude. Pairwise allele-frequency differences and Hudson's FST are also
#' reported between inferred arrangement classes. These summaries describe the
#' regional genetic contrast and help select diagnostic markers. They do not
#' make the arrangement classes populations.
#'
#' Dxy is deliberately not calculated between `AA`, `AB`, and `BB`. In
#' particular, `AB` is an inferred heterokaryotype rather than an independently
#' sampled lineage. Population-level Dxy may be appropriate later when genuine
#' populations or evolutionary lineages are compared.
#'
#' @section Arrangement genotypes and sensitivity datasets:
#' For a three-cluster regional PCA, clusters are ordered along PC1 and labelled
#' `AA`, `AB`, and `BB` by default when quantitative three-cluster support is
#' present. `AA` and `BB` identify the two outer PC1 clusters; neither label
#' establishes the reference, ancestral, derived, or physically inverted
#' arrangement. `AB` is the intermediate cluster and is interpreted as a
#' putative heterokaryotype. Candidates without quantitative three-cluster
#' support retain their numeric cluster IDs but use neutral `Group 1`,
#' `Group 2`, and `Group 3` labels and do not receive arrangement dosages.
#'
#' These are putative arrangement genotypes, not sequence-level breakpoint
#' genotypes. For supported candidates, the individual table includes the
#' arrangement call, its numeric dosage (0, 1, or 2), the original numeric
#' cluster ID, and a relative assignment-confidence score.
#'
#' For every candidate, the function writes whitelists for each arrangement and
#' a homokaryotype-only whitelist containing `AA` and `BB`. It also returns a
#' combined whitelist for individuals classified as an outer arrangement in
#' every candidate region. These datasets support analyses that exclude
#' putative heterokaryotypes. This is useful before haplotype scans: `selscan`,
#' for example, does not accept missing genotype or haplotype data.
#'
#' @section Output and plotting:
#' Following other `radr` detection functions, each call creates a dated
#' `detect_inversions` results folder in the working directory (or below the
#' parent supplied with `path.folder`). It records the function arguments,
#' window and candidate tables, individual PCA scores, cluster summaries, LD
#' summaries, and standard diagnostic plots. PNG and PDF are written by default.
#'
#' @param data A GDS filepath or an open `SeqVarGDSClass` object.
#' @param strata Optional sample metadata supplied as a data frame or a tabular
#'   file. It must contain `INDIVIDUALS`; all other columns are retained. Rows
#'   act as a sample whitelist, and metadata such as `STRATA`, sequencing batch,
#'   library, lane, plate, or caller can be compared with regional PCA and
#'   putative arrangement assignments. The strata file remains sample metadata;
#'   it is not copied into or used to modify the GDS.
#' @param chromosome Optional chromosome or scaffold names to scan. By default,
#'   all chromosomes represented by at least one complete window are scanned.
#' @param chromosome.lengths Optional chromosome-length information supplied as
#'   a named numeric vector, a data frame with `CHROM` and `LENGTH` columns, or
#'   a tabular filepath containing those columns. Explicit values override
#'   lengths stored in the GDS.
#' @param reference.genome Optional reference FASTA or FASTA-index (`.fai`)
#'   filepath used to recover chromosome lengths. When a FASTA is supplied, its
#'   accompanying `FASTA.fai` file must already exist. This is unnecessary when
#'   the GDS retains a VCF contig dictionary with declared sequence lengths.
#' @param window.snps Number of SNPs per window.
#' @param step.snps Number of SNPs between consecutive window starts. Defaults
#'   to `window.snps` (non-overlapping windows).
#' @param window.bp Optional fixed physical window size in base pairs. When
#'   supplied, physical windows are used instead of fixed-SNP windows.
#' @param step.bp Distance in base pairs between physical window starts.
#'   Defaults to `window.bp`.
#' @param window.method Window construction method: `"snps"` for a fixed number
#'   of SNPs, `"bp"` for fixed physical windows, or `"ld"` for experimental
#'   LD-scaled windows. Supplying `window.bp` selects `"bp"`. Default:
#'   \code{window.method = "snps"}.
#' @param ld.window.threshold Adjacent-marker r-squared threshold used to end an
#'   experimental LD-scaled window after the minimum number of markers.
#'   Default: \code{ld.window.threshold = 0.1}.
#' @param ld.window.min.snps Minimum markers in an LD-scaled window. Default:
#'   \code{ld.window.min.snps = 50}.
#' @param ld.window.max.snps Maximum markers examined in an LD-scaled window.
#'   This bounds memory and computation. Default:
#'   \code{ld.window.max.snps = 500}.
#' @param sensitivity.window.snps Additional fixed-SNP window sizes used for a
#'   sensitivity analysis. These runs summarise PC1 variance and LD along the
#'   genome without replacing the primary candidate scan. This additional work
#'   is opt-in; the default `NULL` skips it. A useful focused set is
#'   `c(50, 100, 250, 500)`.
#' @param n.pcs Number of local covariance axes retained per window.
#' @param mds.axes Number of MDS axes used to score unusual windows.
#' @param outlier.quantile Quantile of the robust window score used as the
#'   candidate threshold.
#' @param min.window.snps Minimum number of usable polymorphic SNPs required in
#'   a window after missing-data and variance checks.
#' @param min.call.rate Minimum genotype call rate required for a SNP within a
#'   window.
#' @param min.candidate.windows Minimum number of consecutive candidate windows
#'   required to form a candidate region.
#' @param cluster.k Number of regional PCA clusters. The biological expectation
#'   for a common polymorphic inversion is often three, representing the two
#'   homokaryotypes and their heterokaryotype, but this is diagnostic rather
#'   than proof.
#' @param stability.replicates Number of SNP-resampling replicates used to test
#'   arrangement-assignment stability within every candidate. Use `0` to skip.
#'   Each replicate samples regional SNPs with replacement, repeats PCA and
#'   clustering, and compares assignments after matching cluster order on PC1.
#' @param stability.fraction Fraction of regional SNPs sampled in each stability
#'   replicate. Default: \code{stability.fraction = 0.8}.
#' @param parallel.core Number of independent R workers used for the local-PCA
#'   window scan. The default `1` is memory-conscious. Values greater than one
#'   open a separate read-only GDS connection in each worker, avoiding concurrent
#'   use of one connection. Each worker holds its own sample covariance matrix,
#'   so increase this value gradually on large datasets.
#' @param chromosome.pca Logical indicating whether an independent PCA is
#'   calculated for every chromosome or linkage group and coloured afterward
#'   with each candidate region's arrangement assignments. This diagnoses
#'   whether the same groups separate locally or throughout the genome.
#' @param chromosome.pca.max.snps Maximum number of evenly distributed SNPs
#'   used for each chromosome PCA. This bounds computation while retaining
#'   chromosome-wide coverage. Default: \code{2000}.
#' @param arrangement.labels Three labels, ordered from the lowest to highest
#'   regional PC1 cluster, used when `cluster.k = 3` and quantitative
#'   three-cluster support is present. Default:
#'   \code{arrangement.labels = c("AA", "AB", "BB")}.
#' @param known.regions Optional data frame describing centromeres, regions of
#'   low recombination, assembly gaps, or other annotations. It must contain
#'   `chromosome`, `start`, `end`, and `type` columns. Overlapping annotation
#'   types are reported for each candidate but are not used to select or score
#'   candidate windows.
#' @param ld.max.snps Maximum number of evenly spaced SNPs used for each regional
#'   LD matrix. This bounds memory use without changing the GDS input.
#' @param return.ld Logical indicating whether sampled regional LD matrices are
#'   retained in the result.
#' @param save.plots Logical indicating whether standard PDF and PNG diagnostic
#'   figures are written to the results folder.
#' @param plot.formats One or more of `"png"` and `"pdf"`.
#' @param random.seed Integer seed used for reproducible k-means clustering.
#' @param verbose Logical indicating whether progress messages are printed.
#' @param ... Further arguments for the standard radr workflow. Use
#'   `path.folder` to choose the parent results directory.
#'
#' @return An object of class `detect_inversions` containing:
#'   \itemize{
#'   \item `windows`: coordinates, quality statistics, MDS coordinates, robust
#'     scores, and candidate flags for all analysed windows;
#'   \item `candidates`: one row per contiguous candidate region;
#'   \item `candidate.summaries`: concise narrative summaries of the evidence,
#'     arrangement groups, regional LD, and boundary caution for each candidate;
#'   \item `diagnostics`: regional PCA scores, cluster assignments,
#'     heterozygosity, coverage, allele balance, cluster-number comparison,
#'     assignment stability, sample-metadata audit, marker loadings,
#'     arrangement differentiation, and optional LD matrices;
#'   \item `arrangement.genotypes`: one row per individual and candidate with
#'     putative arrangement genotype and relative assignment confidence;
#'   \item `homokaryotype.whitelist`: candidate-specific `AA` and `BB`
#'     individuals, plus `homokaryotype.all.candidates` for the intersection;
#'   \item `sensitivity`: optional summaries for additional fixed-SNP window
#'     sizes;
#'   \item `chromosome.pca`: independent chromosome or linkage-group PCA
#'     scores and summaries. Candidate-region arrangement labels are joined
#'     only for plotting and never influence these chromosome-wide PCAs;
#'   \item `chromosome.lengths`: declared or estimated chromosome lengths,
#'     observed marker maxima, and the provenance of each length;
#'   \item `path.folder` and `output.files`: locations of written results;
#'   \item `settings`: the effective analysis settings.
#'   }
#'
#' @references Li H, Ralph P (2019). Local PCA shows how the effect of
#' population structure differs along the genome. Genetics, 211, 289-304.
#' \doi{10.1534/genetics.118.301747}.
#'
#' Faria R, Johannesson K, Butlin RK, Westram AM (2019). Evolving inversions.
#' Trends in Ecology & Evolution, 34, 239-248.
#' \doi{10.1016/j.tree.2018.12.005}.
#'
#' Wellenreuther M, Bernatchez L (2018). Eco-evolutionary genomics of
#' chromosomal inversions. Trends in Ecology & Evolution, 33, 427-440.
#' \doi{10.1016/j.tree.2018.04.002}.
#'
#' Hudson RR, Slatkin M, Maddison WP (1992). Estimation of levels of gene flow
#' from DNA sequence data. Genetics, 132, 583-589.
#'
#' Bhatia G, Patterson N, Sankararaman S, Price AL (2013). Estimating and
#' interpreting FST: the impact of rare variants. Genome Research, 23,
#' 1514-1521. \doi{10.1101/gr.154831.113}.
#'
#' @export
#' @author Thierry Gosselin \email{Thierry.Gosselin@@csiro.au}

detect_inversions <- function(
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
) {
  # ---------------------------------------------------------------------------
  # Standard radr setup and reproducibility records
  # ---------------------------------------------------------------------------
  .start <- tgbase::startup(
    package = "radr",
    f.name = "detect_inversions",
    verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)

  rad.dots <- radr_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("path.folder", "internal"),
    deprecated = NULL,
    verbose = FALSE
  )

  path.folder <- tgbase::generate_folder(
    folder = "detect_inversions",
    path.folder = path.folder,
    internal = internal,
    file.date = file.date,
    prefix.int = TRUE,
    verbose = verbose
  )
  tgbase::write_tgbase_tsv(
    data = rad.dots,
    path.folder = path.folder,
    filename = "radr_detect_inversions_args",
    date = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  if (missing(data)) rlang::abort("Argument `data` is required.")

  integer.args <- list(
    window.snps = window.snps,
    step.snps = step.snps,
    n.pcs = n.pcs,
    mds.axes = mds.axes,
    min.window.snps = min.window.snps,
    min.candidate.windows = min.candidate.windows,
    cluster.k = cluster.k,
    parallel.core = parallel.core,
    chromosome.pca.max.snps = chromosome.pca.max.snps,
    ld.window.min.snps = ld.window.min.snps,
    ld.window.max.snps = ld.window.max.snps,
    ld.max.snps = ld.max.snps,
    random.seed = random.seed
  )
  # Validate related integer controls together so the user sees every invalid
  # argument in one message rather than fixing them one at a time.
  bad.integer <- purrr::map_lgl(
    integer.args,
    ~ length(.x) != 1L || !is.numeric(.x) || is.na(.x) ||
      !is.finite(.x) || .x < 1 || .x != as.integer(.x)
  )
  if (any(bad.integer)) {
    rlang::abort(paste0(
      "These arguments must be positive whole numbers: ",
      paste(names(integer.args)[bad.integer], collapse = ", "), "."
    ))
  }
  integer.args <- purrr::map(integer.args, as.integer)
  window.snps <- integer.args$window.snps
  step.snps <- integer.args$step.snps
  n.pcs <- integer.args$n.pcs
  mds.axes <- integer.args$mds.axes
  min.window.snps <- integer.args$min.window.snps
  min.candidate.windows <- integer.args$min.candidate.windows
  cluster.k <- integer.args$cluster.k
  parallel.core <- integer.args$parallel.core
  chromosome.pca.max.snps <- integer.args$chromosome.pca.max.snps
  if (length(stability.replicates) != 1L || !is.numeric(stability.replicates) ||
      is.na(stability.replicates) || !is.finite(stability.replicates) ||
      stability.replicates < 0 ||
      stability.replicates != as.integer(stability.replicates)) {
    rlang::abort("`stability.replicates` must be a non-negative whole number.")
  }
  stability.replicates <- as.integer(stability.replicates)
  ld.window.min.snps <- integer.args$ld.window.min.snps
  ld.window.max.snps <- integer.args$ld.window.max.snps
  ld.max.snps <- integer.args$ld.max.snps
  random.seed <- integer.args$random.seed

  window.method <- match.arg(window.method)
  if (!is.null(window.bp)) window.method <- "bp"
  if (window.method == "bp" && is.null(window.bp)) {
    rlang::abort("`window.bp` is required when `window.method = \"bp\"`.")
  }
  if (!is.null(window.bp)) {
    if (length(window.bp) != 1L || !is.numeric(window.bp) || is.na(window.bp) ||
        !is.finite(window.bp) || window.bp < 1) {
      rlang::abort("`window.bp` must be NULL or one positive number.")
    }
    window.bp <- as.numeric(window.bp)
    if (is.null(step.bp)) step.bp <- window.bp
    if (length(step.bp) != 1L || !is.numeric(step.bp) || is.na(step.bp) ||
        !is.finite(step.bp) || step.bp < 1) {
      rlang::abort("`step.bp` must be one positive number when `window.bp` is used.")
    }
    step.bp <- as.numeric(step.bp)
  }
  if (!is.null(sensitivity.window.snps)) {
    if (!is.numeric(sensitivity.window.snps) || anyNA(sensitivity.window.snps) ||
        any(!is.finite(sensitivity.window.snps)) ||
        any(sensitivity.window.snps < 1) ||
        any(sensitivity.window.snps != as.integer(sensitivity.window.snps))) {
      rlang::abort("`sensitivity.window.snps` must contain positive whole numbers.")
    }
    sensitivity.window.snps <- sort(unique(as.integer(sensitivity.window.snps)))
  }

  if (window.method == "snps" && min.window.snps > window.snps) {
    rlang::abort("`min.window.snps` cannot exceed `window.snps`.")
  }
  .inversion_check_probability(outlier.quantile, "outlier.quantile", open = TRUE)
  .inversion_check_probability(min.call.rate, "min.call.rate", open = FALSE)
  .inversion_check_probability(
    stability.fraction, "stability.fraction", open = TRUE
  )
  .inversion_check_probability(
    ld.window.threshold, "ld.window.threshold", open = FALSE
  )
  if (ld.window.max.snps < ld.window.min.snps) {
    rlang::abort("`ld.window.max.snps` cannot be smaller than `ld.window.min.snps`.")
  }
  if (cluster.k == 3L) {
    if (length(arrangement.labels) != 3L || anyNA(arrangement.labels) ||
        any(!nzchar(as.character(arrangement.labels))) ||
        anyDuplicated(arrangement.labels)) {
      rlang::abort("`arrangement.labels` must contain three unique non-empty labels.")
    }
    arrangement.labels <- as.character(arrangement.labels)
  }
  if (!is.logical(return.ld) || length(return.ld) != 1L || is.na(return.ld)) {
    rlang::abort("`return.ld` must be TRUE or FALSE.")
  }
  if (!is.logical(save.plots) || length(save.plots) != 1L || is.na(save.plots)) {
    rlang::abort("`save.plots` must be TRUE or FALSE.")
  }
  if (!is.logical(chromosome.pca) || length(chromosome.pca) != 1L ||
      is.na(chromosome.pca)) {
    rlang::abort("`chromosome.pca` must be TRUE or FALSE.")
  }
  known.regions <- .inversion_validate_known_regions(known.regions)
  plot.formats <- unique(tolower(as.character(plot.formats)))
  if (!length(plot.formats) || any(!plot.formats %in% c("png", "pdf"))) {
    rlang::abort("`plot.formats` must contain `png`, `pdf`, or both.")
  }
  set.seed(random.seed)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    rlang::abort("`verbose` must be TRUE or FALSE.")
  }

  # ---------------------------------------------------------------------------
  # Open the GDS and construct the ordered marker map
  # ---------------------------------------------------------------------------
  opened.here <- FALSE
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.character(data) && length(data) == 1L && !is.na(data) &&
             file.exists(data)) {
    gds <- SeqArray::seqOpen(data)
    opened.here <- TRUE
  } else {
    rlang::abort("`data` must be a GDS filepath or an open SeqVarGDSClass object.")
  }
  on.exit({
    if (opened.here) try(SeqArray::seqClose(gds), silent = TRUE)
  }, add = TRUE)

  SeqArray::seqFilterPush(gds)
  on.exit(try(SeqArray::seqFilterPop(gds), silent = TRUE), add = TRUE)

  sample.id <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  sample.metadata <- .inversion_read_sample_metadata(strata, sample.id)
  if (!is.null(sample.metadata)) {
    sample.id <- sample.metadata$INDIVIDUALS
    SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
    if (verbose) {
      message(
        "Sample metadata whitelist retained ", length(sample.id),
        " individual(s) and ", ncol(sample.metadata) - 1L,
        " descriptive variable(s)."
      )
    }
  }
  variant.id <- SeqArray::seqGetData(gds, "variant.id")
  chrom <- as.character(SeqArray::seqGetData(gds, "chromosome"))
  position <- suppressWarnings(as.numeric(SeqArray::seqGetData(gds, "position")))

  if (length(sample.id) < max(cluster.k, n.pcs + 1L)) {
    rlang::abort("Too few samples for the requested PCA and clustering settings.")
  }
  if (length(variant.id) != length(chrom) || length(chrom) != length(position)) {
    rlang::abort("GDS variant IDs, chromosomes, and positions have inconsistent lengths.")
  }
  if (anyNA(chrom) || any(!nzchar(chrom)) || anyNA(position)) {
    rlang::abort("All scanned GDS variants require chromosome and position values.")
  }

  if (!is.null(chromosome)) {
    chromosome <- unique(as.character(chromosome))
    unknown <- setdiff(chromosome, unique(chrom))
    if (length(unknown) > 0L) {
      rlang::abort(paste0(
        "Chromosome(s) not found in the active GDS variants: ",
        paste(unknown, collapse = ", "), "."
      ))
    }
    keep <- chrom %in% chromosome
    variant.id <- variant.id[keep]
    chrom <- chrom[keep]
    position <- position[keep]
  }

  marker.order <- order(.inversion_chromosome_order(chrom), position, variant.id)
  marker.table <- data.frame(
    variant_id = variant.id[marker.order],
    chromosome = chrom[marker.order],
    position = position[marker.order],
    stringsAsFactors = FALSE
  )
  chromosome.length.info <- .inversion_chromosome_lengths(
    gds = gds,
    marker.table = marker.table,
    chromosome.lengths = chromosome.lengths,
    reference.genome = reference.genome,
    verbose = verbose
  )
  windows <- if (window.method == "ld") {
    if (verbose) message("Constructing experimental LD-scaled windows...")
    .inversion_make_ld_windows(
      gds = gds,
      marker.table = marker.table,
      sample.id = sample.id,
      min.snps = ld.window.min.snps,
      max.snps = ld.window.max.snps,
      threshold = ld.window.threshold
    )
  } else {
    .inversion_make_windows(
      marker.table = marker.table,
      window.snps = window.snps,
      step.snps = step.snps,
      window.bp = if (window.method == "bp") window.bp else NULL,
      step.bp = if (window.method == "bp") step.bp else NULL
    )
  }
  if (length(windows) < 3L) {
    rlang::abort(
      "At least three complete chromosome-specific windows are required."
    )
  }

  if (verbose) {
    message(
      "Analysing ", length(windows), " windows across ",
      length(unique(purrr::map_chr(windows, "chromosome"))),
      " chromosome(s)..."
    )
  }

  # Each window is read independently. This keeps memory proportional to one
  # window per worker rather than materialising the complete dosage matrix.
  window.results <- .inversion_scan_windows(
    gds = gds,
    windows = windows,
    sample.id = sample.id,
    n.pcs = n.pcs,
    min.call.rate = min.call.rate,
    min.window.snps = min.window.snps,
    parallel.core = parallel.core,
    verbose = verbose
  )

  valid <- purrr::map_lgl(window.results, ~ isTRUE(.x$valid))
  if (sum(valid) < 3L) {
    rlang::abort(
      "Fewer than three windows passed call-rate and polymorphism checks."
    )
  }
  if (verbose && any(!valid)) {
    message(sum(!valid), " window(s) failed quality checks and were not scored.")
  }

  valid.results <- window.results[valid]
  distance.matrix <- .inversion_covariance_distances(valid.results)
  mds.k <- min(mds.axes, nrow(distance.matrix) - 1L)
  mds <- stats::cmdscale(stats::as.dist(distance.matrix), k = mds.k)
  if (is.null(dim(mds))) mds <- matrix(mds, ncol = 1L)
  colnames(mds) <- paste0("MDS", seq_len(ncol(mds)))
  robust.score <- .inversion_robust_score(mds)
  threshold <- as.numeric(stats::quantile(
    robust.score,
    probs = outlier.quantile,
    names = FALSE,
    na.rm = TRUE,
    type = 8
  ))

  window.table <- purrr::map_dfr(window.results, function(x) {
    tibble::tibble(
      window_id = x$window_id,
      chromosome = x$chromosome,
      start = x$start,
      end = x$end,
      n_input_snps = x$n_input_snps,
      n_used_snps = x$n_used_snps,
      mean_call_rate = x$mean_call_rate,
      mean_ld_r2 = x$mean_ld_r2,
      valid = x$valid
    )
  })
  window.table$robust_score <- NA_real_
  window.table$candidate_window <- FALSE
  valid.index <- which(valid)
  window.table$robust_score[valid.index] <- robust.score
  window.table$candidate_window[valid.index] <- robust.score >= threshold
  # MDS can return a run-dependent number of axes; attach each available axis
  # without assuming that MDS1 and MDS2 both exist.
  for (j in seq_len(ncol(mds))) {
    window.table[[colnames(mds)[j]]] <- NA_real_
    window.table[[colnames(mds)[j]]][valid.index] <- mds[, j]
  }

  # ---------------------------------------------------------------------------
  # Join contiguous outlier windows and diagnose each candidate region
  # ---------------------------------------------------------------------------
  candidate.regions <- .inversion_candidate_regions(
    window.table = window.table,
    min.windows = min.candidate.windows
  )
  candidate.regions <- .inversion_annotate_candidates(
    candidate.regions = candidate.regions,
    known.regions = known.regions
  )
  candidate.regions <- .inversion_add_chromosome_context(
    candidate.regions = candidate.regions,
    chromosome.length.info = chromosome.length.info
  )
  # Create diagnostic columns before indexed assignment. This avoids tibble's
  # unknown-column warning while candidate rows are filled one at a time.
  candidate.regions$recommended_cluster_k <- rep(
    NA_integer_, nrow(candidate.regions)
  )
  candidate.regions$assignment_stability <- rep(
    NA_real_, nrow(candidate.regions)
  )

  diagnostics <- vector("list", nrow(candidate.regions))
  if (nrow(candidate.regions) > 0L) {
    if (verbose) {
      message(nrow(candidate.regions), " candidate region(s) selected for diagnostics.")
      message(
        "Candidate diagnostics may take time for large regions. Progress will ",
        "be reported for genotype extraction, regional PCA and clustering, ",
        "LD, resampling, and metadata summaries."
      )
    }
    # Keep this indexed loop because each GDS read, diagnostic list element,
    # and candidate-table row must remain synchronized for the same region.
    for (i in seq_len(nrow(candidate.regions))) {
      # Overlapping scan windows may extend beyond the joined candidate
      # boundaries. Select diagnostic markers by the reported physical interval
      # so PCA, LD, heterozygosity, and exported genotypes all describe exactly
      # the region printed in the candidate table and plot titles.
      region.variant.id <- marker.table$variant_id[
        marker.table$chromosome == candidate.regions$chromosome[i] &
          marker.table$position >= candidate.regions$start[i] &
          marker.table$position <= candidate.regions$end[i]
      ]
      candidate.label <- paste0(
        candidate.regions$candidate_id[i], " | chromosome ",
        candidate.regions$chromosome[i], ": ",
        format(round(candidate.regions$start[i] / 1e6, 3L), nsmall = 3L),
        "-",
        format(round(candidate.regions$end[i] / 1e6, 3L), nsmall = 3L),
        " Mb"
      )
      candidate.start <- proc.time()[["elapsed"]]
      if (verbose) {
        message(
          "\nCandidate ", i, " of ", nrow(candidate.regions), ": ",
          candidate.label
        )
        message(
          "  Reading ", length(region.variant.id),
          " regional SNPs for ", length(sample.id), " individuals..."
        )
      }
      dosage <- .inversion_get_dosage(gds, region.variant.id, sample.id)
      if (verbose) message("  Reading retained depth and allele-count data...")
      coverage <- .inversion_get_coverage(gds, region.variant.id, sample.id)
      diagnostics[[i]] <- .inversion_region_diagnostics(
        dosage = dosage,
        depth = coverage$depth,
        allele.balance = coverage$allele.balance,
        sample.id = sample.id,
        cluster.k = cluster.k,
        min.call.rate = min.call.rate,
        ld.max.snps = ld.max.snps,
        return.ld = return.ld,
        arrangement.labels = arrangement.labels,
        sample.metadata = sample.metadata,
        stability.replicates = stability.replicates,
        stability.fraction = stability.fraction,
        verbose = verbose
      )
      diagnostics[[i]]$coverage_source <- if (length(coverage$source)) {
        paste(coverage$source, collapse = ";")
      } else {
        "none"
      }
      diagnostics[[i]]$technical_summary$coverage_source <-
        diagnostics[[i]]$coverage_source
      diagnostics[[i]]$candidate_id <- candidate.regions$candidate_id[i]
      diagnostics[[i]]$chromosome <- candidate.regions$chromosome[i]
      diagnostics[[i]]$start <- candidate.regions$start[i]
      diagnostics[[i]]$end <- candidate.regions$end[i]
      diagnostics[[i]]$ld_position <- marker.table$position[
        match(diagnostics[[i]]$ld_variant_id, marker.table$variant_id)
      ]
      diagnostics[[i]]$pc1_loadings$chromosome <- candidate.regions$chromosome[i]
      diagnostics[[i]]$pc1_loadings$position <- marker.table$position[
        match(diagnostics[[i]]$pc1_loadings$variant_id, marker.table$variant_id)
      ]
      candidate.regions$n_samples[i] <- nrow(diagnostics[[i]]$scores)
      candidate.regions$n_snps[i] <- diagnostics[[i]]$n_snps
      candidate.regions$regional_mean_ld_r2[i] <- diagnostics[[i]]$mean_ld_r2
      candidate.regions$homokaryotype_mean_ld_r2[i] <-
        diagnostics[[i]]$homokaryotype_mean_ld_r2
      candidate.regions$heterokaryotype_mean_ld_r2[i] <-
        diagnostics[[i]]$heterokaryotype_mean_ld_r2
      candidate.regions$ld_structure_contrast[i] <-
        diagnostics[[i]]$ld_structure_contrast
      candidate.regions$pc1_variance[i] <- diagnostics[[i]]$pc1_variance
      candidate.regions$smallest_cluster_n[i] <- diagnostics[[i]]$smallest_cluster_n
      candidate.regions$smallest_cluster_frequency[i] <-
        diagnostics[[i]]$smallest_cluster_frequency
      candidate.regions$cluster_separation[i] <- diagnostics[[i]]$cluster_separation
      candidate.regions$cluster_compactness[i] <- diagnostics[[i]]$cluster_compactness
      candidate.regions$three_cluster_evidence[i] <-
        diagnostics[[i]]$three_cluster_evidence
      candidate.regions$recommended_cluster_k[i] <-
        diagnostics[[i]]$cluster_models$k[
          which.min(diagnostics[[i]]$cluster_models$bic)
        ]
      candidate.regions$assignment_stability[i] <-
        diagnostics[[i]]$assignment_stability
      candidate.regions$heterozygote_like_middle_cluster[i] <-
        diagnostics[[i]]$heterozygote_like_middle_cluster
      candidate.regions$middle_heterozygosity_excess[i] <-
        diagnostics[[i]]$middle_heterozygosity_excess
      if (verbose) {
        message(
          "  Candidate diagnostics completed in ",
          round(proc.time()[["elapsed"]] - candidate.start, 1L), " sec."
        )
      }
    }
  }
  candidate.regions <- .inversion_evidence_summary(candidate.regions)
  candidate.regions <- .inversion_classify_candidates(candidate.regions)
  candidate.summaries <- .inversion_candidate_summaries(
    candidate.regions = candidate.regions,
    diagnostics = diagnostics,
    arrangement.labels = arrangement.labels
  )

  arrangement.genotypes <- purrr::map_dfr(diagnostics, function(x) {
    if (is.null(x$scores) || !nrow(x$scores)) return(tibble::tibble())
    dplyr::transmute(
      x$scores,
      candidate_id = x$candidate_id,
      chromosome = x$chromosome,
      start = x$start,
      end = x$end,
      individual = .data$individual,
      cluster = .data$cluster,
      arrangement = .data$arrangement,
      arrangement_dosage = .data$arrangement_dosage,
      arrangement_confidence = .data$arrangement_confidence,
      distance_nearest = .data$distance_nearest,
      distance_second = .data$distance_second,
      heterozygosity = .data$heterozygosity,
      call_rate = .data$call_rate,
      mean_depth = .data$mean_depth,
      mean_heterozygote_allele_balance =
        .data$mean_heterozygote_allele_balance,
      assignment_stability = .data$assignment_stability
    )
  })
  if (!nrow(arrangement.genotypes)) {
    arrangement.genotypes <- tibble::tibble(
      candidate_id = character(),
      chromosome = character(),
      start = numeric(),
      end = numeric(),
      individual = character(),
      cluster = integer(),
      arrangement = character(),
      arrangement_dosage = integer(),
      arrangement_confidence = numeric(),
      distance_nearest = numeric(),
      distance_second = numeric(),
      heterozygosity = numeric(),
      call_rate = numeric(),
      mean_depth = numeric(),
      mean_heterozygote_allele_balance = numeric(),
      assignment_stability = numeric()
    )
  }
  homokaryotype.whitelist <- arrangement.genotypes |>
    dplyr::filter(.data$arrangement_dosage %in% c(0L, 2L)) |>
    dplyr::distinct(.data$candidate_id, .data$individual, .data$arrangement)
  homokaryotype.all.candidates <- if (nrow(candidate.regions) > 0L) {
    homokaryotype.whitelist |>
      dplyr::count(.data$individual, name = "n_candidates_homokaryotype") |>
      dplyr::filter(.data$n_candidates_homokaryotype == nrow(candidate.regions)) |>
      dplyr::select(.data$individual)
  } else {
    tibble::tibble(individual = character())
  }

  sensitivity <- .inversion_window_sensitivity(
    gds = gds,
    marker.table = marker.table,
    sample.id = sample.id,
    window.sizes = sensitivity.window.snps,
    min.call.rate = min.call.rate,
    min.window.snps = min.window.snps
  )

  chromosome.pca.results <- if (chromosome.pca) {
    .inversion_chromosome_pca(
      gds = gds,
      marker.table = marker.table,
      sample.id = sample.id,
      max.snps = chromosome.pca.max.snps,
      min.call.rate = min.call.rate,
      parallel.core = parallel.core,
      verbose = verbose
    )
  } else {
    list(scores = tibble::tibble(), summary = tibble::tibble())
  }

  output.files <- .inversion_write_outputs(
    path.folder = path.folder,
    window.table = window.table,
    candidate.regions = candidate.regions,
    diagnostics = diagnostics,
    candidate.summaries = candidate.summaries,
    arrangement.genotypes = arrangement.genotypes,
    homokaryotype.all.candidates = homokaryotype.all.candidates,
    sensitivity = sensitivity,
    chromosome.pca = chromosome.pca.results,
    chromosome.lengths = chromosome.length.info,
    threshold = threshold,
    save.plots = save.plots,
    plot.formats = plot.formats,
    verbose = verbose
  )

  if (!return.ld) {
    diagnostics <- purrr::map(diagnostics, function(x) {
      x$ld_matrices <- NULL
      x
    })
  }

  out <- list(
    windows = tibble::as_tibble(window.table),
    candidates = tibble::as_tibble(candidate.regions),
    diagnostics = diagnostics,
    candidate.summaries = candidate.summaries,
    arrangement.genotypes = arrangement.genotypes,
    homokaryotype.whitelist = homokaryotype.whitelist,
    homokaryotype.all.candidates = homokaryotype.all.candidates,
    sensitivity = sensitivity,
    chromosome.pca = chromosome.pca.results,
    chromosome.lengths = chromosome.length.info,
    path.folder = path.folder,
    output.files = output.files,
    settings = list(
      chromosome = chromosome,
      reference.genome = reference.genome,
      window.snps = window.snps,
      step.snps = step.snps,
      window.bp = window.bp,
      step.bp = step.bp,
      window.method = window.method,
      ld.window.threshold = ld.window.threshold,
      ld.window.min.snps = ld.window.min.snps,
      ld.window.max.snps = ld.window.max.snps,
      sensitivity.window.snps = sensitivity.window.snps,
      stability.replicates = stability.replicates,
      stability.fraction = stability.fraction,
      parallel.core = parallel.core,
      chromosome.pca = chromosome.pca,
      chromosome.pca.max.snps = chromosome.pca.max.snps,
      n.pcs = n.pcs,
      mds.axes = mds.k,
      outlier.quantile = outlier.quantile,
      score.threshold = threshold,
      min.window.snps = min.window.snps,
      min.call.rate = min.call.rate,
      min.candidate.windows = min.candidate.windows,
      cluster.k = cluster.k,
      arrangement.labels = arrangement.labels,
      known.regions = known.regions,
      ld.max.snps = ld.max.snps,
      random.seed = random.seed
    )
  )
  class(out) <- c("detect_inversions", class(out))
  out
}

.inversion_candidate_summaries <- function(
    candidate.regions, diagnostics, arrangement.labels
) {
  if (!nrow(candidate.regions)) return(character())
  summaries <- purrr::map2_chr(
    seq_len(nrow(candidate.regions)),
    diagnostics,
    function(i, diagnostic) {
      candidate <- candidate.regions[i, , drop = FALSE]
      cluster.summary <- diagnostic$cluster_summary
      cluster.labels <- if (
        isTRUE(diagnostic$three_cluster_evidence) &&
          nrow(cluster.summary) == length(arrangement.labels)
      ) {
        arrangement.labels
      } else {
        paste("Group", cluster.summary$cluster)
      }
      cluster.sizes <- paste0(
        cluster.summary$n, " ", cluster.labels,
        collapse = ", "
      )
      heterozygosity <- paste0(
        format(round(cluster.summary$mean_heterozygosity, 3L), nsmall = 3L),
        " ", cluster.labels,
        collapse = ", "
      )
      convincing <- identical(candidate$evidence_strength, "strong")
      paste0(
        candidate$candidate_id, "\n",
        "- Candidate interval: chromosome ", candidate$chromosome, ", ",
        format(candidate$start / 1e6, digits = 4L, nsmall = 3L), "-",
        format(candidate$end / 1e6, digits = 4L, nsmall = 3L), " Mb.\n",
        "- Candidate span: ",
        format(round(candidate$candidate_span_bp / 1e6, 3L), nsmall = 3L),
        " Mb, representing ",
        format(round(candidate$chromosome_percent, 1L), nsmall = 1L),
        "% of a ",
        format(round(candidate$chromosome_length_bp / 1e6, 3L), nsmall = 3L),
        " Mb chromosome",
        if (isTRUE(candidate$chromosome_length_declared)) {
          paste0(" (", candidate$chromosome_length_source, ").\n")
        } else {
          paste0(
            " (estimated from the largest observed marker position; the ",
            "percentage may be overestimated).\n"
          )
        },
        "- Diagnostic SNPs: ", candidate$n_snps,
        ", all inside the reported interval.\n",
        "- Evidence: ", candidate$evidence_score, "/5, ",
        candidate$evidence_strength, ".\n",
        "- PC1 variance: ",
        format(round(100 * candidate$pc1_variance, 1L), nsmall = 1L), "%.", "\n",
        "- Cluster separation: ",
        format(round(candidate$cluster_separation, 2L), nsmall = 2L), ".\n",
        "- Preferred PC1 cluster count by approximate BIC: ",
        candidate$recommended_cluster_k, ".\n",
        if (is.finite(candidate$assignment_stability)) {
          paste0(
            "- SNP-resampling assignment stability: ",
            format(round(100 * candidate$assignment_stability, 1L), nsmall = 1L),
            "%.\n"
          )
        } else {
          "- SNP-resampling assignment stability: not requested.\n"
        },
        "- Cluster sizes: ", cluster.sizes, ".\n",
        "- Mean heterozygosity: ", heterozygosity, ".\n",
        "- Regional mean r2: ",
        format(round(candidate$regional_mean_ld_r2, 5L), nsmall = 5L),
        ", compared with ",
        format(round(candidate$flanking_mean_ld_r2, 5L), nsmall = 5L),
        " in flanking windows.\n",
        if (convincing) {
          "Interpretation: the candidate remains convincing.\n"
        } else {
          "Interpretation: the candidate warrants further evaluation.\n"
        },
        "Caution: this is a threshold-defined candidate core, not a definitive inversion boundary."
      )
    }
  )
  stats::setNames(summaries, candidate.regions$candidate_id)
}

# =============================================================================
# Argument and annotation helpers
# =============================================================================

.inversion_check_probability <- function(x, name, open) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (open) valid <- valid && x > 0 && x < 1
  if (!open) valid <- valid && x >= 0 && x <= 1
  if (!valid) {
    interval <- if (open) "strictly between 0 and 1" else "between 0 and 1"
    rlang::abort(paste0("`", name, "` must be ", interval, "."))
  }
  invisible(x)
}

.inversion_chromosome_key <- function(x) {
  key <- tolower(trimws(as.character(x)))
  key <- sub("^chromosome[_ .-]*", "", key)
  key <- sub("^chr[_ .-]*", "", key)
  numeric.key <- grepl("^[0-9]+$", key)
  key[numeric.key] <- as.character(as.integer(key[numeric.key]))
  key
}

.inversion_read_length_table <- function(x, source) {
  if (is.numeric(x) && !is.null(names(x))) {
    out <- tibble::tibble(
      declared_chromosome = names(x),
      chromosome_length_bp = as.numeric(x)
    )
  } else {
    if (is.character(x) && length(x) == 1L && file.exists(x)) {
      x <- vroom::vroom(x, show_col_types = FALSE, progress = FALSE)
    }
    if (!is.data.frame(x)) {
      rlang::abort(
        "`chromosome.lengths` must be a named numeric vector, data frame, or tabular filepath."
      )
    }
    names(x) <- toupper(names(x))
    chromosome.column <- intersect(
      c("CHROM", "CHROMOSOME", "CONTIG", "SEQUENCE", "ID"), names(x)
    )
    length.column <- intersect(
      c("LENGTH", "LENGTH_BP", "CHROMOSOME_LENGTH", "SEQUENCE_LENGTH"),
      names(x)
    )
    if (!length(chromosome.column) || !length(length.column)) {
      rlang::abort(
        "The chromosome-length table requires chromosome and length columns, such as `CHROM` and `LENGTH`."
      )
    }
    out <- tibble::tibble(
      declared_chromosome = as.character(x[[chromosome.column[1L]]]),
      chromosome_length_bp = suppressWarnings(
        as.numeric(x[[length.column[1L]]])
      )
    )
  }
  out$length_source <- source
  out$declared_length <- TRUE
  invalid <- is.na(out$declared_chromosome) |
    !nzchar(out$declared_chromosome) |
    !is.finite(out$chromosome_length_bp) |
    out$chromosome_length_bp <= 0
  if (any(invalid)) {
    rlang::abort("Chromosome lengths require non-empty names and positive finite values.")
  }
  out |>
    dplyr::distinct(.data$declared_chromosome, .keep_all = TRUE)
}

.inversion_read_fai <- function(reference.genome) {
  if (is.null(reference.genome)) return(NULL)
  if (!is.character(reference.genome) || length(reference.genome) != 1L ||
      is.na(reference.genome) || !nzchar(reference.genome)) {
    rlang::abort("`reference.genome` must be NULL or one FASTA or `.fai` filepath.")
  }
  fai <- if (grepl("\\.fai$", reference.genome, ignore.case = TRUE)) {
    reference.genome
  } else {
    paste0(reference.genome, ".fai")
  }
  if (!file.exists(fai)) {
    rlang::abort(paste0(
      "Reference index not found: ", fai,
      ". Create the FASTA index first or provide `chromosome.lengths`."
    ))
  }
  fields <- readr::read_tsv(
    fai,
    col_names = c("CHROM", "LENGTH", "OFFSET", "LINE_BASES", "LINE_WIDTH"),
    col_types = "cdiii",
    progress = FALSE
  )
  .inversion_read_length_table(
    fields[, c("CHROM", "LENGTH")],
    source = paste0("FASTA index: ", basename(fai))
  )
}

.inversion_gds_contig_lengths <- function(gds) {
  id.node <- gdsfmt::index.gdsn(
    gds, "description/vcf.contig/ID", silent = TRUE
  )
  length.node <- gdsfmt::index.gdsn(
    gds, "description/vcf.contig/length", silent = TRUE
  )
  if (is.null(id.node) || is.null(length.node)) return(NULL)
  ids <- tryCatch(gdsfmt::read.gdsn(id.node), error = function(error) NULL)
  lengths <- tryCatch(
    gdsfmt::read.gdsn(length.node), error = function(error) NULL
  )
  if (is.null(ids) || is.null(lengths) || length(ids) != length(lengths)) {
    return(NULL)
  }
  .inversion_read_length_table(
    data.frame(CHROM = ids, LENGTH = lengths),
    source = "GDS VCF contig dictionary"
  )
}

.inversion_chromosome_lengths <- function(
    gds, marker.table, chromosome.lengths = NULL, reference.genome = NULL,
    verbose = FALSE
) {
  declared <- if (!is.null(chromosome.lengths)) {
    .inversion_read_length_table(
      chromosome.lengths, source = "user-supplied chromosome lengths"
    )
  } else if (!is.null(reference.genome)) {
    .inversion_read_fai(reference.genome)
  } else {
    .inversion_gds_contig_lengths(gds)
  }

  observed <- marker.table |>
    dplyr::group_by(.data$chromosome) |>
    dplyr::summarise(
      observed_marker_max_bp = max(.data$position, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(chromosome = as.character(.data$chromosome))

  if (is.null(declared) || !nrow(declared)) {
    out <- observed |>
      dplyr::mutate(
        chromosome_length_bp = .data$observed_marker_max_bp,
        declared_chromosome = .data$chromosome,
        length_source = "largest observed marker position (underestimate)",
        declared_length = FALSE
      )
    if (verbose) {
      message(
        "Chromosome lengths were unavailable; chromosome percentages use the ",
        "largest observed marker position and are labelled as estimates."
      )
    }
    return(out)
  }

  declared$key <- .inversion_chromosome_key(declared$declared_chromosome)
  observed$key <- .inversion_chromosome_key(observed$chromosome)
  matched <- match(observed$chromosome, declared$declared_chromosome)
  missing <- is.na(matched)
  matched[missing] <- match(observed$key[missing], declared$key)
  out <- observed
  out$chromosome_length_bp <- declared$chromosome_length_bp[matched]
  out$declared_chromosome <- declared$declared_chromosome[matched]
  out$length_source <- declared$length_source[matched]
  out$declared_length <- declared$declared_length[matched]
  unresolved <- is.na(out$chromosome_length_bp)
  out$chromosome_length_bp[unresolved] <- out$observed_marker_max_bp[unresolved]
  out$declared_chromosome[unresolved] <- out$chromosome[unresolved]
  out$length_source[unresolved] <-
    "largest observed marker position (underestimate)"
  out$declared_length[unresolved] <- FALSE
  inconsistent <- out$declared_length &
    out$chromosome_length_bp < out$observed_marker_max_bp
  if (any(inconsistent)) {
    rlang::warn(paste0(
      "Declared chromosome length was shorter than the largest observed marker ",
      "position for: ", paste(out$chromosome[inconsistent], collapse = ", "),
      ". Marker maxima are used for those sequences; verify the reference build."
    ))
    out$chromosome_length_bp[inconsistent] <-
      out$observed_marker_max_bp[inconsistent]
    out$length_source[inconsistent] <- paste0(
      "largest observed marker position; declared length was inconsistent"
    )
    out$declared_length[inconsistent] <- FALSE
  }
  out$key <- NULL
  if (verbose) {
    message(
      "Chromosome lengths: ", sum(out$declared_length), " declared and ",
      sum(!out$declared_length), " estimated from marker positions."
    )
  }
  out
}

.inversion_add_chromosome_context <- function(
    candidate.regions, chromosome.length.info
) {
  if (!nrow(candidate.regions)) {
    candidate.regions$candidate_span_bp <- numeric()
    candidate.regions$chromosome_length_bp <- numeric()
    candidate.regions$chromosome_fraction <- numeric()
    candidate.regions$chromosome_percent <- numeric()
    candidate.regions$left_flank_bp <- numeric()
    candidate.regions$right_flank_bp <- numeric()
    candidate.regions$chromosome_length_source <- character()
    candidate.regions$chromosome_length_declared <- logical()
    return(candidate.regions)
  }
  context <- chromosome.length.info |>
    dplyr::select(
      .data$chromosome, .data$chromosome_length_bp,
      chromosome_length_source = .data$length_source,
      chromosome_length_declared = .data$declared_length
    )
  out <- dplyr::left_join(candidate.regions, context, by = "chromosome")
  out$candidate_span_bp <- pmax(0, out$end - out$start + 1)
  out$chromosome_fraction <- out$candidate_span_bp / out$chromosome_length_bp
  out$chromosome_percent <- 100 * out$chromosome_fraction
  out$left_flank_bp <- pmax(0, out$start - 1)
  out$right_flank_bp <- pmax(0, out$chromosome_length_bp - out$end)
  out
}

.inversion_read_sample_metadata <- function(strata, sample.id) {
  if (is.null(strata)) return(NULL)
  if (is.character(strata) && length(strata) == 1L) {
    if (!file.exists(strata)) rlang::abort("The `strata` file does not exist.")
    strata <- vroom::vroom(
      strata,
      show_col_types = FALSE,
      progress = FALSE
    )
  }
  if (!is.data.frame(strata)) {
    rlang::abort("`strata` must be NULL, a data frame, or a tabular filepath.")
  }
  names(strata) <- stringi::stri_trans_toupper(names(strata))
  if (!"INDIVIDUALS" %in% names(strata)) {
    rlang::abort("`strata` requires an `INDIVIDUALS` column.")
  }
  strata$INDIVIDUALS <- as.character(strata$INDIVIDUALS)
  invalid <- is.na(strata$INDIVIDUALS) | !nzchar(strata$INDIVIDUALS)
  if (any(invalid)) rlang::abort("`strata$INDIVIDUALS` contains missing or empty IDs.")
  if (anyDuplicated(strata$INDIVIDUALS)) {
    rlang::abort("`strata$INDIVIDUALS` must contain unique sample IDs.")
  }
  unknown <- setdiff(strata$INDIVIDUALS, sample.id)
  if (length(unknown)) {
    rlang::inform(paste0(
      length(unknown), " metadata sample(s) were absent from the active GDS and ignored."
    ))
  }
  strata <- strata[strata$INDIVIDUALS %in% sample.id, , drop = FALSE]
  # Preserve the GDS order so every downstream matrix and metadata join remains
  # aligned without relying on a later sort.
  strata <- strata[match(sample.id[sample.id %in% strata$INDIVIDUALS],
                         strata$INDIVIDUALS), , drop = FALSE]
  if (nrow(strata) < 3L) {
    rlang::abort("Fewer than three GDS samples remain after applying `strata`.")
  }
  tibble::as_tibble(strata)
}

.inversion_validate_known_regions <- function(x) {
  if (is.null(x)) {
    return(data.frame(
      chromosome = character(), start = numeric(), end = numeric(),
      type = character(), stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(x)) {
    rlang::abort("`known.regions` must be NULL or a data frame.")
  }
  required <- c("chromosome", "start", "end", "type")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    rlang::abort(paste0(
      "`known.regions` is missing: ", paste(missing, collapse = ", "), "."
    ))
  }
  x <- x |>
    dplyr::select(dplyr::all_of(required)) |>
    dplyr::mutate(
      chromosome = as.character(.data$chromosome),
      type = as.character(.data$type),
      start = suppressWarnings(as.numeric(.data$start)),
      end = suppressWarnings(as.numeric(.data$end))
    )
  bad <- is.na(x$chromosome) | !nzchar(x$chromosome) |
    is.na(x$type) | !nzchar(x$type) |
    !is.finite(x$start) | !is.finite(x$end) | x$start > x$end
  if (any(bad)) {
    rlang::abort(
      "Every `known.regions` row requires a chromosome, type, and finite start <= end."
    )
  }
  x
}

.inversion_annotate_candidates <- function(candidate.regions, known.regions) {
  # Keep annotations descriptive: they flag alternative explanations but never
  # participate in candidate selection or evidence scoring.
  if (!nrow(candidate.regions) || !nrow(known.regions)) {
    return(candidate.regions |>
      dplyr::mutate(
        known_region_overlap = "none",
        n_known_region_overlaps = 0L
      ))
  }

  overlaps <- purrr::pmap(
    candidate.regions[c("chromosome", "start", "end")],
    function(chromosome, start, end) {
      overlap <- known.regions$chromosome == chromosome &
        known.regions$start <= end &
        known.regions$end >= start
      known.regions |>
        dplyr::slice(which(overlap)) |>
        dplyr::pull("type")
    }
  )

  candidate.regions |>
    dplyr::mutate(
      known_region_overlap = purrr::map_chr(overlaps, function(types) {
        types <- sort(unique(types))
        if (length(types)) paste(types, collapse = ";") else "none"
      }),
      n_known_region_overlaps = purrr::map_int(overlaps, length)
    )
}

.inversion_evidence_summary <- function(x) {
  if (!nrow(x)) {
    x$cluster_support <- logical()
    x$heterozygosity_support <- logical()
    x$boundary_support <- logical()
    x$ld_support <- logical()
    x$continuity_support <- logical()
    x$evidence_score <- integer()
    x$evidence_strength <- character()
    return(x)
  }
  cluster.support <- !is.na(x$three_cluster_evidence) & x$three_cluster_evidence
  heterozygosity.support <- is.finite(x$middle_heterozygosity_excess) &
    x$middle_heterozygosity_excess > 0
  boundary.support <- is.finite(x$boundary_contrast) & x$boundary_contrast > 0
  ld.support <- is.finite(x$regional_mean_ld_r2) &
    is.finite(x$flanking_mean_ld_r2) &
    x$regional_mean_ld_r2 > x$flanking_mean_ld_r2
  continuity.support <- x$n_windows >= 2L
  score <- as.integer(cluster.support) + as.integer(heterozygosity.support) +
    as.integer(boundary.support) + as.integer(ld.support) +
    as.integer(continuity.support)
  x |>
    dplyr::mutate(
      cluster_support = cluster.support,
      heterozygosity_support = heterozygosity.support,
      boundary_support = boundary.support,
      ld_support = ld.support,
      continuity_support = continuity.support,
      evidence_score = score,
      evidence_strength = dplyr::case_when(
        score >= 5L ~ "strong",
        score >= 3L ~ "moderate",
        .default = "weak"
      )
    )
}

.inversion_classify_candidates <- function(x) {
  if (!nrow(x)) {
    x$candidate_class <- character()
    x$alternative_explanations <- character()
    return(x)
  }
  annotation <- tolower(x$known_region_overlap)
  technical <- grepl("assembly|gap|repeat|mapping|coverage|missing|batch", annotation)
  recombination <- grepl("centromere|low.recomb|recombination", annotation)
  structural <- grepl("structural|sv|duplication|translocation", annotation)
  x$candidate_class <- dplyr::case_when(
    technical ~ "candidate technical or assembly-associated region",
    recombination ~ "candidate low-recombination or centromeric region",
    structural ~ "candidate structural-variation-associated haploblock",
    x$evidence_strength == "strong" ~ "putative inversion-associated haploblock",
    .default = "unresolved candidate haploblock"
  )
  x$alternative_explanations <- paste(
    "centromere or low recombination; assembly or mapping problem;",
    "introgression; population-specific missingness; other structural variation"
  )
  x
}

.inversion_chromosome_order <- function(x) {
  levels <- .inversion_sequence_layout(x)$levels
  match(as.character(x), levels)
}

.inversion_sequence_layout <- function(x) {
  labels <- unique(as.character(x))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  levels <- labels[stringi::stri_order(labels, numeric = TRUE)]
  lower <- stringi::stri_trans_tolower(labels)
  explicit.scaffold <- any(grepl("scaffold|contig|unitig", lower))
  numeric.labels <- all(grepl("^[0-9]+$", labels))
  chromosome.labels <- all(grepl(
    "^(chr|chromosome|lg|linkage[_. -]?group)[_. -]?[0-9a-z]+$",
    lower
  ))
  type <- dplyr::case_when(
    explicit.scaffold ~ "scaffolds / contigs",
    chromosome.labels ~ "chromosomes / linkage groups",
    numeric.labels && length(labels) < 200L ~
      "numeric sequences compatible with chromosomes / linkage groups",
    length(labels) >= 200L ~ "many sequences, likely scaffolds / contigs",
    .default = "genomic sequences"
  )
  list(levels = levels, type = type, n = length(labels))
}

# =============================================================================
# Window construction and GDS dosage handling
# =============================================================================

.inversion_make_windows <- function(
    marker.table, window.snps, step.snps, window.bp = NULL, step.bp = window.bp
) {
  # Split first: no subsequent operation is allowed to create a window that
  # crosses a chromosome, linkage-group, or scaffold boundary.
  split.markers <- split(marker.table, marker.table$chromosome, drop = TRUE)
  windows <- unlist(purrr::map(split.markers, function(x) {
    x <- x[order(x$position, x$variant_id), , drop = FALSE]
    if (!is.null(window.bp)) {
      first <- floor(min(x$position) / step.bp) * step.bp
      starts <- seq(first, max(x$position), by = step.bp)
      out <- purrr::map(starts, function(start) {
        idx <- which(x$position >= start & x$position < start + window.bp)
        if (!length(idx)) return(NULL)
        list(
          chromosome = x$chromosome[idx[1L]],
          start = start,
          end = start + window.bp - 1,
          variant_id = x$variant_id[idx]
        )
      })
      return(Filter(Negate(is.null), out))
    }
    if (nrow(x) < window.snps) return(list())
    starts <- seq.int(1L, nrow(x) - window.snps + 1L, by = step.snps)
    purrr::map(starts, function(first) {
      idx <- seq.int(first, first + window.snps - 1L)
      list(
        chromosome = x$chromosome[first],
        start = min(x$position[idx]),
        end = max(x$position[idx]),
        variant_id = x$variant_id[idx]
      )
    })
  }), recursive = FALSE)
  windows
}

# Experimental adaptive windows. Starting at the next unused marker, extend a
# window until the median r-squared between adjacent markers drops below the
# requested threshold, or until the maximum size is reached. This is a local
# scaling device, not an LD-block or breakpoint estimator.
.inversion_make_ld_windows <- function(
    gds, marker.table, sample.id, min.snps, max.snps, threshold
) {
  split.markers <- split(marker.table, marker.table$chromosome, drop = TRUE)
  unlist(purrr::map(split.markers, function(x) {
    x <- x[order(x$position, x$variant_id), , drop = FALSE]
    if (nrow(x) < min.snps) return(list())
    windows <- list()
    first <- 1L
    while (first <= nrow(x) - min.snps + 1L) {
      last.available <- min(nrow(x), first + max.snps - 1L)
      idx <- seq.int(first, last.available)
      dosage <- .inversion_get_dosage(gds, x$variant_id[idx], sample.id)
      if (ncol(dosage) > min.snps) {
        adjacent.r2 <- purrr::map_dbl(seq_len(ncol(dosage) - 1L), function(j) {
          value <- suppressWarnings(stats::cor(
            dosage[, j], dosage[, j + 1L], use = "pairwise.complete.obs"
          ))^2
          if (is.finite(value)) value else NA_real_
        })
        rolling <- purrr::map_dbl(seq.int(min.snps, ncol(dosage)), function(n) {
          values <- adjacent.r2[seq_len(n - 1L)]
          if (all(!is.finite(values))) NA_real_ else
            stats::median(values[is.finite(values)])
        })
        stop.at <- which(is.finite(rolling) & rolling < threshold)[1L]
        n.use <- if (is.na(stop.at)) ncol(dosage) else min.snps + stop.at - 1L
      } else {
        n.use <- ncol(dosage)
      }
      use <- idx[seq_len(max(min.snps, n.use))]
      windows[[length(windows) + 1L]] <- list(
        chromosome = x$chromosome[use[1L]],
        start = min(x$position[use]),
        end = max(x$position[use]),
        variant_id = x$variant_id[use]
      )
      first <- max(use) + 1L
    }
    windows
  }), recursive = FALSE)
}

.inversion_get_dosage <- function(gds, variant.id, sample.id) {
  SeqArray::seqFilterPush(gds)
  on.exit(SeqArray::seqFilterPop(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = variant.id, verbose = FALSE)
  dosage <- SeqArray::seqGetData(gds, "$dosage_alt")
  dosage <- as.matrix(dosage)

  if (nrow(dosage) == length(sample.id)) {
    out <- dosage
  } else if (ncol(dosage) == length(sample.id)) {
    out <- t(dosage)
  } else {
    rlang::abort("The GDS dosage matrix dimensions do not match its sample IDs.")
  }
  storage.mode(out) <- "double"
  rownames(out) <- sample.id
  colnames(out) <- as.character(SeqArray::seqGetData(gds, "variant.id"))
  out
}

.inversion_get_coverage <- function(gds, variant.id, sample.id) {
  SeqArray::seqFilterPush(gds)
  on.exit(SeqArray::seqFilterPop(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = variant.id, verbose = FALSE)
  normalise <- function(x) {
    if (!is.matrix(x)) return(NULL)
    if (nrow(x) == length(sample.id)) return(x)
    if (ncol(x) == length(sample.id)) return(t(x))
    NULL
  }
  depth <- tryCatch(
    normalise(SeqArray::seqGetData(gds, "annotation/format/DP")),
    error = function(error) NULL
  )
  ad <- tryCatch(
    SeqArray::seqGetData(gds, "annotation/format/AD"),
    error = function(error) NULL
  )
  allele.balance <- NULL
  if (is.list(ad) && all(c("length", "data") %in% names(ad)) &&
      all(ad$length == 2L)) {
    ad.data <- normalise(ad$data)
    if (!is.null(ad.data) && ncol(ad.data) == 2L * length(variant.id)) {
      ref <- ad.data[, seq.int(1L, ncol(ad.data), by = 2L), drop = FALSE]
      alt <- ad.data[, seq.int(2L, ncol(ad.data), by = 2L), drop = FALSE]
      total <- ref + alt
      allele.balance <- alt / total
      allele.balance[!is.finite(allele.balance)] <- NA_real_
    }
  }
  source <- c(
    if (!is.null(depth)) "annotation/format/DP",
    if (!is.null(allele.balance)) "annotation/format/AD"
  )
  if (is.null(depth) || is.null(allele.balance)) {
    embedded <- .inversion_get_embedded_coverage(gds, variant.id, sample.id)
    if (is.null(depth)) depth <- embedded$depth
    if (is.null(allele.balance)) allele.balance <- embedded$allele.balance
    source <- c(source, embedded$source)
  }
  list(
    depth = depth,
    allele.balance = allele.balance,
    source = unique(source)
  )
}

.inversion_get_embedded_coverage <- function(gds, variant.id, sample.id) {
  empty <- list(depth = NULL, allele.balance = NULL, source = character())
  metadata.node <- NULL
  for (path in c("genometranslator/genotypes.meta", "radiator/genotypes.meta")) {
    candidate <- gdsfmt::index.gdsn(gds$root, path = path, silent = TRUE)
    if (!is.null(candidate)) {
      metadata.node <- candidate
      break
    }
  }
  if (is.null(metadata.node)) return(empty)
  available <- gdsfmt::ls.gdsn(metadata.node)
  want.depth <- "READ_DEPTH" %in% available
  want.alleles <- all(c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH") %in% available)
  if (!want.depth && !want.alleles) return(empty)

  full.sample.id <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(gds$root, "sample.id")
  ))
  full.variant.id <- gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(gds$root, "variant.id")
  )
  sample.index <- match(sample.id, full.sample.id)
  variant.index <- match(variant.id, full.variant.id)
  if (anyNA(sample.index) || anyNA(variant.index)) return(empty)

  # genometranslator writes long genotype metadata marker by marker. Confirm
  # that layout before using efficient contiguous reads; never guess a layout.
  marker.node <- gdsfmt::index.gdsn(metadata.node, "M_SEQ", silent = TRUE)
  if (is.null(marker.node)) return(empty)
  first.block <- gdsfmt::read.gdsn(
    marker.node,
    start = 1L,
    count = min(length(full.sample.id), gdsfmt::objdesp.gdsn(marker.node)$dim)
  )
  if (length(unique(first.block)) != 1L) return(empty)

  read.matrix <- function(node.name) {
    node <- gdsfmt::index.gdsn(metadata.node, node.name, silent = TRUE)
    if (is.null(node)) return(NULL)
    values <- purrr::map(variant.index, function(index) {
      block <- gdsfmt::read.gdsn(
        node,
        start = as.integer((index - 1L) * length(full.sample.id) + 1L),
        count = length(full.sample.id)
      )
      block[sample.index]
    })
    matrix(
      unlist(values, use.names = FALSE),
      nrow = length(sample.id), ncol = length(variant.id)
    )
  }

  depth <- if (want.depth) read.matrix("READ_DEPTH") else NULL
  allele.balance <- NULL
  if (want.alleles) {
    ref <- read.matrix("ALLELE_REF_DEPTH")
    alt <- read.matrix("ALLELE_ALT_DEPTH")
    allele.balance <- alt / (ref + alt)
    allele.balance[!is.finite(allele.balance)] <- NA_real_
    if (is.null(depth)) depth <- ref + alt
  }
  list(
    depth = depth,
    allele.balance = allele.balance,
    source = c(
      if (want.depth) "genotypes.meta/READ_DEPTH",
      if (want.alleles) paste0(
        "genotypes.meta/", c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH")
      )
    )
  )
}

.inversion_prepare_dosage <- function(dosage, min.call.rate) {
  # Marker filtering happens before imputation and is repeated independently
  # for every window or candidate region.
  variant.index <- seq_len(ncol(dosage))
  call.rate <- colMeans(!is.na(dosage))
  keep <- is.finite(call.rate) & call.rate >= min.call.rate
  variant.index <- variant.index[keep]
  dosage <- dosage[, keep, drop = FALSE]
  call.rate <- call.rate[keep]
  if (ncol(dosage) == 0L) {
    return(list(
      dosage = dosage, observed = dosage, call.rate = call.rate,
      variant.index = variant.index
    ))
  }

  observed <- dosage
  means <- colMeans(dosage, na.rm = TRUE)
  missing <- which(is.na(dosage), arr.ind = TRUE)
  if (nrow(missing)) {
    # PCA receives mean-imputed dosages, while `observed` below remains
    # untouched for heterozygosity and pairwise-complete LD calculations.
    dosage[missing] <- means[missing[, 2L]]
  }
  variance <- apply(dosage, 2L, stats::var)
  keep <- is.finite(variance) & variance > sqrt(.Machine$double.eps)
  variant.index <- variant.index[keep]
  list(
    dosage = dosage[, keep, drop = FALSE],
    observed = observed[, keep, drop = FALSE],
    call.rate = call.rate[keep],
    variant.index = variant.index
  )
}

.inversion_mean_ld <- function(dosage) {
  if (ncol(dosage) < 2L) return(NA_real_)
  r <- suppressWarnings(stats::cor(dosage, use = "pairwise.complete.obs"))
  values <- r[upper.tri(r)]^2
  if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
}

# =============================================================================
# Local-PCA window representation and candidate selection
# =============================================================================

.inversion_analyse_window <- function(
    gds, window, window.id, sample.id, n.pcs, min.call.rate,
    min.window.snps
) {
  dosage <- .inversion_get_dosage(gds, window$variant_id, sample.id)
  local <- .inversion_local_covariance(
    dosage = dosage,
    n.pcs = n.pcs,
    min.call.rate = min.call.rate,
    min.window.snps = min.window.snps
  )
  local$window_id <- window.id
  local$chromosome <- window$chromosome
  local$start <- window$start
  local$end <- window$end
  local$n_input_snps <- length(window$variant_id)
  local$variant_id <- window$variant_id
  local
}

.inversion_scan_windows <- function(
    gds, windows, sample.id, n.pcs, min.call.rate, min.window.snps,
    parallel.core, verbose
) {
  n.windows <- length(windows)
  workers <- min(as.integer(parallel.core), n.windows)
  progress.id <- if (verbose) {
    cli::cli_progress_bar(
      name = "Local-PCA windows",
      total = n.windows,
      format = paste0(
        "{cli::pb_name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} ",
        "| ETA {cli::pb_eta}"
      ),
      clear = FALSE
    )
  } else {
    NULL
  }
  on.exit({
    if (!is.null(progress.id)) cli::cli_progress_done(id = progress.id)
  }, add = TRUE)

  if (workers == 1L) {
    return(purrr::map(seq_len(n.windows), function(i) {
      out <- .inversion_analyse_window(
        gds = gds,
        window = windows[[i]],
        window.id = i,
        sample.id = sample.id,
        n.pcs = n.pcs,
        min.call.rate = min.call.rate,
        min.window.snps = min.window.snps
      )
      if (!is.null(progress.id)) cli::cli_progress_update(id = progress.id)
      out
    }))
  }

  if (verbose) {
    cli::cli_alert_info(
      "Using {workers} independent GDS workers for the window scan."
    )
  }
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
  helper.names <- c(
    ".inversion_analyse_window", ".inversion_get_dosage",
    ".inversion_local_covariance", ".inversion_prepare_dosage",
    ".inversion_mean_ld"
  )
  parallel::clusterExport(
    cluster,
    varlist = helper.names,
    envir = environment(.inversion_scan_windows)
  )
  gds.path <- gds$filename
  parallel::clusterCall(
    cluster,
    function(path, ids, window.list, axes, call.rate, minimum.snps) {
      worker.gds <- SeqArray::seqOpen(path)
      SeqArray::seqSetFilter(worker.gds, sample.id = ids, verbose = FALSE)
      assign(".radr_inversion_gds", worker.gds, envir = .GlobalEnv)
      assign(".radr_inversion_windows", window.list, envir = .GlobalEnv)
      assign(".radr_inversion_sample_id", ids, envir = .GlobalEnv)
      assign(".radr_inversion_n_pcs", axes, envir = .GlobalEnv)
      assign(".radr_inversion_min_call_rate", call.rate, envir = .GlobalEnv)
      assign(".radr_inversion_min_window_snps", minimum.snps, envir = .GlobalEnv)
      TRUE
    },
    path = gds.path,
    ids = sample.id,
    window.list = windows,
    axes = n.pcs,
    call.rate = min.call.rate,
    minimum.snps = min.window.snps
  )
  on.exit(try(parallel::clusterCall(cluster, function() {
    if (exists(".radr_inversion_gds", envir = .GlobalEnv)) {
      try(SeqArray::seqClose(get(".radr_inversion_gds", envir = .GlobalEnv)),
          silent = TRUE)
    }
    TRUE
  }), silent = TRUE), add = TRUE)

  # Submit several windows per worker between progress updates. This keeps the
  # progress bar responsive without making scheduler overhead dominate.
  batch.size <- max(workers, workers * 4L)
  batches <- split(
    seq_len(n.windows),
    ceiling(seq_len(n.windows) / batch.size)
  )
  results <- vector("list", n.windows)
  for (batch in batches) {
    batch.results <- parallel::parLapplyLB(cluster, batch, function(i) {
      .inversion_analyse_window(
        gds = get(".radr_inversion_gds", envir = .GlobalEnv),
        window = get(".radr_inversion_windows", envir = .GlobalEnv)[[i]],
        window.id = i,
        sample.id = get(".radr_inversion_sample_id", envir = .GlobalEnv),
        n.pcs = get(".radr_inversion_n_pcs", envir = .GlobalEnv),
        min.call.rate = get(
          ".radr_inversion_min_call_rate", envir = .GlobalEnv
        ),
        min.window.snps = get(
          ".radr_inversion_min_window_snps", envir = .GlobalEnv
        )
      )
    })
    results[batch] <- batch.results
    if (!is.null(progress.id)) {
      cli::cli_progress_update(id = progress.id, inc = length(batch))
    }
  }
  results
}

.inversion_chromosome_pca_one <- function(
    gds, markers, sample.id, max.snps, min.call.rate
) {
  markers <- markers[order(markers$position, markers$variant_id), , drop = FALSE]
  use <- unique(round(seq(
    1L, nrow(markers), length.out = min(nrow(markers), max.snps)
  )))
  dosage <- .inversion_get_dosage(
    gds, markers$variant_id[use], sample.id
  )
  prepared <- .inversion_prepare_dosage(dosage, min.call.rate)
  if (ncol(prepared$dosage) < 2L) return(NULL)
  pca <- stats::prcomp(
    prepared$dosage, center = TRUE, scale. = FALSE, rank. = 2L
  )
  axes <- pca$x[, seq_len(min(2L, ncol(pca$x))), drop = FALSE]
  if (ncol(axes) == 1L) axes <- cbind(axes, PC2 = 0)
  colnames(axes) <- c("PC1", "PC2")
  variance <- pca$sdev^2
  tibble::tibble(
    chromosome = as.character(markers$chromosome[1L]),
    individual = sample.id,
    PC1 = axes[, 1L],
    PC2 = axes[, 2L],
    pc1_variance = variance[1L] / sum(variance),
    pc2_variance = if (length(variance) >= 2L) variance[2L] / sum(variance) else 0,
    n_input_snps = length(use),
    n_used_snps = ncol(prepared$dosage)
  )
}

.inversion_chromosome_pca <- function(
    gds, marker.table, sample.id, max.snps, min.call.rate,
    parallel.core, verbose
) {
  sequence.layout <- .inversion_sequence_layout(marker.table$chromosome)
  chromosome.levels <- sequence.layout$levels
  marker.groups <- split(
    marker.table,
    factor(marker.table$chromosome, levels = chromosome.levels),
    drop = TRUE
  )
  n.groups <- length(marker.groups)
  workers <- min(as.integer(parallel.core), n.groups)
  progress.id <- if (verbose) {
    cli::cli_progress_bar(
      name = "Chromosome PCA",
      total = n.groups,
      format = paste0(
        "{cli::pb_name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} ",
        "| ETA {cli::pb_eta}"
      ),
      clear = FALSE
    )
  } else {
    NULL
  }
  on.exit({
    if (!is.null(progress.id)) cli::cli_progress_done(id = progress.id)
  }, add = TRUE)

  if (workers == 1L) {
    scores <- purrr::map_dfr(marker.groups, function(markers) {
      out <- .inversion_chromosome_pca_one(
        gds, markers, sample.id, max.snps, min.call.rate
      )
      if (!is.null(progress.id)) cli::cli_progress_update(id = progress.id)
      out
    })
  } else {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
    parallel::clusterExport(
      cluster,
      varlist = c(
        ".inversion_chromosome_pca_one", ".inversion_get_dosage",
        ".inversion_prepare_dosage"
      ),
      envir = environment(.inversion_chromosome_pca)
    )
    parallel::clusterCall(
      cluster,
      function(path, ids, groups, maximum, call.rate) {
        worker.gds <- SeqArray::seqOpen(path)
        SeqArray::seqSetFilter(worker.gds, sample.id = ids, verbose = FALSE)
        assign(".radr_chromosome_pca_gds", worker.gds, envir = .GlobalEnv)
        assign(".radr_chromosome_pca_groups", groups, envir = .GlobalEnv)
        assign(".radr_chromosome_pca_ids", ids, envir = .GlobalEnv)
        assign(".radr_chromosome_pca_max", maximum, envir = .GlobalEnv)
        assign(".radr_chromosome_pca_call_rate", call.rate, envir = .GlobalEnv)
        TRUE
      },
      path = gds$filename,
      ids = sample.id,
      groups = marker.groups,
      maximum = max.snps,
      call.rate = min.call.rate
    )
    on.exit(try(parallel::clusterCall(cluster, function() {
      if (exists(".radr_chromosome_pca_gds", envir = .GlobalEnv)) {
        try(SeqArray::seqClose(
          get(".radr_chromosome_pca_gds", envir = .GlobalEnv)
        ), silent = TRUE)
      }
      TRUE
    }), silent = TRUE), add = TRUE)
    indices <- seq_len(n.groups)
    batches <- split(indices, ceiling(indices / max(workers, workers * 2L)))
    results <- vector("list", n.groups)
    for (batch in batches) {
      batch.results <- parallel::parLapplyLB(cluster, batch, function(i) {
        .inversion_chromosome_pca_one(
          gds = get(".radr_chromosome_pca_gds", envir = .GlobalEnv),
          markers = get(".radr_chromosome_pca_groups", envir = .GlobalEnv)[[i]],
          sample.id = get(".radr_chromosome_pca_ids", envir = .GlobalEnv),
          max.snps = get(".radr_chromosome_pca_max", envir = .GlobalEnv),
          min.call.rate = get(
            ".radr_chromosome_pca_call_rate", envir = .GlobalEnv
          )
        )
      })
      results[batch] <- batch.results
      if (!is.null(progress.id)) {
        cli::cli_progress_update(id = progress.id, inc = length(batch))
      }
    }
    scores <- dplyr::bind_rows(results)
  }
  if (!nrow(scores)) {
    return(list(
      scores = tibble::tibble(), summary = tibble::tibble(),
      sequence_layout = sequence.layout
    ))
  }
  scores$chromosome <- factor(scores$chromosome, levels = chromosome.levels)
  summary <- scores |>
    dplyr::distinct(
      .data$chromosome, .data$pc1_variance, .data$pc2_variance,
      .data$n_input_snps, .data$n_used_snps
    ) |>
    dplyr::arrange(.data$chromosome)
  list(
    scores = scores, summary = summary,
    sequence_layout = sequence.layout
  )
}

.inversion_local_covariance <- function(
    dosage, n.pcs, min.call.rate, min.window.snps
) {
  prepared <- .inversion_prepare_dosage(dosage, min.call.rate)
  dosage <- prepared$dosage
  base <- list(
    valid = FALSE,
    n_used_snps = ncol(dosage),
    mean_call_rate = if (length(prepared$call.rate)) {
      mean(prepared$call.rate)
    } else {
      NA_real_
    },
    mean_ld_r2 = .inversion_mean_ld(prepared$observed),
    values = NULL,
    vectors = NULL
  )
  if (ncol(dosage) < min.window.snps) return(base)

  centred <- sweep(dosage, 2L, colMeans(dosage), "-")
  covariance <- tcrossprod(centred) / max(1, ncol(centred) - 1L)
  eig <- eigen(covariance, symmetric = TRUE)
  keep <- which(eig$values > sqrt(.Machine$double.eps))
  keep <- utils::head(keep, n.pcs)
  if (length(keep) < n.pcs) return(base)
  values <- eig$values[keep]
  values <- values / sum(abs(values))

  base$valid <- TRUE
  base$values <- values
  base$vectors <- eig$vectors[, keep, drop = FALSE]
  base
}

.inversion_covariance_distances <- function(results) {
  n <- length(results)
  out <- matrix(0, nrow = n, ncol = n)
  # A symmetric matrix is filled explicitly here. A nested loop avoids
  # allocating every window-pair covariance product at once.
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      cross <- crossprod(results[[i]]$vectors, results[[j]]$vectors)
      cross.term <- sum(
        outer(results[[i]]$values, results[[j]]$values) * cross^2
      )
      squared <- sum(results[[i]]$values^2) +
        sum(results[[j]]$values^2) - 2 * cross.term
      out[i, j] <- out[j, i] <- sqrt(max(0, squared))
    }
  }
  out
}

.inversion_robust_score <- function(x) {
  centre <- apply(x, 2L, stats::median, na.rm = TRUE)
  scale <- apply(x, 2L, stats::mad, constant = 1.4826, na.rm = TRUE)
  fallback <- apply(x, 2L, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <-
    fallback[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)]
  scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <- 1
  z <- sweep(sweep(x, 2L, centre, "-"), 2L, scale, "/")
  sqrt(rowSums(z^2))
}

.inversion_candidate_regions <- function(window.table, min.windows) {
  candidate <- window.table[window.table$candidate_window, , drop = FALSE]
  empty <- data.frame(
    candidate_id = character(),
    chromosome = character(),
    start = numeric(),
    end = numeric(),
    size_bp = numeric(),
    n_windows = integer(),
    max_robust_score = numeric(),
    n_samples = integer(),
    n_snps = integer(),
    regional_mean_ld_r2 = numeric(),
    homokaryotype_mean_ld_r2 = numeric(),
    heterokaryotype_mean_ld_r2 = numeric(),
    ld_structure_contrast = numeric(),
    flanking_mean_ld_r2 = numeric(),
    candidate_median_score = numeric(),
    flanking_max_score = numeric(),
    boundary_contrast = numeric(),
    left_boundary_contrast = numeric(),
    right_boundary_contrast = numeric(),
    internal_transition_max = numeric(),
    pc1_variance = numeric(),
    smallest_cluster_n = integer(),
    smallest_cluster_frequency = numeric(),
    cluster_separation = numeric(),
    cluster_compactness = numeric(),
    three_cluster_evidence = logical(),
    heterozygote_like_middle_cluster = logical(),
    middle_heterozygosity_excess = numeric(),
    stringsAsFactors = FALSE
  )
  if (nrow(candidate) == 0L) return(empty)

  candidate <- candidate[order(
    .inversion_chromosome_order(candidate$chromosome),
    candidate$window_id
  ), , drop = FALSE]
  previous <- c(FALSE,
    candidate$chromosome[-1L] == candidate$chromosome[-nrow(candidate)] &
      candidate$window_id[-1L] == candidate$window_id[-nrow(candidate)] + 1L
  )
  group <- cumsum(!previous)
  groups <- split(candidate, group)
  groups <- purrr::keep(groups, ~ nrow(.x) >= min.windows)
  if (length(groups) == 0L) return(empty)

  purrr::map2_dfr(groups, seq_along(groups), function(x, i) {
    # Immediate noncandidate neighbours define the boundary and LD contrasts.
    # Missing left/right flanks at chromosome ends remain NA rather than being
    # borrowed from another chromosome.
    chromosome.windows <- window.table[
      window.table$chromosome == x$chromosome[1] & window.table$valid,
      , drop = FALSE
    ]
    candidate.index <- match(x$window_id, chromosome.windows$window_id)
    left.index <- min(candidate.index) - 1L
    right.index <- max(candidate.index) + 1L
    left.score <- if (left.index >= 1L) chromosome.windows$robust_score[left.index] else NA_real_
    right.score <- if (right.index <= nrow(chromosome.windows)) chromosome.windows$robust_score[right.index] else NA_real_
    flanking.scores <- c(left.score, right.score)
    flanking.ld <- c(
      if (left.index >= 1L) chromosome.windows$mean_ld_r2[left.index] else NA_real_,
      if (right.index <= nrow(chromosome.windows)) chromosome.windows$mean_ld_r2[right.index] else NA_real_
    )
    candidate.median <- stats::median(x$robust_score, na.rm = TRUE)
    flank.max <- if (all(!is.finite(flanking.scores))) NA_real_ else
      max(flanking.scores, na.rm = TRUE)
    internal.transition <- if (nrow(x) > 1L) {
      max(abs(diff(x$robust_score)), na.rm = TRUE)
    } else {
      NA_real_
    }
    tibble::tibble(
      candidate_id = paste0("INV-CAND-", i),
      chromosome = x$chromosome[1],
      start = min(x$start),
      end = max(x$end),
      size_bp = max(x$end) - min(x$start) + 1,
      n_windows = nrow(x),
      max_robust_score = max(x$robust_score),
      n_samples = NA_integer_,
      n_snps = NA_integer_,
      regional_mean_ld_r2 = NA_real_,
      homokaryotype_mean_ld_r2 = NA_real_,
      heterokaryotype_mean_ld_r2 = NA_real_,
      ld_structure_contrast = NA_real_,
      flanking_mean_ld_r2 = if (all(!is.finite(flanking.ld))) NA_real_ else
        mean(flanking.ld, na.rm = TRUE),
      candidate_median_score = candidate.median,
      flanking_max_score = flank.max,
      boundary_contrast = candidate.median - flank.max,
      left_boundary_contrast = candidate.median - left.score,
      right_boundary_contrast = candidate.median - right.score,
      internal_transition_max = internal.transition,
      pc1_variance = NA_real_,
      smallest_cluster_n = NA_integer_,
      smallest_cluster_frequency = NA_real_,
      cluster_separation = NA_real_,
      cluster_compactness = NA_real_,
      three_cluster_evidence = FALSE,
      heterozygote_like_middle_cluster = FALSE,
      middle_heterozygosity_excess = NA_real_
    )
  })
}

.inversion_region_diagnostics <- function(
    dosage, depth = NULL, allele.balance = NULL, sample.id, cluster.k,
    min.call.rate, ld.max.snps, return.ld,
    arrangement.labels, sample.metadata = NULL, stability.replicates = 0L,
    stability.fraction = 0.8, verbose = FALSE
) {
  # Candidate-level PCA uses every usable SNP in the joined candidate interval,
  # not only the smaller covariance representation used for window scoring.
  if (verbose) message("  Preparing genotypes and applying the call-rate check...")
  prepared <- .inversion_prepare_dosage(dosage, min.call.rate)
  dosage <- prepared$dosage
  if (!is.null(depth)) {
    depth <- depth[, prepared$variant.index, drop = FALSE]
  }
  if (!is.null(allele.balance)) {
    allele.balance <- allele.balance[, prepared$variant.index, drop = FALSE]
  }
  if (ncol(dosage) < 2L) {
    rlang::abort("A candidate region contains fewer than two usable SNPs.")
  }

  if (verbose) {
    message(
      "  Regional PCA and cluster evaluation with ", ncol(dosage),
      " usable SNPs..."
    )
  }
  centred <- sweep(dosage, 2L, colMeans(dosage), "-")
  pca <- stats::prcomp(centred, center = FALSE, scale. = FALSE, rank. = 2L)
  scores <- as.data.frame(pca$x[, seq_len(min(2L, ncol(pca$x))), drop = FALSE])
  names(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores$individual <- sample.id

  pc1 <- pca$x[, 1L]
  pca.variance <- pca$sdev^2
  pc1.variance <- pca.variance[1L] / sum(pca.variance)
  distinct.pc1 <- length(unique(signif(pc1, digits = 12L)))
  cluster.models <- .inversion_cluster_models(pc1, max.k = 3L)
  if (distinct.pc1 >= cluster.k) {
    kmeans.fit <- stats::kmeans(
      pc1, centers = cluster.k, nstart = 50L
    )
    cluster <- kmeans.fit$cluster
  } else {
    cluster <- rep(1L, length(pc1))
  }
  cluster.means <- tapply(pc1, cluster, mean)
  ordered.cluster <- order(cluster.means)
  cluster.labels <- match(cluster, ordered.cluster)
  scores$cluster <- factor(cluster.labels, levels = seq_len(cluster.k))

  centre <- as.numeric(tapply(pc1, scores$cluster, mean))
  distances <- abs(outer(pc1, centre, "-"))
  nearest <- apply(distances, 1L, min)
  second <- apply(distances, 1L, function(z) sort(z, partial = 2L)[2L])
  confidence <- second / pmax(nearest + second, sqrt(.Machine$double.eps))
  scores$distance_nearest <- nearest
  scores$distance_second <- second
  scores$arrangement_confidence <- confidence

  cluster.counts <- as.integer(table(scores$cluster))
  cluster.frequencies <- cluster.counts / length(pc1)
  ordered.centres <- as.numeric(tapply(pc1, scores$cluster, mean))
  within.ss <- sum(purrr::map_dbl(seq_len(cluster.k), function(level) {
    values <- pc1[scores$cluster == level]
    if (!length(values)) return(0)
    sum((values - mean(values))^2)
  }))
  residual.df <- max(1L, length(pc1) - sum(cluster.counts > 0L))
  pooled.within.sd <- sqrt(within.ss / residual.df)
  adjacent.gaps <- diff(ordered.centres)
  cluster.separation <- if (length(adjacent.gaps) && all(is.finite(adjacent.gaps))) {
    min(adjacent.gaps) / max(pooled.within.sd, sqrt(.Machine$double.eps))
  } else {
    NA_real_
  }
  total.ss <- sum((pc1 - mean(pc1))^2)
  cluster.compactness <- if (is.finite(total.ss) && total.ss > 0) {
    1 - within.ss / total.ss
  } else {
    NA_real_
  }
  smallest.cluster.n <- min(cluster.counts)
  smallest.cluster.frequency <- min(cluster.frequencies)
  # Asking k-means for three groups does not establish a three-karyotype
  # pattern. Require usable group sizes and separation relative to residual
  # within-group spread before recording quantitative support.
  three.cluster.evidence <- cluster.k == 3L &&
    all(cluster.counts >= 3L) &&
    smallest.cluster.frequency >= 0.05 &&
    is.finite(cluster.separation) && cluster.separation >= 1

  # AA, AB, and BB communicate the three-arrangement biological hypothesis,
  # so reserve them for candidates with quantitative three-cluster support.
  # Numeric cluster IDs remain available for reproducibility in every case.
  if (three.cluster.evidence) {
    display.labels <- arrangement.labels
    scores$arrangement <- factor(
      display.labels[cluster.labels], levels = display.labels
    )
    scores$arrangement_dosage <- cluster.labels - 1L
  } else {
    display.labels <- paste("Group", seq_len(cluster.k))
    scores$arrangement <- factor(
      display.labels[cluster.labels], levels = display.labels
    )
    scores$arrangement_dosage <- NA_integer_
  }

  observed <- prepared$observed
  heterozygosity <- rowMeans(observed == 1, na.rm = TRUE)
  scores$heterozygosity <- heterozygosity
  scores$call_rate <- rowMeans(!is.na(observed))
  scores$mean_depth <- if (!is.null(depth)) {
    rowMeans(depth, na.rm = TRUE)
  } else {
    NA_real_
  }
  scores$mean_heterozygote_allele_balance <- if (!is.null(allele.balance)) {
    heterozygous <- observed == 1
    values <- allele.balance
    values[!heterozygous] <- NA_real_
    rowMeans(values, na.rm = TRUE)
  } else {
    NA_real_
  }
  scores$mean_depth[!is.finite(scores$mean_depth)] <- NA_real_
  scores$mean_heterozygote_allele_balance[
    !is.finite(scores$mean_heterozygote_allele_balance)
  ] <- NA_real_
  heterozygosity.by.cluster <- tapply(
    heterozygosity, scores$cluster, mean, na.rm = TRUE
  )
  middle <- ceiling(cluster.k / 2)
  heterozygote.like <- cluster.k == 3L &&
    length(heterozygosity.by.cluster) == 3L &&
    is.finite(heterozygosity.by.cluster[middle]) &&
    heterozygosity.by.cluster[middle] == max(heterozygosity.by.cluster)
  outer <- c(1L, cluster.k)
  middle.heterozygosity.excess <- if (
    cluster.k == 3L && all(is.finite(heterozygosity.by.cluster))
  ) {
    heterozygosity.by.cluster[middle] -
      mean(heterozygosity.by.cluster[outer])
  } else {
    NA_real_
  }

  ld.index <- unique(round(seq(1, ncol(dosage), length.out = min(
    ncol(dosage), ld.max.snps
  ))))
  if (verbose) {
    message(
      "  LD diagnostics with ", length(ld.index),
      " evenly distributed SNPs..."
    )
  }
  observed.ld <- observed[, ld.index, drop = FALSE]
  ld.all <- .inversion_ld_matrix(observed.ld)
  mean.ld <- .inversion_ld_summary(ld.all)

  cluster.levels <- levels(scores$cluster)
  ld.by.cluster <- purrr::map(cluster.levels, function(level) {
    .inversion_ld_matrix(observed.ld[scores$cluster == level, , drop = FALSE])
  }) |>
    rlang::set_names(paste0("cluster_", cluster.levels))

  outer.counts <- table(scores$cluster)[outer]
  common.outer <- outer[which.max(outer.counts)]
  common.name <- paste0("cluster_", common.outer)
  ld.common <- ld.by.cluster[[common.name]]
  ld.summary <- tibble::tibble(
    group = c("all", names(ld.by.cluster), "common_homokaryotype"),
    n_samples = c(
      nrow(observed.ld),
      as.integer(table(scores$cluster)),
      as.integer(table(scores$cluster)[common.outer])
    ),
    mean_r2 = c(
      mean.ld,
      purrr::map_dbl(ld.by.cluster, .inversion_ld_summary),
      .inversion_ld_summary(ld.common)
    )
  )
  outer.ld <- purrr::map_dbl(
    ld.by.cluster[outer],
    .inversion_ld_summary
  )
  homokaryotype.mean.ld <- if (all(!is.finite(outer.ld))) NA_real_ else
    stats::weighted.mean(
      outer.ld[is.finite(outer.ld)],
      cluster.counts[outer][is.finite(outer.ld)]
    )
  heterokaryotype.mean.ld <- if (cluster.k == 3L) {
    .inversion_ld_summary(ld.by.cluster[[middle]])
  } else {
    NA_real_
  }
  ld.structure.contrast <- mean.ld - homokaryotype.mean.ld

  if (verbose && stability.replicates > 0L) {
    message(
      "  Assignment-stability resampling: ", stability.replicates,
      " replicates..."
    )
  }
  stability <- .inversion_assignment_stability(
    dosage = dosage,
    reference.cluster = cluster.labels,
    cluster.k = cluster.k,
    replicates = stability.replicates,
    fraction = stability.fraction
  )
  scores$assignment_stability <- stability$individual

  loadings <- tibble::tibble(
    variant_id = colnames(dosage),
    PC1_loading = pca$rotation[, 1L],
    abs_PC1_loading = abs(pca$rotation[, 1L]),
    loading_rank = rank(-abs(pca$rotation[, 1L]), ties.method = "first")
  ) |>
    dplyr::arrange(.data$loading_rank)

  if (verbose) message("  Marker loadings, technical summaries, and metadata audit...")
  metadata.audit <- .inversion_metadata_audit(
    scores = scores,
    sample.metadata = sample.metadata
  )
  arrangement.differentiation <- .inversion_arrangement_differentiation(
    observed = observed,
    cluster = cluster.labels,
    labels = display.labels
  )
  technical.summary <- scores |>
    dplyr::group_by(.data$arrangement) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_call_rate = mean(.data$call_rate, na.rm = TRUE),
      mean_depth = mean(.data$mean_depth, na.rm = TRUE),
      mean_heterozygote_allele_balance = mean(
        .data$mean_heterozygote_allele_balance, na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(
      tidyselect::where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)
    ))

  list(
    n_snps = ncol(dosage),
    scores = tibble::as_tibble(scores),
    cluster_summary = tibble::tibble(
      cluster = factor(seq_len(cluster.k)),
      n = cluster.counts,
      frequency = cluster.frequencies,
      mean_pc1 = ordered.centres,
      sd_pc1 = purrr::map_dbl(seq_len(cluster.k), function(level) {
        values <- pc1[scores$cluster == level]
        if (length(values) < 2L) NA_real_ else stats::sd(values)
      }),
      mean_heterozygosity = as.numeric(heterozygosity.by.cluster)
    ),
    all_clusters_present = all(table(scores$cluster) > 0L),
    heterozygote_like_middle_cluster = heterozygote.like,
    middle_heterozygosity_excess = middle.heterozygosity.excess,
    pc1_variance = pc1.variance,
    smallest_cluster_n = smallest.cluster.n,
    smallest_cluster_frequency = smallest.cluster.frequency,
    cluster_separation = cluster.separation,
    cluster_compactness = cluster.compactness,
    three_cluster_evidence = three.cluster.evidence,
    cluster_models = cluster.models,
    assignment_stability = stability$overall,
    assignment_stability_replicates = stability$replicates,
    metadata_audit = metadata.audit,
    pc1_loadings = loadings,
    arrangement_differentiation = arrangement.differentiation,
    technical_summary = technical.summary,
    mean_ld_r2 = mean.ld,
    homokaryotype_mean_ld_r2 = homokaryotype.mean.ld,
    heterokaryotype_mean_ld_r2 = heterokaryotype.mean.ld,
    ld_structure_contrast = ld.structure.contrast,
    ld_variant_index = ld.index,
    ld_variant_id = colnames(observed)[ld.index],
    common_homokaryotype_cluster = common.outer,
    ld_summary = ld.summary,
    ld_matrices = c(
      list(all = ld.all),
      ld.by.cluster,
      list(common_homokaryotype = ld.common)
    )
  )
}

.inversion_cluster_models <- function(pc1, max.k = 3L) {
  max.k <- min(max.k, length(unique(signif(pc1, 12L))), length(pc1) - 1L)
  purrr::map_dfr(seq_len(max.k), function(k) {
    fit <- if (k == 1L) {
      list(tot.withinss = sum((pc1 - mean(pc1))^2), cluster = rep(1L, length(pc1)))
    } else {
      stats::kmeans(pc1, centers = k, nstart = 50L)
    }
    rss <- max(fit$tot.withinss, .Machine$double.eps)
    # Gaussian equal-variance approximation. It is used only to compare the
    # number of modes in regional PC1, not as formal inversion evidence.
    bic <- length(pc1) * log(rss / length(pc1)) + (2L * k) * log(length(pc1))
    tibble::tibble(k = k, bic = bic, delta_bic = NA_real_)
  }) |>
    dplyr::mutate(delta_bic = .data$bic - min(.data$bic))
}

.inversion_assignment_stability <- function(
    dosage, reference.cluster, cluster.k, replicates, fraction
) {
  empty <- list(
    individual = rep(NA_real_, nrow(dosage)),
    overall = NA_real_,
    replicates = tibble::tibble()
  )
  if (replicates < 1L || ncol(dosage) < 2L) return(empty)
  n.use <- max(2L, round(ncol(dosage) * fraction))
  replicate.results <- purrr::map_dfr(seq_len(replicates), function(replicate) {
    columns <- sample.int(ncol(dosage), n.use, replace = TRUE)
    pca <- stats::prcomp(
      dosage[, columns, drop = FALSE], center = TRUE, scale. = FALSE, rank. = 1L
    )
    pc1 <- pca$x[, 1L]
    if (length(unique(signif(pc1, 12L))) < cluster.k) return(NULL)
    fit <- stats::kmeans(pc1, centers = cluster.k, nstart = 20L)
    ordered <- match(fit$cluster, order(tapply(pc1, fit$cluster, mean)))
    reversed <- cluster.k + 1L - ordered
    if (sum(reversed == reference.cluster) > sum(ordered == reference.cluster)) {
      ordered <- reversed
    }
    tibble::tibble(
      replicate = replicate,
      individual_index = seq_along(ordered),
      stable = ordered == reference.cluster
    )
  })
  if (!nrow(replicate.results)) return(empty)
  individual <- replicate.results |>
    dplyr::group_by(.data$individual_index) |>
    dplyr::summarise(stability = mean(.data$stable), .groups = "drop")
  values <- rep(NA_real_, nrow(dosage))
  values[individual$individual_index] <- individual$stability
  list(
    individual = values,
    overall = mean(values, na.rm = TRUE),
    replicates = replicate.results
  )
}

.inversion_metadata_audit <- function(scores, sample.metadata) {
  empty <- tibble::tibble(
    variable = character(), type = character(), n_levels = integer(),
    pc1_association = numeric(), pc1_p_value = numeric(),
    arrangement_association = numeric(), arrangement_p_value = numeric(),
    interpretation = character()
  )
  if (is.null(sample.metadata) || ncol(sample.metadata) < 2L) return(empty)
  joined <- dplyr::left_join(
    scores,
    sample.metadata,
    by = c("individual" = "INDIVIDUALS")
  )
  variables <- setdiff(names(sample.metadata), "INDIVIDUALS")
  purrr::map_dfr(variables, function(variable) {
    x <- joined[[variable]]
    keep <- !is.na(x) & nzchar(as.character(x)) & is.finite(joined$PC1)
    if (sum(keep) < 3L) return(NULL)
    numeric.x <- if (is.numeric(x)) as.numeric(x[keep]) else rep(NA_real_, sum(keep))
    is.numeric.variable <- is.numeric(x) &&
      length(unique(numeric.x[is.finite(numeric.x)])) > 2L
    if (is.numeric.variable) {
      fit <- stats::lm(joined$PC1[keep] ~ numeric.x)
      arrangement.fit <- stats::anova(
        stats::lm(numeric.x ~ joined$arrangement[keep])
      )
      tibble::tibble(
        variable = variable,
        type = "numeric",
        n_levels = length(unique(numeric.x)),
        pc1_association = summary(fit)$r.squared,
        pc1_p_value = stats::coef(summary(fit))[2L, 4L],
        arrangement_association = arrangement.fit$`Sum Sq`[1L] /
          sum(arrangement.fit$`Sum Sq`),
        arrangement_p_value = arrangement.fit$`Pr(>F)`[1L],
        interpretation = paste(
          "PC1 R-squared; eta-squared for numeric metadata among arrangements"
        )
      )
    } else {
      group <- factor(x[keep])
      if (nlevels(group) < 2L) return(NULL)
      fit <- stats::anova(stats::lm(joined$PC1[keep] ~ group))
      eta2 <- fit$`Sum Sq`[1L] / sum(fit$`Sum Sq`)
      contingency <- table(group, joined$arrangement[keep])
      chi <- suppressWarnings(stats::chisq.test(contingency))
      cramers.v <- sqrt(
        unname(chi$statistic) /
          (sum(contingency) * max(1L, min(dim(contingency)) - 1L))
      )
      tibble::tibble(
        variable = variable,
        type = "categorical",
        n_levels = nlevels(group),
        pc1_association = eta2,
        pc1_p_value = fit$`Pr(>F)`[1L],
        arrangement_association = cramers.v,
        arrangement_p_value = chi$p.value,
        interpretation = paste(
          "PC1 eta-squared; Cramer's V for metadata by arrangement"
        )
      )
    }
  })
}

.inversion_arrangement_differentiation <- function(observed, cluster, labels) {
  pairs <- utils::combn(sort(unique(cluster)), 2L, simplify = FALSE)
  purrr::map_dfr(pairs, function(pair) {
    x1 <- observed[cluster == pair[1L], , drop = FALSE]
    x2 <- observed[cluster == pair[2L], , drop = FALSE]
    n1 <- 2 * colSums(!is.na(x1))
    n2 <- 2 * colSums(!is.na(x2))
    p1 <- colMeans(x1, na.rm = TRUE) / 2
    p2 <- colMeans(x2, na.rm = TRUE) / 2
    valid <- n1 > 1 & n2 > 1 & is.finite(p1) & is.finite(p2)
    numerator <- (p1 - p2)^2 - p1 * (1 - p1) / (n1 - 1) -
      p2 * (1 - p2) / (n2 - 1)
    denominator <- p1 * (1 - p2) + p2 * (1 - p1)
    hudson.fst <- if (any(valid & denominator > 0)) {
      sum(numerator[valid & denominator > 0]) /
        sum(denominator[valid & denominator > 0])
    } else {
      NA_real_
    }
    tibble::tibble(
      arrangement_1 = labels[pair[1L]],
      arrangement_2 = labels[pair[2L]],
      n_1 = nrow(x1),
      n_2 = nrow(x2),
      n_snps = sum(valid),
      mean_absolute_allele_frequency_difference = mean(abs(p1[valid] - p2[valid])),
      hudson_fst = hudson.fst,
      dxy = NA_real_,
      dxy_note = paste(
        "Not calculated: inferred arrangement genotypes are not treated as",
        "independently sampled populations."
      )
    )
  })
}

# =============================================================================
# Linkage-disequilibrium and sensitivity helpers
# =============================================================================

.inversion_ld_matrix <- function(dosage) {
  n.markers <- ncol(dosage)
  empty <- matrix(NA_real_, nrow = n.markers, ncol = n.markers)
  if (nrow(dosage) < 3L || n.markers < 2L) return(empty)
  marker.sd <- apply(dosage, 2L, stats::sd, na.rm = TRUE)
  variable <- is.finite(marker.sd) & marker.sd > sqrt(.Machine$double.eps)
  if (sum(variable) < 2L) return(empty)
  empty[variable, variable] <- suppressWarnings(
    stats::cor(
      dosage[, variable, drop = FALSE],
      use = "pairwise.complete.obs"
    )^2
  )
  empty
}

.inversion_ld_summary <- function(ld) {
  if (ncol(ld) < 2L) return(NA_real_)
  values <- ld[upper.tri(ld)]
  if (!length(values) || all(!is.finite(values))) return(NA_real_)
  mean(values[is.finite(values)])
}

.inversion_window_sensitivity <- function(
    gds, marker.table, sample.id, window.sizes, min.call.rate,
    min.window.snps
) {
  empty <- tibble::tibble(
    window_snps = integer(),
    chromosome = character(),
    start = numeric(),
    end = numeric(),
    n_used_snps = integer(),
    pc1_variance = numeric(),
    mean_ld_r2 = numeric()
  )
  if (is.null(window.sizes) || !length(window.sizes)) return(empty)

  # Sensitivity runs are descriptive and deliberately do not replace the
  # primary candidate calls. Empty/failed windows are dropped from this table.
  rows <- purrr::map_dfr(window.sizes, function(size) {
    windows <- .inversion_make_windows(
      marker.table = marker.table,
      window.snps = size,
      step.snps = size
    )
    if (!length(windows)) return(NULL)
    purrr::map_dfr(windows, function(w) {
      dosage <- .inversion_get_dosage(gds, w$variant_id, sample.id)
      prepared <- .inversion_prepare_dosage(dosage, min.call.rate)
      if (ncol(prepared$dosage) < min(min.window.snps, size)) return(NULL)
      pca <- stats::prcomp(
        prepared$dosage,
        center = TRUE,
        scale. = FALSE,
        rank. = 1L
      )
      variance <- pca$sdev^2
      tibble::tibble(
        window_snps = size,
        chromosome = w$chromosome,
        start = w$start,
        end = w$end,
        n_used_snps = ncol(prepared$dosage),
        pc1_variance = variance[1L] / sum(variance),
        mean_ld_r2 = .inversion_mean_ld(prepared$observed)
      )
    })
  })
  if (!nrow(rows)) return(empty)
  rows
}

.inversion_write_outputs <- function(
    path.folder, window.table, candidate.regions, diagnostics,
    candidate.summaries,
    arrangement.genotypes, homokaryotype.all.candidates, sensitivity,
    chromosome.pca, chromosome.lengths,
    threshold, save.plots, plot.formats, verbose
) {
  files <- character()
  # Centralise table writing so every generated path is returned to the user.
  write_table <- function(x, filename) {
    path <- file.path(path.folder, filename)
    readr::write_tsv(x, path, na = "NA")
    files <<- c(files, path)
  }
  write_table(window.table, "inversion_windows.tsv")
  write_table(candidate.regions, "candidate_inversion_regions.tsv")
  write_table(arrangement.genotypes, "inversion_arrangement_genotypes.tsv")
  write_table(chromosome.lengths, "chromosome_length_context.tsv")
  write_table(
    homokaryotype.all.candidates,
    "homokaryotypes_all_candidates_whitelist.tsv"
  )
  if (nrow(sensitivity) > 0L) {
    write_table(sensitivity, "window_size_sensitivity.tsv")
  }
  if (nrow(chromosome.pca$scores) > 0L) {
    chromosome.scores <- chromosome.pca$scores
    chromosome.scores$chromosome <- as.character(chromosome.scores$chromosome)
    chromosome.summary <- chromosome.pca$summary
    chromosome.summary$chromosome <- as.character(
      chromosome.summary$chromosome
    )
    write_table(chromosome.scores, "chromosome_pca_scores.tsv")
    write_table(chromosome.summary, "chromosome_pca_summary.tsv")
  }
  if (length(candidate.summaries)) {
    summary.path <- file.path(path.folder, "candidate_inversion_summary.txt")
    writeLines(
      paste(unname(candidate.summaries), collapse = "\n\n"),
      con = summary.path
    )
    files <- c(files, summary.path)
  }
  purrr::walk(diagnostics, function(diagnostic) {
    prefix <- diagnostic$candidate_id
    write_table(
      diagnostic$scores,
      paste0(prefix, "_individual_pca.tsv")
    )
    write_table(
      diagnostic$cluster_summary,
      paste0(prefix, "_cluster_summary.tsv")
    )
    write_table(
      diagnostic$ld_summary,
      paste0(prefix, "_ld_summary.tsv")
    )
    write_table(
      diagnostic$cluster_models,
      paste0(prefix, "_cluster_number_models.tsv")
    )
    write_table(
      diagnostic$metadata_audit,
      paste0(prefix, "_sample_metadata_audit.tsv")
    )
    write_table(
      diagnostic$pc1_loadings,
      paste0(prefix, "_regional_PC1_marker_loadings.tsv")
    )
    write_table(
      diagnostic$arrangement_differentiation,
      paste0(prefix, "_arrangement_differentiation.tsv")
    )
    write_table(
      diagnostic$technical_summary,
      paste0(prefix, "_coverage_missingness_allele_balance.tsv")
    )
    if (nrow(diagnostic$assignment_stability_replicates)) {
      write_table(
        diagnostic$assignment_stability_replicates,
        paste0(prefix, "_assignment_stability_replicates.tsv")
      )
    }
    arrangement.table <- diagnostic$scores |>
      dplyr::select(dplyr::all_of(c(
        "individual", "cluster", "arrangement", "arrangement_dosage",
        "arrangement_confidence", "distance_nearest", "distance_second",
        "heterozygosity", "call_rate", "mean_depth",
        "mean_heterozygote_allele_balance", "assignment_stability"
      )))
    write_table(
      arrangement.table,
      paste0(prefix, "_arrangement_genotypes.tsv")
    )
    purrr::walk(unique(as.character(arrangement.table$arrangement)), function(label) {
      write_table(
        dplyr::filter(
          arrangement.table, as.character(.data$arrangement) == label
        ) |>
          dplyr::select(.data$individual),
        paste0(prefix, "_", label, "_individuals_whitelist.tsv")
      )
    })
    write_table(
      dplyr::filter(
        arrangement.table, .data$arrangement_dosage %in% c(0L, 2L)
      ) |>
        dplyr::select(.data$individual, .data$arrangement),
      paste0(prefix, "_homokaryotypes_whitelist.tsv")
    )
  })

  plots <- list()
  if (save.plots) {
    sequence.layout <- .inversion_sequence_layout(window.table$chromosome)
    window.table$chromosome <- factor(
      window.table$chromosome,
      levels = sequence.layout$levels
    )
    if (verbose) {
      full.layout <- chromosome.pca$sequence_layout
      if (!is.null(full.layout) && full.layout$n != sequence.layout$n) {
        omitted <- setdiff(full.layout$levels, sequence.layout$levels)
        omitted.text <- if (length(omitted) <= 10L) {
          paste0("; labels: ", paste(omitted, collapse = ", "))
        } else {
          ""
        }
        message(
          "Genomic layout: ", full.layout$n,
          " sequences are present in the active data; ", sequence.layout$n,
          " formed scan windows and ", length(omitted),
          " did not", omitted.text, ". Figures use natural sequence order."
        )
      } else {
        message(
          "Genomic layout: ", sequence.layout$n, " ", sequence.layout$type,
          "; figures use natural sequence order."
        )
      }
    }
    midpoint <- (window.table$start + window.table$end) / 2
    plot.data <- transform(window.table, midpoint_mb = midpoint / 1e6)
    score.plot.data <- dplyr::filter(
      plot.data,
      is.finite(.data$midpoint_mb), is.finite(.data$robust_score)
    )
    plots$window_scores <- ggplot2::ggplot(
      score.plot.data,
      ggplot2::aes(
        x = midpoint_mb, y = robust_score, colour = candidate_window
      )
    ) +
      ggplot2::geom_hline(yintercept = threshold, linetype = 2L) +
      ggplot2::geom_line(ggplot2::aes(group = chromosome), colour = "grey70") +
      ggplot2::geom_point(size = 1.7) +
      ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
      ggplot2::scale_colour_manual(values = c("FALSE" = "grey35", "TRUE" = "#B2182B")) +
      ggplot2::labs(
        x = "Genomic position (Mb)",
        y = "Robust local-structure score",
        colour = "Candidate window"
      ) +
      ggplot2::theme_bw()

    if (all(c("MDS1", "MDS2") %in% names(window.table))) {
      mds.plot.data <- dplyr::filter(
        window.table,
        is.finite(.data$MDS1), is.finite(.data$MDS2)
      )
      plots$window_mds <- ggplot2::ggplot(
        mds.plot.data,
        ggplot2::aes(
          x = MDS1, y = MDS2, colour = chromosome,
          shape = candidate_window
        )
      ) +
        ggplot2::geom_point(size = 2, alpha = 0.85) +
        ggplot2::labs(
          colour = "Chromosome / LG",
          shape = "Candidate window"
        ) +
        ggplot2::theme_bw()
    }

    ld.plot.data <- dplyr::filter(
      plot.data,
      is.finite(.data$midpoint_mb), is.finite(.data$mean_ld_r2)
    )
    plots$window_ld <- ggplot2::ggplot(
      ld.plot.data,
      ggplot2::aes(x = midpoint_mb, y = mean_ld_r2)
    ) +
      ggplot2::geom_line(colour = "grey50") +
      ggplot2::geom_point(
        ggplot2::aes(colour = candidate_window), size = 1.5
      ) +
      ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
      ggplot2::scale_colour_manual(values = c("FALSE" = "grey35", "TRUE" = "#B2182B")) +
      ggplot2::labs(
        x = "Genomic position (Mb)", y = expression("Mean pairwise " * r^2),
        colour = "Candidate window"
      ) +
      ggplot2::theme_bw()

    call.rate.plot.data <- dplyr::filter(
      plot.data,
      is.finite(.data$midpoint_mb), is.finite(.data$mean_call_rate)
    )
    plots$window_call_rate <- ggplot2::ggplot(
      call.rate.plot.data,
      ggplot2::aes(x = midpoint_mb, y = mean_call_rate)
    ) +
      ggplot2::geom_line(colour = "grey50") +
      ggplot2::geom_point(
        ggplot2::aes(colour = candidate_window), size = 1.5
      ) +
      ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
      ggplot2::scale_colour_manual(values = c("FALSE" = "grey35", "TRUE" = "#B2182B")) +
      ggplot2::labs(
        x = "Genomic position (Mb)", y = "Mean SNP call rate",
        colour = "Candidate window"
      ) +
      ggplot2::theme_bw()

    if (nrow(sensitivity) > 0L) {
      sensitivity.plot <- transform(
        sensitivity,
        midpoint_mb = (start + end) / 2e6,
        window_snps = factor(window_snps),
        chromosome = factor(chromosome, levels = sequence.layout$levels)
      )
      plots$window_sensitivity <- ggplot2::ggplot(
        sensitivity.plot,
        ggplot2::aes(
          x = midpoint_mb, y = pc1_variance, colour = window_snps,
          group = window_snps
        )
      ) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
        ggplot2::labs(
          x = "Genomic position (Mb)", y = "PC1 variance proportion",
          colour = "SNPs per window"
        ) +
        ggplot2::theme_bw()
    }

    # Candidate plots are named lists so users can modify the ggplot objects
    # without rereading any output file.
    candidate.plots <- purrr::map(diagnostics, function(d) {
      out <- list()
      group.title <- if (isTRUE(d$three_cluster_evidence)) {
        "Putative arrangement genotype"
      } else {
        "Algorithmic group"
      }
      pca.scores <- dplyr::filter(
        d$scores,
        is.finite(.data$PC1), is.finite(.data$PC2),
        !is.na(.data$arrangement)
      )
      out[[paste0(d$candidate_id, "_pca")]] <- ggplot2::ggplot(
        pca.scores,
        ggplot2::aes(x = PC1, y = PC2, colour = arrangement)
      ) +
        ggplot2::geom_point(size = 2, alpha = 0.85) +
        ggplot2::labs(
          title = paste0(
            d$candidate_id, " | ", d$chromosome, ":",
            format(d$start, big.mark = ","), "-",
            format(d$end, big.mark = ",")
          ),
          colour = group.title
        ) +
        ggplot2::theme_bw()

      if (nrow(chromosome.pca$scores) > 0L) {
        chromosome.scores <- chromosome.pca$scores |>
          dplyr::left_join(
            d$scores |>
              dplyr::select(.data$individual, .data$arrangement),
            by = "individual"
          )
        chromosome.scores$arrangement <- as.character(
          chromosome.scores$arrangement
        )
        chromosome.scores$arrangement[
          is.na(chromosome.scores$arrangement)
        ] <- "Unassigned"
        chromosome.levels <- levels(chromosome.pca$scores$chromosome)
        chromosome.scores$chromosome <- factor(
          chromosome.scores$chromosome,
          levels = chromosome.levels
        )
        chromosome.scores <- dplyr::filter(
          chromosome.scores,
          is.finite(.data$PC1), is.finite(.data$PC2),
          !is.na(.data$chromosome)
        )
        out[[paste0(d$candidate_id, "_chromosome_pca")]] <-
          ggplot2::ggplot(
            chromosome.scores,
            ggplot2::aes(
              x = .data$PC1, y = .data$PC2,
              colour = .data$arrangement
            )
          ) +
          ggplot2::geom_point(size = 0.9, alpha = 0.6) +
          ggplot2::facet_wrap(~ chromosome, scales = "free") +
          ggplot2::labs(
            title = paste0(
              d$candidate_id,
              " arrangement groups across chromosomes / linkage groups"
            ),
            subtitle = paste0(
              "Colours are inferred only from the candidate interval on ",
              d$chromosome,
              "; each panel is an independent chromosome-wide PCA"
            ),
            x = "Chromosome-specific PC1",
            y = "Chromosome-specific PC2",
            colour = group.title
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            legend.position = "bottom",
            panel.grid.minor = ggplot2::element_blank()
          )
      }

      heterozygosity.scores <- dplyr::filter(
        d$scores,
        is.finite(.data$heterozygosity), !is.na(.data$arrangement)
      )
      out[[paste0(d$candidate_id, "_heterozygosity")]] <-
        ggplot2::ggplot(
          heterozygosity.scores,
          ggplot2::aes(
            x = arrangement, y = heterozygosity, colour = arrangement
          )
        ) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
          ggplot2::labs(
          x = group.title,
          y = "Observed heterozygosity"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")

      loading.scores <- dplyr::filter(
        d$pc1_loadings,
        is.finite(.data$position), is.finite(.data$PC1_loading)
      )
      out[[paste0(d$candidate_id, "_marker_loadings")]] <-
        ggplot2::ggplot(
          loading.scores,
          ggplot2::aes(x = .data$position / 1e6, y = .data$PC1_loading)
        ) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
        ggplot2::geom_point(size = 1.4, alpha = 0.75) +
        ggplot2::labs(
          title = paste0(d$candidate_id, " regional PC1 marker loadings"),
          x = "Genomic position (Mb)", y = "PC1 loading"
        ) +
        ggplot2::theme_bw()

      if (any(is.finite(d$scores$mean_depth))) {
        technical.scores <- dplyr::filter(
          d$scores,
          is.finite(.data$mean_depth), is.finite(.data$call_rate),
          !is.na(.data$arrangement)
        )
        out[[paste0(d$candidate_id, "_coverage_call_rate")]] <-
          ggplot2::ggplot(
            technical.scores,
            ggplot2::aes(
              x = .data$mean_depth, y = .data$call_rate,
              colour = .data$arrangement
            )
          ) +
          ggplot2::geom_point(alpha = 0.65, size = 1.5) +
          ggplot2::labs(
            title = paste0(d$candidate_id, " technical diagnostic"),
            x = "Mean read depth", y = "Regional call rate",
            colour = group.title
          ) +
          ggplot2::theme_bw()
      }

      if (length(d$ld_position) > 1L) {
        ld.all <- d$ld_matrices$all
        ld.common <- d$ld_matrices$common_homokaryotype
        grid <- expand.grid(
          row = seq_along(d$ld_position),
          column = seq_along(d$ld_position)
        )
        # With row mapped to the upward y-axis, row >= column is the visually
        # upper triangle. Keep the displayed triangle consistent with the
        # subtitle and use the lower triangle for the common homokaryotype.
        use.all <- grid$row >= grid$column
        grid$r2 <- ifelse(
          use.all,
          ld.all[cbind(grid$row, grid$column)],
          ld.common[cbind(grid$row, grid$column)]
        )
        grid$group <- ifelse(use.all, "All individuals", "Common homokaryotype")
        axis.breaks <- unique(round(seq(
          1L, length(d$ld_position), length.out = min(6L, length(d$ld_position))
        )))
        axis.labels <- format(
          d$ld_position[axis.breaks] / 1e6,
          digits = 4L,
          trim = TRUE
        )
        out[[paste0(d$candidate_id, "_ld")]] <- ggplot2::ggplot(
          grid,
          ggplot2::aes(x = column, y = row, fill = r2)
        ) +
          ggplot2::geom_raster() +
          ggplot2::coord_equal() +
          ggplot2::scale_x_continuous(
            breaks = axis.breaks,
            labels = axis.labels,
            expand = c(0, 0)
          ) +
          ggplot2::scale_y_continuous(
            breaks = axis.breaks,
            labels = axis.labels,
            expand = c(0, 0)
          ) +
          ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") +
          ggplot2::labs(
            x = "Genomic position (Mb)", y = "Genomic position (Mb)",
            fill = expression(r^2),
            subtitle = paste0(
              "Upper triangle: all individuals; lower triangle: cluster ",
              d$common_homokaryotype_cluster,
              ". White cells: undefined r2, often a fixed marker"
            )
          ) +
          ggplot2::theme_bw()
      }
      out
    })
    plots <- c(plots, unlist(candidate.plots, recursive = FALSE))

    purrr::walk(names(plots), function(name) {
      purrr::walk(plot.formats, function(format) {
        path <- file.path(path.folder, paste0(name, ".", format))
        saved <- tryCatch(
          {
            ggplot2::ggsave(
              filename = path,
              plot = plots[[name]],
              width = if (grepl("_chromosome_pca$", name)) {
                13
              } else if (grepl("_pca$|_heterozygosity$", name)) {
                6
              } else {
                10
              },
              height = if (grepl("_chromosome_pca$", name)) {
                10
              } else if (grepl("_pca$|_heterozygosity$", name)) {
                5
              } else {
                7
              },
              dpi = 300
            )
            TRUE
          },
          error = function(error) {
            warning(
              "Could not write plot `", basename(path), "`: ",
              conditionMessage(error),
              call. = FALSE
            )
            FALSE
          }
        )
        if (saved) files <<- c(files, path)
      })
    })
  }
  if (verbose && length(candidate.summaries)) {
    message(
      "\nCandidate interpretation\n",
      paste(unname(candidate.summaries), collapse = "\n\n")
    )
  }
  if (verbose) {
    message("\nInversion results written to: ", basename(path.folder))
  }
  list(files = normalizePath(files, mustWork = FALSE), plots = plots)
}

# =============================================================================
# User-facing print method
# =============================================================================

#' @export
print.detect_inversions <- function(x, ...) {
  cat("Candidate inversion-associated region scan\n")
  cat("  Windows analysed:", sum(x$windows$valid), "of", nrow(x$windows), "\n")
  cat("  Candidate regions:", nrow(x$candidates), "\n")
  cat("  Results folder:", basename(x$path.folder), "\n")
  if (length(x$candidate.summaries)) {
    cat(
      "\n",
      paste(unname(x$candidate.summaries), collapse = "\n\n"),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}
