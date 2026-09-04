#' Diagnose candidate copy-number variation
#'
#' Screen active GDS loci for relative depth variation. This exploratory
#' diagnostic does not call absolute copy numbers or structural breakpoints.
#'
#' @param data GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional sample metadata data frame or TSV with `INDIVIDUALS`.
#' @param group.column Population column, default `STRATA`.
#' @param batch.columns Optional categorical technical metadata columns.
#' @param risk.annotations Optional locus-level data frame containing
#' `locus.columns` and logical `PARALOG_FLAG`, `LOW_MAPPABILITY`,
#' `REPEAT_OVERLAP`, or `INDEPENDENT_SUPPORT` columns. Use `NA` for unknown.
#' Missing rows or checks are not treated as negative findings. Flags can be
#' prepared from external tools or reviewed `detect_paralogs()` results, but
#' must refer to the same assembly and locus definitions as the GDS.
#' @param locus.columns Marker metadata columns jointly identifying a locus.
#' Defaults to `c("CHROM", "LOCUS")`; valid locus identifiers are required.
#' @param min.samples Minimum samples with usable depth per locus (default 10).
#' @param min.loci Minimum positive-depth loci for normalization (default 100).
#' @param min.ratio Relative depth ratio flagging a gain; its reciprocal flags
#' a loss (default 1.5). Thresholds are exploratory, not significance tests.
#' @param min.carriers Minimum samples exceeding either threshold (default 3).
#' @param metadata.threshold Minimum descriptive eta-squared to flag a metadata
#' association (default 0.5). This is not a p-value or a causal conclusion.
#' @param chunk.size Number of representative markers read per chunk.
#' @param write.files Write tables and a diagnostic plot.
#' @param verbose Display progress.
#' @param ... Common arguments `path.folder` and `internal`.
#'
#' @details
#' \itemize{
#' \item \strong{Locus-level evidence:} One representative variant (smallest
#' VARIANT_ID) is used per locus; depths from SNPs sharing reads are not summed.
#' \item \strong{Depth input:} FORMAT/DP is used, or retained embedded coverage
#' if DP is absent. Missing depth is never changed to zero, and missing
#' genotype calls do not erase available depth.
#' \item \strong{Normalization:} Positive depths are divided by the locus
#' geometric mean. Each sample's median ratio estimates its library factor;
#' scaled depths are then divided by their locus-specific cohort median.
#' \item \strong{Baseline assumptions:} Most baseline loci must have stable
#' copy number. The cohort median is not a known diploid reference. Widespread
#' CNV, global ploidy changes, capture bias and restriction-site dropout can
#' defeat normalization.
#' \item \strong{Allele-balance evidence:} Median absolute deviation from 0.5
#' is calculated among called heterozygotes with positive allele depth at the
#' representative SNP, separately within each stratum. Called counts,
#' heterozygosity and expected heterozygosity accompany this evidence.
#' The pooled summary is descriptive only, not a copy-number call.
#' \item \strong{Safeguards:} Separate evidence and risk tables retain supplied
#' paralogy, mappability, repeat and independent-validation annotations. Missing
#' assessments are labelled NOT_ASSESSED, not passed. These annotations do not
#' automatically remove candidates. No ngsParalog likelihood test, read-mapping
#' audit or normalization-sensitivity benchmark is run by this function.
#' \item \strong{Metadata audit:} Associations use log2(1 + relative depth)
#' and categorical eta-squared. Group-by-batch counts are returned. Small
#' categories and confounded designs require scrutiny; no automatic batch
#' correction is performed.
#' \item \strong{Non-destructive workflow:} No loci or individuals are removed
#' from the GDS, and active filters are restored. Run before permanently
#' discarding paralogous loci.
#' \item \strong{Scope and memory:} Memory scales with samples times loci,
#' not all SNPs. A SNP-only GDS cannot provide a genome-wide WGS CNV scan.
#' \item \strong{Methodological attribution:} McKinney et al. (2017) provide
#' the conceptual basis for using heterozygosity and allelic read ratios to
#' investigate paralogy. This function uses a simpler supporting imbalance
#' summary, not their full HDplot statistic or classification. Its primary
#' candidate flag is based on relative depth. It does not implement the
#' allele-dosage likelihood model of McKinney et al. (2018), Dorant's TMM
#' workflow, or CODEX2.
#' \item \strong{Alternative methods and scientific context:} Comparing
#' alternative approaches is important for testing assumptions and robustness.
#' Dorant et al. demonstrate normalized-depth analysis in marine Rapture data;
#' CODEX2 addresses coverage bias; CNVpytor and Genome STRiP address WGS CNV
#' discovery/genotyping; Dallaire et al. motivate deviant-SNP safeguards.
#' Schrider and Hahn review the evolutionary context. These references identify
#' useful comparisons, not methods implemented or benchmarks completed here.
#' }
#'
#' @return A list with locus statistics, sample normalization factors,
#' relative-depth observations, representative markers, metadata associations,
#' group-by-batch counts, parameters, plot and filter-restoration status.
#' @references
#' McKinney, G. J., Waples, R. K., Seeb, L. W., and Seeb, J. E. (2017).
#' Paralogs are revealed by proportion of heterozygotes and deviations in read
#' ratios in genotyping-by-sequencing data from natural populations.
#' Molecular Ecology Resources, 17, 656-669. \doi{10.1111/1755-0998.12613}.
#'
#' McKinney, G. J., Waples, R. K., Pascal, C. E., Seeb, L. W., and
#' Seeb, J. E. (2018). Resolving allele dosage in duplicated loci using
#' genotyping-by-sequencing data: A path forward for population genetic analysis.
#' Molecular Ecology Resources, 18, 570-579. \doi{10.1111/1755-0998.12763}.
#'
#' Dorant et al. (2020). Copy number variants outperform SNPs to reveal
#' genotype-temperature association in a marine species. Molecular Ecology.
#' \doi{10.1111/mec.15565}.
#'
#' Jiang et al. (2018). CODEX2: full-spectrum copy number variation detection
#' by high-throughput DNA sequencing. Genome Biology, 19, 202.
#' \doi{10.1186/s13059-018-1578-y}.
#'
#' Suvakov et al. (2021). CNVpytor: a tool for copy number variation detection
#' and analysis from read depth and allele imbalance in whole-genome sequencing.
#' GigaScience, 10, giab074. \doi{10.1093/gigascience/giab074}.
#'
#' Handsaker et al. (2015). Large multiallelic copy number variations in humans.
#' Nature Genetics, 47, 296-303. \doi{10.1038/ng.3200}.
#'
#' Dallaire et al. (2023). Widespread Deviant Patterns of Heterozygosity in
#' Whole-Genome Sequencing Due to Autopolyploidy, Repeated Elements, and
#' Duplication. Genome Biology and Evolution, 15, evad229.
#' \doi{10.1093/gbe/evad229}.
#'
#' Schrider, D. R., and Hahn, M. W. (2010). Gene copy-number polymorphism
#' in nature. Proceedings of the Royal Society B, 277, 3213-3221.
#' \doi{10.1098/rspb.2010.1180}.
#' @seealso [detect_paralogs()]
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' cnv <- detect_cnv("study.gds", strata = "samples.tsv",
#'                   batch.columns = c("PLATE", "PROJECT"))
#' }
detect_cnv <- function(
    data, strata = NULL, group.column = "STRATA", batch.columns = character(),
    risk.annotations = NULL,
    locus.columns = c("CHROM", "LOCUS"), min.samples = 10L, min.loci = 100L,
    min.ratio = 1.5, min.carriers = 3L, metadata.threshold = 0.5,
    chunk.size = 2000L, write.files = TRUE, verbose = TRUE, ...
) {
  min.samples <- .paralog_check_count(min.samples, "min.samples", 2L)
  min.loci <- .paralog_check_count(min.loci, "min.loci", 2L)
  min.carriers <- .paralog_check_count(min.carriers, "min.carriers", 1L)
  chunk.size <- .paralog_check_count(chunk.size, "chunk.size", 1L)
  if (length(min.ratio) != 1L || !is.finite(min.ratio) || min.ratio <= 1)
    rlang::abort("`min.ratio` must be finite and greater than one.")
  if (length(metadata.threshold) != 1L || !is.finite(metadata.threshold) ||
      metadata.threshold < 0 || metadata.threshold > 1)
    rlang::abort("`metadata.threshold` must be between zero and one.")
  context <- .analysis_gds_start(data, "detect_cnv", write.files, verbose, ...)
  on.exit(.analysis_gds_finish(context), add = TRUE)
  gds <- context$gds
  if (!genometranslator::detect_biallelic_markers(gds))
    rlang::abort("This candidate diagnostic currently requires biallelic markers.")
  meta <- .analysis_metadata(gds, strata, group.column, TRUE)
  if (!all(batch.columns %in% names(meta$metadata)))
    rlang::abort("Requested batch columns are absent from sample metadata.")
  for (column in batch.columns) {
    x <- as.character(meta$metadata[[column]])
    if (anyNA(x) || any(!nzchar(trimws(x))))
      rlang::abort("Batch metadata must not contain missing or empty values.")
  }
  markers <- .sex_marker_metadata(gds, SeqArray::seqGetData(gds, "variant.id"))
  representatives <- .cnv_representatives(markers, locus.columns)
  risk <- .cnv_risk_annotations(representatives, locus.columns, risk.annotations)
  ids <- representatives$VARIANT_ID
  if (length(ids) < min.loci) rlang::abort("Too few loci for normalization.")
  depth <- matrix(NA_real_, length(meta$sample.id), length(ids))
  imbalance <- rep(NA_real_, length(ids))
  stratum.statistics <- list()
  sources <- character()
  for (index in split(seq_along(ids), ceiling(seq_along(ids) / chunk.size))) {
    SeqArray::seqSetFilter(gds, variant.id = ids[index], verbose = FALSE)
    dp <- tryCatch(SeqArray::seqGetData(gds, "annotation/format/DP"),
                   error = function(e) NULL)
    source <- "FORMAT/DP"
    if (is.list(dp)) dp <- dp$data
    if (is.null(dp)) {
      embedded <- .inversion_get_embedded_coverage(gds, ids[index], meta$sample.id)
      dp <- embedded$depth
      source <- "embedded coverage"
    }
    if (is.null(dp)) rlang::abort("Retained read depth is required for CNV diagnostics.")
    dp <- as.matrix(dp)
    if (!identical(dim(dp), c(length(meta$sample.id), length(index))))
      rlang::abort("Depth dimensions do not match active samples and variants.")
    dp[!is.finite(dp) | dp < 0] <- NA_real_
    depth[, index] <- dp
    sources <- union(sources, source)
    ad <- .paralog_get_allele_depth(gds, ids[index], meta$sample.id)
    gt <- .sex_get_dosage(gds, ids[index], meta$sample.id)
    ratio <- matrix(NA_real_, nrow(gt), ncol(gt))
    if (!is.null(ad$ref)) {
      total <- ad$ref + ad$alt
      ratio <- abs(ad$alt / total - 0.5)
      ratio[is.na(gt) | gt != 1 | !is.finite(total) | total <= 0 |
              ad$ref < 0 | ad$alt < 0] <- NA_real_
      imbalance[index] <- apply(ratio, 2L, .cnv_median)
    }
    for (group in unique(meta$groups)) {
      rows <- meta$groups == group
      x <- gt[rows, , drop = FALSE]
      n <- colSums(!is.na(x))
      frequency <- colSums(x, na.rm = TRUE) / (2 * n)
      h <- colSums(x == 1, na.rm = TRUE) / n
      expected <- 2 * frequency * (1 - frequency)
      h[!is.finite(h)] <- NA_real_
      expected[!is.finite(expected)] <- NA_real_
      stratum.statistics[[length(stratum.statistics) + 1L]] <- tibble::tibble(
        LOCUS_INDEX = index, VARIANT_ID = ids[index], GROUP = group,
        NUMBER_CALLED = n, OBSERVED_HETEROZYGOSITY = h,
        EXPECTED_HETEROZYGOSITY = expected,
        NUMBER_HET_WITH_DEPTH = colSums(is.finite(ratio[rows, , drop = FALSE])),
        HET_ALLELE_IMBALANCE = apply(ratio[rows, , drop = FALSE], 2L, .cnv_median)
      )
    }
    if (verbose) .summary_message("CNV loci read: ", max(index), " / ", length(ids))
  }
  normalized <- .cnv_normalize(depth, min.samples, min.loci)
  relative <- normalized$relative
  statistics <- representatives
  statistics$N_DEPTH <- colSums(is.finite(relative))
  statistics$N_LOW <- colSums(relative <= 1 / min.ratio, na.rm = TRUE)
  statistics$N_HIGH <- colSums(relative >= min.ratio, na.rm = TRUE)
  statistics$HET_ALLELE_IMBALANCE <- imbalance
  statistics$STATUS <- ifelse(statistics$N_DEPTH < min.samples, "UNINFORMATIVE",
    ifelse(pmax(statistics$N_LOW, statistics$N_HIGH) >= min.carriers,
           "CANDIDATE_DEPTH_VARIATION", "WITHIN_THRESHOLD"))
  associations <- list(); batches <- list()
  for (column in unique(c(group.column, batch.columns))) {
    category <- as.character(meta$metadata[[column]])
    associations[[column]] <- tibble::tibble(
      LOCUS_INDEX = seq_along(ids), VARIABLE = column,
      ETA_SQUARED = apply(log2(1 + relative), 2L, .cnv_eta, category = category)
    )
    if (column %in% batch.columns) {
      batches[[column]] <- tibble::as_tibble(as.data.frame(table(
        GROUP = meta$groups, BATCH = category), stringsAsFactors = FALSE))
      batches[[column]]$VARIABLE <- column
    }
  }
  associations <- dplyr::bind_rows(associations)
  associations$FLAG <- is.finite(associations$ETA_SQUARED) &
    associations$ETA_SQUARED >= metadata.threshold
  risk$BATCH_ASSOCIATION <- vapply(seq_along(ids), function(i) {
    a <- associations[associations$LOCUS_INDEX == i &
                        associations$VARIABLE %in% batch.columns, ]
    if (!nrow(a) || !any(is.finite(a$ETA_SQUARED))) return("NOT_ASSESSED")
    if (any(a$FLAG)) "FLAGGED" else "NOT_FLAGGED"
  }, character(1))
  risk$NORMALIZATION_SENSITIVITY <- "NOT_ASSESSED"
  risk$READ_MAPPING_AUDIT <- "NOT_ASSESSED"
  stratum.statistics <- dplyr::bind_rows(stratum.statistics)
  evidence <- statistics[, c("LOCUS_INDEX", "VARIANT_ID", "STATUS", "N_DEPTH",
                             "N_LOW", "N_HIGH"), drop = FALSE]
  names(evidence)[names(evidence) == "STATUS"] <- "DEPTH_EVIDENCE"
  evidence$ALLELE_BALANCE <- vapply(seq_along(ids), function(i) {
    x <- stratum.statistics$HET_ALLELE_IMBALANCE[
      stratum.statistics$LOCUS_INDEX == i]
    if (any(is.finite(x))) "AVAILABLE_BY_STRATUM" else "NOT_ASSESSED"
  }, character(1))
  evidence$ABSOLUTE_COPY_NUMBER <- "NOT_ASSESSED"
  observations <- tibble::tibble(
    INDIVIDUALS = rep(meta$sample.id, length(ids)),
    LOCUS_INDEX = rep(seq_along(ids), each = length(meta$sample.id)),
    RAW_DEPTH = as.vector(depth),
    RELATIVE_DEPTH = as.vector(relative)
  )
  samples <- tibble::tibble(INDIVIDUALS = meta$sample.id,
                            LIBRARY_FACTOR = normalized$factors)
  samples$STATUS <- ifelse(is.finite(samples$LIBRARY_FACTOR),
                            "NORMALIZED", "INSUFFICIENT_BASELINE_DEPTH")
  plot <- ggplot2::ggplot(statistics, ggplot2::aes(
    x = .data$LOCUS_INDEX, y = pmax(.data$N_LOW, .data$N_HIGH), colour = .data$STATUS
  )) + ggplot2::geom_point() + ggplot2::theme_bw() +
    ggplot2::labs(x = "Locus index (not genomic position)",
                  y = "Samples exceeding a relative-depth threshold")
  parameters <- list(min.samples = min.samples, min.loci = min.loci,
    min.ratio = min.ratio, min.carriers = min.carriers,
    metadata.threshold = metadata.threshold, locus.columns = locus.columns,
    batch.columns = batch.columns, group.column = group.column,
    depth.source = sources, normalization = "positive-depth median ratios")
  result <- list(statistics = statistics, samples = samples,
    evidence = evidence, technical.risks = risk,
    stratum.statistics = stratum.statistics,
    sample.metadata = tibble::as_tibble(meta$metadata),
    relative.depth = observations, metadata.associations = associations,
    group.batch.counts = dplyr::bind_rows(batches), parameters = parameters,
    diagnostic.plot = plot, path.folder = context$path.folder)
  if (write.files) {
    for (name in c("statistics", "evidence", "technical.risks", "stratum.statistics",
                   "samples", "sample.metadata", "relative.depth",
                   "metadata.associations", "group.batch.counts")) {
      readr::write_tsv(result[[name]], file.path(context$path.folder,
                                              paste0("cnv_", name, ".tsv")))
    }
    saveRDS(parameters, file.path(context$path.folder, "cnv_parameters.rds"))
    ggplot2::ggsave(file.path(context$path.folder, "cnv_diagnostic.png"),
                    plot, width = 10, height = 6)
  }
  if (verbose) {
    .summary_message("CNV candidates are relative-depth signals, not copy-number calls.")
    if (any(associations$FLAG))
      .summary_message("Strong descriptive metadata associations detected; inspect before interpretation.")
    .summary_message("Review metadata associations and group-by-batch counts.")
  }
  result$active.selection.restored <- .analysis_gds_restore(context)
  class(result) <- c("detect_cnv", "list")
  result
}

.cnv_risk_annotations <- function(representatives, columns, annotations) {
  checks <- c("PARALOG_FLAG", "LOW_MAPPABILITY", "REPEAT_OVERLAP",
              "INDEPENDENT_SUPPORT")
  out <- representatives[, unique(c("LOCUS_INDEX", columns)), drop = FALSE]
  if (!is.null(annotations)) {
    if (!is.data.frame(annotations) || !all(columns %in% names(annotations)))
      rlang::abort("Risk annotations must contain the locus identifier columns.")
    if (anyDuplicated(annotations[, columns, drop = FALSE]))
      rlang::abort("Risk annotations must have exactly one row per locus.")
    for (column in columns) {
      if (anyNA(annotations[[column]]) ||
          any(!nzchar(trimws(as.character(annotations[[column]])))))
        rlang::abort("Risk annotation locus identifiers must not be missing.")
    }
    present <- intersect(checks, names(annotations))
    if (!length(present)) rlang::abort("No supported risk annotation columns supplied.")
    for (check in present) {
      if (!is.logical(annotations[[check]]))
        rlang::abort("Risk annotation flags must be logical TRUE, FALSE or NA.")
    }
    out <- dplyr::left_join(out, annotations[, c(columns, present), drop = FALSE],
                             by = columns)
  }
  for (check in checks) {
    x <- out[[check]]
    if (is.null(x)) x <- rep(NA, nrow(out))
    out[[check]] <- ifelse(is.na(x), "NOT_ASSESSED", ifelse(x,
      if (check == "INDEPENDENT_SUPPORT") "SUPPORTED" else "FLAGGED",
      if (check == "INDEPENDENT_SUPPORT") "NOT_SUPPORTED" else "NOT_FLAGGED"))
  }
  out
}

.cnv_representatives <- function(markers, columns) {
  if (!length(columns) || !all(columns %in% names(markers)))
    rlang::abort("Valid locus metadata columns are required; do not infer them from IDs.")
  for (column in columns) {
    if (anyNA(markers[[column]]) || any(!nzchar(trimws(as.character(markers[[column]])))))
      rlang::abort("Locus metadata contains missing or empty identifiers.")
  }
  markers <- markers[order(markers$VARIANT_ID), , drop = FALSE]
  markers <- markers[!duplicated(markers[, columns, drop = FALSE]), , drop = FALSE]
  markers$LOCUS_INDEX <- seq_len(nrow(markers))
  markers
}

.cnv_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

.cnv_normalize <- function(depth, min.samples, min.loci) {
  positive <- depth; positive[!is.finite(positive) | positive <= 0] <- NA_real_
  usable <- colSums(is.finite(positive)) >= min.samples
  if (sum(usable) < min.loci) rlang::abort("Too few positive-depth baseline loci.")
  base <- exp(colMeans(log(positive[, usable, drop = FALSE]), na.rm = TRUE))
  ratios <- sweep(positive[, usable, drop = FALSE], 2L, base, "/")
  factors <- apply(ratios, 1L, .cnv_median)
  factors[rowSums(is.finite(ratios)) < min.loci] <- NA_real_
  if (sum(is.finite(factors)) < min.samples)
    rlang::abort("Too few samples with adequate depth for normalization.")
  scaled <- sweep(depth, 1L, factors, "/")
  centers <- apply(scaled, 2L, .cnv_median)
  centers[!is.finite(centers) | centers <= 0] <- NA_real_
  list(relative = sweep(scaled, 2L, centers, "/"), factors = factors)
}

.cnv_eta <- function(x, category) {
  good <- is.finite(x) & !is.na(category)
  x <- x[good]; category <- category[good]
  if (length(x) < 3L || length(unique(category)) < 2L) return(NA_real_)
  total <- sum((x - mean(x))^2)
  if (total <= 0) return(NA_real_)
  means <- ave(x, category, FUN = mean)
  sum((means - mean(x))^2) / total
}
