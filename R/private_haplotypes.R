#' Detect private phased haplotypes
#'
#' Construct genuine phased haplotypes across variants sharing a `LOCUS` and
#' identify haplotypes observed in exactly one stratum. Unlike the historical
#' implementation, this function does not label a single-marker diploid
#' genotype string as a haplotype.
#'
#' A sample-locus observation is usable when every variant is called and every
#' heterozygous genotype is phased. Homozygous genotypes do not require phase.
#' Loci with fewer than `min.variants` active variants are not analysed.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional metadata containing `INDIVIDUALS` and `group.column`.
#' @param group.column Metadata column defining strata.
#' @param min.variants Minimum variants required to define a haplotype.
#' @param min.haplotype.count Minimum chromosome count in the private stratum.
#' @param write.files Write result and summary tables.
#' @param verbose Display progress and summary messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#' @return A `private_haplotypes` object containing `private.haplotypes`,
#' `haplotype.statistics`, `stratum.summary`, and output information.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' private <- radr::private_haplotypes(
#'   "phased.gds", strata = "sample_metadata.tsv"
#' )
#' }
private_haplotypes <- function(
    data, strata = NULL, group.column = "STRATA", min.variants = 2L,
    min.haplotype.count = 1L, write.files = TRUE, verbose = TRUE, ...
) {
  min.variants <- .paralog_check_count(min.variants, "min.variants", 1L)
  min.haplotype.count <- .paralog_check_count(
    min.haplotype.count, "min.haplotype.count", 1L
  )
  context <- .analysis_gds_start(
    data, "private_haplotypes", write.files, verbose, ...
  )
  on.exit(.analysis_gds_finish(context), add = TRUE)
  gds <- context$gds
  meta <- .analysis_metadata(gds, strata, group.column, TRUE)
  marker <- .sex_marker_metadata(gds, SeqArray::seqGetData(gds, "variant.id"))
  if (!"LOCUS" %in% names(marker) || all(is.na(marker$LOCUS))) {
    rlang::abort(
      "Private haplotypes require a populated `LOCUS` marker field."
    )
  }
  loci <- split(marker$VARIANT_ID, as.character(marker$LOCUS))
  loci <- loci[lengths(loci) >= min.variants]
  if (!length(loci)) {
    rlang::abort("No active locus contains the required number of variants.")
  }
  observations <- dplyr::bind_rows(lapply(names(loci), function(locus) {
    .private_locus_haplotypes(
      gds, loci[[locus]], locus, meta$sample.id, meta$groups
    )
  }))
  if (!nrow(observations)) {
    rlang::abort(
      "No complete phased sample-locus haplotypes could be constructed."
    )
  }
  statistics <- observations |>
    dplyr::count(
      .data$LOCUS, .data$GROUP, .data$HAPLOTYPE,
      name = "HAPLOTYPE_COUNT"
    ) |>
    dplyr::group_by(.data$LOCUS, .data$GROUP) |>
    dplyr::mutate(
      CHROMOSOMES_CALLED = sum(.data$HAPLOTYPE_COUNT),
      HAPLOTYPE_FREQUENCY = .data$HAPLOTYPE_COUNT /
        .data$CHROMOSOMES_CALLED
    ) |>
    dplyr::ungroup()
  private <- statistics |>
    dplyr::group_by(.data$LOCUS, .data$HAPLOTYPE) |>
    dplyr::mutate(N_STRATA_OBSERVED = dplyr::n_distinct(.data$GROUP)) |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$N_STRATA_OBSERVED == 1L,
      .data$HAPLOTYPE_COUNT >= min.haplotype.count
    ) |>
    dplyr::arrange(.data$GROUP, .data$LOCUS, .data$HAPLOTYPE)
  summary <- private |>
    dplyr::count(.data$GROUP, name = "N_PRIVATE_HAPLOTYPES")
  output.files <- character()
  if (write.files) {
    output.files <- file.path(context$path.folder, c(
      "private_haplotypes.tsv", "haplotype_statistics.tsv",
      "private_haplotypes_by_stratum.tsv"
    ))
    readr::write_tsv(private, output.files[[1L]], na = "NA")
    readr::write_tsv(statistics, output.files[[2L]], na = "NA")
    readr::write_tsv(summary, output.files[[3L]], na = "NA")
  }
  if (verbose) .summary_message("Private phased haplotypes: ", nrow(private), ".")
  restored <- .analysis_gds_restore(context)
  out <- list(
    private.haplotypes = private, haplotype.statistics = statistics,
    stratum.summary = summary, path.folder = context$path.folder,
    output.files = tibble::tibble(files = output.files),
    active.selection.restored = restored
  )
  class(out) <- c("private_haplotypes", class(out))
  out
}

.private_locus_haplotypes <- function(
    gds, variant.id, locus, sample.id, groups
) {
  SeqArray::seqFilterPush(gds)
  on.exit(SeqArray::seqFilterPop(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = variant.id, verbose = FALSE)
  genotype <- SeqArray::seqGetData(gds, "genotype")
  phase <- tryCatch(
    as.matrix(SeqArray::seqGetData(gds, "phase")),
    error = function(error) matrix(0L, length(sample.id), length(variant.id))
  )
  if (nrow(phase) != length(sample.id)) phase <- t(phase)
  allele.definitions <- strsplit(
    as.character(SeqArray::seqGetData(gds, "allele")), ",", fixed = TRUE
  )
  rows <- vector("list", length(sample.id))
  for (s in seq_along(sample.id)) {
    genotype.sample <- genotype[, s, , drop = FALSE]
    dim(genotype.sample) <- c(2L, length(variant.id))
    if (any(genotype.sample < 0L)) next
    heterozygous <- genotype.sample[1L, ] != genotype.sample[2L, ]
    if (any(heterozygous & phase[s, ] != 1L)) next
    haplotypes <- vapply(1:2, function(copy) {
      paste(vapply(seq_along(variant.id), function(v) {
        index <- genotype.sample[copy, v] + 1L
        allele.definitions[[v]][index]
      }, character(1)), collapse = "|")
    }, character(1))
    rows[[s]] <- tibble::tibble(
      LOCUS = as.character(locus), INDIVIDUALS = sample.id[[s]],
      GROUP = groups[[s]], HAPLOTYPE = haplotypes
    )
  }
  dplyr::bind_rows(rows)
}
