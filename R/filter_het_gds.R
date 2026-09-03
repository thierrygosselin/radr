#' Filter heterozygosity outliers
#'
#' Filter active GDS individuals and markers using observed heterozygosity.
#' Individual statistics are calculated first. When individuals are removed,
#' marker statistics are recalculated from the remaining active samples before
#' marker filtering. This avoids applying stale marker heterozygosity estimates.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional metadata containing `INDIVIDUALS` and `group.column`.
#' @param group.column Metadata column defining strata.
#' @param individual.range Permitted individual observed-heterozygosity range.
#'   Use `NULL` to retain all individuals.
#' @param marker.range Permitted marker observed-heterozygosity range within
#'   strata. Use `NULL` to retain all markers.
#' @param strata.threshold Number or proportion of eligible strata that must
#'   violate `marker.range` before marker removal.
#' @param min.individual.call.rate Minimum individual call rate.
#' @param min.marker.call.rate Minimum marker call rate within a stratum.
#' @param chunk.size Number of variants read at a time.
#' @param verbose Display progress and summary messages.
#' @param ... Common argument `path.folder`.
#' @return The filtered open GDS connection with complete filter audit files.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' genome <- radr::filter_het(
#'   genome, strata = "sample_metadata.tsv",
#'   individual.range = c(0.05, 0.35), marker.range = c(0, 0.8)
#' )
#' }
filter_het <- function(
    data, strata = NULL, group.column = "STRATA",
    individual.range = NULL, marker.range = NULL,
    strata.threshold = 1L, min.individual.call.rate = 0,
    min.marker.call.rate = 0.8, chunk.size = 2000L,
    verbose = TRUE, ...
) {
  force(data)
  .het_check_range(individual.range, "individual.range")
  .het_check_range(marker.range, "marker.range")
  .paralog_check_probability(min.individual.call.rate, "min.individual.call.rate")
  .paralog_check_probability(min.marker.call.rate, "min.marker.call.rate")
  if (is.null(individual.range) && is.null(marker.range)) return(data)
  if (!is.numeric(strata.threshold) || length(strata.threshold) != 1L ||
      is.na(strata.threshold) || strata.threshold <= 0) {
    rlang::abort("`strata.threshold` must be positive.")
  }
  .start <- tgbase::startup(package = "radr", f.name = "filter_het", verbose = verbose)
  on.exit(tgbase::teardown(.start), add = TRUE)
  dots <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  unknown <- setdiff(names(dots), "path.folder")
  if (length(unknown)) rlang::abort(paste0("Unknown argument(s): ", paste(unknown, collapse = ", "), "."))
  path.folder <- radr_folder(
    rad.folder = paste0("filter_het_", .start$file.date),
    path.folder = dots$path.folder %||% getwd(), prefix.int = TRUE
  )
  opened <- .filter_gds_open(data)
  gds <- opened$gds
  genometranslator::sync_gds(gds, verbose = FALSE)
  if (!genometranslator::detect_biallelic_markers(gds)) {
    rlang::abort("Heterozygosity filtering requires biallelic markers.")
  }
  selection <- SeqArray::seqGetFilter(gds)
  variant.id <- SeqArray::seqGetData(gds, "variant.id")
  sample.id <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  meta <- .paralog_read_metadata(strata, gds, sample.id, group.column, TRUE)
  missing.meta <- setdiff(sample.id, meta$INDIVIDUALS)
  if (length(missing.meta)) {
    rlang::abort(paste0(
      "Metadata are missing ", length(missing.meta), " active sample(s)."
    ))
  }
  meta <- meta[match(sample.id, meta$INDIVIDUALS), , drop = FALSE]
  chunks <- split(seq_along(variant.id), ceiling(seq_along(variant.id) / chunk.size))
  individual.parts <- lapply(chunks, function(index) {
    dosage <- .sex_get_dosage(gds, variant.id[index], sample.id)
    tibble::tibble(
      INDIVIDUALS = sample.id, CALLED = rowSums(!is.na(dosage)),
      HETEROZYGOUS = rowSums(dosage == 1, na.rm = TRUE)
    )
  })
  individual <- dplyr::bind_rows(individual.parts) |>
    dplyr::group_by(.data$INDIVIDUALS) |>
    dplyr::summarise(
      NUMBER_CALLED = sum(.data$CALLED),
      NUMBER_HETEROZYGOUS = sum(.data$HETEROZYGOUS),
      CALL_RATE = .data$NUMBER_CALLED / length(variant.id),
      OBSERVED_HETEROZYGOSITY = .data$NUMBER_HETEROZYGOUS /
        .data$NUMBER_CALLED, .groups = "drop"
    ) |>
    dplyr::left_join(
      meta[, c("INDIVIDUALS", group.column), drop = FALSE],
      by = "INDIVIDUALS"
    )
  remove.samples <- character()
  if (!is.null(individual.range)) {
    remove.samples <- individual |>
      dplyr::filter(
        .data$CALL_RATE < min.individual.call.rate |
        .data$OBSERVED_HETEROZYGOSITY < individual.range[[1L]] |
        .data$OBSERVED_HETEROZYGOSITY > individual.range[[2L]]
      ) |>
      dplyr::pull(.data$INDIVIDUALS)
  }
  retained.samples <- setdiff(sample.id, remove.samples)
  if (!length(retained.samples)) rlang::abort("Heterozygosity filtering would remove all samples.")
  SeqArray::seqSetFilter(gds, sample.id = retained.samples, verbose = FALSE)
  marker.stats <- summarise_genomic_data(
    gds, strata = meta[meta$INDIVIDUALS %in% retained.samples, , drop = FALSE],
    group.column = group.column, by.strata = TRUE, chunk.size = chunk.size,
    write.files = FALSE, verbose = FALSE, internal = TRUE,
    path.folder = path.folder
  )$stratum.statistics |>
    dplyr::filter(.data$GROUP != "OVERALL") |>
    dplyr::mutate(
      ELIGIBLE = .data$CALL_RATE >= min.marker.call.rate,
      VIOLATION = if (is.null(marker.range)) FALSE else
        .data$ELIGIBLE &
          (.data$OBSERVED_HETEROZYGOSITY < marker.range[[1L]] |
           .data$OBSERVED_HETEROZYGOSITY > marker.range[[2L]])
    )
  evidence <- marker.stats |>
    dplyr::group_by(.data$VARIANT_ID) |>
    dplyr::summarise(
      N_ELIGIBLE_STRATA = sum(.data$ELIGIBLE),
      N_VIOLATING_STRATA = sum(.data$VIOLATION),
      PROPORTION_VIOLATING = dplyr::if_else(
        .data$N_ELIGIBLE_STRATA > 0,
        .data$N_VIOLATING_STRATA / .data$N_ELIGIBLE_STRATA, NA_real_
      ), .groups = "drop"
    ) |>
    dplyr::mutate(REMOVE = if (is.null(marker.range)) FALSE else
      if (strata.threshold < 1) {
        .data$PROPORTION_VIOLATING >= strata.threshold
      } else .data$N_VIOLATING_STRATA >= strata.threshold)
  remove.markers <- evidence$VARIANT_ID[evidence$REMOVE]
  SeqArray::seqResetFilter(gds, verbose = FALSE)
  SeqArray::seqSetFilter(
    gds, sample.id = selection$sample.id,
    variant.id = selection$variant.id, verbose = FALSE
  )
  readr::write_tsv(individual, file.path(path.folder, "individual_heterozygosity.tsv"), na = "NA")
  readr::write_tsv(marker.stats, file.path(path.folder, "marker_stratum_heterozygosity.tsv"), na = "NA")
  readr::write_tsv(evidence, file.path(path.folder, "heterozygosity_filter_evidence.tsv"), na = "NA")
  .filter_gds_apply(
    gds, remove.markers, remove.samples, "filter.het", path.folder,
    .start$file.date, "heterozygosity ranges",
    c(individual.range, marker.range, strata.threshold), verbose
  )
  if (verbose) .summary_message(
    "Heterozygosity filter removed ", length(remove.samples), " samples and ",
    length(remove.markers), " markers."
  )
  gds
}

.het_check_range <- function(x, name) {
  if (is.null(x)) return(invisible(TRUE))
  if (!is.numeric(x) || length(x) != 2L || anyNA(x) ||
      any(x < 0 | x > 1) || x[[1L]] > x[[2L]]) {
    rlang::abort(paste0("`", name, "` must be NULL or an ordered [0, 1] pair."))
  }
  invisible(TRUE)
}
