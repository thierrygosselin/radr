#' Filter markers using FIS
#'
#' Calculate within-stratum FIS from active biallelic GDS genotypes and remove
#' variants, or all variants in their RAD/read locus, when departures exceed
#' the requested limits in enough strata.
#'
#' FIS is `1 - Ho / He` and is undefined for monomorphic marker-stratum
#' combinations. Undefined or low-call-rate combinations are excluded from the
#' evidence count rather than treated as passing or failing.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional metadata containing `INDIVIDUALS` and `group.column`.
#' @param group.column Metadata column defining strata.
#' @param unit Filter individual `"variant"` records or complete `"locus"`
#'   groups defined by marker `LOCUS`.
#' @param fis.min.threshold Minimum permitted FIS.
#' @param fis.max.threshold Maximum permitted FIS.
#' @param strata.threshold Number or proportion of eligible strata that must
#'   violate a limit before removal. Values below one are proportions.
#' @param min.call.rate Minimum marker call rate within a stratum.
#' @param chunk.size Number of variants read at a time.
#' @param verbose Display progress and summary messages.
#' @param ... Common argument `path.folder`.
#' @return The filtered open GDS connection. The GDS, filter history, and full
#' kept/removed audits are updated on disk.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' genome <- radr::filter_fis(
#'   genome, strata = "sample_metadata.tsv",
#'   fis.min.threshold = -0.5, fis.max.threshold = 0.5
#' )
#' }
filter_fis <- function(
    data, strata = NULL, group.column = "STRATA",
    unit = c("variant", "locus"), fis.min.threshold = -Inf,
    fis.max.threshold = Inf, strata.threshold = 1L,
    min.call.rate = 0.8, chunk.size = 2000L, verbose = TRUE, ...
) {
  force(data)
  unit <- match.arg(unit)
  if (!is.numeric(fis.min.threshold) || length(fis.min.threshold) != 1L ||
      !is.numeric(fis.max.threshold) || length(fis.max.threshold) != 1L ||
      fis.min.threshold > fis.max.threshold) {
    rlang::abort("FIS limits must be ordered numeric scalars.")
  }
  if (!is.numeric(strata.threshold) || length(strata.threshold) != 1L ||
      is.na(strata.threshold) || strata.threshold <= 0) {
    rlang::abort("`strata.threshold` must be positive.")
  }
  .paralog_check_probability(min.call.rate, "min.call.rate")
  .start <- tgbase::startup(package = "radr", f.name = "filter_fis", verbose = verbose)
  on.exit(tgbase::teardown(.start), add = TRUE)
  dots <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  unknown <- setdiff(names(dots), "path.folder")
  if (length(unknown)) rlang::abort(paste0("Unknown argument(s): ", paste(unknown, collapse = ", "), "."))
  parent <- dots$path.folder %||% getwd()
  path.folder <- radr_folder(
    rad.folder = paste0("filter_fis_", .start$file.date),
    path.folder = parent, prefix.int = TRUE
  )
  opened <- .filter_gds_open(data)
  gds <- opened$gds
  genometranslator::sync_gds(gds, verbose = FALSE)
  if (!genometranslator::detect_biallelic_markers(gds)) {
    rlang::abort("FIS filtering requires active biallelic markers.")
  }
  summary <- summarise_genomic_data(
    gds, strata = strata, group.column = group.column, by.strata = TRUE,
    chunk.size = chunk.size, write.files = FALSE, verbose = FALSE,
    internal = TRUE, path.folder = path.folder
  )$stratum.statistics |>
    dplyr::filter(.data$GROUP != "OVERALL") |>
    dplyr::mutate(
      ELIGIBLE = .data$CALL_RATE >= min.call.rate & is.finite(.data$FIS),
      VIOLATION = .data$ELIGIBLE &
        (.data$FIS < fis.min.threshold | .data$FIS > fis.max.threshold)
    )
  evidence <- summary |>
    dplyr::group_by(.data$VARIANT_ID) |>
    dplyr::summarise(
      N_ELIGIBLE_STRATA = sum(.data$ELIGIBLE),
      N_VIOLATING_STRATA = sum(.data$VIOLATION),
      PROPORTION_VIOLATING = dplyr::if_else(
        .data$N_ELIGIBLE_STRATA > 0,
        .data$N_VIOLATING_STRATA / .data$N_ELIGIBLE_STRATA, NA_real_
      ),
      MIN_FIS = if (any(.data$ELIGIBLE)) min(.data$FIS[.data$ELIGIBLE]) else NA_real_,
      MAX_FIS = if (any(.data$ELIGIBLE)) max(.data$FIS[.data$ELIGIBLE]) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      REMOVE = if (strata.threshold < 1) {
        .data$PROPORTION_VIOLATING >= strata.threshold
      } else .data$N_VIOLATING_STRATA >= strata.threshold
    )
  marker <- .sex_marker_metadata(gds, SeqArray::seqGetData(gds, "variant.id"))
  evidence <- dplyr::left_join(evidence, marker, by = "VARIANT_ID")
  if (unit == "locus") {
    if (!"LOCUS" %in% names(evidence)) {
      rlang::abort("`unit = \"locus\"` requires marker `LOCUS` metadata.")
    }
    remove.loci <- unique(evidence$LOCUS[evidence$REMOVE])
    evidence$REMOVE <- evidence$LOCUS %in% remove.loci
  }
  remove <- evidence$VARIANT_ID[evidence$REMOVE]
  readr::write_tsv(summary, file.path(path.folder, "fis_stratum_statistics.tsv"), na = "NA")
  readr::write_tsv(evidence, file.path(path.folder, "fis_filter_evidence.tsv"), na = "NA")
  .filter_gds_apply_markers(
    gds, remove, "filter.fis", path.folder, .start$file.date,
    "FIS limits / strata threshold",
    c(fis.min.threshold, fis.max.threshold, strata.threshold), verbose
  )
  if (verbose) .summary_message("FIS-filtered markers: ", length(remove), ".")
  gds
}
