#' Filter Hardy-Weinberg disequilibrium
#'
#' Apply exact mid-p Hardy-Weinberg tests to active biallelic GDS markers within
#' strata and remove markers departing from equilibrium in enough eligible
#' strata. Tests use current active samples and are recalculated on every run.
#'
#' HWE departures are not automatically technical failures. Population
#' mixture, selection, inbreeding, sex linkage, null alleles, paralogy, and
#' genotyping error can all generate departures. This function is a deliberate
#' filter with full evidence tables, not a substitute for diagnosis.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional metadata containing `INDIVIDUALS` and `group.column`.
#' @param group.column Metadata column defining strata.
#' @param p.threshold Mid-p or adjusted-p threshold.
#' @param strata.threshold Number or proportion of eligible strata that must
#'   depart before marker removal.
#' @param adjustment Multiple-testing method: `"BH"` or `"none"`. Adjustment
#'   is performed independently within each stratum.
#' @param min.samples Minimum samples in a stratum.
#' @param min.call.rate Minimum marker call rate within a stratum.
#' @param chunk.size Number of variants read at a time.
#' @param verbose Display progress and summary messages.
#' @param ... Common argument `path.folder`.
#' @return The filtered open GDS connection with evidence and audit files.
#' @references Wigginton JE, Cutler DJ, Abecasis GR (2005). A note on exact
#' tests of Hardy-Weinberg equilibrium. *AJHG*, 76, 887-893.
#' @references Graffelman J, Moreno V (2013). The mid p-value in exact tests for
#' Hardy-Weinberg equilibrium. *SAGMB*, 12, 433-448.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' genome <- radr::filter_hwe(
#'   genome, strata = "sample_metadata.tsv", p.threshold = 1e-4,
#'   strata.threshold = 2
#' )
#' }
filter_hwe <- function(
    data, strata = NULL, group.column = "STRATA", p.threshold = 1e-4,
    strata.threshold = 1L, adjustment = c("BH", "none"),
    min.samples = 10L, min.call.rate = 0.8, chunk.size = 2000L,
    verbose = TRUE, ...
) {
  force(data)
  adjustment <- match.arg(adjustment)
  .paralog_check_probability(p.threshold, "p.threshold")
  .paralog_check_probability(min.call.rate, "min.call.rate")
  min.samples <- .paralog_check_count(min.samples, "min.samples", 2L)
  if (!is.numeric(strata.threshold) || length(strata.threshold) != 1L ||
      is.na(strata.threshold) || strata.threshold <= 0) {
    rlang::abort("`strata.threshold` must be positive.")
  }
  .start <- tgbase::startup(package = "radr", f.name = "filter_hwe", verbose = verbose)
  on.exit(tgbase::teardown(.start), add = TRUE)
  dots <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  unknown <- setdiff(names(dots), "path.folder")
  if (length(unknown)) rlang::abort(paste0("Unknown argument(s): ", paste(unknown, collapse = ", "), "."))
  path.folder <- radr_folder(
    rad.folder = paste0("filter_hwe_", .start$file.date),
    path.folder = dots$path.folder %||% getwd(), prefix.int = TRUE
  )
  opened <- .filter_gds_open(data)
  gds <- opened$gds
  genometranslator::sync_gds(gds, verbose = FALSE)
  if (!genometranslator::detect_biallelic_markers(gds)) {
    rlang::abort("HWE filtering requires active biallelic markers.")
  }
  summary <- summarise_genomic_data(
    gds, strata = strata, group.column = group.column, by.strata = TRUE,
    chunk.size = chunk.size, write.files = FALSE, verbose = FALSE,
    internal = TRUE, path.folder = path.folder
  )$stratum.statistics |>
    dplyr::filter(.data$GROUP != "OVERALL") |>
    dplyr::mutate(
      ELIGIBLE = .data$N_SAMPLES >= min.samples &
        .data$CALL_RATE >= min.call.rate &
        .data$MINOR_ALLELE_COUNT > 0
    )
  tests <- summary |>
    dplyr::group_by(.data$GROUP) |>
    dplyr::group_modify(function(.x, .y) {
      p <- rep(NA_real_, nrow(.x))
      eligible <- which(.x$ELIGIBLE)
      if (length(eligible)) {
        counts <- as.matrix(.x[eligible, c("HOM_REF", "HET", "HOM_ALT")])
        p[eligible] <- HardyWeinberg::HWExactStats(
          X = counts, plinkcode = TRUE, midp = TRUE
        )
      }
      .x$MID_P_VALUE <- p
      .x$ADJUSTED_P_VALUE <- if (adjustment == "none") p else
        stats::p.adjust(p, method = adjustment)
      .x
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      DEPARTURE = .data$ELIGIBLE &
        .data$ADJUSTED_P_VALUE <= p.threshold
    )
  evidence <- tests |>
    dplyr::group_by(.data$VARIANT_ID) |>
    dplyr::summarise(
      N_ELIGIBLE_STRATA = sum(.data$ELIGIBLE),
      N_DEPARTING_STRATA = sum(.data$DEPARTURE, na.rm = TRUE),
      PROPORTION_DEPARTING = dplyr::if_else(
        .data$N_ELIGIBLE_STRATA > 0,
        .data$N_DEPARTING_STRATA / .data$N_ELIGIBLE_STRATA, NA_real_
      ),
      MIN_P_VALUE = if (any(.data$ELIGIBLE))
        min(.data$ADJUSTED_P_VALUE, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(REMOVE = if (strata.threshold < 1) {
      .data$PROPORTION_DEPARTING >= strata.threshold
    } else .data$N_DEPARTING_STRATA >= strata.threshold)
  remove <- evidence$VARIANT_ID[evidence$REMOVE]
  sensitivity <- tidyr::crossing(
    P_THRESHOLD = c(0.05, 0.01, 0.001, 1e-4, 1e-5),
    STRATA_THRESHOLD = seq_len(dplyr::n_distinct(tests$GROUP))
  )
  tests.by.variant <- split(tests, tests$VARIANT_ID)
  sensitivity$N_MARKERS <- mapply(
    FUN = function(p.value, n.strata) {
      sum(vapply(tests.by.variant, function(x) {
        sum(
          x$ELIGIBLE & !is.na(x$ADJUSTED_P_VALUE) &
            x$ADJUSTED_P_VALUE <= p.value
        ) >= n.strata
      }, logical(1)))
    },
    p.value = sensitivity$P_THRESHOLD,
    n.strata = sensitivity$STRATA_THRESHOLD
  )
  readr::write_tsv(tests, file.path(path.folder, "hwe_stratum_tests.tsv"), na = "NA")
  readr::write_tsv(evidence, file.path(path.folder, "hwe_filter_evidence.tsv"), na = "NA")
  readr::write_tsv(sensitivity, file.path(path.folder, "hwe_sensitivity.tsv"), na = "NA")
  .filter_gds_apply_markers(
    gds, remove, "filter.hwe", path.folder, .start$file.date,
    "p / strata / adjustment",
    c(p.threshold, strata.threshold, adjustment), verbose
  )
  if (verbose) .summary_message("HWE-filtered markers: ", length(remove), ".")
  gds
}
