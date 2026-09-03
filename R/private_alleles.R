#' Detect private alleles
#'
#' Identify REF or ALT alleles observed in exactly one stratum among the active
#' samples and biallelic markers in a GDS. Allele counts are calculated from
#' active diploid dosages; absence means a count of zero, not missing data.
#'
#' Results include call rate and allele count so a private allele supported by
#' one chromosome can be distinguished from a well-supported private allele.
#' Use `min.allele.count` and `min.call.rate` to define minimum evidence.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional sample metadata containing `INDIVIDUALS` and
#'   `group.column`. If `NULL`, metadata are read from the GDS.
#' @param group.column Metadata column defining the strata.
#' @param min.allele.count Minimum allele count in its only observed stratum.
#' @param min.call.rate Minimum marker call rate required in every stratum.
#' @param chunk.size Number of variants read at a time.
#' @param write.files Write result tables.
#' @param verbose Display progress and summary messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#' @return A `private_alleles` object with `private.alleles`,
#' `stratum.summary`, `marker.statistics`, and output information.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' private <- radr::private_alleles(
#'   "study.gds", strata = "sample_metadata.tsv"
#' )
#' }
private_alleles <- function(
    data, strata = NULL, group.column = "STRATA",
    min.allele.count = 1L, min.call.rate = 0,
    chunk.size = 2000L, write.files = TRUE, verbose = TRUE, ...
) {
  min.allele.count <- .paralog_check_count(
    min.allele.count, "min.allele.count", 1L
  )
  .paralog_check_probability(min.call.rate, "min.call.rate")
  context <- .analysis_gds_start(
    data, "private_alleles", write.files, verbose, ...
  )
  on.exit(.analysis_gds_finish(context), add = TRUE)
  gds <- context$gds
  if (!genometranslator::detect_biallelic_markers(gds)) {
    rlang::abort("Private-allele analysis requires active biallelic markers.")
  }
  stats <- summarise_genomic_data(
    gds, strata = strata, group.column = group.column, by.strata = TRUE,
    chunk.size = chunk.size, write.files = FALSE, verbose = FALSE,
    internal = TRUE, path.folder = context$path.folder
  )$stratum.statistics |>
    dplyr::filter(.data$GROUP != "OVERALL")
  groups <- unique(stats$GROUP)
  if (length(groups) < 2L) {
    rlang::abort("At least two strata are required to detect private alleles.")
  }
  eligible <- stats |>
    dplyr::group_by(.data$VARIANT_ID) |>
    dplyr::summarise(
      ALL_GROUPS_CALLABLE = all(.data$CALL_RATE >= min.call.rate),
      .groups = "drop"
    )
  long <- dplyr::bind_rows(
    stats |>
      dplyr::transmute(
        VARIANT_ID = .data$VARIANT_ID, GROUP = .data$GROUP,
        ALLELE = "REF", ALLELE_COUNT = .data$REF_ALLELE_COUNT,
        ALLELE_FREQUENCY = .data$REF_FREQUENCY,
        CALL_RATE = .data$CALL_RATE
      ),
    stats |>
      dplyr::transmute(
        VARIANT_ID = .data$VARIANT_ID, GROUP = .data$GROUP,
        ALLELE = "ALT", ALLELE_COUNT = .data$ALT_ALLELE_COUNT,
        ALLELE_FREQUENCY = .data$ALT_FREQUENCY,
        CALL_RATE = .data$CALL_RATE
      )
  )
  private <- long |>
    dplyr::group_by(.data$VARIANT_ID, .data$ALLELE) |>
    dplyr::mutate(N_STRATA_OBSERVED = sum(.data$ALLELE_COUNT > 0)) |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$N_STRATA_OBSERVED == 1L,
      .data$ALLELE_COUNT >= min.allele.count
    ) |>
    dplyr::left_join(eligible, by = "VARIANT_ID") |>
    dplyr::filter(.data$ALL_GROUPS_CALLABLE) |>
    dplyr::arrange(.data$GROUP, .data$VARIANT_ID, .data$ALLELE)
  marker.metadata <- .sex_marker_metadata(
    gds, SeqArray::seqGetData(gds, "variant.id")
  )
  private <- private |>
    dplyr::left_join(marker.metadata, by = "VARIANT_ID") |>
    dplyr::select(
      dplyr::any_of(c("VARIANT_ID", "MARKERS", "CHROM", "POS")),
      dplyr::everything()
    )
  summary <- private |>
    dplyr::count(.data$GROUP, .data$ALLELE, name = "N_PRIVATE_ALLELES") |>
    dplyr::group_by(.data$GROUP) |>
    dplyr::mutate(GROUP_TOTAL = sum(.data$N_PRIVATE_ALLELES)) |>
    dplyr::ungroup()
  output.files <- character()
  if (write.files) {
    output.files <- file.path(context$path.folder, c(
      "private_alleles.tsv", "private_alleles_by_stratum.tsv"
    ))
    readr::write_tsv(private, output.files[[1L]], na = "NA")
    readr::write_tsv(summary, output.files[[2L]], na = "NA")
  }
  if (verbose) {
    .summary_message("Private alleles retained: ", nrow(private), ".")
  }
  restored <- .analysis_gds_restore(context)
  out <- list(
    private.alleles = private, stratum.summary = summary,
    marker.statistics = stats, path.folder = context$path.folder,
    output.files = tibble::tibble(files = output.files),
    active.selection.restored = restored
  )
  class(out) <- c("private_alleles", class(out))
  out
}
