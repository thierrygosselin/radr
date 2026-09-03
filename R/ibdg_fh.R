#' Estimate individual genomic inbreeding with FH
#'
#' Estimate the method-of-moments FH coefficient from active biallelic GDS
#' genotypes. For each individual, observed homozygous genotype counts are
#' compared with expected homozygous counts calculated from allele frequencies
#' at exactly the markers called in that individual.
#'
#' The implemented estimator is `(O - E) / (N - E)`, matching the count-scale
#' form used by PLINK's heterozygosity report. `FH_STRATUM` uses allele
#' frequencies from the individual's stratum; `FH_POOLED` uses frequencies
#' from all active samples. These are genomic method-of-moments estimators, not
#' direct proof that chromosome segments are identical by descent.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional metadata containing `INDIVIDUALS` and `group.column`.
#' @param group.column Metadata column defining strata.
#' @param chunk.size Number of variants read at a time.
#' @param min.call.rate Minimum individual call rate retained in summaries.
#' @param write.files Write result tables and figures.
#' @param verbose Display progress and summary messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#' @return An `ibdg_fh` object with individual and stratum statistics, plots,
#' and output information.
#' @references Keller MC, Visscher PM, Goddard ME (2011). Quantification of
#' inbreeding due to distant ancestors and its detection using dense SNP data.
#' *Genetics*, 189, 237-249.
#' @references Kardos M, Luikart G, Allendorf FW (2015). Measuring individual
#' inbreeding in the age of genomics. *Heredity*, 115, 63-72.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' fh <- radr::ibdg_fh("study.gds", strata = "sample_metadata.tsv")
#' }
ibdg_fh <- function(
    data, strata = NULL, group.column = "STRATA", chunk.size = 2000L,
    min.call.rate = 0, write.files = TRUE, verbose = TRUE, ...
) {
  .paralog_check_probability(min.call.rate, "min.call.rate")
  context <- .analysis_gds_start(data, "ibdg_fh", write.files, verbose, ...)
  on.exit(.analysis_gds_finish(context), add = TRUE)
  gds <- context$gds
  if (!genometranslator::detect_biallelic_markers(gds)) {
    rlang::abort("FH estimation requires active biallelic markers.")
  }
  meta <- .analysis_metadata(gds, strata, group.column, TRUE)
  variant.id <- SeqArray::seqGetData(gds, "variant.id")
  summary <- summarise_genomic_data(
    gds, strata = meta$metadata, group.column = group.column,
    by.strata = TRUE, chunk.size = chunk.size, write.files = FALSE,
    verbose = FALSE, internal = TRUE, path.folder = context$path.folder
  )$stratum.statistics
  frequencies <- summary |>
    dplyr::select(
      "VARIANT_ID", "GROUP", "REF_FREQUENCY", "ALT_FREQUENCY"
    )
  chunks <- split(seq_along(variant.id), ceiling(seq_along(variant.id) / chunk.size))
  accum <- lapply(chunks, function(index) {
    ids <- variant.id[index]
    dosage <- .sex_get_dosage(gds, ids, meta$sample.id)
    called <- !is.na(dosage)
    observed <- rowSums(dosage == 0 | dosage == 2, na.rm = TRUE)
    n.called <- rowSums(called)
    pooled <- frequencies |>
      dplyr::filter(.data$GROUP == "OVERALL", .data$VARIANT_ID %in% ids) |>
      dplyr::arrange(match(.data$VARIANT_ID, ids))
    expected.pooled <- pooled$REF_FREQUENCY^2 + pooled$ALT_FREQUENCY^2
    e.pooled <- rowSums(
      sweep(called, 2L, expected.pooled, `*`), na.rm = TRUE
    )
    e.stratum <- numeric(length(meta$sample.id))
    for (group in unique(meta$groups)) {
      rows <- which(meta$groups == group)
      freq <- frequencies |>
        dplyr::filter(.data$GROUP == group, .data$VARIANT_ID %in% ids) |>
        dplyr::arrange(match(.data$VARIANT_ID, ids))
      expected <- freq$REF_FREQUENCY^2 + freq$ALT_FREQUENCY^2
      e.stratum[rows] <- rowSums(
        sweep(called[rows, , drop = FALSE], 2L, expected, `*`),
        na.rm = TRUE
      )
    }
    tibble::tibble(
      INDIVIDUALS = meta$sample.id, NUMBER_CALLED = n.called,
      OBSERVED_HOMOZYGOTES = observed,
      EXPECTED_HOMOZYGOTES_STRATUM = e.stratum,
      EXPECTED_HOMOZYGOTES_POOLED = e.pooled
    )
  })
  individuals <- dplyr::bind_rows(accum) |>
    dplyr::group_by(.data$INDIVIDUALS) |>
    dplyr::summarise(
      NUMBER_CALLED = sum(.data$NUMBER_CALLED),
      OBSERVED_HOMOZYGOTES = sum(.data$OBSERVED_HOMOZYGOTES),
      EXPECTED_HOMOZYGOTES_STRATUM = sum(
        .data$EXPECTED_HOMOZYGOTES_STRATUM
      ),
      EXPECTED_HOMOZYGOTES_POOLED = sum(
        .data$EXPECTED_HOMOZYGOTES_POOLED
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      CALL_RATE = .data$NUMBER_CALLED / length(variant.id),
      FH_STRATUM = (.data$OBSERVED_HOMOZYGOTES -
        .data$EXPECTED_HOMOZYGOTES_STRATUM) /
        (.data$NUMBER_CALLED - .data$EXPECTED_HOMOZYGOTES_STRATUM),
      FH_POOLED = (.data$OBSERVED_HOMOZYGOTES -
        .data$EXPECTED_HOMOZYGOTES_POOLED) /
        (.data$NUMBER_CALLED - .data$EXPECTED_HOMOZYGOTES_POOLED),
      ELIGIBLE = .data$CALL_RATE >= min.call.rate
    ) |>
    dplyr::left_join(
      meta$metadata[, c("INDIVIDUALS", group.column), drop = FALSE],
      by = "INDIVIDUALS"
    )
  names(individuals)[names(individuals) == group.column] <- "GROUP"
  group.summary <- individuals |>
    dplyr::filter(.data$ELIGIBLE) |>
    dplyr::group_by(.data$GROUP) |>
    dplyr::summarise(
      N_INDIVIDUALS = dplyr::n(),
      MEAN_CALL_RATE = mean(.data$CALL_RATE),
      MEAN_FH_STRATUM = mean(.data$FH_STRATUM, na.rm = TRUE),
      SD_FH_STRATUM = stats::sd(.data$FH_STRATUM, na.rm = TRUE),
      MEAN_FH_POOLED = mean(.data$FH_POOLED, na.rm = TRUE),
      SD_FH_POOLED = stats::sd(.data$FH_POOLED, na.rm = TRUE),
      .groups = "drop"
    )
  plot <- ggplot2::ggplot(
    dplyr::filter(individuals, .data$ELIGIBLE),
    ggplot2::aes(x = .data$GROUP, y = .data$FH_STRATUM)
  ) +
    ggplot2::geom_violin(fill = NA) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.5) +
    ggplot2::coord_flip() + ggplot2::theme_bw() +
    ggplot2::labs(x = group.column, y = "Individual FH")
  output.files <- character()
  if (write.files) {
    output.files <- file.path(context$path.folder, c(
      "fh_individual_statistics.tsv", "fh_stratum_summary.tsv", "fh.png"
    ))
    readr::write_tsv(individuals, output.files[[1L]], na = "NA")
    readr::write_tsv(group.summary, output.files[[2L]], na = "NA")
    ggplot2::ggsave(output.files[[3L]], plot, width = 8, height = 6, dpi = 300)
  }
  if (verbose) .summary_message("FH estimates retained: ", sum(individuals$ELIGIBLE), ".")
  restored <- .analysis_gds_restore(context)
  out <- list(
    individual.statistics = individuals,
    stratum.summary = group.summary, plot = plot,
    path.folder = context$path.folder,
    output.files = tibble::tibble(files = output.files),
    active.selection.restored = restored
  )
  class(out) <- c("ibdg_fh", class(out))
  out
}
