#' Estimate nucleotide diversity
#'
#' Estimate nucleotide diversity (pi), the mean number of pairwise nucleotide
#' differences, from active biallelic variants in a GDS. Estimates are
#' calculated independently for each stratum and for the pooled sample using
#' the finite-sample correction `2n / (2n - 1)`.
#'
#' With `scale = "per_variant"`, pi is averaged across active variant records.
#' This describes diversity among the analysed sites and is not a genome-wide
#' per-base estimate. With `scale = "per_base"`, the sum of pairwise
#' differences is divided by `callable.bases`; this denominator must come from
#' an appropriate callable-genome or callable-target analysis.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional metadata containing `INDIVIDUALS` and `group.column`.
#' @param group.column Metadata column defining strata.
#' @param scale Either `"per_variant"` or `"per_base"`.
#' @param callable.bases Positive callable-base denominator used with
#'   `scale = "per_base"`. A scalar applies to all groups; a named vector can
#'   provide group-specific denominators.
#' @param chunk.size Number of variants read at a time.
#' @param write.files Write result tables and a plot.
#' @param verbose Display progress and summary messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#' @return A `pi` object containing `marker.statistics`,
#' `population.statistics`, `individual.statistics`, `plot`, and output
#' information.
#' @references Nei M, Li WH (1979). Mathematical model for studying genetic
#' variation in terms of restriction endonucleases. *PNAS*, 76, 5269-5273.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' diversity <- radr::pi("study.gds", strata = "sample_metadata.tsv")
#' }
pi <- function(
    data, strata = NULL, group.column = "STRATA",
    scale = c("per_variant", "per_base"), callable.bases = NULL,
    chunk.size = 2000L, write.files = TRUE, verbose = TRUE, ...
) {
  scale <- match.arg(scale)
  context <- .analysis_gds_start(data, "pi", write.files, verbose, ...)
  on.exit(.analysis_gds_finish(context), add = TRUE)
  gds <- context$gds
  if (!genometranslator::detect_biallelic_markers(gds)) {
    rlang::abort("Nucleotide-diversity estimation requires biallelic markers.")
  }
  if (scale == "per_base" && (is.null(callable.bases) ||
      any(!is.finite(callable.bases)) || any(callable.bases <= 0))) {
    rlang::abort(
      "`callable.bases` must be positive when `scale = \"per_base\"`."
    )
  }
  meta <- .analysis_metadata(gds, strata, group.column, TRUE)
  summary <- summarise_genomic_data(
    gds, strata = meta$metadata, group.column = group.column,
    by.strata = TRUE, chunk.size = chunk.size, write.files = FALSE,
    verbose = FALSE, internal = TRUE, path.folder = context$path.folder
  )
  marker <- summary$stratum.statistics |>
    dplyr::mutate(
      PAIRWISE_DIFFERENCE = dplyr::if_else(
        .data$NUMBER_CALLED > 0,
        .data$EXPECTED_HETEROZYGOSITY *
          (2 * .data$NUMBER_CALLED) / (2 * .data$NUMBER_CALLED - 1),
        NA_real_
      )
    )
  population <- marker |>
    dplyr::group_by(.data$GROUP) |>
    dplyr::summarise(
      N_SAMPLES = max(.data$N_SAMPLES),
      N_VARIANTS = dplyr::n(),
      N_CALLED_VARIANTS = sum(.data$NUMBER_CALLED > 0),
      SUM_PAIRWISE_DIFFERENCES = sum(
        .data$PAIRWISE_DIFFERENCE, na.rm = TRUE
      ),
      .groups = "drop"
    )
  if (scale == "per_variant") {
    population <- population |>
      dplyr::mutate(
        DENOMINATOR = .data$N_CALLED_VARIANTS,
        PI = .data$SUM_PAIRWISE_DIFFERENCES / .data$DENOMINATOR,
        SCALE = "per_variant"
      )
  } else {
    denominator <- .pi_denominator(callable.bases, population$GROUP)
    population <- population |>
      dplyr::mutate(
        DENOMINATOR = denominator,
        PI = .data$SUM_PAIRWISE_DIFFERENCES / .data$DENOMINATOR,
        SCALE = "per_base"
      )
  }

  variant.id <- SeqArray::seqGetData(gds, "variant.id")
  chunks <- split(seq_along(variant.id), ceiling(seq_along(variant.id) / chunk.size))
  individual.parts <- lapply(chunks, function(index) {
    dosage <- .sex_get_dosage(gds, variant.id[index], meta$sample.id)
    tibble::tibble(
      INDIVIDUALS = meta$sample.id,
      CALLED = rowSums(!is.na(dosage)),
      HETEROZYGOUS = rowSums(dosage == 1, na.rm = TRUE)
    )
  })
  individual <- dplyr::bind_rows(individual.parts) |>
    dplyr::group_by(.data$INDIVIDUALS) |>
    dplyr::summarise(
      NUMBER_CALLED = sum(.data$CALLED),
      NUMBER_HETEROZYGOUS = sum(.data$HETEROZYGOUS),
      INDIVIDUAL_HETEROZYGOSITY = .data$NUMBER_HETEROZYGOUS /
        .data$NUMBER_CALLED,
      .groups = "drop"
    ) |>
    dplyr::left_join(
      meta$metadata[, c("INDIVIDUALS", group.column), drop = FALSE],
      by = "INDIVIDUALS"
    )
  plot <- ggplot2::ggplot(
    dplyr::filter(population, .data$GROUP != "OVERALL"),
    ggplot2::aes(x = .data$GROUP, y = .data$PI)
  ) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = group.column, y = paste0("Pi (", scale, ")"))
  output.files <- character()
  if (write.files) {
    output.files <- file.path(context$path.folder, c(
      "pi_marker_statistics.tsv", "pi_population_statistics.tsv",
      "individual_heterozygosity.tsv", "pi_by_stratum.png"
    ))
    readr::write_tsv(marker, output.files[[1L]], na = "NA")
    readr::write_tsv(population, output.files[[2L]], na = "NA")
    readr::write_tsv(individual, output.files[[3L]], na = "NA")
    ggplot2::ggsave(output.files[[4L]], plot, width = 8, height = 5, dpi = 300)
  }
  if (verbose) .summary_message("Nucleotide diversity groups: ", nrow(population), ".")
  restored <- .analysis_gds_restore(context)
  out <- list(
    marker.statistics = marker, population.statistics = population,
    individual.statistics = individual, plot = plot,
    path.folder = context$path.folder,
    output.files = tibble::tibble(files = output.files),
    active.selection.restored = restored
  )
  class(out) <- c("pi", class(out))
  out
}

.pi_denominator <- function(callable.bases, groups) {
  if (length(callable.bases) == 1L) return(rep(callable.bases, length(groups)))
  if (is.null(names(callable.bases)) || !all(groups %in% names(callable.bases))) {
    rlang::abort(
      "Group-specific `callable.bases` must be named for every output group."
    )
  }
  as.numeric(callable.bases[groups])
}
