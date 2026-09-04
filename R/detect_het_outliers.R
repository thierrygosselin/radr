#' Diagnose heterozygosity outliers and genotype miscall
#'
#' Compare observed and Hardy-Weinberg expected genotype frequencies for active
#' biallelic GDS markers. The function reports standardized heterozygote
#' residuals by stratum and optionally fits the Anderson heterozygote-miscall
#' model by MCMC to a reproducible marker subset.
#'
#' The miscall model assumes true heterozygotes can be observed as either
#' homozygote with a common probability. Population structure, inbreeding,
#' null alleles, paralogy, allele dropout, and batch effects can violate this
#' model. Its estimate is therefore a diagnostic under explicit assumptions,
#' not a universal sequencing-error rate.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional metadata containing `INDIVIDUALS` and `group.column`.
#' @param group.column Metadata column defining strata.
#' @param z.threshold Absolute heterozygote-residual z-score used to flag
#'   marker-stratum combinations.
#' @param estimate.miscall Fit the heterozygote-miscall MCMC model.
#' @param nreps Number of MCMC iterations.
#' @param burn.in Number of initial MCMC iterations discarded.
#' @param max.mcmc.markers Maximum markers used by the MCMC model.
#' @param min.mcmc.samples Minimum called samples per stratum and marker.
#' @param min.mcmc.markers Minimum polymorphic markers per stratum.
#' @details MCMC runs independently within each metadata stratum. Insufficient
#' groups are skipped, never pooled. Minimums are safeguards, not guarantees
#' of statistical power. A single group assumes population homogeneity.
#' Posterior intervals are conditional on the model, not convergence checks.
#' @param seed Random seed used for marker subsampling and MCMC.
#' @param chunk.size Number of variants read at a time.
#' @param write.files Write statistics and figures.
#' @param verbose Display progress and summary messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#' @return A `detect_het_outliers` object containing marker statistics,
#' outliers, MCMC results when requested, plots, and output information.
#' @author Eric Anderson \email{eric.anderson@noaa.gov} and Thierry Gosselin
#'   \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' outliers <- radr::detect_het_outliers(
#'   "study.gds", strata = "sample_metadata.tsv", burn.in = 500
#' )
#' }
detect_het_outliers <- function(
    data, strata = NULL, group.column = "STRATA", z.threshold = 4,
    estimate.miscall = TRUE, nreps = 2000L, burn.in = 500L,
    max.mcmc.markers = 10000L, min.mcmc.samples = 10L,
    min.mcmc.markers = 10L, seed = 1L, chunk.size = 2000L,
    write.files = TRUE, verbose = TRUE, ...
) {
  .paralog_check_nonnegative(z.threshold, "z.threshold")
  .paralog_check_flag(estimate.miscall, "estimate.miscall")
  nreps <- .paralog_check_count(nreps, "nreps", 2L)
  burn.in <- .paralog_check_count(burn.in, "burn.in", 0L)
  max.mcmc.markers <- .paralog_check_count(
    max.mcmc.markers, "max.mcmc.markers", 1L
  )
  if (burn.in >= nreps) rlang::abort("`burn.in` must be smaller than `nreps`.")
  min.mcmc.samples <- .paralog_check_count(min.mcmc.samples, "min.mcmc.samples", 2L)
  min.mcmc.markers <- .paralog_check_count(min.mcmc.markers, "min.mcmc.markers", 1L)
  if (min.mcmc.markers > max.mcmc.markers) {
    rlang::abort("`min.mcmc.markers` must not exceed `max.mcmc.markers`.")
  }
  context <- .analysis_gds_start(
    data, "detect_het_outliers", write.files, verbose, ...
  )
  on.exit(.analysis_gds_finish(context), add = TRUE)
  gds <- context$gds
  if (!genometranslator::detect_biallelic_markers(gds)) {
    rlang::abort("Heterozygosity diagnostics require biallelic markers.")
  }
  meta <- .analysis_metadata(gds, strata, group.column, TRUE)
  summary <- summarise_genomic_data(
    gds, strata = meta$metadata, group.column = group.column,
    by.strata = TRUE, chunk.size = chunk.size, write.files = FALSE,
    verbose = FALSE, internal = TRUE, path.folder = context$path.folder
  )$stratum.statistics |>
    dplyr::mutate(
      EXPECTED_HET_COUNT = .data$EXPECTED_HETEROZYGOSITY *
        .data$NUMBER_CALLED,
      HET_VARIANCE = .data$NUMBER_CALLED *
        .data$EXPECTED_HETEROZYGOSITY *
        (1 - .data$EXPECTED_HETEROZYGOSITY),
      HET_Z = dplyr::if_else(
        .data$HET_VARIANCE > 0,
        (.data$HET - .data$EXPECTED_HET_COUNT) / sqrt(.data$HET_VARIANCE),
        NA_real_
      ),
      STATUS = dplyr::case_when(
        is.na(.data$HET_Z) ~ "UNINFORMATIVE",
        .data$HET_Z >= z.threshold ~ "HETEROZYGOTE_EXCESS",
        .data$HET_Z <= -z.threshold ~ "HETEROZYGOTE_DEFICIT",
        TRUE ~ "WITHIN_THRESHOLD"
      )
    )
  outliers <- dplyr::filter(
    summary,
    .data$STATUS %in% c("HETEROZYGOTE_EXCESS", "HETEROZYGOTE_DEFICIT")
  )
  trace <- tibble::tibble()
  posterior <- tibble::tibble()
  mcmc.markers <- list()
  if (estimate.miscall) {
    variant.id <- SeqArray::seqGetData(gds, "variant.id")
    groups <- unique(meta$groups)
    if (verbose && length(groups) == 1L) {
      .summary_message("One MCMC group: a homogeneous population is assumed.")
    }
    traces <- posteriors <- list()
    for (group in groups) {
      ids <- meta$sample.id[meta$groups == group]
      eligible <- summary$VARIANT_ID[
        summary$GROUP == group & summary$NUMBER_CALLED >= min.mcmc.samples &
          summary$EXPECTED_HETEROZYGOSITY > 0
      ]
      eligible <- variant.id[variant.id %in% eligible]
      set.seed(seed)
      selected <- if (length(eligible) > max.mcmc.markers) {
        sort(sample(eligible, max.mcmc.markers))
      } else eligible
      mcmc.markers[[group]] <- selected
      status <- if (length(ids) < min.mcmc.samples) {
        "INSUFFICIENT_SAMPLES"
      } else if (length(selected) < min.mcmc.markers) {
        "INSUFFICIENT_MARKERS"
      } else "FITTED"
      row <- tibble::tibble(
        GROUP = group, STATUS = status, N_MARKERS = length(selected),
        N_SAMPLES = length(ids), N_ITERATIONS = if (status == "FITTED") nreps else 0L,
        BURN_IN = burn.in, POSTERIOR_MEAN_MISCALL_RATE = NA_real_,
        POSTERIOR_SD_MISCALL_RATE = NA_real_,
        POSTERIOR_LOWER_95 = NA_real_, POSTERIOR_UPPER_95 = NA_real_,
        ESTIMATED_GENOTYPE_CHANGE_RATE = NA_real_, SEED = seed
      )
      if (verbose) .summary_message("MCMC stratum ", group, ": ", status, ".")
      if (status == "FITTED") {
        dosage <- .het_stratum_dosage(gds, selected, ids)
        set.seed(seed)
        model <- .estimate_m_matrix(dosage, nreps = nreps)
        retained <- model$m[(burn.in + 1L):nreps]
        traces[[group]] <- tibble::tibble(
          GROUP = group, ITERATION = seq_len(nreps), MISCALL_RATE = model$m
        )
        row$POSTERIOR_MEAN_MISCALL_RATE <- mean(retained)
        row$POSTERIOR_SD_MISCALL_RATE <- stats::sd(retained)
        row$POSTERIOR_LOWER_95 <- unname(stats::quantile(retained, 0.025))
        row$POSTERIOR_UPPER_95 <- unname(stats::quantile(retained, 0.975))
        row$ESTIMATED_GENOTYPE_CHANGE_RATE <- model$overall_geno_err_est
      }
      posteriors[[group]] <- row
    }
    trace <- dplyr::bind_rows(traces)
    posterior <- dplyr::bind_rows(posteriors)
  }
  residual.plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = .data$EXPECTED_HETEROZYGOSITY,
      y = .data$OBSERVED_HETEROZYGOSITY,
      colour = .data$STATUS
    )
  ) + ggplot2::geom_point(alpha = 0.45, na.rm = TRUE) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::facet_wrap(~ GROUP) + ggplot2::theme_bw() +
    ggplot2::labs(x = "Expected heterozygosity", y = "Observed heterozygosity")
  trace.plot <- if (nrow(trace) > 0L) {
    ggplot2::ggplot(trace, ggplot2::aes(.data$ITERATION, .data$MISCALL_RATE)) +
      ggplot2::geom_line() + ggplot2::geom_vline(xintercept = burn.in) +
      ggplot2::facet_wrap(~ GROUP) + ggplot2::theme_bw()
  } else NULL
  output.files <- character()
  if (write.files) {
    output.files <- file.path(context$path.folder, c(
      "heterozygosity_marker_statistics.tsv", "heterozygosity_outliers.tsv",
      "heterozygosity_diagnostic.png"
    ))
    readr::write_tsv(summary, output.files[[1L]], na = "NA")
    readr::write_tsv(outliers, output.files[[2L]], na = "NA")
    ggplot2::ggsave(output.files[[3L]], residual.plot, width = 10, height = 7)
    if (estimate.miscall) {
      extra <- file.path(context$path.folder, c(
        "miscall_mcmc_trace.tsv", "miscall_posterior_summary.tsv",
        "miscall_mcmc_trace.png"
      ))
      readr::write_tsv(trace, extra[[1L]])
      readr::write_tsv(posterior, extra[[2L]])
      if (!is.null(trace.plot)) {
        ggplot2::ggsave(extra[[3L]], trace.plot, width = 9, height = 5)
      } else extra <- extra[1:2]
      output.files <- c(output.files, extra)
    }
  }
  if (verbose) .summary_message("Heterozygosity outliers: ", nrow(outliers), ".")
  restored <- .analysis_gds_restore(context)
  out <- list(
    statistics = summary, outliers = outliers, mcmc.trace = trace,
    posterior.summary = posterior, mcmc.variant.id = mcmc.markers,
    diagnostic.plot = residual.plot, trace.plot = trace.plot,
    path.folder = context$path.folder,
    output.files = tibble::tibble(files = output.files),
    active.selection.restored = restored
  )
  class(out) <- c("detect_het_outliers", class(out))
  out
}

.het_stratum_dosage <- function(gds, variant.id, sample.id) {
  SeqArray::seqFilterPush(gds)
  on.exit(SeqArray::seqFilterPop(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
  .sex_get_dosage(gds, variant.id, sample.id)
}

.estimate_m_matrix <- function(dosage, nreps = 2000L, m.init = 0.1,
                               a0 = 0.5, a1 = 0.5, proposal.sd = 0.005) {
  D <- dosage
  D[is.na(D)] <- -1
  N0 <- colSums(D == 0); N1 <- colSums(D == 1); N2 <- colSums(D == 2)
  Z0 <- N0; Z1 <- N1; Z2 <- N2
  m <- rep(NA_real_, nreps); m[[1L]] <- m.init
  p <- rep(0.5, ncol(D))
  for (r in 2:nreps) {
    p <- stats::rbeta(length(Z0), a1 + 2 * Z2 + Z1, a0 + 2 * Z0 + Z1)
    proposed <- m[[r - 1L]] + stats::rnorm(1L, 0, proposal.sd)
    accepted <- FALSE
    if (proposed > 0 && proposed < 1) {
      log.like <- function(rate) sum(
        N0 * log((1 - p)^2 + rate * p * (1 - p)) +
        N1 * log((1 - rate) * 2 * p * (1 - p)) +
        N2 * log(p^2 + rate * p * (1 - p))
      )
      accepted <- log(stats::runif(1L)) <
        log.like(proposed) - log.like(m[[r - 1L]])
    }
    m[[r]] <- if (accepted) proposed else m[[r - 1L]]
    A0 <- stats::rbinom(length(N0), N0, (m[[r]] * p) / (1 - p + m[[r]] * p))
    A2 <- stats::rbinom(length(N2), N2, (m[[r]] * (1 - p)) / (p + m[[r]] * (1 - p)))
    Z0 <- N0 - A0; Z1 <- N1 + A0 + A2; Z2 <- N2 - A2
  }
  simulated <- simulate_genos_from_posterior(D, p, m[[nreps]])
  changed <- simulated != D; changed[D == -1] <- NA
  list(m = m, overall_geno_err_est = mean(changed, na.rm = TRUE))
}
