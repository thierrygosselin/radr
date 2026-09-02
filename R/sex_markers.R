#' Screen a GDS for candidate sex-linked markers
#'
#' `sexy_markers()` compares marker presence, heterozygosity, and, when
#' available, normalized read depth between known females and males. It is a
#' read-only screen: active GDS sample and variant selections are respected and
#' restored, and no filter is applied or written to the GDS.
#'
#' Y-like markers are present mainly in males, W-like markers mainly in females,
#' X-like markers have greater heterozygosity or normalized depth in females,
#' and Z-like markers show the reverse. These labels are candidates, not proof
#' of chromosomal location. Sex-specific dropout, population structure, plate
#' or library effects, paralogy, and mapping bias can produce similar patterns.
#'
#' @section Statistical tests:
#' Presence and heterozygosity use two-sample score tests for proportions.
#' Normalized depth uses Welch tests on `log2(1 + depth)`. P-values are adjusted
#' separately by method with Benjamini-Hochberg FDR. Candidate selection always
#' uses explicit effect-size thresholds; `require.significance = TRUE` also
#' requires the method-specific FDR to pass `fdr.threshold`.
#'
#' @section Read-depth normalization:
#' When depth is available, each sample is divided by its mean positive depth
#' across active markers. This reduces sample-wide sequencing-depth differences
#' but cannot remove marker-by-batch interactions or confounding between sex and
#' plate, lane, library, population, or sampling group.
#'
#' @param data A GDS filepath or an open `SeqVarGDSClass` object.
#' @param strata A data frame or tab-delimited file containing `INDIVIDUALS`
#'   and `sex.column`. If `NULL`, individual metadata are read from the GDS.
#' @param sex.column Metadata column containing sex. `F`, `female`, `M`, and
#'   `male` are recognized without regard to case; other values are unknown.
#' @param presence.threshold Minimum presence in the expected sex for Y/W.
#' @param absence.threshold Maximum presence in the other sex for Y/W.
#' @param min.heterozygosity.difference Minimum absolute female-minus-male
#'   heterozygosity difference for X/Z candidates.
#' @param coverage.ratio.threshold Minimum larger-to-smaller normalized depth
#'   ratio for X/Z coverage candidates.
#' @param coverage.threshold Minimum read depth considered present. Without
#'   depth data, a non-missing genotype is considered present.
#' @param min.samples.per.sex Minimum number of known females and males.
#' @param fdr.threshold Maximum FDR when significance is required.
#' @param require.significance Require effect size and method-specific FDR.
#' @param chunk.size Number of variants read from the GDS at a time.
#' @param save.plots Save diagnostic plots in the result folder.
#' @param plot.formats One or both of `"png"` and `"pdf"`.
#' @param folder.name Optional result-folder stem. Default: `sexy_markers`.
#' @param verbose Display progress and summary messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#'
#' @return A `sexy_markers` object containing marker statistics, candidates,
#'   an assignment-ready Y/W panel, sample summaries, a metadata audit, plots,
#'   and the result-folder path.
#' @export
#' @examples
#' \dontrun{
#' result <- radr::sexy_markers(
#'   data = "study.gds",
#'   strata = "sample_metadata.tsv",
#'   sex.column = "SEX"
#' )
#' result$candidates
#' }
sexy_markers <- function(
    data,
    strata = NULL,
    sex.column = "STRATA",
    presence.threshold = 0.9,
    absence.threshold = 0.1,
    min.heterozygosity.difference = 0.2,
    coverage.ratio.threshold = 1.5,
    coverage.threshold = 1,
    min.samples.per.sex = 5L,
    fdr.threshold = 0.05,
    require.significance = FALSE,
    chunk.size = 2000L,
    save.plots = TRUE,
    plot.formats = c("png", "pdf"),
    folder.name = NULL,
    verbose = TRUE,
    ...
) {
  force(data)
  .start <- tgbase::startup(
    package = "radr", f.name = "sexy_markers", verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)

  dots <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  unknown.dots <- setdiff(names(dots), c("path.folder", "internal"))
  if (length(unknown.dots)) {
    rlang::abort(paste0(
      "Unknown argument(s): ", paste(unknown.dots, collapse = ", "), "."
    ))
  }
  parent.folder <- dots$path.folder %||% getwd()
  internal <- isTRUE(dots$internal)

  .sex_check_probability(presence.threshold, "presence.threshold")
  .sex_check_probability(absence.threshold, "absence.threshold")
  .sex_check_probability(fdr.threshold, "fdr.threshold")
  .sex_check_probability(
    min.heterozygosity.difference,
    "min.heterozygosity.difference"
  )
  if (absence.threshold >= presence.threshold) {
    rlang::abort("`absence.threshold` must be below `presence.threshold`.")
  }
  if (!is.numeric(coverage.ratio.threshold) ||
      length(coverage.ratio.threshold) != 1L ||
      !is.finite(coverage.ratio.threshold) ||
      coverage.ratio.threshold <= 1) {
    rlang::abort("`coverage.ratio.threshold` must be greater than one.")
  }
  .sex_check_nonnegative(coverage.threshold, "coverage.threshold")
  min.samples.per.sex <- .sex_check_count(
    min.samples.per.sex, "min.samples.per.sex", 2L
  )
  chunk.size <- .sex_check_count(chunk.size, "chunk.size", 1L)
  .sex_check_flag(require.significance, "require.significance")
  .sex_check_flag(save.plots, "save.plots")
  .sex_check_flag(verbose, "verbose")
  if (!is.character(sex.column) || length(sex.column) != 1L ||
      is.na(sex.column) || !nzchar(sex.column)) {
    rlang::abort("`sex.column` must be one non-empty column name.")
  }
  plot.formats <- unique(tolower(as.character(plot.formats)))
  if (!length(plot.formats) || any(!plot.formats %in% c("png", "pdf"))) {
    rlang::abort("`plot.formats` must contain `png`, `pdf`, or both.")
  }

  folder.stem <- folder.name %||% "sexy_markers"
  if (!is.character(folder.stem) || length(folder.stem) != 1L ||
      is.na(folder.stem) || !nzchar(folder.stem)) {
    rlang::abort("`folder.name` must be NULL or one non-empty name.")
  }
  dir.create(parent.folder, recursive = TRUE, showWarnings = FALSE)
  parent.folder <- normalizePath(parent.folder, mustWork = TRUE)
  path.folder <- if (internal) {
    parent.folder
  } else {
    radr_folder(
      rad.folder = paste0(folder.stem, "_", file.date),
      path.folder = parent.folder,
      prefix.int = TRUE
    )
  }
  dir.create(path.folder, recursive = TRUE, showWarnings = FALSE)
  if (verbose && !internal) {
    .sex_message("Folder created: ", basename(path.folder))
  }
  arguments <- tibble::tibble(
    argument = c(
      "sex.column", "presence.threshold", "absence.threshold",
      "min.heterozygosity.difference", "coverage.ratio.threshold",
      "coverage.threshold", "min.samples.per.sex", "fdr.threshold",
      "require.significance", "chunk.size", "save.plots", "plot.formats"
    ),
    value = c(
      sex.column, presence.threshold, absence.threshold,
      min.heterozygosity.difference, coverage.ratio.threshold,
      coverage.threshold, min.samples.per.sex, fdr.threshold,
      require.significance, chunk.size, save.plots,
      paste(plot.formats, collapse = ",")
    )
  )
  args.file <- file.path(
    path.folder, paste0("radr_sexy_markers_args_", file.date, ".tsv")
  )
  readr::write_tsv(arguments, args.file)
  if (verbose) {
    .sex_message("Function call and arguments stored in: ", basename(args.file))
  }

  opened.here <- FALSE
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.character(data) && length(data) == 1L && !is.na(data) &&
             file.exists(data) && grepl("\\.gds$", data, ignore.case = TRUE)) {
    gds <- SeqArray::seqOpen(data)
    opened.here <- TRUE
  } else {
    rlang::abort(
      "`data` must be a GDS filepath or an open SeqVarGDSClass object."
    )
  }
  on.exit({
    if (opened.here) try(SeqArray::seqClose(gds), silent = TRUE)
  }, add = TRUE)

  selection.before <- SeqArray::seqGetFilter(gds)
  SeqArray::seqFilterPush(gds)
  filter.pushed <- TRUE
  on.exit({
    if (filter.pushed) try(SeqArray::seqFilterPop(gds), silent = TRUE)
  }, add = TRUE)

  sample.id <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  variant.id <- SeqArray::seqGetData(gds, "variant.id")
  if (!length(sample.id) || !length(variant.id)) {
    rlang::abort("The active GDS selection contains no samples or markers.")
  }
  metadata <- .sex_read_metadata(strata, gds, sample.id, sex.column)
  sample.id <- sample.id[sample.id %in% metadata$INDIVIDUALS]
  metadata <- metadata[match(sample.id, metadata$INDIVIDUALS), , drop = FALSE]
  SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
  sex <- .sex_normalise(metadata[[sex.column]])
  metadata$SEX_RADR <- sex
  female <- sex == "F"
  male <- sex == "M"
  n.female <- sum(female)
  n.male <- sum(male)
  if (n.female < min.samples.per.sex || n.male < min.samples.per.sex) {
    rlang::abort(paste0(
      "At least ", min.samples.per.sex,
      " known females and males are required. Found ", n.female,
      " females and ", n.male, " males."
    ))
  }
  if (verbose) {
    .sex_message(
      "Active GDS selection: ", length(variant.id), " marker(s), ",
      length(sample.id), " sample(s)."
    )
    .sex_message(
      "Known sex: ", n.female, " female(s), ", n.male,
      " male(s), ", sum(sex == "U"), " unknown."
    )
  }
  if (max(n.female, n.male) / min(n.female, n.male) >= 2) {
    warning(
      "Known-sex sample sizes differ by at least two-fold; review effects.",
      call. = FALSE
    )
  }

  marker.metadata <- .sex_marker_metadata(gds, variant.id)
  chunks <- split(
    seq_along(variant.id), ceiling(seq_along(variant.id) / chunk.size)
  )
  if (verbose) {
    .sex_message("Reading marker statistics in ", length(chunks), " chunk(s)...")
  }

  depth.sum <- numeric(length(sample.id))
  depth.n <- integer(length(sample.id))
  depth.available <- FALSE
  for (idx in chunks) {
    depth <- .sex_get_depth(gds, variant.id[idx], sample.id)
    if (!is.null(depth)) {
      depth.available <- TRUE
      positive <- is.finite(depth) & depth > 0
      depth.sum <- depth.sum + rowSums(replace(depth, !positive, 0))
      depth.n <- depth.n + rowSums(positive)
    }
  }
  depth.scale <- depth.sum / depth.n
  depth.scale[!is.finite(depth.scale) | depth.scale <= 0] <- 1
  if (verbose) {
    .sex_message(if (depth.available) {
      "Read depth found; coverage will be normalized within each sample."
    } else {
      "Read depth not found; the coverage method will be skipped."
    })
  }

  pieces <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    idx <- chunks[[i]]
    dosage <- .sex_get_dosage(gds, variant.id[idx], sample.id)
    depth <- if (depth.available) {
      .sex_get_depth(gds, variant.id[idx], sample.id)
    } else NULL
    pieces[[i]] <- .sex_chunk_statistics(
      dosage, depth, depth.scale, female, male, coverage.threshold
    )
    if (verbose && (i == length(chunks) || i %% 10L == 0L)) {
      .sex_message("Processed ", i, " of ", length(chunks), " chunk(s).")
    }
  }
  statistics <- dplyr::bind_cols(
    marker.metadata, dplyr::bind_rows(pieces)
  )
  statistics$presence_fdr <- stats::p.adjust(statistics$presence_p, "BH")
  statistics$heterozygosity_fdr <- stats::p.adjust(
    statistics$heterozygosity_p, "BH"
  )
  statistics$coverage_fdr <- stats::p.adjust(statistics$coverage_p, "BH")

  significant <- function(x) {
    !require.significance | (!is.na(x) & x <= fdr.threshold)
  }
  y.like <- statistics$male_presence >= presence.threshold &
    statistics$female_presence <= absence.threshold &
    significant(statistics$presence_fdr)
  w.like <- statistics$female_presence >= presence.threshold &
    statistics$male_presence <= absence.threshold &
    significant(statistics$presence_fdr)
  x.het <- statistics$heterozygosity_difference >=
    min.heterozygosity.difference &
    significant(statistics$heterozygosity_fdr)
  z.het <- statistics$heterozygosity_difference <=
    -min.heterozygosity.difference &
    significant(statistics$heterozygosity_fdr)
  x.cov <- statistics$coverage_ratio_female_male >=
    coverage.ratio.threshold & significant(statistics$coverage_fdr)
  z.cov <- statistics$coverage_ratio_female_male <=
    1 / coverage.ratio.threshold & significant(statistics$coverage_fdr)

  statistics$candidate_y_like <- y.like
  statistics$candidate_w_like <- w.like
  statistics$candidate_x_heterozygosity <- x.het
  statistics$candidate_z_heterozygosity <- z.het
  statistics$candidate_x_coverage <- x.cov
  statistics$candidate_z_coverage <- z.cov
  statistics$candidate_classes <- vapply(seq_len(nrow(statistics)), function(i) {
    paste(c(
      if (isTRUE(y.like[i])) "Y-like",
      if (isTRUE(w.like[i])) "W-like",
      if (isTRUE(x.het[i])) "X-like:heterozygosity",
      if (isTRUE(z.het[i])) "Z-like:heterozygosity",
      if (isTRUE(x.cov[i])) "X-like:coverage",
      if (isTRUE(z.cov[i])) "Z-like:coverage"
    ), collapse = ";")
  }, character(1))
  candidates <- statistics[nzchar(statistics$candidate_classes), , drop = FALSE]
  assignment.panel <- .sex_assignment_panel(
    statistics = statistics,
    depth.available = depth.available,
    coverage.threshold = coverage.threshold,
    n.female = n.female,
    n.male = n.male
  )

  sample.summary <- tibble::tibble(
    sex = c("F", "M", "U"),
    n = c(n.female, n.male, sum(sex == "U")),
    mean_positive_depth = c(
      mean(depth.scale[female]), mean(depth.scale[male]),
      if (any(sex == "U")) mean(depth.scale[sex == "U"]) else NA_real_
    )
  )
  metadata.audit <- .sex_metadata_audit(metadata, sex.column)
  confounded <- metadata.audit$variable[
    nzchar(metadata.audit$warning)
  ]
  if (length(confounded)) {
    .sex_warning(
      "Recorded sex is strongly associated with metadata: ",
      paste(confounded, collapse = ", "),
      ". Treat marker associations cautiously."
    )
  }
  readr::write_tsv(
    statistics, file.path(path.folder, "sex_marker_statistics.tsv")
  )
  readr::write_tsv(
    candidates, file.path(path.folder, "candidate_sex_markers.tsv")
  )
  readr::write_tsv(
    assignment.panel, file.path(path.folder, "sex_assignment_panel.tsv")
  )
  readr::write_tsv(
    sample.summary, file.path(path.folder, "sex_sample_summary.tsv")
  )
  readr::write_tsv(
    metadata.audit, file.path(path.folder, "sex_metadata_audit.tsv")
  )
  .sex_write_fasta(candidates, path.folder)
  plots <- .sex_make_plots(statistics, depth.available)
  if (save.plots) .sex_save_plots(plots, path.folder, plot.formats)

  SeqArray::seqFilterPop(gds)
  filter.pushed <- FALSE
  selection.after <- SeqArray::seqGetFilter(gds)
  restored <- identical(selection.before$sample.sel, selection.after$sample.sel) &&
    identical(selection.before$variant.sel, selection.after$variant.sel)
  if (!restored) {
    rlang::abort("The active GDS selection could not be restored safely.")
  }
  if (opened.here) {
    SeqArray::seqClose(gds)
    opened.here <- FALSE
  }
  if (verbose) {
    .sex_message("Candidate sex-linked markers: ", nrow(candidates), ".")
    .sex_message("Complete marker statistics: sex_marker_statistics.tsv")
    .sex_message("Candidate table: candidate_sex_markers.tsv")
    .sex_message("Results written to: ", path.folder)
  }
  result <- list(
    statistics = statistics,
    candidates = candidates,
    assignment_panel = assignment.panel,
    sample_summary = sample.summary,
    metadata_audit = metadata.audit,
    plots = plots,
    depth_available = depth.available,
    active_selection_restored = restored,
    path.folder = path.folder
  )
  class(result) <- "sexy_markers"
  result
}

#' @export
print.sexy_markers <- function(x, ...) {
  cat("Candidate sex-linked marker screen\n")
  cat("  Markers tested: ", nrow(x$statistics), "\n", sep = "")
  cat("  Candidates: ", nrow(x$candidates), "\n", sep = "")
  cat("  Depth method: ", if (isTRUE(x$depth_available)) {
    "available"
  } else "not available", "\n", sep = "")
  cat("  Results: ", x$path.folder, "\n", sep = "")
  invisible(x)
}

.sex_message <- function(...) {
  text <- paste0(..., collapse = "")
  message(paste(strwrap(text, width = 80L, exdent = 2L), collapse = "\n"))
}

.sex_warning <- function(...) {
  text <- paste0(..., collapse = "")
  warning(
    paste(strwrap(text, width = 80L, exdent = 2L), collapse = "\n"),
    call. = FALSE
  )
}

.sex_check_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    rlang::abort(paste0("`", name, "` must be one number from zero to one."))
  }
}

.sex_check_nonnegative <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    rlang::abort(paste0("`", name, "` must be one non-negative number."))
  }
}

.sex_check_count <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != as.integer(x) || x < minimum) {
    rlang::abort(paste0(
      "`", name, "` must be a whole number of at least ", minimum, "."
    ))
  }
  as.integer(x)
}

.sex_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    rlang::abort(paste0("`", name, "` must be TRUE or FALSE."))
  }
}

.sex_read_metadata <- function(strata, gds, sample.id, sex.column) {
  if (is.null(strata)) {
    metadata <- tryCatch(
      genometranslator::extract_individuals_metadata(
        gds = gds, whitelist = TRUE
      ), error = function(error) NULL
    )
    if (is.null(metadata)) {
      rlang::abort(
        "No GDS individual metadata found; supply `strata` explicitly."
      )
    }
  } else if (is.data.frame(strata)) {
    metadata <- strata
  } else if (is.character(strata) && length(strata) == 1L &&
             !is.na(strata) && file.exists(strata)) {
    metadata <- readr::read_tsv(
      strata, show_col_types = FALSE, progress = FALSE
    )
  } else {
    rlang::abort("`strata` must be NULL, a data frame, or a TSV file.")
  }
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  missing.columns <- setdiff(c("INDIVIDUALS", sex.column), names(metadata))
  if (length(missing.columns)) {
    rlang::abort(paste0(
      "Sample metadata is missing: ",
      paste(missing.columns, collapse = ", "), "."
    ))
  }
  metadata$INDIVIDUALS <- as.character(metadata$INDIVIDUALS)
  if (anyNA(metadata$INDIVIDUALS) || any(!nzchar(metadata$INDIVIDUALS)) ||
      anyDuplicated(metadata$INDIVIDUALS)) {
    rlang::abort("`INDIVIDUALS` must contain unique, non-missing sample IDs.")
  }
  metadata <- metadata[metadata$INDIVIDUALS %in% sample.id, , drop = FALSE]
  if (!nrow(metadata)) {
    rlang::abort("No metadata sample IDs match the active GDS samples.")
  }
  metadata
}

.sex_normalise <- function(x) {
  value <- toupper(trimws(as.character(x)))
  out <- rep("U", length(value))
  out[value %in% c("F", "FEMALE")] <- "F"
  out[value %in% c("M", "MALE")] <- "M"
  out
}

.sex_marker_metadata <- function(gds, variant.id) {
  fallback <- tibble::tibble(
    VARIANT_ID = variant.id,
    MARKERS = as.character(variant.id),
    CHROM = as.character(SeqArray::seqGetData(gds, "chromosome")),
    POS = suppressWarnings(as.numeric(SeqArray::seqGetData(gds, "position")))
  )
  metadata <- tryCatch(
    genometranslator::extract_markers_metadata(gds = gds, whitelist = TRUE),
    error = function(error) NULL
  )
  if (is.null(metadata) || !nrow(metadata)) return(fallback)
  metadata <- tibble::as_tibble(metadata)
  id.column <- intersect(c("VARIANT_ID", "M_SEQ", "variant.id"), names(metadata))
  if (!length(id.column)) return(fallback)
  matched <- metadata[match(variant.id, metadata[[id.column[[1L]]]]), , drop = FALSE]
  matched <- dplyr::rename_with(matched, toupper)
  if (!"VARIANT_ID" %in% names(matched)) matched$VARIANT_ID <- variant.id
  if (!"MARKERS" %in% names(matched)) matched$MARKERS <- as.character(variant.id)
  leading <- c("VARIANT_ID", "MARKERS", "CHROM", "POS")
  matched[, c(intersect(leading, names(matched)), setdiff(names(matched), leading))]
}

.sex_get_dosage <- function(gds, variant.id, sample.id) {
  SeqArray::seqFilterPush(gds)
  on.exit(SeqArray::seqFilterPop(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = variant.id, verbose = FALSE)
  dosage <- as.matrix(SeqArray::seqGetData(gds, "$dosage_alt"))
  if (nrow(dosage) == length(sample.id)) out <- dosage else
    if (ncol(dosage) == length(sample.id)) out <- t(dosage) else
      rlang::abort("GDS dosage dimensions do not match active sample IDs.")
  storage.mode(out) <- "double"
  out
}

.sex_get_depth <- function(gds, variant.id, sample.id) {
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
  if (is.null(depth)) {
    embedded <- tryCatch(
      .inversion_get_embedded_coverage(gds, variant.id, sample.id),
      error = function(error) NULL
    )
    if (!is.null(embedded)) depth <- embedded$depth
  }
  if (!is.null(depth)) storage.mode(depth) <- "double"
  depth
}

.sex_proportion_test <- function(success.f, total.f, success.m, total.m) {
  pooled <- (success.f + success.m) / (total.f + total.m)
  difference <- success.f / total.f - success.m / total.m
  se <- sqrt(pooled * (1 - pooled) * (1 / total.f + 1 / total.m))
  p <- 2 * stats::pnorm(-abs(difference / se))
  p[!is.finite(p)] <- 1
  p
}

.sex_welch_test <- function(female.matrix, male.matrix) {
  nf <- colSums(is.finite(female.matrix))
  nm <- colSums(is.finite(male.matrix))
  mf <- colMeans(female.matrix, na.rm = TRUE)
  mm <- colMeans(male.matrix, na.rm = TRUE)
  vf <- apply(female.matrix, 2L, stats::var, na.rm = TRUE)
  vm <- apply(male.matrix, 2L, stats::var, na.rm = TRUE)
  se2 <- vf / nf + vm / nm
  statistic <- (mf - mm) / sqrt(se2)
  df <- se2^2 / ((vf / nf)^2 / (nf - 1) + (vm / nm)^2 / (nm - 1))
  p <- 2 * stats::pt(-abs(statistic), df = df)
  p[!is.finite(p) | nf < 2L | nm < 2L] <- NA_real_
  p
}

.sex_chunk_statistics <- function(
    dosage, depth, depth.scale, female, male, coverage.threshold
) {
  called <- is.finite(dosage)
  present <- if (is.null(depth)) called else
    is.finite(depth) & depth >= coverage.threshold
  f.present <- colSums(present[female, , drop = FALSE])
  m.present <- colSums(present[male, , drop = FALSE])
  f.total <- rep(sum(female), ncol(dosage))
  m.total <- rep(sum(male), ncol(dosage))
  f.called <- colSums(called[female, , drop = FALSE])
  m.called <- colSums(called[male, , drop = FALSE])
  heterozygous <- called & abs(dosage - 1) < 1e-8
  f.het <- colSums(heterozygous[female, , drop = FALSE])
  m.het <- colSums(heterozygous[male, , drop = FALSE])
  f.het.rate <- f.het / f.called
  m.het.rate <- m.het / m.called
  if (!is.null(depth)) {
    normalized <- sweep(depth, 1L, depth.scale, "/")
    normalized[!is.finite(normalized)] <- NA_real_
    f.depth <- colMeans(normalized[female, , drop = FALSE], na.rm = TRUE)
    m.depth <- colMeans(normalized[male, , drop = FALSE], na.rm = TRUE)
    ratio <- f.depth / m.depth
    ratio[!is.finite(ratio)] <- NA_real_
    coverage.p <- .sex_welch_test(
      log2(1 + normalized[female, , drop = FALSE]),
      log2(1 + normalized[male, , drop = FALSE])
    )
  } else {
    f.depth <- m.depth <- ratio <- coverage.p <- rep(NA_real_, ncol(dosage))
  }
  tibble::tibble(
    female_presence = f.present / f.total,
    male_presence = m.present / m.total,
    presence_difference = female_presence - male_presence,
    presence_p = .sex_proportion_test(f.present, f.total, m.present, m.total),
    female_call_rate = f.called / f.total,
    male_call_rate = m.called / m.total,
    female_heterozygosity = f.het.rate,
    male_heterozygosity = m.het.rate,
    heterozygosity_difference = f.het.rate - m.het.rate,
    heterozygosity_p = .sex_proportion_test(f.het, f.called, m.het, m.called),
    female_normalized_depth = f.depth,
    male_normalized_depth = m.depth,
    coverage_ratio_female_male = ratio,
    coverage_p = coverage.p
  )
}

.sex_metadata_audit <- function(metadata, sex.column) {
  other <- setdiff(names(metadata), c("INDIVIDUALS", sex.column, "SEX_RADR"))
  if (!length(other)) return(tibble::tibble(
    variable = character(), known_n = integer(), levels = integer(),
    cramers_v = numeric(), warning = character()
  ))
  dplyr::bind_rows(lapply(other, function(variable) {
    value <- metadata[[variable]]
    keep <- metadata$SEX_RADR %in% c("F", "M") & !is.na(value)
    tab <- table(metadata$SEX_RADR[keep], as.character(value[keep]))
    n <- sum(tab)
    v <- NA_real_
    if (n > 0L && nrow(tab) > 1L && ncol(tab) > 1L) {
      chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
      v <- sqrt(unname(chi$statistic) /
        (n * min(nrow(tab) - 1L, ncol(tab) - 1L)))
    }
    tibble::tibble(
      variable = variable, known_n = n, levels = ncol(tab), cramers_v = v,
      warning = ifelse(is.finite(v) && v >= 0.5,
        "strong association with recorded sex", "")
    )
  }))
}

.sex_make_plots <- function(statistics, depth.available) {
  theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  make <- function(x, y, title, xlab, ylab) {
    ggplot2::ggplot(statistics, ggplot2::aes(x = {{ x }}, y = {{ y }})) +
      ggplot2::geom_point(alpha = 0.45, size = 1.2) +
      ggplot2::labs(x = xlab, y = ylab, title = title) + theme
  }
  presence <- make(
    presence_difference, -log10(pmax(presence_fdr, 1e-300)),
    "Sex-biased marker presence", "Female - male marker presence",
    "-log10 presence FDR"
  )
  heterozygosity <- make(
    heterozygosity_difference,
    -log10(pmax(heterozygosity_fdr, 1e-300)),
    "Sex-biased marker heterozygosity", "Female - male heterozygosity",
    "-log10 heterozygosity FDR"
  )
  coverage <- NULL
  if (depth.available) coverage <- make(
    log2(coverage_ratio_female_male),
    -log10(pmax(coverage_fdr, 1e-300)),
    "Sex-biased normalized depth", "log2 female / male normalized depth",
    "-log10 coverage FDR"
  ) + ggplot2::geom_point(na.rm = TRUE)
  list(presence = presence, heterozygosity = heterozygosity, coverage = coverage)
}

.sex_save_plots <- function(plots, path.folder, formats) {
  for (name in names(plots)) {
    if (is.null(plots[[name]])) next
    for (format in formats) ggplot2::ggsave(
      file.path(path.folder, paste0("sex_markers_", name, ".", format)),
      plots[[name]], width = 8, height = 5, units = "in", dpi = 300
    )
  }
}

.sex_write_fasta <- function(candidates, path.folder) {
  sequence.column <- intersect(
    c("SEQUENCE", "SEQ", "ALLELE_SEQUENCE"), names(candidates)
  )
  if (!nrow(candidates) || !length(sequence.column)) return(invisible(NULL))
  sequence <- as.character(candidates[[sequence.column[[1L]]]])
  keep <- !is.na(sequence) & nzchar(sequence)
  if (!any(keep)) return(invisible(NULL))
  labels <- gsub("[^A-Za-z0-9_.:-]", "_", candidates$MARKERS[keep])
  readr::write_lines(
    as.vector(rbind(paste0(">", labels), sequence[keep])),
    file.path(path.folder, "candidate_sex_markers.fasta")
  )
}

.sex_assignment_panel <- function(
    statistics, depth.available, coverage.threshold, n.female, n.male
) {
  keep <- statistics$candidate_y_like | statistics$candidate_w_like
  panel <- statistics[keep, , drop = FALSE]
  if (!nrow(panel)) {
    return(tibble::tibble(
      VARIANT_ID = statistics$VARIANT_ID[0],
      MARKERS = character(),
      ASSIGNMENT_DIRECTION = character(),
      EXPECTED_PRESENT_SEX = character(),
      EXPECTED_ABSENT_SEX = character(),
      PRESENCE_SOURCE = character(),
      COVERAGE_THRESHOLD = numeric(),
      PRESENCE_EFFECT = numeric(),
      PRESENCE_FDR = numeric(),
      DISCOVERY_FEMALES = integer(),
      DISCOVERY_MALES = integer()
    ))
  }
  direction <- ifelse(panel$candidate_y_like, "Y-like", "W-like")
  tibble::tibble(
    VARIANT_ID = panel$VARIANT_ID,
    MARKERS = panel$MARKERS,
    ASSIGNMENT_DIRECTION = direction,
    EXPECTED_PRESENT_SEX = ifelse(direction == "Y-like", "M", "F"),
    EXPECTED_ABSENT_SEX = ifelse(direction == "Y-like", "F", "M"),
    PRESENCE_SOURCE = if (depth.available) "read_depth" else "genotype_call",
    COVERAGE_THRESHOLD = if (depth.available) coverage.threshold else NA_real_,
    PRESENCE_EFFECT = panel$presence_difference,
    PRESENCE_FDR = panel$presence_fdr,
    DISCOVERY_FEMALES = n.female,
    DISCOVERY_MALES = n.male
  )
}
