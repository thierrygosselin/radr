#' Detect candidate paralogous markers
#'
#' `detect_paralogs()` implements a conservative, GDS-based version of the
#' HDplot diagnostic described by McKinney et al. (2017). It combines the
#' proportion of called samples that are heterozygous (`H_CALLED`) with the
#' deviation of pooled REF and ALT reads from a 1:1 ratio among heterozygous
#' genotypes (`D`). Candidate labels prioritize review; they are not proof of
#' paralogy.
#'
#' Excess heterozygosity or allele imbalance can also be caused by population
#' mixture, sex linkage, copy-number variation, contamination, allele dropout,
#' mapping bias, collapsed repeats, or genotype-calling errors. Confirm
#' candidates with population consistency, sequence or reference placement,
#' segregation, copy-number, and assembly evidence whenever possible.
#'
#' @section HD statistics:
#' The function reports both heterozygosity definitions from the published
#' implementation. `H_ALL` divides the number of heterozygotes by all samples
#' and is sensitive to missingness. `H_CALLED` divides by called samples and is
#' the default HDplot x-axis. `D` is the signed binomial z-score for departure
#' from equal REF and ALT read counts among heterozygotes; `ABS_D` is used for
#' symmetric candidate screening. MAF is calculated from called diploid
#' dosages in the active sample set.
#'
#' @section Candidate categories:
#' Each pooled or stratum-specific marker receives one status:
#' \itemize{
#'   \item `INSUFFICIENT_DATA`: one or more evidence requirements failed;
#'   \item `SINGLETON_LIKE`: neither screening threshold was crossed;
#'   \item `EXCESS_HET_ONLY`: only the heterozygosity threshold was crossed;
#'   \item `READ_IMBALANCE_ONLY`: only the read-deviation threshold was crossed;
#'   \item `PARALOG_CANDIDATE`: both thresholds were crossed.
#' }
#' These empirical thresholds are screening controls, not universal paralog
#' boundaries. Inspect the full distribution, MAF, missingness, depth, strata,
#' and technical metadata before removing markers.
#'
#' @section Strata:
#' With `by.strata = TRUE`, diagnostics are calculated for the pooled active
#' samples and independently for each value of `group.column`. A signal that
#' occurs in only one population or technical group deserves additional review:
#' it can reflect real lineage-specific duplication, but also population
#' structure, batch effects, or small group size. Pooled and stratified results
#' are therefore both retained.
#'
#' @param data A GDS filepath or an open `SeqVarGDSClass` object.
#' @param strata Optional data frame or tab-delimited metadata file containing
#'   `INDIVIDUALS` and `group.column`. If `NULL`, metadata are read from
#'   the GDS.
#' @param group.column Metadata column used for stratified diagnostics.
#' @param by.strata Calculate diagnostics independently within each group in
#'   addition to the pooled analysis.
#' @param min.call.rate Minimum called-genotype proportion required.
#' @param min.samples.per.group Minimum active samples required in a diagnostic
#'   group.
#' @param min.heterozygotes Minimum heterozygous genotypes with usable
#'   allele-depth values required for read-ratio inference.
#' @param min.heterozygote.depth Minimum pooled REF plus ALT depth across
#'   heterozygotes.
#' @param heterozygosity.threshold Minimum `H_CALLED` for excess
#'   heterozygosity evidence.
#' @param deviation.threshold Minimum `ABS_D` for read-imbalance evidence.
#' @param chunk.size Number of active variants read from the GDS at a time.
#' @param save.plots Save HD diagnostic plots.
#' @param plot.formats One or both of `"png"` and `"pdf"`.
#' @param verbose Display progress and summary messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#'
#' @return A `detect_paralogs` object containing:
#' \itemize{
#'   \item `statistics`: pooled marker statistics and cross-stratum evidence;
#'   \item `strata.statistics`: pooled and group-specific HD statistics;
#'   \item `candidates`: markers flagged in the pooled data or any stratum;
#'   \item `thresholds`: screening and minimum-evidence settings;
#'   \item `plots`: generated ggplot objects;
#'   \item `path.folder` and `output.files`;
#'   \item `active.selection.restored`: filter-restoration check.
#' }
#' The function is read-only and does not change the persistent GDS filters.
#'
#' @references
#' McKinney GJ, Waples RK, Seeb LW, Seeb JE (2017). Paralogs are revealed by
#' proportion of heterozygotes and deviations in read ratios in
#' genotyping-by-sequencing data from natural populations. *Molecular Ecology
#' Resources*, 17, 656-669. https://doi.org/10.1111/1755-0998.12613
#'
#' @seealso
#' \code{\link{detect_het_outliers}}, \code{\link{detect_mixed_genomes}},
#' \code{\link{sexy_markers}}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' paralogs <- radr::detect_paralogs(
#'   data = "study.gds",
#'   strata = "sample_metadata.tsv",
#'   group.column = "STRATA"
#' )
#' paralogs$candidates
#' }
detect_paralogs <- function(
    data,
    strata = NULL,
    group.column = "STRATA",
    by.strata = TRUE,
    min.call.rate = 0.8,
    min.samples.per.group = 10L,
    min.heterozygotes = 5L,
    min.heterozygote.depth = 50,
    heterozygosity.threshold = 0.5,
    deviation.threshold = 4,
    chunk.size = 2000L,
    save.plots = TRUE,
    plot.formats = c("png", "pdf"),
    verbose = TRUE,
    ...
) {
  force(data)
  .start <- tgbase::startup(
    package = "radr", f.name = "detect_paralogs", verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)

  dots <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  unknown <- setdiff(names(dots), c("path.folder", "internal"))
  if (length(unknown)) {
    rlang::abort(paste0(
      "Unknown argument(s): ", paste(unknown, collapse = ", "), "."
    ))
  }
  parent.folder <- dots$path.folder %||% getwd()
  internal <- isTRUE(dots$internal)

  .paralog_check_probability(min.call.rate, "min.call.rate")
  .paralog_check_probability(
    heterozygosity.threshold, "heterozygosity.threshold"
  )
  min.samples.per.group <- .paralog_check_count(
    min.samples.per.group, "min.samples.per.group", 2L
  )
  min.heterozygotes <- .paralog_check_count(
    min.heterozygotes, "min.heterozygotes", 1L
  )
  chunk.size <- .paralog_check_count(chunk.size, "chunk.size", 1L)
  .paralog_check_nonnegative(
    min.heterozygote.depth, "min.heterozygote.depth"
  )
  .paralog_check_nonnegative(deviation.threshold, "deviation.threshold")
  .paralog_check_flag(by.strata, "by.strata")
  .paralog_check_flag(save.plots, "save.plots")
  .paralog_check_flag(verbose, "verbose")
  if (!is.character(group.column) || length(group.column) != 1L ||
      is.na(group.column) || !nzchar(group.column)) {
    rlang::abort("`group.column` must be one non-empty column name.")
  }
  plot.formats <- unique(tolower(as.character(plot.formats)))
  if (!length(plot.formats) || any(!plot.formats %in% c("png", "pdf"))) {
    rlang::abort("`plot.formats` must contain `png`, `pdf`, or both.")
  }

  dir.create(parent.folder, recursive = TRUE, showWarnings = FALSE)
  parent.folder <- normalizePath(parent.folder, mustWork = TRUE)
  path.folder <- if (internal) {
    parent.folder
  } else {
    radr_folder(
      rad.folder = paste0("detect_paralogs_", file.date),
      path.folder = parent.folder,
      prefix.int = TRUE
    )
  }
  dir.create(path.folder, recursive = TRUE, showWarnings = FALSE)
  if (verbose && !internal) {
    .paralog_message("Folder created: ", basename(path.folder))
  }

  arguments <- tibble::tibble(
    argument = c(
      "group.column", "by.strata", "min.call.rate",
      "min.samples.per.group", "min.heterozygotes",
      "min.heterozygote.depth", "heterozygosity.threshold",
      "deviation.threshold", "chunk.size", "save.plots", "plot.formats"
    ),
    value = as.character(c(
      group.column, by.strata, min.call.rate, min.samples.per.group,
      min.heterozygotes, min.heterozygote.depth,
      heterozygosity.threshold, deviation.threshold, chunk.size,
      save.plots, paste(plot.formats, collapse = ",")
    ))
  )
  args.file <- file.path(
    path.folder, paste0("radr_detect_paralogs_args_", file.date, ".tsv")
  )
  readr::write_tsv(arguments, args.file)

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
  if (!genometranslator::detect_biallelic_markers(gds)) {
    rlang::abort("Paralog diagnostics require active biallelic markers.")
  }

  metadata <- .paralog_read_metadata(
    strata, gds, sample.id, group.column, by.strata
  )
  sample.id <- sample.id[sample.id %in% metadata$INDIVIDUALS]
  metadata <- metadata[
    match(sample.id, metadata$INDIVIDUALS), , drop = FALSE
  ]
  SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
  groups <- if (by.strata) as.character(metadata[[group.column]]) else
    rep("POOLED", length(sample.id))
  if (anyNA(groups) || any(!nzchar(trimws(groups)))) {
    rlang::abort(paste0(
      "`", group.column, "` contains missing or empty values."
    ))
  }
  if (any(groups == "OVERALL")) {
    rlang::abort("`OVERALL` is reserved for pooled diagnostics.")
  }
  if (by.strata && identical(unique(groups), "POOLED")) {
    by.strata <- FALSE
  }

  marker.metadata <- .sex_marker_metadata(gds, variant.id)
  chunks <- split(
    seq_along(variant.id), ceiling(seq_along(variant.id) / chunk.size)
  )
  if (verbose) {
    .paralog_message(
      "Active GDS selection: ", length(variant.id), " marker(s), ",
      length(sample.id), " sample(s)."
    )
    .paralog_message(
      "Calculating HD statistics in ", length(chunks), " chunk(s)."
    )
  }

  pieces <- vector("list", length(chunks))
  depth.sources <- character()
  for (i in seq_along(chunks)) {
    index <- chunks[[i]]
    dosage <- .sex_get_dosage(gds, variant.id[index], sample.id)
    allele.depth <- .paralog_get_allele_depth(
      gds, variant.id[index], sample.id
    )
    if (is.null(allele.depth$ref) || is.null(allele.depth$alt)) {
      rlang::abort(paste0(
        "Allele-specific depth is required. Import AD or DArT REF/ALT ",
        "counts into the GDS."
      ))
    }
    depth.sources <- c(depth.sources, allele.depth$source)
    pieces[[i]] <- .paralog_chunk_statistics(
      variant.id = variant.id[index],
      dosage = dosage,
      ref.depth = allele.depth$ref,
      alt.depth = allele.depth$alt,
      groups = groups,
      by.strata = by.strata,
      min.call.rate = min.call.rate,
      min.samples.per.group = min.samples.per.group,
      min.heterozygotes = min.heterozygotes,
      min.heterozygote.depth = min.heterozygote.depth,
      heterozygosity.threshold = heterozygosity.threshold,
      deviation.threshold = deviation.threshold
    )
    if (verbose && (i == length(chunks) || i %% 10L == 0L)) {
      .paralog_message(
        "Processed ", i, " of ", length(chunks), " chunk(s)."
      )
    }
  }
  strata.statistics <- dplyr::bind_rows(pieces)
  pooled <- strata.statistics[
    strata.statistics$GROUP == "OVERALL", , drop = FALSE
  ]
  group.statistics <- strata.statistics[
    strata.statistics$GROUP != "OVERALL", , drop = FALSE
  ]
  consistency <- if (nrow(group.statistics)) {
    group.statistics |>
      dplyr::group_by(.data$VARIANT_ID) |>
      dplyr::summarise(
        ELIGIBLE_STRATA = sum(.data$ELIGIBLE),
        CANDIDATE_STRATA = sum(
          .data$STATUS == "PARALOG_CANDIDATE"
        ),
        EXCESS_HET_STRATA = sum(.data$EXCESS_HETEROZYGOSITY),
        READ_IMBALANCE_STRATA = sum(.data$READ_IMBALANCE),
        .groups = "drop"
      )
  } else {
    tibble::tibble(
      VARIANT_ID = variant.id,
      ELIGIBLE_STRATA = 0L,
      CANDIDATE_STRATA = 0L,
      EXCESS_HET_STRATA = 0L,
      READ_IMBALANCE_STRATA = 0L
    )
  }
  statistics <- marker.metadata |>
    dplyr::left_join(pooled, by = "VARIANT_ID") |>
    dplyr::left_join(consistency, by = "VARIANT_ID") |>
    dplyr::mutate(
      ANY_PARALOG_CANDIDATE =
        .data$STATUS == "PARALOG_CANDIDATE" |
        .data$CANDIDATE_STRATA > 0L,
      STRATUM_SPECIFIC_CANDIDATE =
        .data$STATUS != "PARALOG_CANDIDATE" &
        .data$CANDIDATE_STRATA > 0L,
      STRATUM_INCONSISTENT =
        .data$CANDIDATE_STRATA > 0L &
        .data$CANDIDATE_STRATA < .data$ELIGIBLE_STRATA
    )
  candidates <- statistics[
    statistics$ANY_PARALOG_CANDIDATE %in% TRUE, , drop = FALSE
  ]

  thresholds <- tibble::tibble(
    parameter = c(
      "min.call.rate", "min.samples.per.group", "min.heterozygotes",
      "min.heterozygote.depth", "heterozygosity.threshold",
      "deviation.threshold"
    ),
    value = c(
      min.call.rate, min.samples.per.group, min.heterozygotes,
      min.heterozygote.depth, heterozygosity.threshold,
      deviation.threshold
    )
  )
  plots <- .paralog_make_plots(
    strata.statistics, heterozygosity.threshold, deviation.threshold
  )
  output.files <- c(
    args.file,
    file.path(path.folder, "paralog_marker_statistics.tsv"),
    file.path(path.folder, "paralog_stratum_statistics.tsv"),
    file.path(path.folder, "candidate_paralog_markers.tsv"),
    file.path(path.folder, "paralog_thresholds.tsv")
  )
  readr::write_tsv(statistics, output.files[[2L]])
  readr::write_tsv(strata.statistics, output.files[[3L]])
  readr::write_tsv(candidates, output.files[[4L]])
  readr::write_tsv(thresholds, output.files[[5L]])
  if (save.plots) {
    plot.files <- .paralog_save_plots(
      plots, path.folder, plot.formats
    )
    output.files <- c(output.files, plot.files)
  }

  SeqArray::seqFilterPop(gds)
  filter.pushed <- FALSE
  selection.after <- SeqArray::seqGetFilter(gds)
  restored <- identical(
    selection.before$sample.sel, selection.after$sample.sel
  ) && identical(
    selection.before$variant.sel, selection.after$variant.sel
  )
  if (!restored) {
    rlang::abort("The active GDS selection was not restored.")
  }

  if (verbose) {
    .paralog_message(
      "Candidate markers: ", nrow(candidates), " of ",
      nrow(statistics), "."
    )
    .paralog_message(
      "Allele-depth source: ", paste(unique(depth.sources), collapse = ", "),
      "."
    )
    .paralog_message("Paralog diagnostics written to: ", path.folder)
  }

  out <- list(
    statistics = tibble::as_tibble(statistics),
    strata.statistics = tibble::as_tibble(strata.statistics),
    candidates = tibble::as_tibble(candidates),
    thresholds = thresholds,
    plots = plots,
    allele.depth.source = unique(depth.sources),
    path.folder = path.folder,
    output.files = tibble::tibble(files = output.files),
    active.selection.restored = restored
  )
  class(out) <- c("detect_paralogs", class(out))
  out
}

.paralog_chunk_statistics <- function(
    variant.id, dosage, ref.depth, alt.depth, groups, by.strata,
    min.call.rate, min.samples.per.group, min.heterozygotes,
    min.heterozygote.depth, heterozygosity.threshold, deviation.threshold
) {
  if (!identical(dim(dosage), dim(ref.depth)) ||
      !identical(dim(dosage), dim(alt.depth))) {
    rlang::abort("Dosage and allele-depth matrices have different dimensions.")
  }
  group.index <- list(OVERALL = seq_len(nrow(dosage)))
  if (by.strata) {
    group.index <- c(
      group.index,
      split(seq_along(groups), groups, drop = TRUE)
    )
  }
  dplyr::bind_rows(lapply(names(group.index), function(group.name) {
    index <- group.index[[group.name]]
    x <- dosage[index, , drop = FALSE]
    ref <- ref.depth[index, , drop = FALSE]
    alt <- alt.depth[index, , drop = FALSE]
    called <- is.finite(x)
    heterozygous <- called & abs(x - 1) < 1e-8
    depth.complete <- heterozygous & is.finite(ref) & is.finite(alt) &
      ref >= 0 & alt >= 0
    number.called <- colSums(called)
    number.het <- colSums(heterozygous)
    number.het.depth <- colSums(depth.complete)
    n.samples <- nrow(x)
    allele.frequency <- colSums(x, na.rm = TRUE) / (2 * number.called)
    maf <- pmin(allele.frequency, 1 - allele.frequency)
    ref.sum <- colSums(replace(ref, !depth.complete, 0), na.rm = TRUE)
    alt.sum <- colSums(replace(alt, !depth.complete, 0), na.rm = TRUE)
    total.depth <- ref.sum + alt.sum
    ratio <- ref.sum / total.depth
    d <- (ref.sum - total.depth / 2) / sqrt(total.depth * 0.25)
    ratio[!is.finite(ratio)] <- NA_real_
    d[!is.finite(d)] <- NA_real_
    call.rate <- number.called / n.samples
    h.called <- number.het / number.called
    h.called[!is.finite(h.called)] <- NA_real_
    h.all <- number.het / n.samples
    eligible <- n.samples >= min.samples.per.group &
      call.rate >= min.call.rate &
      number.het.depth >= min.heterozygotes &
      total.depth >= min.heterozygote.depth &
      is.finite(h.called) & is.finite(d)
    excess <- eligible & h.called >= heterozygosity.threshold
    imbalance <- eligible & abs(d) >= deviation.threshold
    status <- rep("INSUFFICIENT_DATA", ncol(x))
    status[eligible] <- "SINGLETON_LIKE"
    status[eligible & excess & !imbalance] <- "EXCESS_HET_ONLY"
    status[eligible & !excess & imbalance] <- "READ_IMBALANCE_ONLY"
    status[eligible & excess & imbalance] <- "PARALOG_CANDIDATE"
    reasons <- vapply(seq_len(ncol(x)), function(j) {
      if (eligible[j]) return("")
      paste(c(
        if (n.samples < min.samples.per.group) "too_few_samples",
        if (!is.finite(call.rate[j]) || call.rate[j] < min.call.rate)
          "low_call_rate",
        if (number.het.depth[j] < min.heterozygotes)
          "too_few_heterozygotes_with_depth",
        if (!is.finite(total.depth[j]) ||
            total.depth[j] < min.heterozygote.depth)
          "low_heterozygote_depth",
        if (!is.finite(h.called[j])) "undefined_heterozygosity",
        if (!is.finite(d[j])) "undefined_read_deviation"
      ), collapse = ";")
    }, character(1))
    tibble::tibble(
      VARIANT_ID = variant.id,
      GROUP = group.name,
      N_SAMPLES = n.samples,
      NUMBER_CALLED = number.called,
      NUMBER_MISSING = n.samples - number.called,
      CALL_RATE = call.rate,
      MISSING_PROP = 1 - call.rate,
      NUMBER_HET = number.het,
      NUMBER_HET_WITH_DEPTH = number.het.depth,
      H_ALL = unname(h.all),
      H_CALLED = unname(h.called),
      ALT_FREQUENCY = allele.frequency,
      MAF = maf,
      DEPTH_REF_HET = ref.sum,
      DEPTH_ALT_HET = alt.sum,
      HETEROZYGOTE_TOTAL_DEPTH = total.depth,
      REF_RATIO_HET = ratio,
      D = d,
      ABS_D = abs(d),
      ELIGIBLE = eligible,
      EXCESS_HETEROZYGOSITY = excess,
      READ_IMBALANCE = imbalance,
      STATUS = status,
      INSUFFICIENT_REASON = reasons
    )
  }))
}

.paralog_get_allele_depth <- function(gds, variant.id, sample.id) {
  SeqArray::seqFilterPush(gds)
  on.exit(SeqArray::seqFilterPop(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = variant.id, verbose = FALSE)
  normalise <- function(x) {
    if (!is.matrix(x)) return(NULL)
    if (nrow(x) == length(sample.id)) return(x)
    if (ncol(x) == length(sample.id)) return(t(x))
    NULL
  }
  ad <- tryCatch(
    SeqArray::seqGetData(gds, "annotation/format/AD"),
    error = function(error) NULL
  )
  if (is.list(ad) && all(c("length", "data") %in% names(ad)) &&
      all(ad$length == 2L)) {
    ad.data <- normalise(ad$data)
    if (!is.null(ad.data) &&
        ncol(ad.data) == 2L * length(variant.id)) {
      ref <- ad.data[, seq.int(1L, ncol(ad.data), 2L), drop = FALSE]
      alt <- ad.data[, seq.int(2L, ncol(ad.data), 2L), drop = FALSE]
      storage.mode(ref) <- storage.mode(alt) <- "double"
      return(list(
        ref = ref, alt = alt, source = "annotation/format/AD"
      ))
    }
  }
  embedded <- tryCatch(
    .inversion_get_embedded_coverage(gds, variant.id, sample.id),
    error = function(error) NULL
  )
  if (!is.null(embedded) && !is.null(embedded$depth) &&
      !is.null(embedded$allele.balance)) {
    ref <- embedded$depth * (1 - embedded$allele.balance)
    alt <- embedded$depth * embedded$allele.balance
    storage.mode(ref) <- storage.mode(alt) <- "double"
    return(list(
      ref = ref, alt = alt,
      source = paste(embedded$source, collapse = "+")
    ))
  }
  list(ref = NULL, alt = NULL, source = character())
}

.paralog_read_metadata <- function(
    strata, gds, sample.id, group.column, by.strata
) {
  if (is.null(strata)) {
    metadata <- tryCatch(
      genometranslator::extract_individuals_metadata(
        gds = gds, whitelist = TRUE
      ),
      error = function(error) NULL
    )
    if (is.null(metadata)) {
      metadata <- data.frame(INDIVIDUALS = sample.id)
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
  if (!"INDIVIDUALS" %in% names(metadata)) {
    rlang::abort("Sample metadata must contain `INDIVIDUALS`.")
  }
  metadata$INDIVIDUALS <- as.character(metadata$INDIVIDUALS)
  if (anyNA(metadata$INDIVIDUALS) || any(!nzchar(metadata$INDIVIDUALS)) ||
      anyDuplicated(metadata$INDIVIDUALS)) {
    rlang::abort("`INDIVIDUALS` must contain unique, non-missing sample IDs.")
  }
  metadata <- metadata[
    metadata$INDIVIDUALS %in% sample.id, , drop = FALSE
  ]
  if (!nrow(metadata)) {
    rlang::abort("No metadata sample IDs match the active GDS samples.")
  }
  if (by.strata && !group.column %in% names(metadata)) {
    warning(
      paste0(
        "`", group.column,
        "` is unavailable; only pooled diagnostics will be calculated."
      ),
      call. = FALSE
    )
    metadata[[group.column]] <- "POOLED"
  }
  if (!group.column %in% names(metadata)) {
    metadata[[group.column]] <- "POOLED"
  }
  metadata
}

.paralog_make_plots <- function(
    statistics, heterozygosity.threshold, deviation.threshold
) {
  theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )
  hd <- ggplot2::ggplot(
    statistics,
    ggplot2::aes(
      x = .data$H_CALLED, y = .data$ABS_D,
      colour = .data$STATUS, size = .data$MISSING_PROP
    )
  ) +
    ggplot2::geom_vline(
      xintercept = heterozygosity.threshold,
      linetype = "dashed", colour = "grey45"
    ) +
    ggplot2::geom_hline(
      yintercept = deviation.threshold,
      linetype = "dashed", colour = "grey45"
    ) +
    ggplot2::geom_point(alpha = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~ GROUP) +
    ggplot2::scale_size_area(
      name = "Missing proportion", max_size = 4
    ) +
    ggplot2::labs(
      x = "Heterozygotes among called samples (H)",
      y = "Absolute read-ratio deviation (|D|)",
      colour = "Diagnostic status",
      title = "Candidate paralog HD diagnostic"
    ) +
    theme
  ratio <- ggplot2::ggplot(
    statistics,
    ggplot2::aes(
      x = .data$H_CALLED, y = .data$REF_RATIO_HET,
      colour = .data$STATUS, size = .data$MISSING_PROP
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0.5, linetype = "dashed", colour = "grey45"
    ) +
    ggplot2::geom_point(alpha = 0.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~ GROUP) +
    ggplot2::scale_size_area(
      name = "Missing proportion", max_size = 4
    ) +
    ggplot2::labs(
      x = "Heterozygotes among called samples (H)",
      y = "REF read proportion among heterozygotes",
      colour = "Diagnostic status",
      title = "Heterozygosity and heterozygote allele balance"
    ) +
    theme
  list(hd = hd, ratio = ratio)
}

.paralog_save_plots <- function(plots, path.folder, formats) {
  unlist(lapply(names(plots), function(name) {
    vapply(formats, function(format) {
      file <- file.path(
        path.folder, paste0("paralog_", name, "_plot.", format)
      )
      arguments <- list(
        filename = file, plot = plots[[name]], width = 24, height = 16,
        units = "cm", dpi = 300, limitsize = FALSE
      )
      if (format == "pdf") arguments$useDingbats <- FALSE
      do.call(ggplot2::ggsave, arguments)
      file
    }, character(1))
  }), use.names = FALSE)
}

.paralog_check_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0 || x > 1) {
    rlang::abort(paste0("`", name, "` must be between zero and one."))
  }
}

.paralog_check_nonnegative <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0) {
    rlang::abort(paste0("`", name, "` must be a non-negative number."))
  }
}

.paralog_check_count <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < minimum || x != as.integer(x)) {
    rlang::abort(paste0(
      "`", name, "` must be a whole number of at least ", minimum, "."
    ))
  }
  as.integer(x)
}

.paralog_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    rlang::abort(paste0("`", name, "` must be TRUE or FALSE."))
  }
}

.paralog_message <- function(...) {
  text <- paste0(...)
  wrapped <- strwrap(text, width = 80L, exdent = 2L)
  message(paste(wrapped, collapse = "\n"))
}
