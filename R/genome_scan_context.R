#' Summarise genomic context around genome-scan signals
#'
#' Build a window-level context table for interpreting regional genome-scan
#' peaks. The table places marker density, call rate, heterozygosity, minor
#' allele frequency, depth when available, local LD, candidate inversion or
#' structural-region annotations, and user-supplied scan statistics beside one
#' another. It does not decide whether a region is under selection.
#'
#' @description A regional peak is easier to interpret when the analyst can ask
#' whether it overlaps low recombination, unusual marker density, increased
#' missingness, a sequencing or mapping batch, a candidate inversion, or a
#' signal supported by a method based on a different genomic signature.
#'
#' Use this function after basic quality control but before interpreting genome
#' scans. When candidate inversion regions exist, compare at least three views:
#' the complete genome, a collinear sensitivity analysis excluding candidate
#' regions, and an inversion-specific analysis. Candidate regions are annotated
#' but never silently excluded.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param window.snps Number of SNPs per fixed-SNP window. Default:
#'   \code{window.snps = 250}.
#' @param step.snps Number of SNPs between window starts. Default:
#'   \code{step.snps = window.snps}.
#' @param window.bp Optional physical window size in base pairs. Supplying it
#'   replaces fixed-SNP windows. Default: \code{window.bp = NULL}.
#' @param step.bp Physical distance between window starts. Default:
#'   \code{step.bp = window.bp}.
#' @param inversion.regions Optional table, or the result from
#'   [detect_inversions()]. A table requires `chromosome`, `start`, and `end`;
#'   `candidate_id` and `candidate_class` are retained when present. Default:
#'   \code{inversion.regions = NULL}.
#' @param scan.statistics Optional marker-level table with chromosome, position,
#'   and one or more scan statistics. Column names are matched without regard to
#'   case. Numeric statistics are summarised by window mean, median, and maximum.
#'   Default: \code{scan.statistics = NULL}.
#' @param ld.max.snps Maximum evenly spaced SNPs used for LD in each window.
#'   Default: \code{ld.max.snps = 250}.
#' @param filename Output filename stem. Default:
#'   \code{filename = "genome_scan_context"}.
#' @param verbose Logical. Display progress messages. Default:
#'   \code{verbose = TRUE}.
#' @param ... Standard `radr` workflow arguments, including `path.folder`.
#'
#' @return A list containing the window context table, normalised region table,
#' plots, data-source information, and output paths.
#'
#' @details The context table is descriptive. Overlap with an inversion-like
#' haploblock, centromere, low-recombination region, assembly gap, or technical
#' anomaly changes the interpretation of a peak but does not by itself validate
#' or invalidate selection. Repeat analyses after excluding putative
#' heterokaryotypes when arrangement calls are available, and seek support from
#' methods based on different summaries.
#'
#' @references Booker TR, Yeaman S, Whitlock MC (2020). Variation in
#' recombination rate affects detection of outliers in genome scans under
#' neutrality. Molecular Ecology, 29, 4274-4279.
#'
#' Faria R, Johannesson K, Butlin RK, Westram AM (2019). Evolving inversions.
#' Trends in Ecology & Evolution, 34, 239-248.
#'
#' @export
genome_scan_context <- function(
    data,
    window.snps = 250L,
    step.snps = window.snps,
    window.bp = NULL,
    step.bp = window.bp,
    inversion.regions = NULL,
    scan.statistics = NULL,
    ld.max.snps = 250L,
    filename = "genome_scan_context",
    verbose = TRUE,
    ...
) {
  .start <- tgbase::startup(
    package = "radr", f.name = "genome_scan_context", verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)
  dots <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  path.parent <- if (is.null(dots$path.folder)) getwd() else dots$path.folder
  path.folder <- tgbase::generate_folder(
    folder = "genome_scan_context", path.folder = path.parent,
    internal = isTRUE(dots$internal), file.date = .start$file.date,
    prefix.int = TRUE, verbose = verbose
  )
  if (missing(data)) rlang::abort("Argument `data` is required.")
  window.snps <- as.integer(window.snps)
  step.snps <- as.integer(step.snps)
  ld.max.snps <- as.integer(ld.max.snps)
  if (any(!is.finite(c(window.snps, step.snps, ld.max.snps))) ||
      any(c(window.snps, step.snps, ld.max.snps) < 2L)) {
    rlang::abort("Window and LD marker counts must be whole numbers of at least 2.")
  }

  opened.here <- FALSE
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.character(data) && length(data) == 1L && file.exists(data)) {
    gds <- SeqArray::seqOpen(data)
    opened.here <- TRUE
  } else {
    rlang::abort("`data` must be a GDS filepath or open SeqVarGDSClass object.")
  }
  on.exit(if (opened.here) try(SeqArray::seqClose(gds), silent = TRUE), add = TRUE)
  SeqArray::seqFilterPush(gds)
  on.exit(try(SeqArray::seqFilterPop(gds), silent = TRUE), add = TRUE)

  sample.id <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  markers <- tibble::tibble(
    variant_id = SeqArray::seqGetData(gds, "variant.id"),
    chromosome = as.character(SeqArray::seqGetData(gds, "chromosome")),
    position = as.numeric(SeqArray::seqGetData(gds, "position"))
  ) |>
    dplyr::arrange(
      .inversion_chromosome_order(.data$chromosome),
      .data$position, .data$variant_id
    )
  windows <- .inversion_make_windows(
    marker.table = markers,
    window.snps = window.snps,
    step.snps = step.snps,
    window.bp = window.bp,
    step.bp = step.bp
  )
  if (!length(windows)) rlang::abort("No complete genomic windows were produced.")

  if (verbose) message("Summarising ", length(windows), " genomic windows...")
  context <- purrr::map2_dfr(windows, seq_along(windows), function(w, i) {
    dosage <- .inversion_get_dosage(gds, w$variant_id, sample.id)
    mean.depth <- .genome_context_mean_depth(gds, w$variant_id)
    called <- !is.na(dosage)
    alt.count <- colSums(dosage, na.rm = TRUE)
    chromosomes.called <- 2 * colSums(called)
    af <- alt.count / chromosomes.called
    maf <- pmin(af, 1 - af)
    select <- unique(round(seq(
      1, ncol(dosage), length.out = min(ncol(dosage), ld.max.snps)
    )))
    tibble::tibble(
      window_id = i,
      chromosome = w$chromosome,
      start = w$start,
      end = w$end,
      midpoint = (w$start + w$end) / 2,
      width_bp = w$end - w$start + 1,
      n_markers = ncol(dosage),
      markers_per_mb = ncol(dosage) / ((w$end - w$start + 1) / 1e6),
      mean_call_rate = mean(called),
      mean_missing_rate = 1 - mean(called),
      mean_heterozygosity = mean(dosage == 1, na.rm = TRUE),
      mean_depth = mean.depth,
      mean_maf = mean(maf, na.rm = TRUE),
      median_maf = stats::median(maf, na.rm = TRUE),
      mean_ld_r2 = .inversion_mean_ld(dosage[, select, drop = FALSE])
    )
  })

  regions <- .genome_context_regions(inversion.regions)
  if (nrow(regions)) {
    overlap <- purrr::pmap(
      context[c("chromosome", "start", "end")],
      function(chromosome, start, end) which(
        regions$chromosome == chromosome & regions$start <= end & regions$end >= start
      )
    )
    context$region_overlap <- purrr::map_chr(overlap, function(i) {
      if (!length(i)) "none" else paste(regions$candidate_id[i], collapse = ";")
    })
    context$region_class <- purrr::map_chr(overlap, function(i) {
      if (!length(i)) "none" else paste(unique(regions$candidate_class[i]), collapse = ";")
    })
  } else {
    context$region_overlap <- "none"
    context$region_class <- "none"
  }
  context <- .genome_context_join_statistics(context, scan.statistics)

  plot <- ggplot2::ggplot(
    context,
    ggplot2::aes(
      x = .data$midpoint / 1e6, y = .data$mean_missing_rate,
      colour = .data$region_overlap != "none"
    )
  ) +
    ggplot2::geom_line(ggplot2::aes(group = chromosome), colour = "grey70") +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(~ chromosome, scales = "free_x") +
    ggplot2::scale_colour_manual(values = c("FALSE" = "grey35", "TRUE" = "#B2182B")) +
    ggplot2::labs(
      x = "Genomic position (Mb)", y = "Mean missing-genotype rate",
      colour = "Candidate region"
    ) +
    ggplot2::theme_bw()

  table.path <- file.path(path.folder, paste0(filename, ".tsv"))
  plot.path <- file.path(path.folder, paste0(filename, "_missingness.png"))
  readr::write_tsv(context, table.path, na = "NA")
  ggplot2::ggsave(plot.path, plot, width = 10, height = 7, dpi = 300)
  if (verbose) message("Genome-scan context written: ", basename(table.path))
  source <- tryCatch(
    genometranslator::extract_data_source(gds), error = function(e) NA_character_
  )
  structure(list(
    context = context,
    regions = regions,
    plot = plot,
    data.source = source,
    path.folder = path.folder,
    output.files = c(table.path, plot.path)
  ), class = c("genome_scan_context", "list"))
}

.genome_context_mean_depth <- function(gds, variant.id) {
  SeqArray::seqFilterPush(gds)
  on.exit(SeqArray::seqFilterPop(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = variant.id, verbose = FALSE)
  depth <- tryCatch(
    SeqArray::seqGetData(gds, "annotation/format/DP"),
    error = function(e) NULL
  )
  if (is.null(depth)) return(NA_real_)
  depth <- suppressWarnings(as.numeric(depth))
  if (!any(is.finite(depth))) NA_real_ else mean(depth, na.rm = TRUE)
}

.genome_context_regions <- function(x) {
  if (is.null(x)) return(tibble::tibble(
    candidate_id = character(), chromosome = character(), start = numeric(),
    end = numeric(), candidate_class = character()
  ))
  if (inherits(x, "detect_inversions")) x <- x$candidates
  x <- tibble::as_tibble(x)
  names(x) <- tolower(names(x))
  if (!all(c("chromosome", "start", "end") %in% names(x))) {
    rlang::abort("`inversion.regions` requires chromosome, start, and end columns.")
  }
  if (!"candidate_id" %in% names(x)) x$candidate_id <- paste0("REGION-", seq_len(nrow(x)))
  if (!"candidate_class" %in% names(x)) x$candidate_class <- "candidate genomic region"
  dplyr::transmute(
    x, candidate_id = as.character(.data$candidate_id),
    chromosome = as.character(.data$chromosome), start = as.numeric(.data$start),
    end = as.numeric(.data$end),
    candidate_class = as.character(.data$candidate_class)
  )
}

.genome_context_join_statistics <- function(context, statistics) {
  if (is.null(statistics)) return(context)
  statistics <- tibble::as_tibble(statistics)
  names(statistics) <- tolower(names(statistics))
  chrom.name <- intersect(c("chromosome", "chrom", "chr"), names(statistics))[1L]
  position.name <- intersect(c("position", "pos", "bp"), names(statistics))[1L]
  if (is.na(chrom.name) || is.na(position.name)) {
    rlang::abort("`scan.statistics` requires chromosome and position columns.")
  }
  statistics$chromosome <- as.character(statistics[[chrom.name]])
  statistics$position <- as.numeric(statistics[[position.name]])
  numeric.names <- setdiff(
    names(statistics)[vapply(statistics, is.numeric, logical(1))], position.name
  )
  if (!length(numeric.names)) return(context)
  summaries <- purrr::map_dfr(seq_len(nrow(context)), function(i) {
    use <- statistics$chromosome == context$chromosome[i] &
      statistics$position >= context$start[i] & statistics$position <= context$end[i]
    values <- statistics[use, numeric.names, drop = FALSE]
    out <- purrr::map_dfc(numeric.names, function(name) {
      value <- values[[name]]
      values <- c(
        if (length(value)) mean(value, na.rm = TRUE) else NA_real_,
        if (length(value)) stats::median(value, na.rm = TRUE) else NA_real_,
        if (any(is.finite(value))) max(value, na.rm = TRUE) else NA_real_
      )
      tibble::as_tibble(stats::setNames(
        as.list(values), paste0(name, c("_mean", "_median", "_max"))
      ))
    })
    dplyr::mutate(out, window_id = context$window_id[i], .before = 1L)
  })
  dplyr::left_join(context, summaries, by = "window_id")
}

#' @export
print.genome_scan_context <- function(x, ...) {
  cat("Genome-scan context\n")
  cat("  Windows:", nrow(x$context), "\n")
  cat("  Annotated candidate regions:", nrow(x$regions), "\n")
  cat("  Results folder:", x$path.folder, "\n")
  invisible(x)
}
