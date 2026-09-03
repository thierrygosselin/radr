#' Evaluate read support for rare alleles
#'
#' `detect_rare_alleles()` identifies rare alleles and evaluates their
#' allele-specific read support. It does not diagnose null alleles or allele
#' dropout directly, and it does not remove markers. The observed minor allele
#' is recalculated from the active samples, so results do not depend on whether
#' that allele is encoded as REF or ALT.
#'
#' @section Evidence categories:
#' \itemize{
#'   \item `MONOMORPHIC`: no minor allele is present;
#'   \item `NOT_RARE`: minor-allele count exceeds the requested maximum;
#'   \item `INSUFFICIENT_DATA`: sample size or call rate is inadequate;
#'   \item `SUPPORTED_RARE`: the rare allele passes the requested checks;
#'   \item `LOW_READ_SUPPORT`: carrier depth is incomplete or weak;
#'   \item `ALLELE_BALANCE_WARNING`: heterozygote balance is outside bounds;
#'   \item `DEPTH_AND_BALANCE_WARNING`: both warnings are present.
#' }
#' These are review categories, not universal biological boundaries.
#' Sequencing error, mapping bias, dropout, contamination, paralogy, batch
#' effects, and genuine rare variation can produce overlapping patterns.
#'
#' @section Stratified diagnostics:
#' With `by.strata = TRUE`, statistics are calculated in the pooled samples and
#' independently within each value of `group.column`. Minor-allele identity is
#' inferred within each group and can switch between REF and ALT.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object.
#' @param strata Optional sample metadata data frame or tab-delimited file.
#' @param group.column Metadata column used for stratified diagnostics.
#' @param by.strata Calculate pooled and group-specific statistics.
#' @param max.minor.allele.count Largest minor-allele count classified as rare.
#' @param min.call.rate Minimum called-genotype proportion.
#' @param min.samples.per.group Minimum samples in a diagnostic group.
#' @param min.depth.completeness Minimum proportion of carriers with REF and
#'   ALT depth.
#' @param min.minor.read.depth Minimum minor-allele depth in one carrier.
#' @param min.total.minor.depth Minimum summed minor depth across carriers.
#' @param max.low.support.fraction Maximum proportion of weakly supported
#'   carriers.
#' @param min.heterozygote.balance Lower acceptable minor-read fraction in
#'   heterozygotes; the upper bound is one minus this value.
#' @param min.heterozygotes.balance Minimum heterozygotes with depth required
#'   for the balance check.
#' @param chunk.size Number of active variants processed per GDS chunk.
#' @param save.plots Save diagnostic plots.
#' @param plot.formats One or both of `"png"` and `"pdf"`.
#' @param verbose Display concise progress messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#'
#' @return A `detect_rare_alleles` object with marker and stratum
#' statistics, review candidates, thresholds, plots, output paths, depth
#' provenance, and confirmation that the incoming GDS filter was restored.
#'
#' @seealso \code{\link{detect_paralogs}},
#'   \code{\link{detect_biallelic_problems}}, \code{\link{filter_ma}}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' support <- radr::detect_rare_alleles(
#'   data = "study.gds", strata = "sample_metadata.tsv"
#' )
#' support$candidates
#' }
detect_rare_alleles <- function(
    data, strata = NULL, group.column = "STRATA", by.strata = TRUE,
    max.minor.allele.count = 3L, min.call.rate = 0.8,
    min.samples.per.group = 10L, min.depth.completeness = 0.8,
    min.minor.read.depth = 3, min.total.minor.depth = 6,
    max.low.support.fraction = 0.5, min.heterozygote.balance = 0.2,
    min.heterozygotes.balance = 1L, chunk.size = 2000L,
    save.plots = TRUE, plot.formats = c("png", "pdf"), verbose = TRUE, ...
) {
  force(data)
  .start <- tgbase::startup(
    package = "radr", f.name = "detect_rare_alleles",
    verbose = verbose
  )
  stamp <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)
  dots <- rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)
  unknown <- setdiff(names(dots), c("path.folder", "internal"))
  if (length(unknown)) rlang::abort(paste0(
    "Unknown argument(s): ", paste(unknown, collapse = ", "), "."
  ))
  parent <- dots$path.folder %||% getwd()
  internal <- isTRUE(dots$internal)
  max.minor.allele.count <- .paralog_check_count(
    max.minor.allele.count, "max.minor.allele.count", 1L
  )
  min.samples.per.group <- .paralog_check_count(
    min.samples.per.group, "min.samples.per.group", 2L
  )
  min.heterozygotes.balance <- .paralog_check_count(
    min.heterozygotes.balance, "min.heterozygotes.balance", 1L
  )
  chunk.size <- .paralog_check_count(chunk.size, "chunk.size", 1L)
  for (nm in c("min.call.rate", "min.depth.completeness",
               "max.low.support.fraction", "min.heterozygote.balance")) {
    .paralog_check_probability(get(nm), nm)
  }
  if (min.heterozygote.balance >= 0.5) rlang::abort(
    "`min.heterozygote.balance` must be smaller than 0.5."
  )
  .paralog_check_nonnegative(min.minor.read.depth, "min.minor.read.depth")
  .paralog_check_nonnegative(min.total.minor.depth, "min.total.minor.depth")
  .paralog_check_flag(by.strata, "by.strata")
  .paralog_check_flag(save.plots, "save.plots")
  .paralog_check_flag(verbose, "verbose")
  if (!is.character(group.column) || length(group.column) != 1L ||
      is.na(group.column) || !nzchar(group.column)) rlang::abort(
    "`group.column` must be one non-empty column name."
  )
  plot.formats <- unique(tolower(as.character(plot.formats)))
  if (!length(plot.formats) || any(!plot.formats %in% c("png", "pdf"))) {
    rlang::abort("`plot.formats` must contain `png`, `pdf`, or both.")
  }
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  parent <- normalizePath(parent, mustWork = TRUE)
  folder <- if (internal) parent else radr_folder(
    rad.folder = paste0("detect_rare_alleles_", stamp),
    path.folder = parent, prefix.int = TRUE
  )
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  if (verbose && !internal) .rare_message("Folder created: ", basename(folder))
  thresholds <- tibble::tibble(
    parameter = c(
      "max.minor.allele.count", "min.call.rate", "min.samples.per.group",
      "min.depth.completeness", "min.minor.read.depth",
      "min.total.minor.depth", "max.low.support.fraction",
      "min.heterozygote.balance", "min.heterozygotes.balance"
    ),
    value = c(
      max.minor.allele.count, min.call.rate, min.samples.per.group,
      min.depth.completeness, min.minor.read.depth, min.total.minor.depth,
      max.low.support.fraction, min.heterozygote.balance,
      min.heterozygotes.balance
    )
  )
  args <- dplyr::bind_rows(
    dplyr::mutate(thresholds, value = as.character(.data$value)),
    tibble::tibble(
    parameter = c("group.column", "by.strata", "chunk.size", "save.plots",
                  "plot.formats"),
    value = c(group.column, by.strata, chunk.size, save.plots,
              paste(plot.formats, collapse = ","))
    )
  )
  args.file <- file.path(
    folder, paste0("radr_detect_rare_alleles_args_", stamp, ".tsv")
  )
  readr::write_tsv(args, args.file)
  opened <- FALSE
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.character(data) && length(data) == 1L && file.exists(data) &&
             grepl("\\.gds$", data, ignore.case = TRUE)) {
    gds <- SeqArray::seqOpen(data)
    opened <- TRUE
  } else rlang::abort(
    "`data` must be a GDS filepath or an open SeqVarGDSClass object."
  )
  on.exit(if (opened) try(SeqArray::seqClose(gds), silent = TRUE), add = TRUE)
  before <- SeqArray::seqGetFilter(gds)
  SeqArray::seqFilterPush(gds)
  pushed <- TRUE
  on.exit(if (pushed) try(SeqArray::seqFilterPop(gds), silent = TRUE), add = TRUE)
  samples <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  variants <- SeqArray::seqGetData(gds, "variant.id")
  if (!length(samples) || !length(variants)) rlang::abort(
    "The active GDS selection contains no samples or markers."
  )
  if (!genometranslator::detect_biallelic_markers(gds)) rlang::abort(
    "Rare-allele diagnostics require active biallelic markers."
  )
  metadata <- .paralog_read_metadata(
    strata, gds, samples, group.column, by.strata
  )
  samples <- samples[samples %in% metadata$INDIVIDUALS]
  metadata <- metadata[match(samples, metadata$INDIVIDUALS), , drop = FALSE]
  SeqArray::seqSetFilter(gds, sample.id = samples, verbose = FALSE)
  groups <- if (by.strata) as.character(metadata[[group.column]]) else
    rep("POOLED", length(samples))
  if (anyNA(groups) || any(!nzchar(trimws(groups)))) rlang::abort(paste0(
    "`", group.column, "` contains missing or empty values."
  ))
  if (any(groups == "OVERALL")) rlang::abort(
    "`OVERALL` is reserved for pooled diagnostics."
  )
  if (by.strata && identical(unique(groups), "POOLED")) by.strata <- FALSE
  marker.metadata <- .sex_marker_metadata(gds, variants)
  chunks <- split(seq_along(variants), ceiling(seq_along(variants) / chunk.size))
  if (verbose) {
    .rare_message("Active selection: ", length(variants), " marker(s), ",
                  length(samples), " sample(s).")
    .rare_message("Evaluating ", length(chunks), " marker chunk(s).")
  }
  pieces <- vector("list", length(chunks))
  sources <- character()
  for (i in seq_along(chunks)) {
    idx <- chunks[[i]]
    dosage <- .sex_get_dosage(gds, variants[idx], samples)
    depth <- .paralog_get_allele_depth(gds, variants[idx], samples)
    if (is.null(depth$ref) || is.null(depth$alt)) rlang::abort(paste0(
      "Allele-specific depth is required. Import AD or DArT REF/ALT ",
      "counts into the GDS."
    ))
    sources <- c(sources, depth$source)
    pieces[[i]] <- .rare_chunk(
      variants[idx], dosage, depth$ref, depth$alt, groups, by.strata,
      max.minor.allele.count, min.call.rate, min.samples.per.group,
      min.depth.completeness, min.minor.read.depth, min.total.minor.depth,
      max.low.support.fraction, min.heterozygote.balance,
      min.heterozygotes.balance
    )
    if (verbose && (i == length(chunks) || i %% 10L == 0L)) {
      .rare_message("Processed ", i, " of ", length(chunks), " chunk(s).")
    }
  }
  strata.stats <- dplyr::bind_rows(pieces)
  pooled <- strata.stats[strata.stats$GROUP == "OVERALL", , drop = FALSE]
  grouped <- strata.stats[strata.stats$GROUP != "OVERALL", , drop = FALSE]
  consistency <- if (nrow(grouped)) grouped |>
    dplyr::group_by(.data$VARIANT_ID) |>
    dplyr::summarise(
      RARE_STRATA = sum(.data$IS_RARE),
      REVIEW_STRATA = sum(.data$REVIEW_CANDIDATE),
      LOW_SUPPORT_STRATA = sum(.data$LOW_READ_SUPPORT),
      BALANCE_WARNING_STRATA = sum(.data$ALLELE_BALANCE_WARNING),
      .groups = "drop"
    ) else tibble::tibble(
      VARIANT_ID = variants, RARE_STRATA = 0L, REVIEW_STRATA = 0L,
      LOW_SUPPORT_STRATA = 0L, BALANCE_WARNING_STRATA = 0L
    )
  statistics <- marker.metadata |>
    dplyr::left_join(pooled, by = "VARIANT_ID") |>
    dplyr::left_join(consistency, by = "VARIANT_ID") |>
    dplyr::mutate(
      ANY_REVIEW_CANDIDATE = .data$REVIEW_CANDIDATE |
        .data$REVIEW_STRATA > 0L,
      STRATUM_SPECIFIC_WARNING = !.data$REVIEW_CANDIDATE &
        .data$REVIEW_STRATA > 0L
    )
  candidates <- statistics[statistics$ANY_REVIEW_CANDIDATE %in% TRUE,
                           , drop = FALSE]
  plots <- .rare_plots(strata.stats, max.minor.allele.count,
                       min.minor.read.depth, min.heterozygote.balance)
  files <- c(
    args.file, file.path(folder, "rare_allele_marker_statistics.tsv"),
    file.path(folder, "rare_allele_stratum_statistics.tsv"),
    file.path(folder, "rare_allele_review_candidates.tsv"),
    file.path(folder, "rare_allele_thresholds.tsv")
  )
  readr::write_tsv(statistics, files[[2L]])
  readr::write_tsv(strata.stats, files[[3L]])
  readr::write_tsv(candidates, files[[4L]])
  readr::write_tsv(thresholds, files[[5L]])
  if (save.plots) files <- c(files, .rare_save_plots(
    plots, folder, plot.formats
  ))
  SeqArray::seqFilterPop(gds)
  pushed <- FALSE
  after <- SeqArray::seqGetFilter(gds)
  restored <- identical(before$sample.sel, after$sample.sel) &&
    identical(before$variant.sel, after$variant.sel)
  if (!restored) rlang::abort("The active GDS selection was not restored.")
  if (verbose) {
    .rare_message("Markers requiring review: ", nrow(candidates), ".")
    .rare_message("Allele-depth source: ",
                  paste(unique(sources), collapse = ", "), ".")
    .rare_message("Rare-allele diagnostics written to: ", folder)
  }
  out <- list(
    statistics = tibble::as_tibble(statistics),
    strata.statistics = tibble::as_tibble(strata.stats),
    candidates = tibble::as_tibble(candidates), thresholds = thresholds,
    plots = plots, allele.depth.source = unique(sources), path.folder = folder,
    output.files = tibble::tibble(files = files),
    active.selection.restored = restored
  )
  class(out) <- c("detect_rare_alleles", class(out))
  out
}

.rare_chunk <- function(
    ids, dosage, ref, alt, groups, by.strata, max.mac, min.call.rate,
    min.samples, min.completeness, min.depth, min.total.depth,
    max.low.fraction, min.balance, min.hets
) {
  if (!identical(dim(dosage), dim(ref)) || !identical(dim(dosage), dim(alt))) {
    rlang::abort("Dosage and allele-depth matrices have different dimensions.")
  }
  sets <- list(OVERALL = seq_len(nrow(dosage)))
  if (by.strata) sets <- c(sets, split(seq_along(groups), groups, drop = TRUE))
  dplyr::bind_rows(lapply(names(sets), function(group) {
    ii <- sets[[group]]
    x <- dosage[ii, , drop = FALSE]
    rd <- ref[ii, , drop = FALSE]
    ad <- alt[ii, , drop = FALSE]
    called <- is.finite(x)
    n <- nrow(x)
    nc <- colSums(called)
    ac <- colSums(replace(x, !called, 0), na.rm = TRUE)
    rc <- 2 * nc - ac
    minor.alt <- ac <= rc
    mac <- pmin(ac, rc)
    maf <- mac / (2 * nc)
    maf[!is.finite(maf)] <- NA_real_
    call.rate <- nc / n
    carrier <- matrix(FALSE, nrow(x), ncol(x))
    if (any(minor.alt)) carrier[, minor.alt] <-
      x[, minor.alt, drop = FALSE] > 0
    if (any(!minor.alt)) carrier[, !minor.alt] <-
      x[, !minor.alt, drop = FALSE] < 2
    carrier[!called] <- FALSE
    complete <- is.finite(rd) & is.finite(ad) & rd >= 0 & ad >= 0
    md <- rd
    od <- ad
    if (any(minor.alt)) {
      md[, minor.alt] <- ad[, minor.alt, drop = FALSE]
      od[, minor.alt] <- rd[, minor.alt, drop = FALSE]
    }
    carrier.complete <- carrier & complete
    carriers <- colSums(carrier)
    carriers.depth <- colSums(carrier.complete)
    completeness <- carriers.depth / carriers
    completeness[!is.finite(completeness)] <- NA_real_
    supported <- replace(md, !carrier.complete, NA_real_)
    total.depth <- colSums(supported, na.rm = TRUE)
    median.depth <- apply(supported, 2L, stats::median, na.rm = TRUE)
    median.depth[!is.finite(median.depth)] <- NA_real_
    low.n <- colSums(carrier.complete & md < min.depth)
    low.fraction <- low.n / carriers.depth
    low.fraction[!is.finite(low.fraction)] <- NA_real_
    het.complete <- called & abs(x - 1) < 1e-8 & complete
    nhet <- colSums(het.complete)
    balance <- md / (md + od)
    balance[!het.complete] <- NA_real_
    median.balance <- apply(balance, 2L, stats::median, na.rm = TRUE)
    median.balance[!is.finite(median.balance)] <- NA_real_
    eligible <- n >= min.samples & call.rate >= min.call.rate
    rare <- eligible & mac >= 1 & mac <= max.mac
    low <- rare & (!is.finite(completeness) | completeness < min.completeness |
      total.depth < min.total.depth | !is.finite(low.fraction) |
      low.fraction > max.low.fraction)
    imbalance <- rare & nhet >= min.hets & is.finite(median.balance) &
      (median.balance < min.balance | median.balance > 1 - min.balance)
    review <- rare & (low | imbalance)
    status <- rep("INSUFFICIENT_DATA", ncol(x))
    status[eligible & mac == 0] <- "MONOMORPHIC"
    status[eligible & mac > max.mac] <- "NOT_RARE"
    status[rare] <- "SUPPORTED_RARE"
    status[rare & low] <- "LOW_READ_SUPPORT"
    status[rare & imbalance] <- "ALLELE_BALANCE_WARNING"
    status[rare & low & imbalance] <- "DEPTH_AND_BALANCE_WARNING"
    reason <- vapply(seq_len(ncol(x)), function(j) {
      if (eligible[j]) return("")
      paste(c(if (n < min.samples) "too_few_samples",
              if (call.rate[j] < min.call.rate) "low_call_rate"),
            collapse = ";")
    }, character(1))
    tibble::tibble(
      VARIANT_ID = ids, GROUP = group, N_SAMPLES = n,
      NUMBER_CALLED = unname(nc), CALL_RATE = unname(call.rate),
      ALT_ALLELE_COUNT = unname(ac), REF_ALLELE_COUNT = unname(rc),
      MINOR_ALLELE = unname(ifelse(minor.alt, "ALT", "REF")),
      MINOR_ALLELE_COUNT = unname(mac),
      MINOR_ALLELE_FREQUENCY = unname(maf),
      NUMBER_CARRIERS = unname(carriers),
      CARRIERS_WITH_DEPTH = unname(carriers.depth),
      DEPTH_COMPLETENESS = unname(completeness),
      TOTAL_MINOR_READ_DEPTH = unname(total.depth),
      MEDIAN_MINOR_READ_DEPTH = unname(median.depth),
      LOW_SUPPORT_CARRIERS = unname(low.n),
      LOW_SUPPORT_FRACTION = unname(low.fraction),
      HETEROZYGOTES_WITH_DEPTH = unname(nhet),
      MEDIAN_MINOR_READ_FRACTION_HET = unname(median.balance),
      ELIGIBLE = unname(eligible), IS_RARE = unname(rare),
      LOW_READ_SUPPORT = unname(low),
      ALLELE_BALANCE_WARNING = unname(imbalance),
      REVIEW_CANDIDATE = unname(review), STATUS = status,
      INELIGIBLE_REASON = reason
    )
  }))
}

.rare_plots <- function(x, max.mac, min.depth, min.balance) {
  pooled <- x[x$GROUP == "OVERALL", , drop = FALSE]
  rare <- x[x$IS_RARE %in% TRUE, , drop = FALSE]
  list(
    minor.allele.count = ggplot2::ggplot(
      pooled, ggplot2::aes(x = .data$MINOR_ALLELE_COUNT)
    ) + ggplot2::geom_histogram(binwidth = 1, boundary = -0.5) +
      ggplot2::geom_vline(xintercept = max.mac + 0.5, linetype = 2) +
      ggplot2::labs(x = "Minor-allele count", y = "Markers") +
      ggplot2::theme_bw(),
    minor.allele.depth = ggplot2::ggplot(
      rare, ggplot2::aes(.data$MINOR_ALLELE_COUNT,
                         .data$MEDIAN_MINOR_READ_DEPTH,
                         colour = .data$STATUS)
    ) + ggplot2::geom_hline(yintercept = min.depth, linetype = 2) +
      ggplot2::geom_point(alpha = 0.65, na.rm = TRUE) +
      ggplot2::facet_wrap(ggplot2::vars(.data$GROUP)) +
      ggplot2::labs(x = "Minor-allele count", y = "Median supporting depth",
                    colour = "Diagnostic") + ggplot2::theme_bw(),
    heterozygote.balance = ggplot2::ggplot(
      rare, ggplot2::aes(.data$MINOR_ALLELE_COUNT,
                         .data$MEDIAN_MINOR_READ_FRACTION_HET,
                         colour = .data$STATUS)
    ) + ggplot2::geom_hline(
      yintercept = c(min.balance, 1 - min.balance), linetype = 2
    ) + ggplot2::geom_point(alpha = 0.65, na.rm = TRUE) +
      ggplot2::facet_wrap(ggplot2::vars(.data$GROUP)) +
      ggplot2::labs(x = "Minor-allele count", y = "Minor-read fraction",
                    colour = "Diagnostic") + ggplot2::theme_bw()
  )
}

.rare_save_plots <- function(plots, folder, formats) {
  files <- character()
  for (nm in names(plots)) for (format in formats) {
    file <- file.path(folder, paste0(nm, ".", format))
    args <- list(filename = file, plot = plots[[nm]], width = 9, height = 6,
                 units = "in", dpi = 300)
    if (format == "pdf") args$useDingbats <- FALSE
    do.call(ggplot2::ggsave, args)
    files <- c(files, file)
  }
  files
}

.rare_message <- function(...) {
  text <- paste0(...)
  cat(paste(strwrap(text, width = 80), collapse = "\n"), "\n", sep = "")
  invisible(text)
}
