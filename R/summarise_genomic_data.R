#' Summarise genomic data
#'
#' Calculate allele-frequency, missingness, and heterozygosity statistics from
#' active biallelic markers and samples in a GDS. Genotypes are read in chunks.
#'
#' `ALT_FREQUENCY` is the VCF ALT-allele frequency. In contrast,
#' `MINOR_ALLELE_FREQUENCY` is the smaller REF or ALT frequency and is always
#' between zero and 0.5. In stratified results,
#' `POOLED_MINOR_ALLELE_FREQUENCY` tracks the allele that is minor in the
#' pooled sample even where that allele is locally common or fixed.
#'
#' @section FIS:
#' `FIS` is `1 - Ho / He`. It is `NA` for monomorphic markers because expected
#' heterozygosity is zero and the coefficient is undefined.
#'
#' @param data A GDS filepath or an open `SeqVarGDSClass` object.
#' @param strata Optional data frame or TSV containing `INDIVIDUALS` and
#'   `group.column`. If `NULL`, metadata are read from the GDS.
#' @param group.column Metadata column used for stratified summaries.
#' @param by.strata Calculate group summaries in addition to pooled summaries.
#' @param chunk.size Number of active variants read at a time.
#' @param digits Decimal places used in returned and written tables.
#' @param write.files Write summary tables to disk.
#' @param verbose Display progress and summary messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#'
#' @return A `summarise_genomic_data` object containing `marker.statistics`,
#' `stratum.statistics`, `group.summary`, `path.folder`, `output.files`, and
#' `active.selection.restored`. The function does not alter persistent filters.
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' x <- radr::summarise_genomic_data(
#'   "study.gds", strata = "sample_metadata.tsv"
#' )
#' x$marker.statistics
#' }
summarise_genomic_data <- function(
    data, strata = NULL, group.column = "STRATA", by.strata = TRUE,
    chunk.size = 2000L, digits = 6L, write.files = TRUE,
    verbose = TRUE, ...
) {
  force(data)
  .start <- tgbase::startup(
    package = "radr", f.name = "summarise_genomic_data", verbose = verbose
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
  chunk.size <- .paralog_check_count(chunk.size, "chunk.size", 1L)
  digits <- .paralog_check_count(digits, "digits", 0L)
  .paralog_check_flag(by.strata, "by.strata")
  .paralog_check_flag(write.files, "write.files")
  .paralog_check_flag(verbose, "verbose")
  if (!is.character(group.column) || length(group.column) != 1L ||
      is.na(group.column) || !nzchar(group.column)) {
    rlang::abort("`group.column` must be one non-empty column name.")
  }

  dir.create(parent.folder, recursive = TRUE, showWarnings = FALSE)
  parent.folder <- normalizePath(parent.folder, mustWork = TRUE)
  path.folder <- if (internal) parent.folder else radr_folder(
    rad.folder = paste0("summarise_genomic_data_", file.date),
    path.folder = parent.folder, prefix.int = TRUE
  )
  dir.create(path.folder, recursive = TRUE, showWarnings = FALSE)
  if (verbose && !internal) {
    .summary_message("Folder created: ", basename(path.folder))
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
  if (!genometranslator::detect_biallelic_markers(gds)) {
    rlang::abort("Genomic summaries require active biallelic markers.")
  }

  metadata <- .paralog_read_metadata(
    strata, gds, sample.id, group.column, by.strata
  )
  sample.id <- sample.id[sample.id %in% metadata$INDIVIDUALS]
  metadata <- metadata[
    match(sample.id, metadata$INDIVIDUALS), , drop = FALSE
  ]
  SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
  groups <- as.character(metadata[[group.column]])
  if (anyNA(groups) || any(!nzchar(trimws(groups)))) {
    rlang::abort(paste0("`", group.column, "` contains missing values."))
  }
  if (by.strata && "OVERALL" %in% groups) {
    rlang::abort("`OVERALL` is reserved for the pooled GDS summary.")
  }

  marker.metadata <- .sex_marker_metadata(gds, variant.id)
  if (verbose) {
    .summary_message(
      "Summarising ", length(variant.id), " markers for ",
      length(sample.id), " samples."
    )
  }
  chunks <- split(
    seq_along(variant.id), ceiling(seq_along(variant.id) / chunk.size)
  )
  statistics <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    index <- chunks[[i]]
    statistics[[i]] <- .summary_chunk_statistics(
      variant.id[index],
      .sex_get_dosage(gds, variant.id[index], sample.id),
      groups, by.strata
    )
    if (verbose && (i == length(chunks) || i %% 10L == 0L)) {
      .summary_message("Processed marker chunk ", i, " of ", length(chunks), ".")
    }
  }
  stratum.statistics <- dplyr::bind_rows(statistics)
  pooled <- stratum.statistics |>
    dplyr::filter(.data$GROUP == "OVERALL") |>
    dplyr::select(
      "VARIANT_ID",
      POOLED_MINOR_ALLELE = "MINOR_ALLELE"
    )
  stratum.statistics <- stratum.statistics |>
    dplyr::left_join(pooled, by = "VARIANT_ID") |>
    dplyr::mutate(
      POOLED_MINOR_ALLELE_FREQUENCY = dplyr::if_else(
        .data$POOLED_MINOR_ALLELE == "REF",
        .data$REF_FREQUENCY, .data$ALT_FREQUENCY
      )
    )
  marker.statistics <- stratum.statistics |>
    dplyr::filter(.data$GROUP == "OVERALL") |>
    dplyr::select(-"GROUP") |>
    dplyr::left_join(marker.metadata, by = "VARIANT_ID") |>
    dplyr::select(
      dplyr::any_of(c("VARIANT_ID", "MARKERS", "CHROM", "POS")),
      dplyr::everything()
    )
  group.summary <- .summary_group_statistics(stratum.statistics)

  count.columns <- c(
    "VARIANT_ID", "N_SAMPLES", "NUMBER_CALLED", "NUMBER_MISSING",
    "HOM_REF", "HET", "HOM_ALT", "REF_ALLELE_COUNT",
    "ALT_ALLELE_COUNT", "MINOR_ALLELE_COUNT"
  )
  decimals <- setdiff(
    names(stratum.statistics)[
      vapply(stratum.statistics, is.numeric, logical(1))
    ], count.columns
  )
  stratum.statistics[decimals] <- lapply(
    stratum.statistics[decimals], round, digits = digits
  )
  marker.decimals <- intersect(decimals, names(marker.statistics))
  marker.statistics[marker.decimals] <- lapply(
    marker.statistics[marker.decimals], round, digits = digits
  )
  group.decimals <- setdiff(
    names(group.summary)[vapply(group.summary, is.numeric, logical(1))],
    c("N_SAMPLES", "N_MARKERS_TOTAL", "N_MARKERS_CALLED")
  )
  group.summary[group.decimals] <- lapply(
    group.summary[group.decimals], round, digits = digits
  )

  output.files <- character()
  if (write.files) {
    names.out <- c(
      "genomic_marker_statistics.tsv",
      "genomic_stratum_statistics.tsv",
      "genomic_group_summary.tsv"
    )
    tables <- list(marker.statistics, stratum.statistics, group.summary)
    output.files <- file.path(path.folder, names.out)
    for (i in seq_along(output.files)) {
      readr::write_tsv(tables[[i]], output.files[[i]], na = "NA")
    }
    args.file <- file.path(
      path.folder,
      paste0("radr_summarise_genomic_data_args_", file.date, ".tsv")
    )
    readr::write_tsv(tibble::tibble(
      argument = c(
        "group.column", "by.strata", "chunk.size", "digits", "write.files"
      ),
      value = as.character(c(
        group.column, by.strata, chunk.size, digits, write.files
      ))
    ), args.file)
    output.files <- c(output.files, args.file)
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
    .summary_message("Genomic summaries written to: ", path.folder)
  }

  out <- list(
    marker.statistics = tibble::as_tibble(marker.statistics),
    stratum.statistics = tibble::as_tibble(stratum.statistics),
    group.summary = tibble::as_tibble(group.summary),
    path.folder = path.folder,
    output.files = tibble::tibble(files = output.files),
    active.selection.restored = restored
  )
  class(out) <- c("summarise_genomic_data", class(out))
  out
}

.summary_chunk_statistics <- function(
    variant.id, dosage, groups, by.strata = TRUE
) {
  group.index <- list(OVERALL = seq_len(nrow(dosage)))
  if (by.strata) {
    group.index <- c(group.index, split(seq_along(groups), groups))
  }
  dplyr::bind_rows(lapply(names(group.index), function(group) {
    x <- dosage[group.index[[group]], , drop = FALSE]
    called <- colSums(!is.na(x))
    hom.ref <- colSums(x == 0, na.rm = TRUE)
    het <- colSums(x == 1, na.rm = TRUE)
    hom.alt <- colSums(x == 2, na.rm = TRUE)
    alt.count <- colSums(x, na.rm = TRUE)
    ref.count <- 2 * called - alt.count
    alt.frequency <- ifelse(called > 0, alt.count / (2 * called), NA_real_)
    ref.frequency <- ifelse(called > 0, ref.count / (2 * called), NA_real_)
    maf <- pmin(ref.frequency, alt.frequency)
    ho <- ifelse(called > 0, het / called, NA_real_)
    he <- 2 * ref.frequency * alt.frequency
    fis <- ifelse(is.finite(he) & he > 0, 1 - ho / he, NA_real_)
    minor <- ifelse(
      is.na(maf), NA_character_,
      ifelse(alt.frequency <= ref.frequency, "ALT", "REF")
    )
    tibble::tibble(
      VARIANT_ID = variant.id, GROUP = group, N_SAMPLES = nrow(x),
      NUMBER_CALLED = called, NUMBER_MISSING = nrow(x) - called,
      CALL_RATE = called / nrow(x), HOM_REF = hom.ref, HET = het,
      HOM_ALT = hom.alt, REF_ALLELE_COUNT = ref.count,
      ALT_ALLELE_COUNT = alt.count, REF_FREQUENCY = ref.frequency,
      ALT_FREQUENCY = alt.frequency, MINOR_ALLELE = minor,
      MINOR_ALLELE_COUNT = ifelse(
        minor == "ALT", alt.count,
        ifelse(minor == "REF", ref.count, NA_real_)
      ),
      MINOR_ALLELE_FREQUENCY = maf,
      OBSERVED_HETEROZYGOSITY = ho,
      EXPECTED_HETEROZYGOSITY = he, FIS = fis
    )
  }))
}

.summary_group_statistics <- function(statistics) {
  safe.mean <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }
  statistics |>
    dplyr::group_by(.data$GROUP) |>
    dplyr::summarise(
      N_SAMPLES = max(.data$N_SAMPLES),
      N_MARKERS_TOTAL = dplyr::n(),
      N_MARKERS_CALLED = sum(.data$NUMBER_CALLED > 0),
      MEAN_CALL_RATE = safe.mean(.data$CALL_RATE),
      MEAN_ALT_FREQUENCY = safe.mean(.data$ALT_FREQUENCY),
      MEAN_MINOR_ALLELE_FREQUENCY = safe.mean(
        .data$MINOR_ALLELE_FREQUENCY
      ),
      MEAN_OBSERVED_HETEROZYGOSITY = safe.mean(
        .data$OBSERVED_HETEROZYGOSITY
      ),
      MEAN_EXPECTED_HETEROZYGOSITY = safe.mean(
        .data$EXPECTED_HETEROZYGOSITY
      ),
      MEAN_MARKER_FIS = safe.mean(.data$FIS),
      WEIGHTED_FIS = {
        expected <- sum(
          .data$EXPECTED_HETEROZYGOSITY * .data$NUMBER_CALLED,
          na.rm = TRUE
        )
        observed <- sum(.data$HET, na.rm = TRUE)
        if (expected > 0) 1 - observed / expected else NA_real_
      },
      .groups = "drop"
    )
}

.summary_message <- function(...) {
  message(paste(strwrap(paste0(...), width = 80, exdent = 2), collapse = "\n"))
}
