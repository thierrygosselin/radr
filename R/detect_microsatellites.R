#' Detect candidate microsatellite motifs
#'
#' Scans marker sequences for perfect tandem repeats, also called simple
#' sequence repeats (SSRs). A detected repeat is a sequence motif consistent
#' with a microsatellite; it is not evidence that the locus is polymorphic or
#' suitable as a microsatellite marker.
#'
#' The scanner is implemented internally and requires no downloaded software.
#' Repeat thresholds are explicit and saved with every analysis. By default,
#' motif lengths one through six require 10, 6, 5, 4, 4, and 4 consecutive
#' copies, respectively. Imperfect or interrupted repeats are not detected.
#'
#' @section Interpretation:
#' Markers are assigned one of three states:
#' \itemize{
#'   \item `SSR_DETECTED`: at least one qualifying perfect repeat was found;
#'   \item `NO_SSR_DETECTED`: a sequence was scanned but no repeat passed;
#'   \item `SEQUENCE_UNAVAILABLE`: the marker could not be evaluated.
#' }
#' Markers without sequence are never placed among markers without detected
#' repeats. Removing an SSR-containing marker is a separate filtering decision.
#'
#' @param data A GDS filepath or open `SeqVarGDSClass` object. Marker metadata
#'   must contain `SEQUENCE`.
#' @param motif.lengths Integer motif lengths to scan.
#' @param min.repeats Minimum consecutive copies for each motif length. Supply
#'   one value per `motif.lengths`, or one named value for every length.
#' @param canonicalize Group rotational and reverse-complement equivalents under
#'   one canonical motif label.
#' @param chunk.size Number of marker sequences scanned per chunk.
#' @param save.plots Save summary plots.
#' @param plot.formats One or both of `"png"` and `"pdf"`.
#' @param verbose Display concise progress messages.
#' @param ... Common arguments: `path.folder` and `internal`.
#'
#' @return A `detect_microsatellites` object containing marker status, repeat
#' hits, candidate markers, markers with no detected repeats, markers lacking
#' sequence, thresholds, plots, output paths, and confirmation that the active
#' GDS filter was unchanged.
#'
#' @references
#' Wang X, Lu P, Luo Z (2013). GMATo: A novel tool for the identification and
#' analysis of microsatellites in large genomes. *Bioinformation*, 9, 541-544.
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @note Thanks to Peter Grewe for suggesting this diagnostic for radr.
#' @export
#' @examples
#' \dontrun{
#' ssr <- radr::detect_microsatellites(data = "study.gds")
#' ssr$candidates
#' }
detect_microsatellites <- function(
    data,
    motif.lengths = 1:6,
    min.repeats = c(10, 6, 5, 4, 4, 4),
    canonicalize = TRUE,
    chunk.size = 50000L,
    save.plots = TRUE,
    plot.formats = c("png", "pdf"),
    verbose = TRUE,
    ...
) {
  force(data)
  .start <- tgbase::startup(
    package = "radr", f.name = "detect_microsatellites", verbose = verbose
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
  motif.lengths <- .ssr_validate_lengths(motif.lengths)
  min.repeats <- .ssr_validate_repeats(min.repeats, motif.lengths)
  chunk.size <- .paralog_check_count(chunk.size, "chunk.size", 1L)
  .paralog_check_flag(canonicalize, "canonicalize")
  .paralog_check_flag(save.plots, "save.plots")
  .paralog_check_flag(verbose, "verbose")
  plot.formats <- unique(tolower(as.character(plot.formats)))
  if (!length(plot.formats) || any(!plot.formats %in% c("png", "pdf"))) {
    rlang::abort("`plot.formats` must contain `png`, `pdf`, or both.")
  }

  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  parent <- normalizePath(parent, mustWork = TRUE)
  folder <- if (internal) parent else radr_folder(
    rad.folder = paste0("detect_microsatellites_", stamp),
    path.folder = parent, prefix.int = TRUE
  )
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  if (verbose && !internal) .ssr_message("Folder created: ", basename(folder))

  thresholds <- tibble::tibble(
    MOTIF_LENGTH = motif.lengths,
    MINIMUM_REPEATS = unname(min.repeats),
    MINIMUM_ARRAY_LENGTH = motif.lengths * unname(min.repeats)
  )
  arguments <- tibble::tibble(
    argument = c("motif.lengths", "min.repeats", "canonicalize",
                 "chunk.size", "save.plots", "plot.formats"),
    value = c(paste(motif.lengths, collapse = ","),
              paste(min.repeats, collapse = ","), canonicalize, chunk.size,
              save.plots, paste(plot.formats, collapse = ","))
  )
  args.file <- file.path(
    folder, paste0("radr_detect_microsatellites_args_", stamp, ".tsv")
  )
  readr::write_tsv(arguments, args.file)

  opened <- FALSE
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.character(data) && length(data) == 1L && !is.na(data) &&
             file.exists(data) && grepl("\\.gds$", data, ignore.case = TRUE)) {
    gds <- SeqArray::seqOpen(data)
    opened <- TRUE
  } else rlang::abort(
    "`data` must be a GDS filepath or an open SeqVarGDSClass object."
  )
  on.exit(if (opened) try(SeqArray::seqClose(gds), silent = TRUE), add = TRUE)
  selection.before <- SeqArray::seqGetFilter(gds)
  active.id <- SeqArray::seqGetData(gds, "variant.id")
  if (!length(active.id)) rlang::abort(
    "The active GDS selection contains no markers."
  )

  metadata <- genometranslator::extract_markers_metadata(
    gds = gds,
    markers.meta.select = c(
      "VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "SEQUENCE"
    ),
    whitelist = FALSE,
    verbose = FALSE
  )
  if (!"VARIANT_ID" %in% names(metadata)) rlang::abort(
    "Marker metadata must contain `VARIANT_ID`."
  )
  if (!"SEQUENCE" %in% names(metadata)) rlang::abort(paste0(
    "Marker metadata does not contain `SEQUENCE`. Microsatellite motifs ",
    "cannot be reconstructed from genotypes."
  ))
  metadata <- metadata[metadata$VARIANT_ID %in% active.id, , drop = FALSE]
  metadata <- metadata[match(active.id, metadata$VARIANT_ID), , drop = FALSE]
  if (anyNA(metadata$VARIANT_ID)) rlang::abort(
    "Marker metadata are incomplete for the active GDS selection."
  )
  metadata <- metadata |>
    dplyr::distinct(.data$VARIANT_ID, .keep_all = TRUE) |>
    dplyr::mutate(
      SEQUENCE = toupper(trimws(as.character(.data$SEQUENCE))),
      SEQUENCE_LENGTH = stringi::stri_length(.data$SEQUENCE),
      SEQUENCE_AVAILABLE = !is.na(.data$SEQUENCE) & nzchar(.data$SEQUENCE)
    )
  available <- which(metadata$SEQUENCE_AVAILABLE)
  chunks <- split(available, ceiling(seq_along(available) / chunk.size))
  if (!length(chunks)) chunks <- list()
  if (verbose) {
    .ssr_message(
      "Active markers: ", nrow(metadata), "; sequences available: ",
      length(available), "."
    )
    .ssr_message("Scanning perfect tandem repeats in ", length(chunks),
                 " chunk(s).")
  }
  hits <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    idx <- chunks[[i]]
    hits[[i]] <- .ssr_scan_sequences(
      metadata[idx, , drop = FALSE], motif.lengths, min.repeats,
      canonicalize
    )
    if (verbose && (i == length(chunks) || i %% 10L == 0L)) {
      .ssr_message("Processed ", i, " of ", length(chunks), " chunk(s).")
    }
  }
  hits <- dplyr::bind_rows(hits)
  if (!nrow(hits)) hits <- .ssr_empty_hits(metadata)
  hit.ids <- unique(hits$VARIANT_ID)
  marker.status <- metadata |>
    dplyr::mutate(
      SSR_STATUS = dplyr::case_when(
        !.data$SEQUENCE_AVAILABLE ~ "SEQUENCE_UNAVAILABLE",
        .data$VARIANT_ID %in% hit.ids ~ "SSR_DETECTED",
        TRUE ~ "NO_SSR_DETECTED"
      )
    ) |>
    dplyr::left_join(
      if (nrow(hits)) hits |>
        dplyr::count(.data$VARIANT_ID, name = "NUMBER_OF_REPEATS") else
        tibble::tibble(VARIANT_ID = integer(), NUMBER_OF_REPEATS = integer()),
      by = "VARIANT_ID"
    ) |>
    dplyr::mutate(NUMBER_OF_REPEATS = dplyr::coalesce(
      .data$NUMBER_OF_REPEATS, 0L
    ))
  candidates <- marker.status[
    marker.status$SSR_STATUS == "SSR_DETECTED", , drop = FALSE
  ]
  without <- marker.status[
    marker.status$SSR_STATUS == "NO_SSR_DETECTED", , drop = FALSE
  ]
  unavailable <- marker.status[
    marker.status$SSR_STATUS == "SEQUENCE_UNAVAILABLE", , drop = FALSE
  ]
  plots <- .ssr_make_plots(hits)
  files <- c(
    args.file,
    file.path(folder, "microsatellite_repeat_thresholds.tsv"),
    file.path(folder, "microsatellite_hits.tsv"),
    file.path(folder, "microsatellite_marker_status.tsv"),
    file.path(folder, "candidate_microsatellite_markers.tsv"),
    file.path(folder, "markers_without_detected_microsatellites.tsv"),
    file.path(folder, "markers_without_sequence.tsv")
  )
  readr::write_tsv(thresholds, files[[2L]])
  readr::write_tsv(hits, files[[3L]])
  readr::write_tsv(marker.status, files[[4L]])
  readr::write_tsv(candidates, files[[5L]])
  readr::write_tsv(without, files[[6L]])
  readr::write_tsv(unavailable, files[[7L]])
  if (save.plots) files <- c(files, .ssr_save_plots(
    plots, folder, plot.formats
  ))
  selection.after <- SeqArray::seqGetFilter(gds)
  restored <- identical(selection.before$sample.sel, selection.after$sample.sel) &&
    identical(selection.before$variant.sel, selection.after$variant.sel)
  if (!restored) rlang::abort("The active GDS selection changed unexpectedly.")
  if (verbose) {
    .ssr_message("Repeat arrays detected: ", nrow(hits), ".")
    .ssr_message("Markers with detected repeats: ", nrow(candidates), ".")
    .ssr_message("Markers without sequence: ", nrow(unavailable), ".")
    .ssr_message("Microsatellite diagnostics written to: ", folder)
  }
  out <- list(
    marker.status = tibble::as_tibble(marker.status),
    hits = tibble::as_tibble(hits),
    candidates = tibble::as_tibble(candidates),
    without.detected.repeats = tibble::as_tibble(without),
    sequence.unavailable = tibble::as_tibble(unavailable),
    thresholds = thresholds, plots = plots, path.folder = folder,
    output.files = tibble::tibble(files = files),
    active.selection.restored = restored
  )
  class(out) <- c("detect_microsatellites", class(out))
  out
}

.ssr_scan_sequences <- function(metadata, lengths, repeats, canonicalize) {
  sequences <- metadata$SEQUENCE
  pieces <- vector("list", length(lengths))
  for (z in seq_along(lengths)) {
    k <- lengths[[z]]
    minimum <- repeats[[z]]
    pattern <- paste0(
      "([ACGT]{", k, "})", intToUtf8(92), "1{", minimum - 1L, ",}"
    )
    matches <- stringi::stri_match_all_regex(
      sequences, pattern, omit_no_match = TRUE
    )
    locations <- stringi::stri_locate_all_regex(
      sequences, pattern, omit_no_match = TRUE
    )
    rows <- lapply(seq_along(matches), function(i) {
      if (!nrow(matches[[i]])) return(NULL)
      tibble::tibble(
        ROW = i,
        START = unname(locations[[i]][, "start"]),
        END = unname(locations[[i]][, "end"]),
        REPEAT_SEQUENCE = matches[[i]][, 1L],
        MOTIF = matches[[i]][, 2L]
      )
    })
    found <- dplyr::bind_rows(rows)
    if (!nrow(found)) next
    found <- found |>
      dplyr::mutate(
        MOTIF_LENGTH = k,
        REPEAT_COUNT = stringi::stri_length(.data$REPEAT_SEQUENCE) / k,
        ARRAY_LENGTH = stringi::stri_length(.data$REPEAT_SEQUENCE),
        PRIMITIVE_MOTIF = unname(vapply(
          .data$MOTIF, .ssr_primitive_motif, character(1)
        ))
      ) |>
      dplyr::filter(stringi::stri_length(.data$PRIMITIVE_MOTIF) == k) |>
      dplyr::mutate(
        CANONICAL_MOTIF = if (canonicalize) unname(vapply(
          .data$MOTIF, .ssr_canonical_motif, character(1)
        )) else .data$MOTIF
      )
    if (!nrow(found)) next
    pieces[[z]] <- dplyr::bind_cols(
      metadata[found$ROW, setdiff(names(metadata),
                                  c("SEQUENCE", "SEQUENCE_AVAILABLE")),
               drop = FALSE],
      dplyr::select(found, -"ROW")
    )
  }
  out <- dplyr::bind_rows(pieces)
  if (!nrow(out)) return(out)
  dplyr::arrange(out, .data$VARIANT_ID, .data$START, .data$MOTIF_LENGTH)
}

.ssr_primitive_motif <- function(motif) {
  n <- nchar(motif)
  divisors <- which(n %% seq_len(n) == 0L)
  for (k in divisors) {
    unit <- substr(motif, 1L, k)
    if (paste(rep(unit, n / k), collapse = "") == motif) return(unit)
  }
  motif
}

.ssr_canonical_motif <- function(motif) {
  chars <- strsplit(motif, "", fixed = TRUE)[[1L]]
  rotations <- vapply(seq_along(chars), function(i) paste0(
    c(chars[i:length(chars)], chars[seq_len(i - 1L)]), collapse = ""
  ), character(1))
  complement <- chartr("ACGT", "TGCA", rev(chars))
  reverse.rotations <- vapply(seq_along(complement), function(i) paste0(
    c(complement[i:length(complement)], complement[seq_len(i - 1L)]),
    collapse = ""
  ), character(1))
  sort(c(rotations, reverse.rotations))[[1L]]
}

.ssr_validate_lengths <- function(x) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x < 1) || any(x != as.integer(x)) || anyDuplicated(x)) {
    rlang::abort("`motif.lengths` must contain unique positive integers.")
  }
  as.integer(x)
}

.ssr_validate_repeats <- function(x, lengths) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x < 2) || any(x != as.integer(x))) {
    rlang::abort("`min.repeats` must contain whole numbers of at least two.")
  }
  if (!is.null(names(x)) && all(as.character(lengths) %in% names(x))) {
    x <- x[as.character(lengths)]
  } else if (length(x) != length(lengths)) {
    rlang::abort("Supply one `min.repeats` value per motif length.")
  }
  stats::setNames(as.integer(x), lengths)
}

.ssr_empty_hits <- function(metadata) {
  columns <- metadata[0, setdiff(names(metadata),
                                 c("SEQUENCE", "SEQUENCE_AVAILABLE")),
                      drop = FALSE]
  dplyr::bind_cols(columns, tibble::tibble(
    START = integer(), END = integer(), REPEAT_SEQUENCE = character(),
    MOTIF = character(), MOTIF_LENGTH = integer(), REPEAT_COUNT = numeric(),
    ARRAY_LENGTH = integer(), PRIMITIVE_MOTIF = character(),
    CANONICAL_MOTIF = character()
  ))
}

.ssr_make_plots <- function(hits) {
  if (!nrow(hits)) return(list())
  list(
    motif.length = ggplot2::ggplot(
      hits, ggplot2::aes(x = factor(.data$MOTIF_LENGTH))
    ) + ggplot2::geom_bar() +
      ggplot2::labs(x = "Motif length", y = "Repeat arrays") +
      ggplot2::theme_bw(),
    repeat.count = ggplot2::ggplot(
      hits, ggplot2::aes(x = .data$REPEAT_COUNT,
                         fill = factor(.data$MOTIF_LENGTH))
    ) + ggplot2::geom_histogram(binwidth = 1, boundary = 0) +
      ggplot2::labs(x = "Repeat count", y = "Repeat arrays",
                    fill = "Motif length") + ggplot2::theme_bw()
  )
}

.ssr_save_plots <- function(plots, folder, formats) {
  files <- character()
  for (nm in names(plots)) for (format in formats) {
    file <- file.path(folder, paste0("microsatellite_", nm, ".", format))
    args <- list(filename = file, plot = plots[[nm]], width = 8, height = 5,
                 units = "in", dpi = 300)
    if (format == "pdf") args$useDingbats <- FALSE
    do.call(ggplot2::ggsave, args)
    files <- c(files, file)
  }
  files
}

.ssr_message <- function(...) {
  text <- paste0(...)
  cat(paste(strwrap(text, width = 80), collapse = "\n"), "\n", sep = "")
  invisible(text)
}
