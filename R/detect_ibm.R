#' @rdname detect_ibm
#' @title Detect identity-by-missingness structure
#' @description Computes a binary identity-by-missingness matrix directly from a
#' GDS object and returns a heatmap together with marker- and individual-level
#' missingness summaries. The heatmap reveals whether missing genotypes are
#' concentrated in particular individuals, markers, or groups, helping users
#' decide whether individual- or marker-level filtering should be investigated
#' first.
#'
#' @param gds A \code{SeqVarGDSClass} object or a path to a GDS file.
#'
#' @param strata (optional) Path to a strata file with a minimum of 2 columns:
#' \code{INDIVIDUALS} and the grouping column selected in \code{strata.select}.
#' Default: \code{strata = NULL}.
#'
#' @param strata.select (character) Column used to facet individuals in the
#' heatmap.
#' Default: \code{strata.select = "STRATA"}.
#'
#' @param sort.individuals (character) Sorting method for individuals.
#' Choices are: \code{"input"}, \code{"missingness"}, \code{"strata"}.
#' Keeping the input order is recommended for the initial diagnostic plot
#' because it can reveal sequencing, plate, library-preparation, or sample-
#' processing structure.
#' Default: \code{sort.individuals = "input"}.
#'
#' @param sort.markers (character) Sorting method for markers.
#' Choices are: \code{"input"}, \code{"missingness"}, \code{"position"}.
#' Keeping the input order is recommended for the initial diagnostic plot.
#' Default: \code{sort.markers = "input"}.
#'
#' @param sample.max (optional, integer) Maximum number of individuals plotted.
#' Default: \code{sample.max = NULL}.
#'
#' @param marker.max (optional, integer) Maximum number of markers plotted.
#' Default: \code{marker.max = 50000}.
#'
#' @param filename Optional name for writing the heatmap directly to a PNG file
#' inside the function results folder.
#' When no \code{.png} extension is supplied, it is added automatically.
#' The plot is not displayed by this argument; it remains available in the
#' returned \code{heatmap} component. Default: \code{filename = NULL}.
#'
#' @param image.width Width of the PNG in pixels.
#' Default: \code{image.width = 1800}.
#'
#' @param image.height Height of the PNG in pixels.
#' Default: \code{image.height = 2400}.
#'
#' @param image.res PNG resolution in pixels per inch.
#' Default: \code{image.res = 150}.
#'
#' @param facet (logical) Should the heatmap be faceted by \code{strata.select}?
#' Default: \code{facet = TRUE}.
#'
#' @param parallel.core (optional, integer) Number of cores.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.
#' @inheritParams radr_common_arguments
#'
#' @param ... (optional) Further arguments passed to advanced mode.
#' Use \code{path.folder} to select the parent results directory.
#'
#' @details Missing genotypes are shown in black and observed genotypes in grey.
#' Broad vertical bands indicate individuals with elevated missingness, whereas
#' broad horizontal bands indicate markers with elevated missingness. Structure
#' restricted to a stratum can indicate a group- or batch-specific issue. These
#' patterns can guide whether individual or marker filtering should be examined
#' first; they are diagnostic evidence rather than an automatic filtering rule.
#' The heatmap is intended to reveal dataset-level patterns and guide filtering
#' strategy. It is not intended to identify specific markers or individuals;
#' consequently, individual and marker labels are deliberately omitted.
#'
#' For the first inspection, the defaults preserve the current GDS order for
#' both individuals and markers. This is intentional: ordering inherited from
#' sequencing, plates, libraries, or sample processing can expose technical
#' structure that would be obscured by sorting on missingness. The sorted modes
#' are best used as complementary views after examining the input-order plot.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{heatmap}: a \code{ggplot2} object containing the raster
#'   \item \code{ibm.matrix}: integer matrix, 0 = missing, 1 = genotyped
#'   \item \code{raster}: raster representation used by the heatmap
#'   \item \code{image.file}: path to the PNG, or \code{NULL} when no file was written
#'   \item \code{individuals.missingness}: tibble with per-individual missingness
#'   \item \code{markers.missingness}: tibble with per-marker missingness
#'   \item \code{plot.data}: \code{NULL}; retained for compatibility because the
#'   raster renderer does not generate a long plotting table
#' }
#'
#' @export
detect_ibm <- function(
    gds,
    strata = NULL,
    strata.select = "STRATA",
    sort.individuals = "input",
    sort.markers = "input",
    sample.max = NULL,
    marker.max = 50000,
    filename = NULL,
    image.width = 1800L,
    image.height = 2400L,
    image.res = 150,
    facet = TRUE,
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE,
    ...
) {

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "radr",
    f.name = "detect_ibm",
    verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)

  if (missing(gds)) rlang::abort("Input gds missing")

  rad.dots <- radr_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("path.folder", "internal"),
    deprecated = NULL,
    verbose = FALSE
  )

  # Results folder -------------------------------------------------------------
  path.folder <- tgbase::generate_folder(
    folder = "detect_ibm",
    path.folder = path.folder,
    internal = internal,
    file.date = file.date,
    prefix.int = TRUE,
    verbose = verbose
  )

  tgbase::write_tgbase_tsv(
    data = rad.dots,
    path.folder = path.folder,
    filename = "radr_detect_ibm_args",
    date = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  valid.sort.ind <- c("input", "missingness", "strata")
  valid.sort.mark <- c("input", "missingness", "position")

  if (!is.logical(facet) || length(facet) != 1L || is.na(facet)) {
    rlang::abort("facet must be TRUE or FALSE")
  }

  if (!is.null(sample.max) &&
      (!is.numeric(sample.max) || length(sample.max) != 1L ||
       is.na(sample.max) || sample.max < 1)) {
    rlang::abort("sample.max must be NULL or a positive number")
  }

  if (!is.null(marker.max) &&
      (!is.numeric(marker.max) || length(marker.max) != 1L ||
       is.na(marker.max) || marker.max < 1)) {
    rlang::abort("marker.max must be NULL or a positive number")
  }

  if (!is.null(sample.max)) sample.max <- as.integer(sample.max)
  if (!is.null(marker.max)) marker.max <- as.integer(marker.max)

  if (!is.null(filename) &&
      (!is.character(filename) || length(filename) != 1L ||
       is.na(filename) || !nzchar(filename))) {
    rlang::abort("filename must be NULL or a single non-empty path")
  }

  image.arguments <- list(
    image.width = image.width,
    image.height = image.height,
    image.res = image.res
  )
  invalid.image.argument <- vapply(
    image.arguments,
    function(x) {
      !is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x <= 0
    },
    logical(1)
  )
  if (any(invalid.image.argument)) {
    rlang::abort(paste0(
      names(image.arguments)[which(invalid.image.argument)[1L]],
      " must be a single positive number"
    ))
  }

  image.width <- as.integer(image.width)
  image.height <- as.integer(image.height)

  if (!sort.individuals %in% valid.sort.ind) {
    rlang::abort(
      message = stringi::stri_join(
        "sort.individuals must be one of: ",
        stringi::stri_join(valid.sort.ind, collapse = ", ")
      )
    )
  }

  if (!sort.markers %in% valid.sort.mark) {
    rlang::abort(
      message = stringi::stri_join(
        "sort.markers must be one of: ",
        stringi::stri_join(valid.sort.mark, collapse = ", ")
      )
    )
  }

  gds.opened.here <- FALSE

  if (inherits(gds, "character")) {
    gds <- SeqArray::seqOpen(gds)
    gds.opened.here <- TRUE
  }

  if (!inherits(gds, "SeqVarGDSClass")) {
    rlang::abort("gds must be a SeqVarGDSClass object or a path to a GDS file")
  }

  SeqArray::seqFilterPush(gds)
  on.exit({
    try(SeqArray::seqFilterPop(gds), silent = TRUE)
    if (gds.opened.here) {
      try(SeqArray::seqClose(gds), silent = TRUE)
    }
  }, add = TRUE)

  sample.ids <- SeqArray::seqGetData(gds, "sample.id")
  marker.ids <- SeqArray::seqGetData(gds, "variant.id")

  n.samples <- length(sample.ids)
  n.markers <- length(marker.ids)

  chrom <- rep(NA_character_, n.markers)
  pos <- rep(NA_integer_, n.markers)

  if (sort.markers == "position") {
    chrom <- tryCatch(
      as.character(SeqArray::seqGetData(gds, "chromosome")),
      error = function(e) rep(NA_character_, n.markers)
    )
    pos <- tryCatch(
      as.integer(SeqArray::seqGetData(gds, "position")),
      error = function(e) rep(NA_integer_, n.markers)
    )

    if (length(chrom) != n.markers) {
      chrom <- rep(NA_character_, n.markers)
    }
    if (length(pos) != n.markers) {
      pos <- rep(NA_integer_, n.markers)
    }
  }

  individuals.missing.count <- integer(n.samples)

  markers.missing.count <- SeqArray::seqApply(
    gds,
    "genotype",
    FUN = function(x) {
      missing <- colSums(is.na(x)) == nrow(x)
      individuals.missing.count <<- individuals.missing.count +
        as.integer(missing)
      sum(missing)
    },
    margin = "by.variant",
    as.is = "integer"
  )

  individuals.missingness <- tibble::tibble(
    INDIVIDUALS = sample.ids,
    MISSING_GENOTYPE = individuals.missing.count,
    MARKER_NUMBER = n.markers,
    MISSING_GENOTYPE_PROP = MISSING_GENOTYPE / MARKER_NUMBER,
    PERCENT = round(MISSING_GENOTYPE_PROP * 100, 2)
  )

  markers.missingness <- tibble::tibble(
    MARKERS = marker.ids,
    CHROM = chrom,
    POS = pos,
    MISSING_GENOTYPE = markers.missing.count,
    INDIVIDUALS_NUMBER = n.samples,
    MISSING_GENOTYPE_PROP = MISSING_GENOTYPE / INDIVIDUALS_NUMBER,
    PERCENT = round(MISSING_GENOTYPE_PROP * 100, 2)
  )

  if (is.null(strata)) {
    strata.df <- tibble::tibble(
      INDIVIDUALS = sample.ids,
      STRATA = "all"
    )
    if (strata.select != "STRATA") strata.df[[strata.select]] <- "all"
  } else {
    strata.df <- genometranslator::read_strata(
      strata = strata,
      verbose = FALSE
    )$strata

    if (!"INDIVIDUALS" %in% names(strata.df)) {
      rlang::abort("strata must contain an INDIVIDUALS column")
    }

    if (!strata.select %in% names(strata.df)) {
      rlang::abort(
        message = stringi::stri_join(
          "Column '", strata.select, "' not found in strata"
        )
      )
    }

    strata.df <- dplyr::semi_join(
      strata.df,
      tibble::tibble(INDIVIDUALS = sample.ids),
      by = "INDIVIDUALS"
    )

    missing.strata <- setdiff(sample.ids, strata.df$INDIVIDUALS)
    if (length(missing.strata) > 0L) {
      rlang::abort(
        message = paste0(
          "strata is missing ", length(missing.strata),
          " individual(s) present in the GDS"
        )
      )
    }
  }

  sample.order <- sample.ids
  marker.order <- marker.ids

  if (sort.individuals == "missingness") {
    sample.order <- individuals.missingness %>%
      dplyr::arrange(MISSING_GENOTYPE_PROP, INDIVIDUALS) %>%
      dplyr::pull(INDIVIDUALS)
  }

  if (sort.individuals == "strata") {
    sample.order <- strata.df %>%
      dplyr::select(INDIVIDUALS, dplyr::all_of(strata.select)) %>%
      dplyr::distinct() %>%
      dplyr::arrange(.data[[strata.select]], INDIVIDUALS) %>%
      dplyr::pull(INDIVIDUALS)
  }

  if (sort.markers == "missingness") {
    marker.order <- markers.missingness %>%
      dplyr::arrange(MISSING_GENOTYPE_PROP, MARKERS) %>%
      dplyr::pull(MARKERS)
  }

  if (sort.markers == "position" &&
      (all(is.na(chrom)) || all(is.na(pos)))) {
    rlang::abort(
      "sort.markers = 'position' requires chromosome and position information"
    )
  }

  if (sort.markers == "position") {
    marker.order <- markers.missingness %>%
      dplyr::arrange(CHROM, POS, MARKERS) %>%
      dplyr::pull(MARKERS)
  }

  if (!is.null(sample.max)) {
    sample.order <- sample.order[seq_len(min(sample.max, length(sample.order)))]
  }

  if (!is.null(marker.max)) {
    marker.order <- marker.order[seq_len(min(marker.max, length(marker.order)))]
  }

  SeqArray::seqSetFilter(
    gds,
    sample.id = sample.order,
    variant.id = marker.order,
    verbose = FALSE
  )

  ibm.blocks <- SeqArray::seqApply(
    gds,
    "genotype",
    FUN = function(x) {
      as.integer(!(colSums(is.na(x)) == nrow(x)))
    },
    margin = "by.variant",
    as.is = "list"
  )

  ibm.matrix <- do.call(cbind, ibm.blocks)
  rownames(ibm.matrix) <- SeqArray::seqGetData(gds, "sample.id")
  colnames(ibm.matrix) <- SeqArray::seqGetData(gds, "variant.id")

  # seqSetFilter() selects IDs but retains their physical GDS order. Reorder the
  # matrix explicitly so the requested diagnostic order is preserved.
  ibm.matrix <- ibm.matrix[
    match(as.character(sample.order), rownames(ibm.matrix)),
    match(as.character(marker.order), colnames(ibm.matrix)),
    drop = FALSE
  ]

  # Render the binary matrix directly. The first raster row is drawn at the top,
  # so marker rows are reversed to retain the orientation of the former ggplot
  # tile heatmap, where the first marker was at the bottom.
  raster.values <- t(ibm.matrix)
  raster.values <- raster.values[nrow(raster.values):1L, , drop = FALSE]
  raster.image <- grDevices::as.raster(
    ifelse(raster.values == 0L, "black", "grey")
  )

  n.samples.plot <- nrow(ibm.matrix)
  n.markers.plot <- ncol(ibm.matrix)
  legend.data <- tibble::tibble(
    x = 1,
    y = 1,
    Missingness = factor(
      c("genotyped", "missing"),
      levels = c("genotyped", "missing")
    )
  )

  heatmap <- ggplot2::ggplot(
    legend.data,
    ggplot2::aes(x = x, y = y, fill = Missingness)
  ) +
    ggplot2::annotation_raster(
      raster = raster.image,
      xmin = 0.5,
      xmax = n.samples.plot + 0.5,
      ymin = 0.5,
      ymax = n.markers.plot + 0.5,
      interpolate = FALSE
    ) +
    ggplot2::geom_point(shape = 22, alpha = 0) +
    ggplot2::scale_fill_manual(values = c("grey", "black")) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(alpha = 1, size = 5)
      )
    ) +
    ggplot2::labs(
      x = "Individuals",
      y = "Markers",
      fill = NULL
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0.5, n.samples.plot + 0.5),
      ylim = c(0.5, n.markers.plot + 0.5),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        size = 10, family = "Helvetica", face = "bold"
      ),
      axis.title.y = ggplot2::element_text(
        size = 10, family = "Helvetica", face = "bold"
      ),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = if (facet) 24 else 5.5)
    )

  if (facet) {
    strata.plot <- strata.df[
      match(rownames(ibm.matrix), strata.df$INDIVIDUALS),
      strata.select,
      drop = TRUE
    ]
    strata.runs <- rle(as.character(strata.plot))
    strata.ends <- cumsum(strata.runs$lengths)
    strata.starts <- c(1L, head(strata.ends, -1L) + 1L)
    strata.labels <- tibble::tibble(
      x = (strata.starts + strata.ends) / 2,
      y = n.markers.plot + 0.5,
      label = strata.runs$values
    )

    heatmap <- heatmap +
      ggplot2::geom_vline(
        xintercept = head(strata.ends, -1L) + 0.5,
        colour = "white",
        linewidth = 0.25,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = strata.labels,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        vjust = -0.6,
        size = 3.5,
        family = "Helvetica",
        fontface = "bold"
      )
  }

  image.file <- NULL
  if (!is.null(filename)) {
    tgbase::check_package(package = "ragg")

    if (!grepl("[.]png$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".png")
    }
    filename <- file.path(path.folder, basename(filename))

    ragg::agg_png(
      filename = filename,
      width = image.width,
      height = image.height,
      units = "px",
      res = image.res,
      background = "white"
    )
    tryCatch(
      print(heatmap),
      finally = grDevices::dev.off()
    )

    image.file <- normalizePath(filename, mustWork = TRUE)
    if (verbose) message("Heatmap written: ", basename(image.file))
  }

  return(list(
    heatmap = heatmap,
    ibm.matrix = ibm.matrix,
    raster = raster.image,
    image.file = image.file,
    individuals.missingness = individuals.missingness,
    markers.missingness = markers.missingness,
    plot.data = NULL
  ))
}#END detect_ibm
