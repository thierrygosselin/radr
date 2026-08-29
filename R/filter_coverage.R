# Filter markers coverage
#' @name filter_coverage
#' @title Filter markers mean coverage
#' @description This function is designed to remove/blacklist markers
#' based on mean coverage information.
#'
#' **Filter target:** Markers.
#'
#' \strong{Statistics}: Mean marker coverage. Genotype read depth is averaged
#' across active individuals for each marker.

#' @param filter.coverage (optional, string) 2 options:
#' \itemize{
#' \item character string \code{filter.coverage = "outliers"} will use as
#' thresholds the lower and higher outlier values in the box plot.
#' \item numeric vector \code{filter.coverage = c(10, 200)} for the
#' marker mean-coverage lower and upper bounds.
#' }
#' Default: \code{filter.coverage = NULL}.

#' @param filename (optional, character)
#' Default: \code{filename = NULL}.

#' @inheritParams radr_common_arguments

#' @section Advanced mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#' \item \code{filter.common.markers} (optional, logical).
#' Default: \code{filter.common.markers = FALSE},
#' Documented in \code{\link{filter_common_markers}}.
#' \item \code{filter.monomorphic} (logical, optional) Should the monomorphic
#' markers present in the dataset be filtered out ?
#' Default: \code{filter.monomorphic = TRUE}.
#' Documented in \code{\link{filter_monomorphic}}.
#' \item \code{path.folder}: to write ouput in a specific path
#' (used internally in radr).
#' Default: \code{path.folder = getwd()}.
#' If the supplied directory doesn't exist, it's created.
#' }


#' @section Interactive version:
#'
#' The interactive mode first calculates marker coverage, writes and displays
#' the coverage distribution and helper plots, and then asks:
#' \enumerate{
#' \item \code{"Choose the min mean coverage threshold (e.g. 7 or 10):"}
#' \item \code{"Choose the max mean coverage threshold (e.g. 100 or 300):"}
#' }
#' Markers outside the inclusive interval are blacklisted. Use
#' \code{interactive.filter = FALSE} with an explicit two-value
#' \code{filter.coverage} for a reproducible analysis.


#' @rdname filter_coverage
#' @export

#' @return The filtered data in the same representation as the input. GDS
#' marker metadata and active variants are updated in place, so the underlying
#' GDS file is modified. Coverage tables, figures, marker lists, and filtering
#' parameters are written to the function output folder when applicable.

#' @examples
#' \dontrun{
#' genome <- genometranslator::read_genome("my_genome.gds")
#'
#' # Inspect marker coverage and choose lower and upper limits interactively.
#' genome <- radr::filter_coverage(data = genome)
#'
#' # Alternatively, start from a separate unfiltered GDS for a scripted run.
#' scripted_genome <- genometranslator::read_genome("my_genome_scripted.gds")
#' scripted_genome <- radr::filter_coverage(
#'   data = scripted_genome,
#'   interactive.filter = FALSE,
#'   filter.coverage = c(10, 200)
#' )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

filter_coverage <- function(
    data,
    interactive.filter = TRUE,
    filter.coverage = NULL,
    filename = NULL,
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE,
    ...
) {
  # Force piped input before printing this function's startup banner.
  force(data)

  # interactive.filter = TRUE
  # data = gds
  # filter.coverage = NULL
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE
  # path.folder = NULL
  # parameters <- NULL
  # force.stats <- TRUE
  # internal = FALSE
  if (is.null(filter.coverage) && !interactive.filter) return(data)
  if (interactive.filter) verbose <- TRUE

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "radr",
    f.name = "filter_coverage",
    verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)
  # on.exit(rm(list = setdiff(ls(envir = sys.frame(-1L)), obj.keeper), envir = sys.frame(-1L)))

  # Function call and dotslist -----------------------------------------------
  rad.dots <- radr_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("path.folder", "parameters", "force.stats", "internal"),
    verbose = FALSE
  )

  # Checking for missing and/or default arguments ----------------------------
  if (missing(data)) rlang::abort("data is missing")

  # Folders---------------------------------------------------------------------
  path.folder <- generate_folder(
    rad.folder = "filter_coverage",
    path.folder = path.folder,
    internal = internal,
    file.date = file.date,
    verbose = verbose)

  # write the dots file
  tgbase::write_tgbase_tsv(
    data = rad.dots,
    path.folder = path.folder,
    filename = "radr_filter_coverage_args",
    date = TRUE,
    internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  # Message about steps taken during the process -----------------------------
  if (interactive.filter) {
    message("Interactive mode: on\n")
    message("Step 1. Visualization and helper table")
    message("Step 2. Filtering markers based on mean coverage\n\n")
  }


  # Import data ---------------------------------------------------------------
  if (verbose) cli::cli_progress_step("Importing data ")
  # Detect format
  tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)
  data.type <- genometranslator::detect_genomic_format(data)
  if (!data.type %in% c("SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }
  if (data.type == "gds.file") {
    data <- genometranslator::read_genome(data, verbose = verbose)
    data.type <- "SeqVarGDSClass"
  }
  cli::cli_progress_done()

  # Filter parameter file: initiate ------------------------------------------
  filters.parameters <- filter_parameters(
    generate = TRUE,
    initiate = TRUE,
    update = FALSE,
    parameter.obj = parameters,
    data = data,
    path.folder = path.folder,
    file.date = file.date,
    internal = internal,
    verbose = verbose)

  # Verify that coverage information is present in the data...
  depth.info <- check_coverage(gds = data, genotypes.metadata.check = TRUE, stacks.haplo.check = TRUE, dart.check = TRUE)
  if (is.null(depth.info)) {
    message("\n\nCoverate information is not available for this dataset, returning GDS...")
    return(data)
  }


  # Step 1. Visuals ----------------------------------------------------------
  if (interactive.filter) message("\nStep 1. Coverage visualization and helper table\n")

  # Generate coverage stats-----------------------------------------------------
  info <- generate_stats(
    gds = data,
    individuals = FALSE,
    markers = TRUE,
    missing = FALSE,
    coverage = TRUE,
    allele.coverage = FALSE,
    mac = FALSE,
    heterozygosity = FALSE,
    snp.per.locus = FALSE,
    snp.position.read = FALSE,
    force.stats = force.stats,
    plot = FALSE,
    path.folder = path.folder,
    file.date = file.date,
    parallel.core = parallel.core
  )

  stats <- info$m.stats
  info <- info$m.info


  # identify outliers: low and high -----------------------------------------
  out.low <- floor(stats$OUTLIERS_LOW[stats$GROUP == "mean coverage"]*1000)/1000
  out.high <- floor(stats$OUTLIERS_HIGH[stats$GROUP == "mean coverage"]*1000)/1000
  if (verbose) message("Generating coverage statistics: without outliers")

  # Helper table -------------------------------------------------------------
  # filter.coverage <- "outliers"
  # filter.coverage <- c(10, 100)
  if (!is.null(filter.coverage)) {
    if (length(filter.coverage) > 1) {
      combined.info <- TRUE
      cov.low <- filter.coverage[1]
      cov.high <- filter.coverage[2]
    } else {
      if (is.character(filter.coverage)) {
        combined.info <- FALSE
        cov.low <- out.low
        cov.high <- out.high
        filter.coverage <- c(cov.low, cov.high)
      } else {
        rlang::abort("Unknown mean coverage thresholds used")
      }
    }
  } else {
    combined.info <- FALSE
  }

  # combined.info <- TRUE
  # cov.low <- 10
  # cov.high <- 100
  min.c <- floor(stats$MIN[stats$GROUP == "mean coverage"]*1000)/1000
  med.c <- floor(stats$MEDIAN[stats$GROUP == "mean coverage"]*1000)/1000
  max.c <- floor(stats$MAX[stats$GROUP == "mean coverage"]*1000)/1000

  if (combined.info) {
    # include outliers
    max.out <- floor(stats$OUTLIERS_HIGH[stats$GROUP == "mean coverage"]*1000)/1000

    if (out.low != cov.low) {
      cl.range <- min(out.low,cov.low):max(out.low,cov.low)
    } else {
      cl.range <- out.low:med.c
    }

    if (out.high != cov.high) {
      ch.range <- min(out.high,cov.high):max(out.high,cov.high)
    } else {
      ch.range <- out.high
    }

  } else {
    # use min and max values with outliers
    #MIN
    if (min.c != out.low) {
      cl.range <- out.low:min.c
    } else {
      cl.range <- min.c:med.c
    }

    #MAX
    if (max.c != out.high) {
      ch.range <- out.high:max.c
    } else {
      ch.range <- max.c
    }
  }

  # reduce size of points...
  if (length(cl.range) > 100) cl.range <- seq(min(cl.range), max(cl.range), by = 10)
  if (length(ch.range) > 100) ch.range <- seq(min(ch.range), max(ch.range), by = 10)


  if (verbose) message("Generating mean coverage helper table...")
  how_many_markers <- function(threshold, x, low = FALSE) {
    if (low) {
      nrow(dplyr::filter(x, COVERAGE_MEAN >= threshold))
    } else {
      nrow(dplyr::filter(x, COVERAGE_MEAN <= threshold))
    }
  }#End how_many_markers

  n.markers <- nrow(info)
  helper.table.low <- tibble::tibble(COVERAGE_LOW = cl.range) %>%
    dplyr::mutate(
      WHITELISTED_MARKERS = purrr::map_int(.x = cl.range, .f = how_many_markers, x = info, low = TRUE),
      BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
    ) %>%
    readr::write_tsv(
      x = .,
      file = file.path(path.folder, "coverage.low.helper.table.tsv"))

  if (nrow(helper.table.low) > 1) {
    markers.plot.low <- radr_helper_plot(
      data = tgbase::trans_long(
        x = helper.table.low,
        cols = "COVERAGE_LOW",
        names_to = "LIST",
        values_to = "MARKERS"
      ),
      stats = "COVERAGE_LOW",
      x.axis.title =  "Minimum mean coverage allowed",
      x.breaks = cl.range,
      plot.filename = file.path(path.folder, "coverage.low.helper.plot")
    )
    print(markers.plot.low)
  }

  helper.table.high <- tibble::tibble(COVERAGE_HIGH = ch.range) %>%
    dplyr::mutate(
      WHITELISTED_MARKERS = purrr::map_int(.x = ch.range, .f = how_many_markers, x = info),
      BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
    ) %>%
    readr::write_tsv(
      x = .,
      file = file.path(path.folder, "coverage.high.helper.table.tsv"))

  if (nrow(helper.table.high) > 1) {
    markers.plot.high <- radr_helper_plot(
      data = tgbase::trans_long(
        x = helper.table.high,
        cols = "COVERAGE_HIGH",
        names_to = "LIST",
        values_to = "MARKERS"
      ),
      stats = "COVERAGE_HIGH",
      x.axis.title =  "Maximum mean coverage allowed",
      x.breaks = floor(ch.range),
      text.orientation = "vertical",
      plot.filename = file.path(path.folder, "coverage.high.helper.plot")
    )
    print(markers.plot.high)
  }
  if (verbose) message("Files written: helper tables and plots")


  # Step 2. Thresholds selection ---------------------------------------------
  if (interactive.filter) {
    filter.coverage <- c(min.c, max.c)
    message("\nStep 2. Filtering markers based on mean coverage\n")
    filter.coverage[1] <- tgbase::question(
      x = "Choose the min mean coverage threshold(e.g. 7 or 10): ", minmax = c(1, 10000))
  }
  if (interactive.filter) {
    filter.coverage[2] <- tgbase::question(
      x = "Choose the max mean coverage threshold (e.g. 100 or 300): ", minmax = c(1, 10000))
  }

  # Whitelist and blacklist --------------------------------------------------
  # want <- c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS")
  bl <- info %>%
    dplyr::filter(
      COVERAGE_MEAN < filter.coverage[1] |
        COVERAGE_MEAN > filter.coverage[2]
    ) %$% VARIANT_ID
  markers.meta <- extract_markers_metadata(gds = data, whitelist = FALSE) %>%
    dplyr::mutate(
      FILTERS = dplyr::if_else(
        VARIANT_ID %in% bl, "filter.mean.coverage", FILTERS
      )
    )
  # Update GDS
  genometranslator::update_genome_gds(
    gds = data,
    node.name = "markers.meta",
    value = markers.meta,
    sync = TRUE
  )


  tgbase::write_tgbase_tsv(
    data = markers.meta %>% dplyr::filter(FILTERS == "filter.mean.coverage"),
    path.folder = path.folder,
    filename = "blacklist.markers.coverage",
    date = TRUE,
    internal = internal,
    write.message = "standard",
    verbose = verbose
  )

  tgbase::write_tgbase_tsv(
    data = markers.meta %>% dplyr::filter(FILTERS == "whitelist"),
    path.folder = path.folder,
    filename = "whitelist.markers.coverage",
    date = TRUE,
    internal = internal,
    write.message = "standard",
    verbose = verbose
  )

  # Update parameters --------------------------------------------------------
  filters.parameters <- filter_parameters(
    generate = FALSE,
    initiate = FALSE,
    update = TRUE,
    parameter.obj = filters.parameters,
    data = data,
    filter.name = "Filter coverage min / max",
    param.name = "filter.coverage",
    values = paste(filter.coverage, collapse = " / "),
    path.folder = path.folder,
    file.date = file.date,
    internal = internal,
    verbose = verbose)

  # results ------------------------------------------------------------------
  radr_results_message(
    rad.message = stringi::stri_join("\nFilter mean coverage thresholds: ",
                                     paste(filter.coverage, collapse = " / ")),
    filters.parameters,
    internal,
    verbose
  )

  return(data)
}#End filter_coverage
