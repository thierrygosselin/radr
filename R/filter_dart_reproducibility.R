# SNP number per haplotype
#' @name filter_dart_reproducibility
#' @title Filter data based on DArT reproducibility statistics
#' @description This filter removes markers below a certain threshold.
#' Based on the repoducibility column found in DArT files.
#'
#' **Filter target:** Markers.
#'
#' \strong{Statistics}: Reproducibility (established by DArT)


# Most arguments are inherited from tidy_genome
#' @inheritParams radr_common_arguments

#' @param filter.reproducibility (double, character) This is best decided after viewing the figures.
#' Usually values higher than 0.95 are not uncommon.
#' The value can also be character: \code{filter.reproducibility = "outliers"}.
#' Using this, will remove outlier markers using the lower outlier statistics.
#' Default: \code{filter.reproducibility = NULL}.


#' @rdname filter_dart_reproducibility
#' @export

#' @details
#' \strong{Interactive version}
#'
#' There are 2 steps in the interactive version to visualize and filter
#' the data based on the reproducibility value:
#'
#' Step 1. Visualization using a box plot
#'
#' Step 2. Choose the filtering threshold
#'
#'
#' @return A list in the global environment with 6 objects:
#' \enumerate{
#' \item $whitelist.markers
#' \item $blacklist.markers
#' \item $filters.parameters
#' }
#'
#' The object can be isolated in separate object outside the list by
#' following the example below.

#' @examples
#' \dontrun{
#' spotted.cod <- genometranslator::read_dart(
#'     data = "Combined_1514and1614_SNP_80Callrate.csv",
#'     strata = "strata.dart.spotted.cod.tsv"
#' )
#' turtle.filtered <- radr::filter_dart_reproducibility(
#'     data = spotted.cod,
#'     filter.reproducibility = 0.97
#' )
#' }

filter_dart_reproducibility <- function(
  data,
  interactive.filter = TRUE,
  filter.reproducibility = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  # interactive.filter <- TRUE
  # # data <- gds
  # path.folder <- NULL
  # parameters <- NULL
  # filename <- NULL
  # filter.reproducibility <- NULL
  # parallel.core <- parallel::detectCores() - 1
  # verbose = TRUE
  # internal <- FALSE


  # obj.keeper <- c(ls(envir = globalenv()), "data")

  if (interactive.filter || !is.null(filter.reproducibility)) {
    if (interactive.filter) verbose <- TRUE
    # Common startup -----------------------------------------------------------
    .start <- tgbase::startup(
      package = "radr",
      f.name = "filter_dart_reproducibility",
      verbose = verbose
    )
    file.date <- .start$file.date
    on.exit(tgbase::teardown(.start), add = TRUE)
    # on.exit(rm(list = setdiff(ls(envir = sys.frame(-1L)), obj.keeper), envir = sys.frame(-1L)))
    res <- list()

    # Function call and dotslist -----------------------------------------------
    rad.dots <- radr_dots(
      func.name = as.list(sys.call())[[1]],
      fd = rlang::fn_fmls_names(),
      args.list = as.list(environment()),
      dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
      keepers = c("path.folder", "parameters", "internal"),
      verbose = FALSE
    )

    # Checking for missing and/or default arguments ------------------------------
    if (missing(data)) rlang::abort("data is missing")

    # Folders---------------------------------------------------------------------
    path.folder <- generate_folder(
      rad.folder = "filter_dart_reproducibility",
      path.folder = path.folder,
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    # write the dots file
    tgbase::write_tgbase_tsv(
      data = rad.dots,
      path.folder = path.folder,
      filename = "radr_filter_dart_reproducibility_args",
      date = TRUE,
      internal = internal,
      write.message = "Function call and arguments stored in: ",
      verbose = verbose
    )

    # interactive.filter steps -------------------------------------------------
    if (interactive.filter) {
      message("\nInteractive mode: on")
      message("2 steps to visualize and filter the data based on reproducibility:")
      message("Step 1. Visualization")
      message("Step 2. Choose the filtering threshold\n\n")
    }

    # File type detection----------------------------------------------------------
    data.type <- genometranslator::detect_genomic_format(data)
    if (!data.type %in% c("tbl_df", "SeqVarGDSClass", "gds.file")) {
      rlang::abort("Input not supported for this function: read function documentation")
    }


    if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
      tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)

      if (data.type == "gds.file") {
        data <- genometranslator::read_genome(data, verbose = verbose)
        data.type <- "SeqVarGDSClass"
      }
    } else {
      if (is.vector(data)) data <- genometranslator::read_genome(data = data, import.metadata = TRUE)
      data.type <- "tbl_df"
    }

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

    # Whitelist and blacklist --------------------------------------------------
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "COL", "REP_AVG")
    if (data.type == "SeqVarGDSClass") {
      markers.meta <- bl <- extract_markers_metadata(gds = data, whitelist = FALSE)
    } else {
      markers.meta <- bl <- dplyr::select(data, tidyselect::any_of(want))
    }
    # Check that required info is present in data:
    if (!tibble::has_name(markers.meta, "REP_AVG")) {
      message("This filter requires REP_AVG info, skipping filtering...")
      return(data)
    }

    # Generating statistics ----------------------------------------------------
    rep.stats <- tgbase::tibble_stats(x = markers.meta$REP_AVG, group = "reproducibility")
    readr::write_tsv(x = rep.stats, file = file.path(path.folder, "dart_reproducibility_stats.tsv"))
    if (verbose) message("File written: dart_reproducibility_stats.tsv")

    # Generate box plot ---------------------------------------------------------
    # Previously using IQR and I think it's an error that was introduced.
    # pretty sure it should be the outlier low threshold..
    # outlier.rep <- ceiling(rep.stats[[8]] * 1000) / 1000
    outlier.rep <- ceiling(rep.stats[["OUTLIERS_LOW"]] * 1000) / 1000
    bp.filename <- stringi::stri_join("dart_reproducibility_boxplot_", file.date, ".pdf")
    rep.fig <- tgbase::boxplot_stats(
      data = rep.stats,
      title =  "DArT marker's reproducibility",
      subtitle = stringi::stri_join("\nLower outlier: ", outlier.rep),
      x.axis.title = NULL,
      y.axis.title = "Reproducibility (proportion)",
      bp.filename = bp.filename,
      path.folder = path.folder)
    if (verbose) message("File written: ", bp.filename)

    # Helper table -------------------------------------------------------------
    if (verbose) message("Generating helper table...")
    how_many_markers <- function(threshold, x) {
      nrow(dplyr::filter(x, REP_AVG >= threshold))
    }#End how_many_markers

    snp.range <- seq(0.9, 1, 0.005)
    n.markers <- nrow(markers.meta)
    helper.table <- tibble::tibble(REPRODUCIBILITY = snp.range) %>%
      dplyr::mutate(
        WHITELISTED_MARKERS = purrr::map_int(.x = snp.range, .f = how_many_markers, x = markers.meta),
        BLACKLISTED_MARKERS = n.markers - WHITELISTED_MARKERS
      ) %>%
      readr::write_tsv(
        x = .,
        file = file.path(path.folder, "dart.reproducibility.helper.table.tsv"))

    # figures
    markers.plot <- ggplot2::ggplot(
      data = helper.table  %<>%
        tidyr::pivot_longer(
          data = .,
          cols = -REPRODUCIBILITY,
          names_to = "LIST",
          values_to = "MARKERS"
        ),
      ggplot2::aes(x = REPRODUCIBILITY, y = MARKERS)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
      ggplot2::geom_vline(ggplot2::aes(xintercept = as.numeric(outlier.rep)), color = "yellow") +
      ggplot2::scale_x_continuous(name = "Reproducibility (proportion)", breaks = snp.range) +
      ggplot2::scale_y_continuous(name = "Number of markers") +
      ggplot2::theme_bw()+
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", angle = 90, hjust = 1, vjust = 0.5)
      ) +
      ggplot2::facet_grid(LIST ~. , scales = "free", space = "free")
    print(markers.plot)

    # save
    ggplot2::ggsave(
      filename = file.path(path.folder, "dart.reproducibility.helper.plot.pdf"),
      plot = markers.plot,
      width = 20,
      height = 15,
      dpi = 300,
      units = "cm",
      limitsize = FALSE,
      useDingbats = FALSE
      )

    helper.table <- markers.plot <- NULL
    if (verbose) message("Files written: helper tables and plots")


    # Step 2. Thresholds selection ---------------------------------------------
    if (interactive.filter) {
      message("\nStep 2. Filtering markers based on markers reproducibility\n")

      filter.reproducibility <- tgbase::question(
        x = "Do you still want to blacklist markers? (y/n):",
        answer.opt = c("y", "n"))

      if (filter.reproducibility == "y") {
        message("2 options to blacklist markers based on reproducibility:")
        message("1. use the outlier statistic")
        message("2. enter your own threshold")
        outlier.stats <- tgbase::question(
          x = "Enter the option (1 or 2): ", minmax = c(1, 2))
        if (outlier.stats == 1) {
          filter.reproducibility <- "outliers"
        } else {
          filter.reproducibility <- tgbase::question(
            x = "Enter the proportion threshold (0-1)
The minimum reproducibility tolerated:", minmax = c(0, 1))
        }
        outlier.stats <- NULL
      } else {
        filter.reproducibility <- NULL
      }
    }


    # Filtering ----------------------------------------------------------------
    if (!is.null(filter.reproducibility)) {
      if (!purrr::is_double(filter.reproducibility)) {
        message("\nRemoving outlier markers based on reproducibility statistic: ", outlier.rep)
        filter.reproducibility <- outlier.rep
      } else {
        message("\nRemoving markers based on reproducibility statistic: ", filter.reproducibility)
      }



      if (data.type == "SeqVarGDSClass") {
        markers.meta %<>%
          dplyr::mutate(
            FILTERS = dplyr::if_else(REP_AVG < filter.reproducibility,
                                     "filter.dart.reproducibility", FILTERS)
          )
        wl <- markers.meta %>% dplyr::filter(FILTERS == "whitelist")
        # Update GDS
        genometranslator::update_genome_gds(
          gds = data,
          node.name = "markers.meta",
          value = markers.meta,
          sync = TRUE
        )

      } else {
        # Apply the filter to the tidy data
        wl <- dplyr::filter(markers.meta, REP_AVG >= filter.reproducibility)
        data  %<>% dplyr::filter(MARKERS %in% wl$MARKERS)

      }
      # Whitelist and Blacklist of markers
      readr::write_tsv(
        x = wl,
        file = file.path(path.folder, "whitelist.dart.reproducibility.tsv"),
        append = FALSE, col_names = TRUE)

      bl %<>% dplyr::setdiff(wl) %>%
        readr::write_tsv(
          x = .,
          file = file.path(path.folder, "blacklist.dart.reproducibility.tsv"),
          append = FALSE, col_names = TRUE)
      # saving whitelist and blacklist
      if (verbose) message("File written: whitelist.dart.reproducibility.tsv")
      if (verbose) message("File written: blacklist.dart.reproducibility.tsv")



      # Update parameters --------------------------------------------------------
      filters.parameters <- filter_parameters(
        generate = FALSE,
        initiate = FALSE,
        update = TRUE,
        parameter.obj = filters.parameters,
        data = data,
        filter.name = "Filter DArT reproducibility",
        param.name = "filter.reproducibility",
        values = filter.reproducibility,
        path.folder = path.folder,
        file.date = file.date,
        internal = internal,
        verbose = verbose
      )

      # Return -----------------------------------------------------------------------
      radr_results_message(
        rad.message = stringi::stri_join("\nFilter DArT reproducibility: ",
                                         filter.reproducibility),
        filters.parameters,
        internal,
        verbose
      )
    }
    return(data)
  } else {
    return(data)
  }
} #End filter_dart_reproducibility
