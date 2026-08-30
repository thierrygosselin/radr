# radr package startup message ---------------------------------------------
# .onAttach <- function(libname, pkgname) {
#   radr.version <- utils::packageDescription("radr", fields = "Version")
#   radr.build <- utils::packageDescription("radr", fields = "Built")
#   startup.message <- stringi::stri_join(
#     "******************************* IMPORTANT NOTICE *******************************\n",
#     "radr v.", radr.version, " was modified heavily.\n",
#     "Read functions documentation and available vignettes.\n\n",
#     "For reproducibility:\n",
#     "    radr version: ", radr.version,"\n",
#     "    radr build date: ", radr.build,"\n",
#     "    Keep zenodo DOI.\n",
#     "********************************************************************************",
#     sep = "")
#   packageStartupMessage(startup.message)
# }

# radr_function_header -----------------------------------------------------
#' @title radr_function_header
#' @description Generate function header
#' @rdname radr_function_header
#' @keywords internal
#' @export
radr_function_header <- function(f.name = NULL, start = TRUE, verbose = TRUE) {
  if (is.null(f.name)) invisible(NULL)
  if (start) {
    if (verbose) {
      message(stringi::stri_pad_both(str = "", width = 80L, pad = "#"))
      message(stringi::stri_pad_both(str = paste0(" radr::", f.name, " "), width = 80L, pad = "#"))
      message(stringi::stri_pad_both(str = "", width = 80L, pad = "#"), "\n")
    }
  } else {
    if (verbose) {
      message(stringi::stri_pad_both(str = paste0(" completed ", f.name, " "), width = 80L, pad = "#"), "\n")
    }
  }
}# End radr_function_header

# radr_startup--------------------------------------------------------------
#' @title Common startup helper for radr functions
#'
#' @description
#' Perform standard radr initialisation steps inside a function:
#' \itemize{
#'   \item print a header with \code{tgbase::function_header()};
#'   \item record execution date/time;
#'   \item record and temporarily modify global options
#'     (e.g. \code{width}, \code{future.globals.maxSize});
#'   \item start a timing object with \code{tgbase::tic()}.
#' }
#'
#' The function returns a small list with everything needed for a matching
#' teardown helper. You should typically pair this with
#' \code{radr_teardown()} inside an \code{on.exit()} call in the calling
#' function.
#'
#' @param f.name (character) Name of the calling function
#'   (e.g. \code{"read_vcf"}). Used only for logging in
#'   \code{tgbase::function_header()}.
#'
#' @param verbose (logical) When \code{TRUE}, the helper prints the execution
#'   date/time and passes verbosity to the teardown helper.
#'   Default: \code{verbose = TRUE}.
#'
#' @param width (integer) Temporary value for \code{getOption("width")}.
#'   Default: \code{width = 70}.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{file.date} – character timestamp \code{"YYYYMMDD@HHMM"};
#'   \item \code{old.dir} – working directory to restore;
#'   \item \code{opt.width} – original \code{options("width")};
#'   \item \code{opt.future} – original \code{options("future.globals.maxSize")};
#'   \item \code{timing} – timer object from \code{tgbase::tic()};
#'   \item \code{f.name} – function name;
#'   \item \code{verbose} – verbosity flag.
#' }
#'
#' @keywords internal
#' @export
radr_startup <- function(f.name, verbose = TRUE, width = 70L) {

  # Standard header
  tgbase::function_header(f.name = f.name, verbose = verbose)

  # Execution timestamp
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (isTRUE(verbose)) {
    message("Execution date@time: ", file.date)
  }

  # Save current state
  old.dir    <- getwd()
  opt.width  <- getOption("width")
  opt.future <- getOption("future.globals.maxSize")

  # Temporary options
  options(future.globals.maxSize = Inf)
  options(width = width)

  # Start timer
  timing <- tgbase::tic()

  # Return everything needed for teardown
  list(
    file.date  = file.date,
    old.dir    = old.dir,
    opt.width  = opt.width,
    opt.future = opt.future,
    timing     = timing,
    f.name     = f.name,
    verbose    = verbose
  )
} # End radr_startup


# radr_teardown-------------------------------------------------------------
#' @title Common teardown helper for radr functions
#'
#' @description
#' Companion to \code{radr_startup()}. Restores working directory and
#' options, stops the timer, and prints a closing header.
#'
#' Intended to be called from \code{on.exit()} in the calling function.
#'
#' @param start.obj (list) The object returned by \code{radr_startup()}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects.
#'
#' @keywords internal
#' @export
radr_teardown <- function(start.obj) {

  # Be defensive about missing elements
  if (!is.null(start.obj$old.dir)) {
    try(setwd(start.obj$old.dir), silent = TRUE)
  }

  if (!is.null(start.obj$opt.width)) {
    options(width = start.obj$opt.width)
  }

  if (!is.null(start.obj$opt.future)) {
    options(future.globals.maxSize = start.obj$opt.future)
  }

  # Stop timer (if present)
  if (!is.null(start.obj$timing)) {
    tgbase::toc(start.obj$timing, verbose = isTRUE(start.obj$verbose))
  }

  # Closing header
  tgbase::function_header(
    f.name  = start.obj$f.name %||% "radr_function",
    start   = FALSE,
    verbose = isTRUE(start.obj$verbose)
  )

  invisible(NULL)
} # End radr_teardown



# radr_question ---------------------------------------------------------
#' @title radr_question
#' @description Ask to enter a word or number
#' @rdname radr_question
#' @keywords internal
#' @export
radr_question <- function(x, answer.opt = NULL, minmax = NULL) {
  .Deprecated("tgbase::question", package = "radr")
  tgbase::question(x = x, answer.opt = answer.opt, minmax = minmax)
}#End radr_question


# filter_parameters -----------------------------------------------------------
#' Track radr filtering parameters
#'
#' A package-facing wrapper around
#' [genometranslator::genome_parameters()] that maintains the common genomic
#' operation history during radr filtering. For GDS input, every completed
#' operation also writes an identifier-level audit beside the usual parameter
#' file: complete marker and individual tables for records newly removed and
#' records still active, plus `filter_audit_manifest.tsv`. Consequently this
#' audit also covers state-changing functions whose names do not start with
#' `filter_`, including [detect_duplicate_genomes()] and
#' [detect_mixed_genomes()].
#'
#' @param ... Arguments passed to [genometranslator::genome_parameters()].
#' @rdname filter_parameters
#' @export
#' @keywords internal
filter_parameters <- function(...) {
  args <- list(...)
  parameter.obj <- args$parameter.obj
  before <- if (!is.null(parameter.obj)) parameter.obj$.radr_filter_audit else NULL

  res <- do.call(genometranslator::genome_parameters, args)
  if (is.null(res)) return(NULL)

  data <- args$data
  is.gds <- inherits(data, "SeqVarGDSClass")
  if (is.gds && isTRUE(args$update) && !is.null(before)) {
    after <- .capture_filter_audit(data)
    .write_filter_audit(
      before = before,
      after = after,
      path.folder = args$path.folder,
      filter.name = args$filter.name,
      param.name = args$param.name,
      values = args$values,
      verbose = isTRUE(args$verbose)
    )
    res$.radr_filter_audit <- after
  } else if (is.gds && (isTRUE(args$generate) || isTRUE(args$initiate))) {
    res$.radr_filter_audit <- .capture_filter_audit(data)
  } else if (!is.null(before)) {
    res$.radr_filter_audit <- before
  }
  res
}

# Capture both the complete metadata and the currently active identifiers. This
# is deliberately kept in the transient parameter object rather than the GDS:
# the persistent GDS history stores compact operation summaries, while the
# potentially large identifier-level audit is written beside each operation.
.capture_filter_audit <- function(data) {
  markers <- genometranslator::extract_markers_metadata(
    gds = data,
    whitelist = FALSE
  )
  individuals <- genometranslator::extract_individuals_metadata(
    gds = data,
    whitelist = FALSE
  )
  marker.fallback <- if ("MARKERS" %in% names(markers)) {
    as.character(markers$MARKERS)
  } else {
    key.fields <- intersect(c("CHROM", "LOCUS", "POS"), names(markers))
    if (!length(key.fields)) as.character(seq_len(nrow(markers))) else
      do.call(paste, c(markers[key.fields], sep = ":"))
  }
  marker.key <- if (
    "VARIANT_ID" %in% names(markers) &&
      any(!is.na(markers$VARIANT_ID) & nzchar(as.character(markers$VARIANT_ID)))
  ) {
    variant.id <- as.character(markers$VARIANT_ID)
    valid <- !is.na(variant.id) & nzchar(variant.id)
    ifelse(valid, paste0("VARIANT_ID:", variant.id), paste0("MARKER:", marker.fallback))
  } else {
    paste0("MARKER:", marker.fallback)
  }
  active.markers <- if ("FILTERS" %in% names(markers)) {
    marker.key[is.na(markers$FILTERS) | markers$FILTERS == "whitelist"]
  } else marker.key
  active.individuals <- if ("FILTERS" %in% names(individuals)) {
    as.character(individuals$INDIVIDUALS[
      is.na(individuals$FILTERS) | individuals$FILTERS == "whitelist"
    ])
  } else as.character(individuals$INDIVIDUALS)

  list(
    markers = markers,
    marker.key = marker.key,
    active.markers = unique(active.markers),
    individuals = individuals,
    active.individuals = unique(active.individuals)
  )
}

.filter_audit_slug <- function(x) {
  x <- tolower(if (is.null(x) || !length(x) || is.na(x[[1L]])) "operation" else x[[1L]])
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

.write_filter_audit <- function(
    before,
    after,
    path.folder = NULL,
    filter.name = NULL,
    param.name = NULL,
    values = NULL,
    verbose = TRUE
) {
  if (is.null(path.folder)) path.folder <- getwd()
  dir.create(path.folder, recursive = TRUE, showWarnings = FALSE)
  manifest.path <- file.path(path.folder, "filter_audit_manifest.tsv")
  audit.index <- if (file.exists(manifest.path)) {
    nrow(readr::read_tsv(manifest.path, show_col_types = FALSE, progress = FALSE)) + 1L
  } else 1L
  stem <- sprintf("audit_%02d_%s", audit.index, .filter_audit_slug(filter.name))

  marker.removed <- setdiff(before$active.markers, after$active.markers)
  individual.removed <- setdiff(before$active.individuals, after$active.individuals)
  marker.kept.rows <- after$markers[after$marker.key %in% after$active.markers, , drop = FALSE]
  marker.removed.rows <- before$markers[before$marker.key %in% marker.removed, , drop = FALSE]
  individual.kept.rows <- after$individuals[
    as.character(after$individuals$INDIVIDUALS) %in% after$active.individuals,
    , drop = FALSE
  ]
  individual.removed.rows <- before$individuals[
    as.character(before$individuals$INDIVIDUALS) %in% individual.removed,
    , drop = FALSE
  ]

  files <- c(
    markers.kept = paste0(stem, "_markers_kept.tsv"),
    markers.removed = paste0(stem, "_markers_removed.tsv"),
    individuals.kept = paste0(stem, "_individuals_kept.tsv"),
    individuals.removed = paste0(stem, "_individuals_removed.tsv")
  )
  readr::write_tsv(marker.kept.rows, file.path(path.folder, files[["markers.kept"]]))
  readr::write_tsv(marker.removed.rows, file.path(path.folder, files[["markers.removed"]]))
  readr::write_tsv(individual.kept.rows, file.path(path.folder, files[["individuals.kept"]]))
  readr::write_tsv(individual.removed.rows, file.path(path.folder, files[["individuals.removed"]]))

  manifest <- tibble::tibble(
    AUDIT_INDEX = audit.index,
    FILTER = if (is.null(filter.name)) "" else as.character(filter.name),
    PARAMETER = if (is.null(param.name)) "" else as.character(param.name),
    VALUE = if (is.null(values)) "" else paste(values, collapse = " / "),
    MARKERS_BEFORE = length(before$active.markers),
    MARKERS_REMOVED = length(marker.removed),
    MARKERS_KEPT = length(after$active.markers),
    INDIVIDUALS_BEFORE = length(before$active.individuals),
    INDIVIDUALS_REMOVED = length(individual.removed),
    INDIVIDUALS_KEPT = length(after$active.individuals),
    MARKERS_KEPT_FILE = files[["markers.kept"]],
    MARKERS_REMOVED_FILE = files[["markers.removed"]],
    INDIVIDUALS_KEPT_FILE = files[["individuals.kept"]],
    INDIVIDUALS_REMOVED_FILE = files[["individuals.removed"]]
  )
  readr::write_tsv(
    manifest,
    manifest.path,
    append = file.exists(manifest.path),
    col_names = !file.exists(manifest.path)
  )
  if (verbose) {
    message(
      "Filter audit written: ", length(marker.removed), " marker(s) and ",
      length(individual.removed), " individual(s) newly removed; ",
      length(after$active.markers), " marker(s) and ",
      length(after$active.individuals), " individual(s) kept."
    )
  }
  invisible(manifest)
}

# radr_results_message ------------------------------------------------------------
#' @title radr_results_message
#' @description Message printed at the end of most radr functions
#' @keywords internal
#' @export
radr_results_message <- function(
  rad.message = NULL,
  filters.parameters,
  internal = FALSE,
  verbose = TRUE
) {
  if (!internal) {
    # if (verbose) cat("################################### RESULTS ####################################\n")
    if (!is.null(rad.message)) message(rad.message)
    message("Number of individuals / strata / chrom / locus / SNP:")
    if (verbose) message("    Before: ", filters.parameters$filters.parameters$BEFORE)
    message("    Blacklisted: ", filters.parameters$filters.parameters$BLACKLIST)
    if (verbose) message("    After: ", filters.parameters$filters.parameters$AFTER)
  }
}#End radr_results_message

# radr_folder--------------------------------------------------------------------
#' @title radr_folder
#' @description Generate the rad folders
#' @param path.folder path of the folder
#' @param prefix.int Use an integer prefix padded left with 0.
#' Default: \code{prefix.int = TRUE}.
#' @keywords internal
#' @export
#' @rdname radr_folder
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

radr_folder <- function(rad.folder, path.folder = NULL, prefix.int = TRUE) {
  if (is.null(path.folder)) path.folder <- getwd()
  if (prefix.int) {
    existing.dirs <- list.dirs(path = path.folder, full.names = FALSE, recursive = FALSE)
    if (length(existing.dirs) > 0) {
      check <- existing.dirs %>%
        stringi::stri_extract_first_regex(str = ., pattern = "^[0-9]*_") %>%
        stringi::stri_remove_na(x = .)

      if (length(check) > 0L) {
        check %<>%
          stringi::stri_extract_first_regex(
            str = .,
            pattern = "\\d{2}"
          ) %>%
          stringi::stri_replace_first_regex(
            str = .,
            pattern = "^[0]",
            # pattern = "0",
            replacement = ""
          ) %>%
          as.integer(x = .) %>%
          sort

        last.num <- utils::tail(x = check, 1)
        check.num <- length(check)
        if (identical(check.num, last.num)) {
          select.last <- max(check.num, last.num, na.rm = TRUE)
        } else {
          select.last <- check.num <- last.num
        }
        existing.dirs <- select.last + 1L
      } else {
        existing.dirs <- 1L
      }
    } else {
      existing.dirs <- 1L
    }

    rad.folder <- existing.dirs %>%
      as.character %>%
      stringi::stri_pad_left(str = ., width = 2, pad = 0) %>%
      stringi::stri_join(., "_", rad.folder)
  }
  radr.folder.full.path <- file.path(path.folder, rad.folder)
  return(radr.folder.full.path)
}#End radr_folder


# generate_squeleton_folders----------------------------------------------------
#' @title generate_squeleton_folders
#' @description Generate squeleton folders
#' @keywords internal
#' @export
generate_squeleton_folders <- function(
  fp = 0L,
  path.folder = NULL,
  interactive.filter = TRUE,
  ...
) {

  # test
  # fp = 0L
  # file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  # interactive.filter = TRUE

  if (is.null(path.folder)) path.folder <- getwd()
  folders.labels <- c(
    "filter_dart_reproducibility",
    "filter_individuals", "filter_individuals", "filter_individuals",
    "filter_common_markers",
    "filter_ma",
    "filter_coverage",
    "filter_genotyping",
    "filter_snp_position_read",
    "filter_snp_number",
    "filter_ld", "filter_ld",
    "detect_mixed_genomes",
    "detect_duplicate_genomes",
    "filter_hwe")

  if (!interactive.filter) {
    get.filters <- ls(envir = as.environment(1))
    need <- c(
      "filter.reproducibility",
      "filter.individuals.missing",
      "filter.individuals.heterozygosity",
      "filter.individuals.coverage.total",
      "filter.common.markers",
      "filter.ma",
      "filter.coverage",
      "filter.genotyping",
      "filter.snp.position.read",
      "filter.snp.number",
      "filter.short.ld",
      "filter.long.ld",
      "detect.mixed.genomes",
      "detect.duplicate.genomes",
      "filter.hwe")
    folders <- purrr::keep(.x = get.filters, .p = get.filters %in% need)
    wanted_filters <- function(x) {
      !is.null(rlang::eval_tidy(rlang::parse_expr(x)))
    }
    folders <- purrr::keep(.x = folders, .p = wanted_filters)
    folders <- factor(
      x = folders,
      levels = need,
      labels = folders.labels,
      ordered = TRUE
    ) %>%
      droplevels(.) %>%
      unique %>%
      sort %>%
      as.character
  } else {
    folders <- unique(folders.labels)
  }

  folders <- c("radr", folders)

  res <- list()
  fp.loop <- fp
  temp <- NULL
  for (f in folders) {
    # message("Processing: ", f)
    temp <- folder_prefix(
      prefix.int = fp.loop,
      prefix.name = f,
      path.folder = path.folder)
    res[[f]] <- temp$folder.prefix
    fp.loop <- temp$prefix.int
  }
  return(res)
}#End generate_squeleton_folders
# generate_filename-------------------------------------------------------------
#' @title Filename radr
#' @description Generate a filename object
#' @name generate_filename
#' @rdname generate_filename
#' @keywords internal
#' @export
generate_filename <- function(
  name.shortcut = NULL,
  path.folder = NULL,
  date = TRUE,
  extension = c(
    "tsv", "gds.rad", "rad", "gds", "gen", "dat",
    "genind", "genlight", "gtypes", "vcf", "colony",
    "bayescan", "gsisim", "hierfstat", "hzar", "ldna",
    "pcadapt", "related", "stockr", "structure", "arlequin",
    "arrow.parquet"
  )
) {

  if (is.null(path.folder)) path.folder <- getwd()

  # date and time-
  if (is.character(date)) {
    file.date <- stringi::stri_join("_", date)
  } else if (date) {
    file.date <- stringi::stri_join("_", format(Sys.time(), "%Y%m%d@%H%M"))
  } else {
    file.date <- ""
  }

  # path.folder
  if (!dir.exists(path.folder)) dir.create(path.folder)

  # Extension
  want <- c("tsv", "gds.rad", "rad", "gds", "gen", "dat", "genind", "genlight", "gtypes",
            "vcf", "colony", "bayescan", "gsisim", "hierfstat", "hzar", "ldna",
            "pcadapt", "plink", "related", "stockr", "structure", "arlequin",
            "arrow.parquet")
  extension <- match.arg(extension, want)

  # note to myself: currently excluded output :
  # "fineradstructure", "maverick", "plink", "betadiv"


  # with same extension
  # extension <- "tsv"
  if (extension %in% c("tsv", "gds.rad", "rad", "gds", "vcf", "colony", "ldna",
                       "arrow.parquet")) {
    extension <- stringi::stri_join(file.date, ".", extension)
  }


  # Radr saveRDS
  # extension <- "genind"
  if (extension %in% c("genind", "genlight", "gtypes", "stockr")) {
    extension <- stringi::stri_join("_", extension, file.date, ".RData")
  }

  # Radr tsv
  if (extension %in% c("tsv")) {
    extension <- stringi::stri_join("_", extension, file.date, ".tsv")
  }

  # Radr txt
  if (extension %in% c("bayescan", "pcadapt", "related")) {
    extension <- stringi::stri_join("_", extension, file.date, ".txt")
  }

  # Radr csv
  if (extension %in% c("hzar", "arlequin")) {
    extension <- stringi::stri_join("_", extension, file.date, ".csv")
  }

  # custom
  if (extension == "gen") extension <- stringi::stri_join("_genepop", file.date, ".gen")
  if (extension == "dat") extension <- stringi::stri_join("_fstat", file.date, ".dat")
  if (extension == "hierfstat") extension <- stringi::stri_join("_hierfstat", file.date, ".dat")
  if (extension == "structure") extension <- stringi::stri_join("_structure", file.date, ".str")

  # Filename
  if (is.null(name.shortcut)) {
    filename <- stringi::stri_join("radr", extension)
  } else {
    filename.problem <- file.exists(stringi::stri_join(name.shortcut, extension))
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_radr", extension)
    } else {
      filename <- stringi::stri_join(name.shortcut, extension)
    }
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join("duplicated_", filename)
    }
  }


  # Include path.folder in returned object
  return(res = list(filename.short = filename, filename = file.path(path.folder, filename)))
}#End generate_filename


# generate_folder---------------------------------------------------------------

#' @title generate_folder
#' @description Generate a folder based on ...
#' @name generate_folder
#' @param rad.folder Name of the rad folder
#' @param internal (optional, logical) Is the function internal or not
#' @param append.date Include the date and time with the folder.
#' Default: \code{append.date = TRUE}.
#' @param file.date The file date included in as argument/value or with
#' default \code{file.date = NULL}, generated by the fucntion.
#' @inheritParams radr_folder
#' @inheritParams radr_common_arguments
#' @keywords internal
#' @export
#' @rdname generate_folder
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

generate_folder <- function(
    rad.folder = NULL,
    path.folder = NULL,
    internal = FALSE,
    append.date = TRUE,
    file.date = NULL,
    prefix.int = TRUE,
    verbose = FALSE
) {

  if (internal) {
    rad.folder <- NULL
    f.temp <- radr.folder.full.path <- path.folder
  }

  if (!is.null(rad.folder)) {
    f.temp <- radr.folder.full.path <- radr_folder(
      rad.folder = rad.folder,
      path.folder = path.folder,
      prefix.int = prefix.int
      )
  }



  if (is.null(file.date)) {
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")# Date and time
  }

  if (is.null(radr.folder.full.path)) {
    radr.folder.full.path <- getwd()
  } else {
    #working directory in the path?
    wd.present <- TRUE %in% unique(stringi::stri_detect_fixed(str = radr.folder.full.path, pattern = c(getwd(), paste0(getwd(), "/"))))
    date.present <- TRUE %in% unique(stringi::stri_detect_fixed(str = radr.folder.full.path, pattern = "@"))
    if (!date.present && append.date) radr.folder.full.path <- stringi::stri_join(radr.folder.full.path, file.date, sep = "_")
    if (!wd.present) radr.folder.full.path <- file.path(getwd(), radr.folder.full.path)
    if (verbose && !identical(f.temp, radr.folder.full.path)) message("Folder created: ", basename(radr.folder.full.path))
  }
  if (!dir.exists(radr.folder.full.path)) dir.create(radr.folder.full.path)
  return(radr.folder.full.path)
}#End generate_folder


# folder_prefix-----------------------------------------------------------------
#' @title folder_prefix
#' @description Generate a seq and folder prefix
#' @name folder_prefix
#' @rdname folder_prefix
#' @keywords internal
#' @export
folder_prefix <- function(
  prefix.int = NULL,
  prefix.name = NULL,
  path.folder = NULL
) {
  if (is.null(path.folder)) {
    path.folder <- getwd()
  } else {
    if (stringi::stri_sub(str = path.folder, from = -1, length = 1) == "/") {
      path.folder <- stringi::stri_replace_last_regex(
        str = path.folder,
        pattern = "[/$]",
        replacement = "")
    }
  }

  if (is.null(prefix.int)) {
    prefix.int <- 0L
  } else {
    if (is.list(prefix.int)) {
      prefix.int <- as.integer(prefix.int$prefix.int) + 1L
    } else {
      prefix.int <- as.integer(prefix.int) + 1L
    }
  }

  if (is.null(prefix.name)) {
    folder.prefix <- stringi::stri_join(
      stringi::stri_pad_left(
        str = prefix.int, width = 2, pad = 0
      ), "_"
    )
  } else {
    folder.prefix <- stringi::stri_join(
      stringi::stri_pad_left(
        str = prefix.int, width = 2, pad = 0
      ),
      prefix.name,
      sep = "_"
    )
  }
  folder.prefix <- file.path(path.folder, folder.prefix)
  res = list(prefix.int = prefix.int, folder.prefix = folder.prefix)
}#End folder_prefix


# radr_snakecase------------------------------------------------------------
#' @title radr_snakecase
#' @description Transform CamelCase to SCREAMING snake_cases
#' @name radr_snakecase
#' @rdname radr_snakecase
#' @keywords internal
#' @export
radr_snakecase <- function(x) {
  x |>
    stringi::stri_replace_all_regex("([A-Za-z])([A-Z])([a-z])", "$1_$2$3") |>
    stringi::stri_replace_all_fixed(".", "_") |>
    stringi::stri_replace_all_regex("([a-z])([A-Z])", "$1_$2") |>
    stringi::stri_trans_toupper()


}#End radr_snakecase


# radr_packages_dep---------------------------------------------------------
#' @title radr_packages_dep
#' @description Verify required packages
#' @rdname radr_packages_dep
#' @keywords internal
#' @export
radr_packages_dep <- function(package, cran = TRUE, bioc = FALSE) {
  .Deprecated("tgbase::check_package", package = "radr")
  tgbase::check_package(package = package, cran = cran, bioc = bioc)
}#End radr_packages_dep

# tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)
# requireNamespace
# installed.packages


# radr_clock ---------------------------------------------------------------
#' @title radr_tic
#' @description radr tictoc function
#' @rdname radr_tic
#' @keywords internal
#' @export
radr_tic <- function(timing = proc.time()) {
    invisible(timing)
}# End radr_tic

#' @title radr_toc
#' @description radr tictoc function
#' @rdname radr_toc
#' @keywords internal
#' @export
radr_toc <- function(
  timing,
  end.message = "Computation time, overall:",
  verbose = TRUE
) {
  if (verbose) message("\n", end.message, " ", round((proc.time() - timing)[[3]]), " sec")
}# End radr_toc
