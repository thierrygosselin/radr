#' Filter individuals and markers using ID lists
#'
#' Apply any combination of individual and marker lists to a prepared genome
#' GDS. A whitelist keeps matching records from the current selection. A
#' blacklist then excludes matching records. Earlier filter reasons are retained.
#'
#' @md
#' @param data A writable genometranslator GDS connection or path to a GDS file.
#' @param whitelist.individuals,blacklist.individuals Character vector of sample
#'   IDs or a TSV with an `INDIVIDUALS` column. IDs are matched exactly,
#'   including leading zeros.
#' @param whitelist.markers,blacklist.markers Character vector of `MARKERS`
#'   names, numeric vector of GDS `VARIANT_ID` values, or a TSV with either
#'   column. When both columns are supplied they must identify the same rows.
#'   Variant IDs refer to the GDS from which they were exported.
#' @details All lists default to `NULL`. Duplicate IDs are counted once.
#'   Unknown or ambiguous IDs and selections that leave no individuals or
#'   markers raise an error before filtering. Only active, previously
#'   whitelisted records can acquire the label `filter.lists`. Genotypes are
#'   not rescanned. Successive calls narrow the current selection.
#'
#'   Each run writes its input lists, the cumulative filter history and the
#'   standard radr audit of retained and newly excluded metadata records to a
#'   numbered, dated folder in the working directory. To compare independent
#'   scenarios, restore all list exclusions with
#'   `genometranslator::reset_filters(gds, filter.lists = TRUE)` between runs.
#' @return The open GDS, with its metadata and active selection updated in
#'   place. Close the returned connection when finished. With four NULL lists,
#'   returns the supplied object or path as is.
#' @examples
#' \dontrun{
#' gds <- filter_lists("swordfish.gds", whitelist.individuals = adult_ids,
#'                     blacklist.markers = "inversion_markers.tsv")
#' genometranslator::reset_filters(gds, filter.lists = TRUE)
#' gds <- filter_lists(gds, whitelist.individuals = adult_ids)
#' }
#' @seealso [filter_whitelist()], [genometranslator::list_filters()]
#' @export
filter_lists <- function(
    data = NULL,
    whitelist.individuals = NULL,
    blacklist.individuals = NULL,
    whitelist.markers = NULL,
    blacklist.markers = NULL
) {
  if (is.null(data)) rlang::abort("Supply a genometranslator GDS in `data`.")
  inputs <- list(
    whitelist.individuals = whitelist.individuals,
    blacklist.individuals = blacklist.individuals,
    whitelist.markers = whitelist.markers,
    blacklist.markers = blacklist.markers
  )
  inputs <- purrr::discard(inputs, is.null)
  if (!length(inputs)) return(data)

  .start <- tgbase::startup(package = "radr", f.name = "filter_lists")
  on.exit(tgbase::teardown(.start), add = TRUE)
  context <- .filter_gds_open(data)
  gds <- context$gds
  completed <- FALSE
  on.exit(if (context$opened && !completed) SeqArray::seqClose(gds), add = TRUE)
  if (isTRUE(gds$readonly)) rlang::abort("Open the GDS with write access.")
  root <- if (!is.null(gdsfmt::index.gdsn(gds, "genometranslator", silent = TRUE))) {
    "genometranslator"
  } else "radiator"
  for (field in c("individuals.meta/INDIVIDUALS", "markers.meta/VARIANT_ID")) {
    if (is.null(gdsfmt::index.gdsn(gds, paste(root, field, sep = "/"), silent = TRUE))) {
      rlang::abort(paste("Missing genome metadata field", field))
    }
  }

  # Resolve all lists and validate both selections before creating outputs.
  inputs <- purrr::imap(inputs, .read_filter_list)
  individuals <- genometranslator::extract_individuals_metadata(
    gds, ind.field.select = c("INDIVIDUALS", "FILTERS")
  )
  markers <- genometranslator::extract_markers_metadata(
    gds, markers.meta.select = c("VARIANT_ID", "MARKERS", "FILTERS")
  )
  active.i <- SeqArray::seqGetData(gds, "sample.id")
  active.m <- SeqArray::seqGetData(gds, "variant.id")
  samples <- .plan_filter_list(individuals, "INDIVIDUALS", active.i,
                               inputs$whitelist.individuals, inputs$blacklist.individuals)
  variants <- .plan_filter_list(markers, "VARIANT_ID", active.m,
                                inputs$whitelist.markers, inputs$blacklist.markers)

  path.folder <- tgbase::generate_folder("filter_lists", prefix.int = TRUE,
                                         file.date = .start$file.date, verbose = TRUE)
  purrr::iwalk(inputs, function(ids, name) {
    tgbase::write_tgbase_tsv(ids, path.folder = path.folder, filename = name, date = FALSE)
  })
  parameters <- filter_parameters(generate = TRUE, initiate = TRUE, update = FALSE,
                                  data = gds, path.folder = path.folder, file.date = .start$file.date)
  # Keep both metadata tables, the active selection and history consistent
  # if writing a node or the identifier-level audit fails.
  touched <- character()
  previous <- list("individuals.meta/FILTERS" = as.character(individuals$FILTERS),
                   "markers.meta/FILTERS" = as.character(markers$FILTERS))
  tryCatch({
    for (plan in list(samples, variants)) {
      if (any(plan$remove)) {
        touched <- c(touched, plan$node)
        labels <- previous[[plan$node]]
        labels[plan$remove] <- "filter.lists"
        genometranslator::update_genome_gds(gds, node.name = plan$node, value = labels)
      }
    }
    genometranslator::sync_gds(gds, samples = samples$keep, variant.id = variants$keep)
    result <- filter_parameters(
      update = TRUE, parameter.obj = parameters, data = gds, path.folder = path.folder,
      filter.name = "filter_lists", param.name = paste(names(inputs), collapse = " / "),
      values = paste(vapply(inputs, nrow, integer(1)), collapse = " / "),
      comments = basename(path.folder)
    )
    gdsfmt::sync.gds(gds)
  }, error = function(e) {
    for (node in touched) {
      genometranslator::update_genome_gds(gds, node.name = node, value = previous[[node]])
    }
    genometranslator::sync_gds(gds, samples = active.i, variant.id = active.m)
    genometranslator::update_genome_gds(gds, node.name = "filter_history",
                                        value = parameters$filter.history)
    unlink(path.folder, recursive = TRUE)
    stop(e)
  })
  completed <- TRUE
  radr_results_message("Filter individuals and markers by list", result)
  gds
}

.read_filter_list <- function(x, name) {
  keys <- if (endsWith(name, ".individuals")) "INDIVIDUALS" else c("VARIANT_ID", "MARKERS")
  if (is.character(x) && length(x) == 1L && !is.na(x) && file.exists(x)) {
    ids <- readr::read_tsv(x, col_types = readr::cols(.default = "c"), na = "",
                           trim_ws = FALSE, name_repair = "check_unique", progress = FALSE)
    readr::stop_for_problems(ids)
    ids <- dplyr::select(ids, dplyr::any_of(keys))
    if (!ncol(ids)) rlang::abort(paste(name, "requires a column named", paste(keys, collapse = " or ")))
  } else {
    valid <- is.character(x) || ("VARIANT_ID" %in% keys && is.numeric(x))
    if (!valid || !is.null(dim(x))) rlang::abort(paste("Invalid vector or TSV path for", name))
    key <- if (is.character(x)) utils::tail(keys, 1L) else "VARIANT_ID"
    ids <- tibble::tibble(!!key := as.character(x))
  }
  if (anyNA(ids) || any(vapply(ids, function(id) any(!nzchar(trimws(id))), logical(1)))) {
    rlang::abort(paste("Missing or empty ID in", name))
  }
  dplyr::distinct(ids)
}

.plan_filter_list <- function(metadata, key, active, whitelist, blacklist) {
  match.rows <- function(ids) {
    if (is.null(ids)) return(integer())
    matches <- purrr::imap(ids, function(id, field) {
      if (!field %in% names(metadata)) rlang::abort(paste("Missing metadata column", field))
      values <- as.character(metadata[[field]])
      rows <- match(id, values)
      if (anyNA(rows)) {
        rlang::abort(paste("Unknown", field, "in list:", paste(utils::head(id[is.na(rows)]), collapse = ", ")))
      }
      if (any(id %in% values[duplicated(values)])) rlang::abort(paste("Ambiguous", field, "in list."))
      rows
    })
    if (length(matches) == 2L && !identical(matches[[1L]], matches[[2L]])) {
      rlang::abort("VARIANT_ID and MARKERS identify different markers in the TSV.")
    }
    matches[[1L]]
  }
  include <- match.rows(whitelist)
  exclude <- match.rows(blacklist)
  eligible <- metadata$FILTERS %in% "whitelist" & metadata[[key]] %in% active
  rows <- seq_len(nrow(metadata))
  remove <- eligible & ((!is.null(whitelist) & !rows %in% include) | rows %in% exclude)
  keep <- metadata[[key]][eligible & !remove]
  if (!length(keep)) rlang::abort(paste("No active", key, "would remain."))
  node <- if (key == "INDIVIDUALS") "individuals.meta/FILTERS" else "markers.meta/FILTERS"
  list(keep = keep, remove = remove, node = node)
}
